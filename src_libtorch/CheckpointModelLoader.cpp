#include "CheckpointModelLoader.h"

#include <algorithm>
#include <array>
#include <cctype>
#include <fstream>
#include <iomanip>
#include <sstream>
#include <stdexcept>

namespace {

uint16_t read_u16(const char* data) {
  const auto* bytes = reinterpret_cast<const unsigned char*>(data);
  return static_cast<uint16_t>(bytes[0]) |
         (static_cast<uint16_t>(bytes[1]) << 8);
}

uint32_t read_u32(const char* data) {
  const auto* bytes = reinterpret_cast<const unsigned char*>(data);
  return static_cast<uint32_t>(bytes[0]) |
         (static_cast<uint32_t>(bytes[1]) << 8) |
         (static_cast<uint32_t>(bytes[2]) << 16) |
         (static_cast<uint32_t>(bytes[3]) << 24);
}

std::string format_bytes(uint64_t bytes) {
  constexpr std::array<const char*, 5> units = {"B", "KiB", "MiB", "GiB",
                                                "TiB"};
  auto value = static_cast<double>(bytes);
  std::size_t unit = 0;
  while (value >= 1024.0 && unit + 1 < units.size()) {
    value /= 1024.0;
    ++unit;
  }

  std::ostringstream out;
  out << std::fixed << std::setprecision(unit == 0 ? 0 : 2) << value << ' '
      << units[unit];
  return out.str();
}

std::string format_shape(const std::vector<int64_t>& sizes) {
  std::ostringstream out;
  out << '[';
  for (std::size_t i = 0; i < sizes.size(); ++i) {
    if (i > 0) {
      out << ", ";
    }
    out << sizes[i];
  }
  out << ']';
  return out.str();
}

CheckpointModelLoader::TensorInfo make_tensor_info(const std::string& name,
                                                   const torch::Tensor& tensor) {
  CheckpointModelLoader::TensorInfo info;
  info.name = name;
  info.dtype = c10::toString(tensor.scalar_type());
  info.device = tensor.device().str();
  info.sizes = tensor.sizes().vec();
  info.elements = tensor.numel();
  info.requires_grad = tensor.requires_grad();
  return info;
}

void keep_largest_entries(
    std::vector<CheckpointModelLoader::ArchiveEntry>& entries,
    CheckpointModelLoader::ArchiveEntry entry) {
  entries.push_back(std::move(entry));
  std::sort(entries.begin(), entries.end(), [](const auto& lhs,
                                               const auto& rhs) {
    return lhs.uncompressed_size > rhs.uncompressed_size;
  });
  if (entries.size() > 8) {
    entries.pop_back();
  }
}

bool looks_like_tensor_name(const std::string& value) {
  return value.find('.') != std::string::npos &&
         (value.find("weight") != std::string::npos ||
          value.find("bias") != std::string::npos ||
          value.find("encoder") != std::string::npos ||
          value.find("decoder") != std::string::npos ||
          value.find("transformer") != std::string::npos ||
          value.find("norm") != std::string::npos);
}

bool is_readable_text(const std::string& value) {
  return std::all_of(value.begin(), value.end(), [](const char ch) {
    const auto byte = static_cast<unsigned char>(ch);
    return std::isprint(byte) != 0;
  });
}

std::vector<std::string> extract_pickle_strings(const std::vector<char>& data) {
  std::vector<std::string> values;
  const auto append = [&](std::size_t offset, std::size_t length) {
    if (length == 0 || length > 4096 || offset + length > data.size()) {
      return;
    }
    std::string value(data.data() + offset, data.data() + offset + length);
    if (is_readable_text(value)) {
      values.push_back(std::move(value));
    }
  };

  for (std::size_t i = 0; i < data.size(); ++i) {
    const auto opcode = static_cast<unsigned char>(data[i]);
    if (opcode == 'X' && i + 5 <= data.size()) {
      const auto length = read_u32(data.data() + i + 1);
      append(i + 5, length);
      i += 4;
    } else if ((opcode == 'U' || opcode == 0x8c) && i + 2 <= data.size()) {
      const auto length = static_cast<unsigned char>(data[i + 1]);
      append(i + 2, length);
      i += 1;
    } else if (opcode == 'T' && i + 5 <= data.size()) {
      const auto length = read_u32(data.data() + i + 1);
      append(i + 5, length);
      i += 4;
    }
  }

  return values;
}

bool is_checkpoint_key(const std::string& value) {
  static const std::array<std::string, 9> keys = {
      "epoch",     "global_step", "pytorch-lightning_version",
      "state_dict", "optimizer_states", "lr_schedulers",
      "callbacks", "hyper_parameters", "loops"};
  return std::find(keys.begin(), keys.end(), value) != keys.end();
}

void append_unique_limited(std::vector<std::string>& values, std::string value,
                           std::size_t limit) {
  if (values.size() >= limit) {
    return;
  }
  if (std::find(values.begin(), values.end(), value) == values.end()) {
    values.push_back(std::move(value));
  }
}

}  // namespace

bool CheckpointModelLoader::load(const std::filesystem::path& checkpoint_path,
                                 const c10::Device& device) {
  reset();
  checkpoint_path_ = checkpoint_path;

  std::error_code error;
  const auto file_size = std::filesystem::file_size(checkpoint_path_, error);
  if (error) {
    last_error_ = "Cannot read checkpoint file size: " + error.message();
    return false;
  }
  file_size_ = file_size;

  try {
    torch::NoGradGuard no_grad;
    module_.emplace(torch::jit::load(checkpoint_path_.string(), device));
    module_->eval();
    collect_module_summary();
    inspect_zip_archive(checkpoint_path_);
    return true;
  } catch (const c10::Error& error) {
    last_error_ = error.what_without_backtrace();
  } catch (const std::exception& error) {
    last_error_ = error.what();
  }

  inspect_zip_archive(checkpoint_path_);
  if (archive_.has_pickle && !archive_.has_torchscript_code) {
    last_error_ =
        "not a TorchScript export; found Python checkpoint archive/data.pkl";
  }
  return archive_.is_zip;
}

bool CheckpointModelLoader::loaded() const {
  return module_.has_value() || archive_.is_zip;
}

bool CheckpointModelLoader::has_torchscript_module() const {
  return module_.has_value();
}

const std::string& CheckpointModelLoader::last_error() const {
  return last_error_;
}

const torch::jit::script::Module& CheckpointModelLoader::module() const {
  if (!module_) {
    throw std::logic_error("No TorchScript module has been loaded.");
  }
  return *module_;
}

bool CheckpointModelLoader::is_python_checkpoint() const {
  return archive_.is_zip && archive_.has_pickle && !archive_.has_torchscript_code;
}

const std::vector<std::string>& CheckpointModelLoader::checkpoint_keys() const {
  return archive_.checkpoint_keys;
}

const std::vector<std::string>& CheckpointModelLoader::sample_state_dict_names()
    const {
  return archive_.sample_tensor_names;
}

const std::optional<std::string>&
CheckpointModelLoader::pytorch_lightning_version() const {
  return archive_.pytorch_lightning_version;
}

std::string CheckpointModelLoader::summary(const std::size_t max_items) const {
  std::ostringstream out;
  out << "Checkpoint: " << checkpoint_path_.string() << '\n';
  out << "File size: " << format_bytes(file_size_) << '\n';

  if (module_) {
    out << "Format: TorchScript module\n";
    out << "Methods: ";
    if (method_names_.empty()) {
      out << "(none reported)";
    } else {
      for (std::size_t i = 0; i < method_names_.size(); ++i) {
        if (i > 0) {
          out << ", ";
        }
        out << method_names_[i];
      }
    }
    out << '\n';
    out << "Parameters: " << parameters_.size() << " tensors, "
        << parameter_count_ << " values";
    if (trainable_parameter_count_ != parameter_count_) {
      out << " (" << trainable_parameter_count_ << " trainable)";
    }
    out << '\n';
    out << "Buffers: " << buffers_.size() << " tensors, " << buffer_count_
        << " values\n";

    const auto print_tensors = [&](const char* label,
                                   const std::vector<TensorInfo>& tensors) {
      out << label << ":\n";
      if (tensors.empty()) {
        out << "  (none)\n";
        return;
      }
      const auto count = std::min(max_items, tensors.size());
      for (std::size_t i = 0; i < count; ++i) {
        const auto& tensor = tensors[i];
        out << "  " << tensor.name << " " << format_shape(tensor.sizes) << ' '
            << tensor.dtype << ' ' << tensor.device << " elements="
            << tensor.elements;
        if (tensor.requires_grad) {
          out << " trainable";
        }
        out << '\n';
      }
      if (tensors.size() > count) {
        out << "  ... " << (tensors.size() - count) << " more\n";
      }
    };

    print_tensors("Parameter tensors", parameters_);
    print_tensors("Buffer tensors", buffers_);
    return out.str();
  }

  if (archive_.is_zip) {
    out << "Format: PyTorch zip checkpoint";
    if (archive_.has_pickle) {
      out << " with archive/data.pkl";
    }
    out << '\n';
    out << "TorchScript code: "
        << (archive_.has_torchscript_code ? "present" : "not found") << '\n';
    out << "Archive entries: " << archive_.entry_count << '\n';
    out << "Tensor storage entries: " << archive_.tensor_storage_count << '\n';
    out << "Compressed payload: " << format_bytes(archive_.compressed_bytes)
        << '\n';
    out << "Uncompressed payload: "
        << format_bytes(archive_.uncompressed_bytes) << '\n';
    if (!last_error_.empty()) {
      out << "LibTorch module load: failed (" << last_error_ << ")\n";
    }
    if (!archive_.checkpoint_keys.empty()) {
      out << "Checkpoint keys:";
      for (const auto& key : archive_.checkpoint_keys) {
        out << ' ' << key;
      }
      out << '\n';
    }
    if (archive_.pytorch_lightning_version) {
      out << "PyTorch-Lightning version: "
          << *archive_.pytorch_lightning_version << '\n';
    }
    if (!archive_.sample_tensor_names.empty()) {
      out << "Sample state_dict names:\n";
      for (const auto& name : archive_.sample_tensor_names) {
        out << "  " << name << '\n';
      }
    }
    out << "Largest archive entries:\n";
    for (const auto& entry : archive_.largest_entries) {
      out << "  " << entry.name << " "
          << format_bytes(entry.uncompressed_size) << '\n';
    }
    out << "Note: Python torch.save or Lightning .ckpt files need the original "
           "Python model class to restore weights. Export with torch.jit.save "
           "to run the model directly from LibTorch.\n";
    return out.str();
  }

  out << "Format: unknown\n";
  if (!last_error_.empty()) {
    out << "Error: " << last_error_ << '\n';
  }
  return out.str();
}

void CheckpointModelLoader::reset() {
  checkpoint_path_.clear();
  module_.reset();
  parameters_.clear();
  buffers_.clear();
  method_names_.clear();
  archive_ = ArchiveSummary{};
  last_error_.clear();
  file_size_ = 0;
  parameter_count_ = 0;
  trainable_parameter_count_ = 0;
  buffer_count_ = 0;
}

void CheckpointModelLoader::collect_module_summary() {
  if (!module_) {
    return;
  }

  for (const auto& parameter : module_->named_parameters(true)) {
    auto info = make_tensor_info(parameter.name, parameter.value);
    parameter_count_ += info.elements;
    if (info.requires_grad) {
      trainable_parameter_count_ += info.elements;
    }
    parameters_.push_back(std::move(info));
  }

  for (const auto& buffer : module_->named_buffers(true)) {
    auto info = make_tensor_info(buffer.name, buffer.value);
    buffer_count_ += info.elements;
    buffers_.push_back(std::move(info));
  }

  for (const auto& method : module_->get_methods()) {
    method_names_.push_back(method.name());
  }
}

bool CheckpointModelLoader::inspect_zip_archive(
    const std::filesystem::path& checkpoint_path) {
  std::ifstream input(checkpoint_path, std::ios::binary);
  if (!input) {
    return false;
  }

  input.seekg(0, std::ios::end);
  const auto file_size = static_cast<uint64_t>(input.tellg());
  const auto tail_size = static_cast<std::size_t>(
      std::min<uint64_t>(file_size, 65557));
  std::vector<char> tail(tail_size);
  input.seekg(static_cast<std::streamoff>(file_size - tail_size),
              std::ios::beg);
  input.read(tail.data(), static_cast<std::streamsize>(tail.size()));

  int64_t eocd_offset = -1;
  for (int64_t i = static_cast<int64_t>(tail.size()) - 22; i >= 0; --i) {
    if (read_u32(tail.data() + i) == 0x06054b50) {
      eocd_offset = i;
      break;
    }
  }

  if (eocd_offset < 0) {
    return false;
  }

  const char* eocd = tail.data() + eocd_offset;
  const auto entry_count = read_u16(eocd + 10);
  const auto central_directory_size = read_u32(eocd + 12);
  const auto central_directory_offset = read_u32(eocd + 16);
  if (entry_count == 0 || central_directory_size == 0) {
    archive_.is_zip = true;
    return true;
  }

  input.seekg(static_cast<std::streamoff>(central_directory_offset),
              std::ios::beg);
  std::vector<char> directory(central_directory_size);
  input.read(directory.data(), static_cast<std::streamsize>(directory.size()));
  if (!input) {
    return false;
  }

  archive_ = ArchiveSummary{};
  archive_.is_zip = true;

  std::optional<ArchiveEntry> pickle_entry;
  std::size_t cursor = 0;
  for (uint16_t i = 0; i < entry_count && cursor + 46 <= directory.size();
       ++i) {
    const char* header = directory.data() + cursor;
    if (read_u32(header) != 0x02014b50) {
      break;
    }

    const auto compressed_size = read_u32(header + 20);
    const auto uncompressed_size = read_u32(header + 24);
    const auto name_length = read_u16(header + 28);
    const auto extra_length = read_u16(header + 30);
    const auto comment_length = read_u16(header + 32);
    const auto compression_method = read_u16(header + 10);
    const auto local_header_offset = read_u32(header + 42);
    const auto entry_size =
        static_cast<std::size_t>(46 + name_length + extra_length +
                                 comment_length);
    if (cursor + entry_size > directory.size()) {
      break;
    }

    std::string name(header + 46, header + 46 + name_length);
    ++archive_.entry_count;
    archive_.compressed_bytes += compressed_size;
    archive_.uncompressed_bytes += uncompressed_size;
    if (name == "archive/data.pkl") {
      archive_.has_pickle = true;
    }
    if (name.rfind("archive/code/", 0) == 0) {
      archive_.has_torchscript_code = true;
    }
    if (name.rfind("archive/data/", 0) == 0) {
      ++archive_.tensor_storage_count;
    }
    ArchiveEntry entry{name, compressed_size, uncompressed_size,
                       local_header_offset, compression_method};
    if (name == "archive/data.pkl") {
      pickle_entry = entry;
    }
    keep_largest_entries(archive_.largest_entries, entry);

    cursor += entry_size;
  }

  if (pickle_entry && pickle_entry->compression_method == 0) {
    input.seekg(static_cast<std::streamoff>(pickle_entry->local_header_offset),
                std::ios::beg);
    std::array<char, 30> local_header{};
    input.read(local_header.data(),
               static_cast<std::streamsize>(local_header.size()));
    if (input && read_u32(local_header.data()) == 0x04034b50) {
      const auto local_name_length = read_u16(local_header.data() + 26);
      const auto local_extra_length = read_u16(local_header.data() + 28);
      input.seekg(static_cast<std::streamoff>(local_name_length +
                                             local_extra_length),
                  std::ios::cur);
      std::vector<char> pickle_data(pickle_entry->uncompressed_size);
      input.read(pickle_data.data(),
                 static_cast<std::streamsize>(pickle_data.size()));
      if (input) {
        const auto strings = extract_pickle_strings(pickle_data);
        for (std::size_t i = 0; i < strings.size(); ++i) {
          const auto& value = strings[i];
          if (is_checkpoint_key(value)) {
            append_unique_limited(archive_.checkpoint_keys, value, 16);
          }
          if (value == "pytorch-lightning_version" &&
              i + 1 < strings.size()) {
            archive_.pytorch_lightning_version = strings[i + 1];
          }
          if (looks_like_tensor_name(value)) {
            append_unique_limited(archive_.sample_tensor_names, value, 16);
          }
        }
      }
    }
  }

  return true;
}
