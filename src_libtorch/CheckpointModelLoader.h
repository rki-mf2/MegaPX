#pragma once

#include <torch/script.h>

#include <cstdint>
#include <filesystem>
#include <optional>
#include <string>
#include <vector>

class CheckpointModelLoader {
 public:
  struct TensorInfo {
    std::string name;
    std::string dtype;
    std::string device;
    std::vector<int64_t> sizes;
    int64_t elements = 0;
    bool requires_grad = false;
  };

  struct ArchiveEntry {
    std::string name;
    uint64_t compressed_size = 0;
    uint64_t uncompressed_size = 0;
    uint64_t local_header_offset = 0;
    uint16_t compression_method = 0;
  };

  bool load(const std::filesystem::path& checkpoint_path,
            const c10::Device& device = c10::Device(c10::kCPU));

  bool loaded() const;
  bool has_torchscript_module() const;
  const std::string& last_error() const;
  std::string summary(std::size_t max_items = 12) const;

  const torch::jit::script::Module& module() const;
  bool is_python_checkpoint() const;
  const std::vector<std::string>& checkpoint_keys() const;
  const std::vector<std::string>& sample_state_dict_names() const;
  const std::optional<std::string>& pytorch_lightning_version() const;

 private:
  struct ArchiveSummary {
    bool is_zip = false;
    bool has_pickle = false;
    bool has_torchscript_code = false;
    uint64_t entry_count = 0;
    uint64_t tensor_storage_count = 0;
    uint64_t compressed_bytes = 0;
    uint64_t uncompressed_bytes = 0;
    std::vector<ArchiveEntry> largest_entries;
    std::vector<std::string> checkpoint_keys;
    std::vector<std::string> sample_tensor_names;
    std::optional<std::string> pytorch_lightning_version;
  };

  void reset();
  void collect_module_summary();
  bool inspect_zip_archive(const std::filesystem::path& checkpoint_path);

  std::filesystem::path checkpoint_path_;
  std::optional<torch::jit::script::Module> module_;
  std::vector<TensorInfo> parameters_;
  std::vector<TensorInfo> buffers_;
  std::vector<std::string> method_names_;
  ArchiveSummary archive_;
  std::string last_error_;
  uint64_t file_size_ = 0;
  int64_t parameter_count_ = 0;
  int64_t trainable_parameter_count_ = 0;
  int64_t buffer_count_ = 0;
};
