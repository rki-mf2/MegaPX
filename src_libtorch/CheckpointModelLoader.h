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

/*
* @fn load
* @brief Loads a TorchScript module or inspects a PyTorch checkpoint archive.
* @signature bool load(const std::filesystem::path& checkpoint_path, const c10::Device& device = c10::Device(c10::kCPU));
* @param checkpoint_path: path to the checkpoint or TorchScript model file.
* @param device: LibTorch device used when loading TorchScript.
* @throws None.
* @return True when a TorchScript module was loaded or a PyTorch zip checkpoint was recognized.
*/
  bool load(const std::filesystem::path& checkpoint_path,
            const c10::Device& device = c10::Device(c10::kCPU));

/*
* @fn loaded
* @brief Reports whether the loader has a module or recognized checkpoint archive.
* @signature bool loaded() const;
* @throws None.
* @return True when load() recognized usable model or archive information.
*/
  bool loaded() const;

/*
* @fn has_torchscript_module
* @brief Reports whether the loaded file is directly runnable as TorchScript.
* @signature bool has_torchscript_module() const;
* @throws None.
* @return True when a TorchScript module is loaded.
*/
  bool has_torchscript_module() const;

/*
* @fn last_error
* @brief Returns the last module-loading or archive-inspection error message.
* @signature const std::string& last_error() const;
* @throws None.
* @return Last error string.
*/
  const std::string& last_error() const;

/*
* @fn summary
* @brief Formats checkpoint, module, tensor, and archive metadata.
* @signature std::string summary(std::size_t max_items = 12) const;
* @param max_items: maximum number of tensors to print per tensor section.
* @throws None.
* @return Checkpoint summary text.
*/
  std::string summary(std::size_t max_items = 12) const;

/*
* @fn module
* @brief Returns the loaded TorchScript module.
* @signature const torch::jit::script::Module& module() const;
* @throws std::logic_error when no TorchScript module has been loaded.
* @return Reference to the loaded TorchScript module.
*/
  const torch::jit::script::Module& module() const;

/*
* @fn is_python_checkpoint
* @brief Reports whether the archive looks like a Python/PyTorch-Lightning checkpoint.
* @signature bool is_python_checkpoint() const;
* @throws None.
* @return True when the file is a Python checkpoint rather than TorchScript.
*/
  bool is_python_checkpoint() const;

/*
* @fn checkpoint_keys
* @brief Returns discovered top-level checkpoint keys from archive/data.pkl.
* @signature const std::vector<std::string>& checkpoint_keys() const;
* @throws None.
* @return Vector of checkpoint key names.
*/
  const std::vector<std::string>& checkpoint_keys() const;

/*
* @fn sample_state_dict_names
* @brief Returns representative state_dict tensor names discovered in the checkpoint.
* @signature const std::vector<std::string>& sample_state_dict_names() const;
* @throws None.
* @return Vector of sampled tensor names.
*/
  const std::vector<std::string>& sample_state_dict_names() const;

/*
* @fn pytorch_lightning_version
* @brief Returns the PyTorch-Lightning version recorded in the checkpoint, when present.
* @signature const std::optional<std::string>& pytorch_lightning_version() const;
* @throws None.
* @return Optional PyTorch-Lightning version string.
*/
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

/*
* @fn reset
* @brief Clears all loaded module state, archive metadata, and summary counters.
* @signature void reset();
* @throws None.
* @return None.
*/
  void reset();

/*
* @fn collect_module_summary
* @brief Collects TorchScript method, parameter, and buffer metadata.
* @signature void collect_module_summary();
* @throws None.
* @return None.
*/
  void collect_module_summary();

/*
* @fn inspect_zip_archive
* @brief Inspects a PyTorch zip archive without requiring Python model code.
* @signature bool inspect_zip_archive(const std::filesystem::path& checkpoint_path);
* @param checkpoint_path: path to the zip checkpoint to inspect.
* @throws None.
* @return True when the file is recognized as a zip archive.
*/
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
