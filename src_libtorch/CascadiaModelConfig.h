#pragma once

#include <filesystem>
#include <string>

struct CascadiaModelConfig {
  std::filesystem::path model_path = "model/cascadia.pt";
  std::filesystem::path spectrum_path;
  std::string device = "cpu";
  std::string tokenizer = "massivekb";
  std::filesystem::path modifications_path;

  int d_model = 512;
  int n_layers = 9;
  int n_head = 8;
  int dim_feedforward = 1024;
  double dropout = 0.0;
  double rt_width = 2.0;
  int max_charge = 10;

  int batch_size = 32;
  int augmentation_width = 2;
  int candidate_max_charge = 4;
  int scan_width = 1;
  int top_n_peaks = 150;
  int max_sequence_length = 64;
  double score_threshold = 0.8;

/*
* @fn from_toml
* @brief Loads Cascadia model and sequence inference settings from a TOML file.
* @signature static CascadiaModelConfig from_toml(const std::filesystem::path& config_path);
* @param config_path: path to the TOML configuration file.
* @throws std::runtime_error when the config file cannot be opened or contains invalid assignments.
* @return Parsed CascadiaModelConfig.
*/
  static CascadiaModelConfig from_toml(
      const std::filesystem::path& config_path);

/*
* @fn summary
* @brief Builds a human-readable summary of the configured model and sequence inference settings.
* @signature std::string summary() const;
* @throws None.
* @return Configuration summary text.
*/
  std::string summary() const;
};
