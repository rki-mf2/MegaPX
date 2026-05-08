#pragma once

#include <torch/torch.h>

#include <filesystem>
#include <string>
#include <vector>

class CascadiaMzmlReader {
 public:
  struct Spectrum {
    int ms_level = 0;
    std::string id;
    double retention_time_seconds = 0.0;
    double selected_ion_mz = 0.0;
    int precursor_charge = 0;
    double isolation_target_mz = 0.0;
    double isolation_lower_offset = 0.0;
    double isolation_upper_offset = 0.0;
    std::vector<double> mz;
    std::vector<double> intensity;
  };

  struct AugmentedSpectrum {
    double precursor_mz = 0.0;
    int charge = 0;
    double retention_time_seconds = 0.0;
    std::vector<std::array<float, 4>> peaks;
  };

  struct TensorBatch {
    torch::Tensor spectra;
    torch::Tensor precursors;
    std::vector<double> retention_times;
    std::vector<double> precursor_mz;
    std::vector<int> charges;
  };

  struct Options {
    std::size_t top_n = 150;
    int scan_width = 1;
    int max_charge = 4;
  };

  std::vector<Spectrum> read(const std::filesystem::path& mzml_path) const;

  std::vector<AugmentedSpectrum> build_augmented_spectra(
      const std::vector<Spectrum>& spectra,
      const Options& options) const;

  TensorBatch to_tensors(const std::vector<AugmentedSpectrum>& spectra) const;
};
