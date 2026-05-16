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

/*
* @fn read
* @brief Reads spectra and metadata from an mzML file.
* @signature std::vector<Spectrum> read(const std::filesystem::path& mzml_path) const;
* @param mzml_path: path to the mzML input file.
* @throws std::runtime_error when the mzML file cannot be opened or contains unsupported binary precision.
* @return Vector of parsed spectra with m/z and intensity arrays.
*/
  std::vector<Spectrum> read(const std::filesystem::path& mzml_path) const;

/*
* @fn build_augmented_spectra
* @brief Builds charge-candidate augmented spectra from parsed MS1/MS2 spectra.
* @signature std::vector<AugmentedSpectrum> build_augmented_spectra(const std::vector<Spectrum>& spectra, const Options& options) const;
* @param spectra: parsed mzML spectra.
* @param options: peak, scan-width, and charge-candidate settings.
* @throws None.
* @return Vector of augmented spectra ready for tensor conversion.
*/
  std::vector<AugmentedSpectrum> build_augmented_spectra(
      const std::vector<Spectrum>& spectra,
      const Options& options) const;

/*
* @fn to_tensors
* @brief Converts augmented spectra into LibTorch tensors and matching candidate metadata.
* @signature TensorBatch to_tensors(const std::vector<AugmentedSpectrum>& spectra) const;
* @param spectra: augmented spectra to batch.
* @throws None.
* @return TensorBatch containing spectra, precursor tensors, retention times, precursor m/z values, and charges.
*/
  TensorBatch to_tensors(const std::vector<AugmentedSpectrum>& spectra) const;
};
