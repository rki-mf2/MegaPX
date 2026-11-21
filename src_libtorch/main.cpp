#include "CascadiaModelConfig.h"
#include "CascadiaMzmlReader.h"
#include "CascadiaPeptideTokenizer.h"
#include "CascadiaSequenceForward.h"
#include "CheckpointModelLoader.h"
#include "InferenceParameterInspector.h"

#include <torch/torch.h>

#include <algorithm>
#include <cmath>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <regex>
#include <stdexcept>

namespace {

/*
* @fn print_tensor_shape
* @brief Prints a tensor name, shape, and dtype to standard output.
* @signature void print_tensor_shape(const char* name, const torch::Tensor& tensor);
* @param name: label printed before the tensor shape.
* @param tensor: tensor whose rank, dimensions, and dtype are printed.
* @throws None.
* @return None.
*/
void print_tensor_shape(const char* name, const torch::Tensor& tensor) {
  std::cout << name << ": [";
  for (int64_t i = 0; i < tensor.dim(); ++i) {
    if (i > 0) {
      std::cout << ", ";
    }
    std::cout << tensor.size(i);
  }
  std::cout << "] " << tensor.dtype() << '\n';
}

/*
* @fn ssl_output_path
* @brief Resolves the SSL output path and appends .ssl when no extension is provided.
* @signature std::filesystem::path ssl_output_path(const std::filesystem::path& output_path);
* @param output_path: configured output path or Cascadia-style output prefix.
* @throws None.
* @return Output path ending in .ssl.
*/
std::filesystem::path ssl_output_path(const std::filesystem::path& output_path) {
  if (output_path.empty()) {
    return {};
  }
  if (output_path.extension() == ".ssl") {
    return output_path;
  }
  auto ssl_path = output_path;
  ssl_path += ".ssl";
  return ssl_path;
}

/*
* @fn normalized_sequence_for_ssl
* @brief Applies Cascadia SSL sequence filtering rules for terminal modification tokens.
* @signature std::string normalized_sequence_for_ssl(const std::string& sequence);
* @param sequence: predicted MassiveKB-style peptide sequence.
* @throws None.
* @return Sequence when it is valid for SSL output, or an empty string when it should be skipped.
*/
std::string normalized_sequence_for_ssl(const std::string& sequence) {
  if (sequence.empty()) {
    return {};
  }

  const auto unmodified = std::regex_replace(sequence, std::regex("\\[.*?\\]"),
                                             "");
  if (sequence.find('-') != std::string::npos &&
      unmodified.rfind('-') != 0) {
    return {};
  }
  return sequence;
}

/*
* @fn run_cascadia_inference
* @brief Runs greedy Cascadia sequence inference over tensor batches and prints sample candidates.
* @signature void run_cascadia_inference(const CascadiaModelConfig& config, const CascadiaMzmlReader::TensorBatch& tensors);
* @param config: Cascadia runtime and decoding settings.
* @param tensors: batched spectra, precursors, and candidate metadata.
* @throws c10::Error when the TorchScript model cannot be loaded or executed.
* @throws std::exception when decoding or tensor access fails.
* @return None.
*/
void run_cascadia_inference(const CascadiaModelConfig& config,
                            const CascadiaMzmlReader::TensorBatch& tensors) {
  CascadiaSequenceForward forward(config);
  forward.load();

  const auto tokenizer = CascadiaPeptideTokenizer::massivekb();
  const int64_t total = tensors.spectra.size(0);
  const int64_t batch_size = std::max<int64_t>(1, config.batch_size);

  std::cout << "Running greedy Cascadia sequence inference\n";
  std::cout << "Candidates: " << total << ", batch_size: " << batch_size
            << ", max_sequence_length: " << config.max_sequence_length << '\n';

  const auto output_path = ssl_output_path(config.output_path);
  std::ofstream ssl_output;
  if (!output_path.empty()) {
    ssl_output.open(output_path);
    if (!ssl_output) {
      throw std::runtime_error("Cannot open SSL output file: " +
                               output_path.string());
    }
    ssl_output << "file\tscan\tcharge\tsequence\tscore-type\tscore"
               << "\tretention-time\tstart-time\tend-time\n";
    ssl_output.flush();
    std::cout << "Writing SSL results to: " << output_path.string() << '\n';
  }

  double max_retention_time = 0.0;
  for (const auto retention_time : tensors.retention_times) {
    max_retention_time = std::max(max_retention_time, retention_time);
  }
  const bool write_minutes = max_retention_time > 500.0;
  const double retention_scale = write_minutes ? 60.0 : 1.0;
  const double time_width =
      (config.augmentation_width * config.rt_width) / retention_scale;

  int64_t printed = 0;
  int64_t written = 0;
  for (int64_t begin = 0; begin < total; begin += batch_size) {
    const int64_t end = std::min(begin + batch_size, total);
    const auto rows = torch::indexing::Slice(begin, end);
    const auto spectra_batch = tensors.spectra.index({rows});
    const auto precursor_batch = tensors.precursors.index({rows});
    const auto result = forward.greedy_decode(
        spectra_batch, precursor_batch, tokenizer, config.max_sequence_length);

    const auto scores = result.peptide_log_scores.contiguous();
    const auto score_acc = scores.accessor<float, 1>();
    for (std::size_t i = 0; i < result.sequences.size(); ++i) {
      const auto global_index = begin + static_cast<int64_t>(i);
      if (printed < 10) {
        std::cout << "  candidate " << global_index << " charge="
                  << tensors.charges[static_cast<std::size_t>(global_index)]
                  << " precursor_mz="
                  << tensors.precursor_mz[static_cast<std::size_t>(global_index)]
                  << " rt="
                  << tensors.retention_times[static_cast<std::size_t>(global_index)]
                  << " sequence=" << result.sequences[i]
                  << " log_score=" << std::fixed << std::setprecision(4)
                  << score_acc[static_cast<int64_t>(i)] << '\n';
        ++printed;
      }

      if (ssl_output) {
        const auto sequence = normalized_sequence_for_ssl(result.sequences[i]);
        const auto confidence =
            std::exp(static_cast<double>(score_acc[static_cast<int64_t>(i)]));
        if (!sequence.empty() && confidence > config.score_threshold) {
          const auto metadata_index = static_cast<std::size_t>(global_index);
          const double retention_time =
              tensors.retention_times[metadata_index] / retention_scale;
          ssl_output << config.spectrum_path.string() << '\t'
                     << global_index << '\t'
                     << tensors.charges[metadata_index] << '\t'
                     << sequence << "\tUNKNOWN\t"
                     << confidence << '\t'
                     << retention_time << '\t'
                     << (retention_time - time_width) << '\t'
                     << (retention_time + time_width) << '\n';
          ++written;
        }
      }
    }
  }

  if (ssl_output) {
    std::cout << "SSL rows written: " << written << '\n';
  }

  if (total > printed) {
    std::cout << "  ... " << (total - printed)
              << " more candidates decoded\n";
  }
}

}  // namespace

/*
* @fn main
* @brief Command-line entry point for checkpoint inspection and optional mzML inference.
* @signature int main(int argc, char** argv);
* @param argc: command-line argument count.
* @param argv: command-line argument values.
* @throws None.
* @return Process exit code.
*/
int main(int argc, char** argv) {
  if (argc > 1) {
    const std::filesystem::path input_path(argv[1]);
    std::filesystem::path model_path = input_path;
    CascadiaModelConfig config;
    bool has_config = false;

    if (input_path.extension() == ".toml") {
      try {
        config = CascadiaModelConfig::from_toml(input_path);
        has_config = true;
        std::cout << config.summary() << std::endl;
        model_path = config.model_path;
      } catch (const std::exception& error) {
        std::cerr << "Config error: " << error.what() << std::endl;
        return 1;
      }
    }

    CheckpointModelLoader loader;
    const bool loaded = loader.load(model_path);
    std::cout << loader.summary() << std::endl;
    InferenceParameterInspector inspector(loader);
    std::cout << inspector.summary() << std::endl;

    if (has_config && !config.spectrum_path.empty()) {
      try {
        CascadiaMzmlReader reader;
        const auto spectra = reader.read(config.spectrum_path);
        CascadiaMzmlReader::Options options;
        options.top_n = static_cast<std::size_t>(config.top_n_peaks);
        options.scan_width = config.scan_width;
        options.max_charge = config.candidate_max_charge;
        const auto augmented = reader.build_augmented_spectra(spectra, options);
        const auto tensors = reader.to_tensors(augmented);

        std::cout << "mzML input\n";
        std::cout << "Raw spectra: " << spectra.size() << '\n';
        std::cout << "Augmented candidate spectra: " << augmented.size()
                  << '\n';
        print_tensor_shape("spectra tensor", tensors.spectra);
        print_tensor_shape("precursors tensor", tensors.precursors);
        if (loader.has_torchscript_module()) {
          run_cascadia_inference(config, tensors);
        } else {
          std::cout << "Inference skipped: model.path is not a TorchScript "
                       "export yet. Export the Lightning checkpoint to .pt and "
                       "put that path in the TOML model.path field.\n";
        }
      } catch (const std::exception& error) {
        std::cerr << "mzML error: " << error.what() << std::endl;
        return 1;
      }
    }

    return loaded ? 0 : 1;
  }

  std::cout << "\nPass a checkpoint/model path or TOML config:\n"
            << "  torch_loader model/cascadia.ckpt\n"
            << "  torch_loader src_libtorch/cascadia_sequence.example.toml\n";
  return 0;
}
