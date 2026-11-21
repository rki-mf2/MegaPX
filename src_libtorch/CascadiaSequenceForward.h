#pragma once

#include "CascadiaModelConfig.h"
#include "CascadiaPeptideTokenizer.h"

#include <torch/script.h>

#include <filesystem>
#include <string>

class CascadiaSequenceForward {
 public:
  struct Output {
    torch::Tensor token_logits;
    torch::Tensor precursor_prediction;
    torch::Tensor fragment_logits;
  };

  struct GreedyResult {
    torch::Tensor token_ids;
    torch::Tensor amino_acid_confidence;
    torch::Tensor peptide_log_scores;
    std::vector<std::string> sequences;
  };

/*
* @fn CascadiaSequenceForward
* @brief Constructs a Cascadia TorchScript forward runner from model configuration.
* @signature explicit CascadiaSequenceForward(CascadiaModelConfig config);
* @param config: model path, runtime device, and inference settings.
* @throws None.
* @return None.
*/
  explicit CascadiaSequenceForward(CascadiaModelConfig config);

/*
* @fn load
* @brief Loads the configured TorchScript model onto the configured device.
* @signature void load();
* @throws c10::Error when the TorchScript model cannot be loaded.
* @return None.
*/
  void load();

/*
* @fn loaded
* @brief Reports whether a TorchScript model has been loaded.
* @signature bool loaded() const;
* @throws None.
* @return True when the model is loaded.
*/
  bool loaded() const;

/*
* @fn forward
* @brief Runs the Cascadia TorchScript forward method.
* @signature Output forward(const torch::Tensor& spectra, const torch::Tensor& precursors, const torch::Tensor& partial_sequence_tokens);
* @param spectra: spectrum tensor with shape [batch, peaks, 4].
* @param precursors: precursor tensor with shape [batch, 2].
* @param partial_sequence_tokens: partial peptide token tensor with shape [batch, sequence_length].
* @throws std::logic_error when load() has not been called.
* @throws std::invalid_argument when input tensor shapes are invalid.
* @throws std::runtime_error when the TorchScript output has an unsupported structure.
* @return Token logits, precursor prediction, and fragment logits.
*/
  Output forward(const torch::Tensor& spectra,
                 const torch::Tensor& precursors,
                 const torch::Tensor& partial_sequence_tokens);

/*
* @fn next_token_logits
* @brief Returns logits for the next token position from a forward pass.
* @signature torch::Tensor next_token_logits(const torch::Tensor& spectra, const torch::Tensor& precursors, const torch::Tensor& partial_sequence_tokens);
* @param spectra: spectrum tensor with shape [batch, peaks, 4].
* @param precursors: precursor tensor with shape [batch, 2].
* @param partial_sequence_tokens: partial peptide token tensor with shape [batch, sequence_length].
* @throws std::runtime_error when token logits have fewer than two dimensions.
* @return Tensor of next-token logits.
*/
  torch::Tensor next_token_logits(const torch::Tensor& spectra,
                                  const torch::Tensor& precursors,
                                  const torch::Tensor& partial_sequence_tokens);

/*
* @fn greedy_decode
* @brief Greedily decodes peptide token sequences from spectra and precursor tensors.
* @signature GreedyResult greedy_decode(const torch::Tensor& spectra, const torch::Tensor& precursors, const CascadiaPeptideTokenizer& tokenizer, int64_t max_length);
* @param spectra: spectrum tensor with shape [batch, peaks, 4].
* @param precursors: precursor tensor with shape [batch, 2].
* @param tokenizer: tokenizer used to detect stop tokens and detokenize ids.
* @param max_length: maximum decoded peptide length.
* @throws std::invalid_argument when max_length is not positive.
* @return Greedy decoding tensors and peptide sequence strings.
*/
  GreedyResult greedy_decode(const torch::Tensor& spectra,
                             const torch::Tensor& precursors,
                             const CascadiaPeptideTokenizer& tokenizer,
                             int64_t max_length);

/*
* @fn config
* @brief Returns the model configuration used by this forward runner.
* @signature const CascadiaModelConfig& config() const;
* @throws None.
* @return Reference to the current CascadiaModelConfig.
*/
  const CascadiaModelConfig& config() const;

 private:
/*
* @fn device
* @brief Resolves the configured runtime device.
* @signature c10::Device device() const;
* @throws std::runtime_error when CUDA is requested but unavailable.
* @return CPU or CUDA LibTorch device.
*/
  c10::Device device() const;

/*
* @fn validate_inputs
* @brief Validates tensor ranks, feature sizes, and batch-size compatibility.
* @signature void validate_inputs(const torch::Tensor& spectra, const torch::Tensor& precursors, const torch::Tensor& partial_sequence_tokens) const;
* @param spectra: spectrum tensor to validate.
* @param precursors: precursor tensor to validate.
* @param partial_sequence_tokens: partial sequence token tensor to validate.
* @throws std::invalid_argument when tensor shapes are incompatible with Cascadia inference.
* @return None.
*/
  void validate_inputs(const torch::Tensor& spectra,
                       const torch::Tensor& precursors,
                       const torch::Tensor& partial_sequence_tokens) const;

  CascadiaModelConfig config_;
  torch::jit::script::Module module_;
  bool loaded_ = false;
};
