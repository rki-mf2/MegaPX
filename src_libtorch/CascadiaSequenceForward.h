#pragma once

#include "CascadiaModelConfig.h"

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

  explicit CascadiaSequenceForward(CascadiaModelConfig config);

  void load();
  bool loaded() const;

  Output forward(const torch::Tensor& spectra,
                 const torch::Tensor& precursors,
                 const torch::Tensor& partial_sequence_tokens);

  torch::Tensor next_token_logits(const torch::Tensor& spectra,
                                  const torch::Tensor& precursors,
                                  const torch::Tensor& partial_sequence_tokens);

  const CascadiaModelConfig& config() const;

 private:
  c10::Device device() const;
  void validate_inputs(const torch::Tensor& spectra,
                       const torch::Tensor& precursors,
                       const torch::Tensor& partial_sequence_tokens) const;

  CascadiaModelConfig config_;
  torch::jit::script::Module module_;
  bool loaded_ = false;
};
