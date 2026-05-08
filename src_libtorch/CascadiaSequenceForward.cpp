#include "CascadiaSequenceForward.h"

#include <torch/torch.h>

#include <algorithm>
#include <cctype>
#include <stdexcept>
#include <utility>
#include <vector>

namespace {

torch::Tensor require_tensor(const c10::IValue& value,
                             const std::string& name) {
  if (!value.isTensor()) {
    throw std::runtime_error("TorchScript forward output '" + name +
                             "' is not a Tensor.");
  }
  return value.toTensor();
}

}  // namespace

CascadiaSequenceForward::CascadiaSequenceForward(CascadiaModelConfig config)
    : config_(std::move(config)) {}

void CascadiaSequenceForward::load() {
  module_ = torch::jit::load(config_.model_path.string(), device());
  module_.eval();
  loaded_ = true;
}

bool CascadiaSequenceForward::loaded() const {
  return loaded_;
}

CascadiaSequenceForward::Output CascadiaSequenceForward::forward(
    const torch::Tensor& spectra,
    const torch::Tensor& precursors,
    const torch::Tensor& partial_sequence_tokens) {
  if (!loaded_) {
    throw std::logic_error("CascadiaSequenceForward::load() must be called first.");
  }

  validate_inputs(spectra, precursors, partial_sequence_tokens);

  torch::NoGradGuard no_grad;
  std::vector<c10::IValue> inputs{
      spectra.to(device()),
      precursors.to(device()),
      partial_sequence_tokens.to(device(), torch::kInt64),
  };

  const auto result = module_.forward(std::move(inputs));
  if (result.isTensor()) {
    return {result.toTensor(), torch::Tensor(), torch::Tensor()};
  }

  if (!result.isTuple()) {
    throw std::runtime_error(
        "Cascadia TorchScript forward must return logits Tensor or "
        "(token_logits, precursor_prediction, fragment_logits).");
  }

  const auto tuple = result.toTuple();
  const auto& elements = tuple->elements();
  if (elements.empty()) {
    throw std::runtime_error("Cascadia TorchScript forward returned an empty tuple.");
  }

  Output output;
  output.token_logits = require_tensor(elements[0], "token_logits");
  if (elements.size() > 1 && !elements[1].isNone()) {
    output.precursor_prediction =
        require_tensor(elements[1], "precursor_prediction");
  }
  if (elements.size() > 2 && !elements[2].isNone()) {
    output.fragment_logits = require_tensor(elements[2], "fragment_logits");
  }
  return output;
}

torch::Tensor CascadiaSequenceForward::next_token_logits(
    const torch::Tensor& spectra,
    const torch::Tensor& precursors,
    const torch::Tensor& partial_sequence_tokens) {
  auto output = forward(spectra, precursors, partial_sequence_tokens);
  if (output.token_logits.dim() < 2) {
    throw std::runtime_error("token_logits must have at least 2 dimensions.");
  }

  const auto last_index = output.token_logits.size(1) - 1;
  return output.token_logits.index({torch::indexing::Slice(), last_index,
                                    torch::indexing::Slice()});
}

const CascadiaModelConfig& CascadiaSequenceForward::config() const {
  return config_;
}

c10::Device CascadiaSequenceForward::device() const {
  auto requested = config_.device;
  std::transform(requested.begin(), requested.end(), requested.begin(),
                 [](unsigned char ch) {
                   return static_cast<char>(std::tolower(ch));
                 });

  if (requested == "cuda" || requested == "gpu") {
    if (!torch::cuda::is_available()) {
      throw std::runtime_error("Config requested CUDA, but CUDA is not available.");
    }
    return c10::Device(c10::kCUDA);
  }
  return c10::Device(c10::kCPU);
}

void CascadiaSequenceForward::validate_inputs(
    const torch::Tensor& spectra,
    const torch::Tensor& precursors,
    const torch::Tensor& partial_sequence_tokens) const {
  if (spectra.dim() != 3 || spectra.size(2) != 4) {
    throw std::invalid_argument(
        "spectra must have shape [batch, peaks, 4] for "
        "(mz, intensity, retention_time, ms_level).");
  }
  if (precursors.dim() != 2 || precursors.size(1) != 2) {
    throw std::invalid_argument(
        "precursors must have shape [batch, 2] for "
        "(neutral_precursor_mass, charge).");
  }
  if (partial_sequence_tokens.dim() != 2) {
    throw std::invalid_argument(
        "partial_sequence_tokens must have shape [batch, sequence_length].");
  }
  if (spectra.size(0) != precursors.size(0) ||
      spectra.size(0) != partial_sequence_tokens.size(0)) {
    throw std::invalid_argument(
        "spectra, precursors, and partial_sequence_tokens must have the same "
        "batch size.");
  }
}
