#include "CascadiaSequenceForward.h"

#include <torch/torch.h>
#include <torch/csrc/jit/runtime/graph_executor.h>
#include <torch/csrc/jit/python/update_graph_executor_opt.h>

#include <algorithm>
#include <cctype>
#include <stdexcept>
#include <utility>
#include <vector>

namespace {

/*
* @fn require_tensor
* @brief Converts a TorchScript IValue output to a tensor and validates its type.
* @signature torch::Tensor require_tensor(const c10::IValue& value, const std::string& name);
* @param value: TorchScript output value.
* @param name: logical output name used in error messages.
* @throws std::runtime_error when value is not a tensor.
* @return Tensor stored in the IValue.
*/
torch::Tensor require_tensor(const c10::IValue& value,
                             const std::string& name) {
  if (!value.isTensor()) {
    throw std::runtime_error("TorchScript forward output '" + name +
                             "' is not a Tensor.");
  }
  return value.toTensor();
}

}  // namespace

/*
* @fn CascadiaSequenceForward
* @brief Constructs a Cascadia TorchScript forward runner from model configuration.
* @signature CascadiaSequenceForward::CascadiaSequenceForward(CascadiaModelConfig config);
* @param config: model path, runtime device, and inference settings.
* @throws None.
* @return None.
*/
CascadiaSequenceForward::CascadiaSequenceForward(CascadiaModelConfig config)
    : config_(std::move(config)) {}

/*
* @fn load
* @brief Loads the configured TorchScript model onto the configured device.
* @signature void CascadiaSequenceForward::load();
* @throws c10::Error when the TorchScript model cannot be loaded.
* @return None.
*/
void CascadiaSequenceForward::load() {
  // This CLI owns its Torch runtime. Disable plan profiling/optimization for
  // memory-constrained inference when requested; model weights are unchanged.
  torch::jit::setGraphExecutorOptimize(config_.jit_optimize);
  torch::jit::getProfilingMode().store(config_.jit_optimize);
  torch::jit::getExecutorMode().store(config_.jit_optimize);
  torch::jit::ExtraFilesMap extra{{"cascadia.vocab.tsv", ""}, {"cascadia.metadata.json", ""}};
  module_ = torch::jit::load(config_.model_path.string(), device(), extra);
  tokenizer_ = CascadiaPeptideTokenizer::from_bundle(extra.at("cascadia.vocab.tsv"));
  metadata_ = extra.at("cascadia.metadata.json");
  if (!module_.find_method("encode_spectra") || !module_.find_method("decode_tokens"))
    throw std::runtime_error("Model bundle lacks cached encoder methods; re-export the checkpoint");
  module_.eval();
  loaded_ = true;
}

/*
* @fn loaded
* @brief Reports whether a TorchScript model has been loaded.
* @signature bool CascadiaSequenceForward::loaded() const;
* @throws None.
* @return True when the model is loaded.
*/
bool CascadiaSequenceForward::loaded() const {
  return loaded_;
}

/*
* @fn forward
* @brief Runs the Cascadia TorchScript forward method.
* @signature CascadiaSequenceForward::Output CascadiaSequenceForward::forward(const torch::Tensor& spectra, const torch::Tensor& precursors, const torch::Tensor& partial_sequence_tokens);
* @param spectra: spectrum tensor with shape [batch, peaks, 4].
* @param precursors: precursor tensor with shape [batch, 2].
* @param partial_sequence_tokens: partial peptide token tensor with shape [batch, sequence_length].
* @throws std::logic_error when load() has not been called.
* @throws std::invalid_argument when input tensor shapes are invalid.
* @throws std::runtime_error when the TorchScript output has an unsupported structure.
* @return Token logits, precursor prediction, and fragment logits.
*/
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

/*
* @fn next_token_logits
* @brief Returns logits for the next token position from a forward pass.
* @signature torch::Tensor CascadiaSequenceForward::next_token_logits(const torch::Tensor& spectra, const torch::Tensor& precursors, const torch::Tensor& partial_sequence_tokens);
* @param spectra: spectrum tensor with shape [batch, peaks, 4].
* @param precursors: precursor tensor with shape [batch, 2].
* @param partial_sequence_tokens: partial peptide token tensor with shape [batch, sequence_length].
* @throws std::runtime_error when token logits have fewer than two dimensions.
* @return Tensor of next-token logits.
*/
torch::Tensor CascadiaSequenceForward::next_token_logits(
    const torch::Tensor& spectra,
    const torch::Tensor& precursors,
    const torch::Tensor& partial_sequence_tokens) {
  auto output = forward(spectra, precursors, partial_sequence_tokens);
  if (output.token_logits.dim() < 2) {
    throw std::runtime_error("token_logits must have at least 2 dimensions.");
  }
  if (output.token_logits.dim() == 2) {
    return output.token_logits;
  }

  const auto last_index = output.token_logits.size(1) - 1;
  return output.token_logits.index({torch::indexing::Slice(), last_index,
                                    torch::indexing::Slice()});
}

/*
* @fn greedy_decode
* @brief Greedily decodes peptide token sequences from spectra and precursor tensors.
* @signature CascadiaSequenceForward::GreedyResult CascadiaSequenceForward::greedy_decode(const torch::Tensor& spectra, const torch::Tensor& precursors, const CascadiaPeptideTokenizer& tokenizer, int64_t max_length);
* @param spectra: spectrum tensor with shape [batch, peaks, 4].
* @param precursors: precursor tensor with shape [batch, 2].
* @param tokenizer: tokenizer used to detect stop tokens and detokenize ids.
* @param max_length: maximum decoded peptide length.
* @throws std::invalid_argument when max_length is not positive.
* @return Greedy decoding tensors and peptide sequence strings.
*/
CascadiaSequenceForward::GreedyResult CascadiaSequenceForward::greedy_decode(
    const torch::Tensor& spectra,
    const torch::Tensor& precursors,
    const CascadiaPeptideTokenizer& tokenizer,
    const int64_t max_length) {
  if (max_length <= 0) {
    throw std::invalid_argument("max_length must be positive.");
  }

  torch::InferenceMode inference;
  auto spectrum_device = spectra.to(device());
  auto precursor_device = precursors.to(device());
  auto encoded = module_.get_method("encode_spectra")({spectrum_device}).toTuple();
  const auto memory = encoded->elements()[0];
  const auto memory_mask = encoded->elements()[1];
  const auto batch_size = spectra.size(0);
  auto tokens = torch::empty({batch_size, 0},
                             torch::TensorOptions().dtype(torch::kInt64)
                                 .device(device()));
  auto aa_confidence = torch::empty({batch_size, 0},
                                    torch::TensorOptions().dtype(torch::kFloat32)
                                        .device(device()));
  auto finished = torch::zeros({batch_size},
                               torch::TensorOptions().dtype(torch::kBool)
                                   .device(device()));

  for (int64_t step = 0; step < max_length; ++step) {
    auto all_logits = module_.get_method("decode_tokens")({memory, memory_mask, precursor_device, tokens}).toTensor();
    if (all_logits.size(2) != tokenizer.vocab_size_with_padding())
      throw std::runtime_error("Model logits disagree with embedded vocabulary");
    auto logits = all_logits.index({torch::indexing::Slice(), -1, torch::indexing::Slice()});
    auto scores = torch::softmax(logits, 1);
    auto next_tokens = std::get<1>(scores.max(1)).to(torch::kInt64);
    auto next_confidence = scores.gather(1, next_tokens.unsqueeze(1)).squeeze(1);

    if (config_.output_mode != "python") {
      next_tokens = torch::where(finished, torch::full_like(next_tokens, tokenizer.stop_token_id()), next_tokens);
      next_confidence = torch::where(finished, torch::ones_like(next_confidence), next_confidence);
    }

    tokens = torch::cat({tokens, next_tokens.unsqueeze(1)}, 1);
    aa_confidence = torch::cat({aa_confidence, next_confidence.unsqueeze(1)}, 1);
    finished = finished.logical_or(next_tokens == tokenizer.stop_token_id());

    if (config_.output_mode != "python" && finished.all().item<bool>()) {
      break;
    }
  }

  GreedyResult result;
  result.token_ids = tokens.to(torch::kCPU);
  result.amino_acid_confidence = aa_confidence.to(torch::kCPU);
  result.sequences = tokenizer.detokenize(result.token_ids);
  if (config_.output_mode == "python") {
    // Cascadia 0.0.7 uses string length, including modification labels, and
    // continues decoding after the first stop. Preserve this default exactly.
    std::vector<torch::Tensor> log_scores;
    for (int64_t row = 0; row < batch_size; ++row) {
      int64_t length = 0;
      for (unsigned char ch : result.sequences[row]) if ((ch & 0xc0) != 0x80) ++length;
      auto count = std::min<int64_t>(length + 1, aa_confidence.size(1));
      log_scores.push_back(torch::log(aa_confidence[row].slice(0, 0, count)).sum());
    }
    result.peptide_log_scores = torch::stack(log_scores).to(torch::kCPU);
  } else {
    result.peptide_log_scores = torch::log(aa_confidence.clamp_min(1.0e-12)).sum(1).to(torch::kCPU);
  }
  return result;
}

/*
* @fn config
* @brief Returns the model configuration used by this forward runner.
* @signature const CascadiaModelConfig& CascadiaSequenceForward::config() const;
* @throws None.
* @return Reference to the current CascadiaModelConfig.
*/
const CascadiaModelConfig& CascadiaSequenceForward::config() const {
  return config_;
}

/*
* @fn device
* @brief Resolves the configured runtime device.
* @signature c10::Device CascadiaSequenceForward::device() const;
* @throws std::runtime_error when CUDA is requested but unavailable.
* @return CPU or CUDA LibTorch device.
*/
c10::Device CascadiaSequenceForward::device() const {
  auto requested = config_.device;
  std::transform(requested.begin(), requested.end(), requested.begin(),
                 [](unsigned char ch) {
                   return static_cast<char>(std::tolower(ch));
                 });

  if (requested != "auto" && requested != "cpu" && requested != "cuda" && requested != "gpu")
    throw std::runtime_error("Supported devices: auto, cpu, cuda");
  if (requested == "cuda" || requested == "gpu" || (requested == "auto" && torch::cuda::is_available())) {
    if (!torch::cuda::is_available()) {
      throw std::runtime_error("Config requested CUDA, but CUDA is not available.");
    }
    return c10::Device(c10::kCUDA);
  }
  return c10::Device(c10::kCPU);
}

/*
* @fn validate_inputs
* @brief Validates tensor ranks, feature sizes, and batch-size compatibility.
* @signature void CascadiaSequenceForward::validate_inputs(const torch::Tensor& spectra, const torch::Tensor& precursors, const torch::Tensor& partial_sequence_tokens) const;
* @param spectra: spectrum tensor to validate.
* @param precursors: precursor tensor to validate.
* @param partial_sequence_tokens: partial sequence token tensor to validate.
* @throws std::invalid_argument when tensor shapes are incompatible with Cascadia inference.
* @return None.
*/
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
