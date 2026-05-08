#include "CascadiaPeptideTokenizer.h"

#include <algorithm>
#include <stdexcept>

CascadiaPeptideTokenizer CascadiaPeptideTokenizer::massivekb() {
  std::vector<std::string> tokens = {
      "$",
      "[+25.980265]-",
      "[Acetyl]-",
      "[Ammonia-loss]-",
      "[Carbamyl]-",
      "A",
      "C",
      "C[Carbamidomethyl]",
      "D",
      "E",
      "F",
      "G",
      "H",
      "K",
      "L",
      "M",
      "M[Oxidation]",
      "N",
      "N[Deamidated]",
      "P",
      "Q",
      "Q[Deamidated]",
      "R",
      "S",
      "T",
      "V",
      "W",
      "Y",
  };
  std::sort(tokens.begin(), tokens.end());
  return CascadiaPeptideTokenizer(std::move(tokens));
}

CascadiaPeptideTokenizer::CascadiaPeptideTokenizer(
    std::vector<std::string> tokens) {
  reverse_index_.reserve(tokens.size() + 1);
  reverse_index_.push_back({});
  for (const auto& token : tokens) {
    if (token == "$") {
      stop_token_id_ = static_cast<int64_t>(reverse_index_.size());
    }
    reverse_index_.push_back(token);
  }
  if (stop_token_id_ == 0) {
    throw std::invalid_argument("Tokenizer is missing the stop token.");
  }
}

int64_t CascadiaPeptideTokenizer::stop_token_id() const {
  return stop_token_id_;
}

int64_t CascadiaPeptideTokenizer::vocab_size_with_padding() const {
  return static_cast<int64_t>(reverse_index_.size());
}

std::vector<std::string> CascadiaPeptideTokenizer::detokenize(
    const torch::Tensor& token_ids) const {
  const auto tokens = token_ids.to(torch::kCPU).to(torch::kInt64).contiguous();
  if (tokens.dim() != 2) {
    throw std::invalid_argument("token_ids must have shape [batch, length].");
  }

  std::vector<std::string> sequences;
  const auto accessor = tokens.accessor<int64_t, 2>();
  for (int64_t row = 0; row < tokens.size(0); ++row) {
    std::string sequence;
    for (int64_t col = 0; col < tokens.size(1); ++col) {
      const auto id = accessor[row][col];
      if (id == 0) {
        continue;
      }
      if (id == stop_token_id_) {
        break;
      }
      if (id < 0 || id >= static_cast<int64_t>(reverse_index_.size())) {
        sequence += "<UNK:" + std::to_string(id) + ">";
        continue;
      }
      sequence += reverse_index_[static_cast<std::size_t>(id)];
    }
    sequences.push_back(std::move(sequence));
  }
  return sequences;
}
