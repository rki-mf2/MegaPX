#include "CascadiaPeptideTokenizer.h"

#include <algorithm>
#include <stdexcept>

/*
* @fn massivekb
* @brief Creates the MassiveKB peptide tokenizer vocabulary used by Cascadia.
* @signature CascadiaPeptideTokenizer CascadiaPeptideTokenizer::massivekb();
* @throws std::invalid_argument when the tokenizer vocabulary is missing the stop token.
* @return CascadiaPeptideTokenizer initialized with the MassiveKB token set.
*/
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

/*
* @fn CascadiaPeptideTokenizer
* @brief Constructs a tokenizer from peptide tokens and records the stop token id.
* @signature CascadiaPeptideTokenizer::CascadiaPeptideTokenizer(std::vector<std::string> tokens);
* @param tokens: peptide token vocabulary.
* @throws std::invalid_argument when the stop token is missing.
* @return None.
*/
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

/*
* @fn stop_token_id
* @brief Returns the token id that marks the end of a peptide sequence.
* @signature int64_t CascadiaPeptideTokenizer::stop_token_id() const;
* @throws None.
* @return Stop token id.
*/
int64_t CascadiaPeptideTokenizer::stop_token_id() const {
  return stop_token_id_;
}

/*
* @fn vocab_size_with_padding
* @brief Returns the vocabulary size including the zero padding token.
* @signature int64_t CascadiaPeptideTokenizer::vocab_size_with_padding() const;
* @throws None.
* @return Vocabulary size including padding.
*/
int64_t CascadiaPeptideTokenizer::vocab_size_with_padding() const {
  return static_cast<int64_t>(reverse_index_.size());
}

/*
* @fn detokenize
* @brief Converts token id tensors into peptide sequence strings.
* @signature std::vector<std::string> CascadiaPeptideTokenizer::detokenize(const torch::Tensor& token_ids) const;
* @param token_ids: integer tensor with shape [batch, length].
* @throws std::invalid_argument when token_ids does not have shape [batch, length].
* @return Peptide sequence strings for each batch row.
*/
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
