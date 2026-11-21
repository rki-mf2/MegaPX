#pragma once

#include <torch/torch.h>

#include <string>
#include <vector>

class CascadiaPeptideTokenizer {
 public:
/*
* @fn massivekb
* @brief Creates the MassiveKB peptide tokenizer vocabulary used by Cascadia.
* @signature static CascadiaPeptideTokenizer massivekb();
* @throws std::invalid_argument when the tokenizer vocabulary is missing the stop token.
* @return CascadiaPeptideTokenizer initialized with the MassiveKB token set.
*/
  static CascadiaPeptideTokenizer massivekb();

/*
* @fn stop_token_id
* @brief Returns the token id that marks the end of a peptide sequence.
* @signature int64_t stop_token_id() const;
* @throws None.
* @return Stop token id.
*/
  int64_t stop_token_id() const;

/*
* @fn vocab_size_with_padding
* @brief Returns the vocabulary size including the zero padding token.
* @signature int64_t vocab_size_with_padding() const;
* @throws None.
* @return Vocabulary size including padding.
*/
  int64_t vocab_size_with_padding() const;

/*
* @fn detokenize
* @brief Converts token id tensors into peptide sequence strings.
* @signature std::vector<std::string> detokenize(const torch::Tensor& token_ids) const;
* @param token_ids: integer tensor with shape [batch, length].
* @throws std::invalid_argument when token_ids does not have shape [batch, length].
* @return Peptide sequence strings for each batch row.
*/
  std::vector<std::string> detokenize(const torch::Tensor& token_ids) const;

 private:
/*
* @fn CascadiaPeptideTokenizer
* @brief Constructs a tokenizer from sorted peptide tokens and records the stop token id.
* @signature explicit CascadiaPeptideTokenizer(std::vector<std::string> tokens);
* @param tokens: peptide token vocabulary.
* @throws std::invalid_argument when the stop token is missing.
* @return None.
*/
  explicit CascadiaPeptideTokenizer(std::vector<std::string> tokens);

  std::vector<std::string> reverse_index_;
  int64_t stop_token_id_ = 0;
};
