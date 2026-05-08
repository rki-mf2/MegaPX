#pragma once

#include <torch/torch.h>

#include <string>
#include <vector>

class CascadiaPeptideTokenizer {
 public:
  static CascadiaPeptideTokenizer massivekb();

  int64_t stop_token_id() const;
  int64_t vocab_size_with_padding() const;
  std::vector<std::string> detokenize(const torch::Tensor& token_ids) const;

 private:
  explicit CascadiaPeptideTokenizer(std::vector<std::string> tokens);

  std::vector<std::string> reverse_index_;
  int64_t stop_token_id_ = 0;
};
