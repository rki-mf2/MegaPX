#include "CascadiaPeptideTokenizer.h"

#include <algorithm>
#include <cctype>
#include <stdexcept>
#include <sstream>
#include <set>
#include <cmath>
#include <iomanip>

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

CascadiaPeptideTokenizer CascadiaPeptideTokenizer::from_bundle(const std::string& tsv) {
  std::istringstream input(tsv);
  std::string line;
  std::getline(input, line);
  if (line != "#cascadia-bundle-v1")
    throw std::runtime_error("Missing versioned vocabulary. Re-export the checkpoint with export_cascadia_torchscript.py.");
  std::vector<std::string> tokens;
  std::vector<double> masses{0.0};
  std::set<std::string> seen;
  while (std::getline(input, line)) {
    const auto tab = line.find('\t');
    if (tab == std::string::npos) throw std::runtime_error("Invalid embedded vocabulary row");
    auto token = line.substr(0, tab);
    auto mass = std::stod(line.substr(tab + 1));
    if (token.empty() || !seen.insert(token).second || !std::isfinite(mass))
      throw std::runtime_error("Invalid embedded token/mass");
    tokens.push_back(token);
    masses.push_back(mass);
  }
  CascadiaPeptideTokenizer result(std::move(tokens));
  result.masses_ = std::move(masses);
  return result;
}

double CascadiaPeptideTokenizer::peptide_mz(const torch::Tensor& ids, int charge) const {
  if (charge <= 0 || masses_.size() != reverse_index_.size())
    throw std::runtime_error("Peptide mass requires bundled vocabulary masses and positive charge");
  auto values = ids.to(torch::kCPU).to(torch::kInt64).contiguous();
  const auto* data = values.data_ptr<int64_t>();
  double mass = 18.0105647;
  for (int64_t i = 0; i < values.numel(); ++i) {
    auto id = data[i];
    if (id == stop_token_id_) break;
    if (id < 0 || id >= static_cast<int64_t>(masses_.size())) throw std::runtime_error("Token outside vocabulary");
    mass += masses_[id];
  }
  return mass / charge + 1.00727646688;
}

std::string CascadiaPeptideTokenizer::ssl_sequence(const torch::Tensor& ids) const {
  auto values = ids.to(torch::kCPU).to(torch::kInt64).contiguous();
  std::ostringstream out;
  for (int64_t i = 0; i < values.numel(); ++i) {
    auto id = values.data_ptr<int64_t>()[i];
    if (id == stop_token_id_) break;
    if (id == 0) continue;
    const auto& token = reverse_index_.at(id);
    auto left = token.find('['), right = token.rfind(']');
    if (left == std::string::npos) { out << token; continue; }
    const auto label = token.substr(left + 1, right - left - 1);
    // Python leaves numeric mass labels untouched when Unimod has no name match.
    if (!label.empty() && (label[0] == '+' || label[0] == '-' || std::isdigit(static_cast<unsigned char>(label[0])))) {
      out << token;
      continue;
    }
    double shift = masses_.at(id);
    auto prefix = token.substr(0, left);
    for (char c : prefix) {
      if (c < 'A' || c > 'Z') continue;
      auto it = std::find(reverse_index_.begin(), reverse_index_.end(), std::string(1, c));
      if (it == reverse_index_.end()) throw std::runtime_error("Modified residue lacks base mass");
      shift -= masses_[std::distance(reverse_index_.begin(), it)];
    }
    std::ostringstream formatted;
    formatted << std::showpos << std::fixed << std::setprecision(2) << shift;
    auto shift_text = formatted.str();
    while (shift_text.size() > 2 && shift_text.back() == '0' && shift_text[shift_text.size()-2] != '.') shift_text.pop_back();
    out << prefix << '[' << shift_text << ']' << token.substr(right+1);
  }
  return out.str();
}
