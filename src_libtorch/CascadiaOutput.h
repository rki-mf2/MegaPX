#pragma once
#include <cctype>
#include <iomanip>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

struct CascadiaPeptideDetails {
  std::string sequence;
  std::string modifications;
  bool valid = true;
};

// Positions count residues, starting at 1; terminal modifications are explicit.
inline CascadiaPeptideDetails cascadia_peptide_details(const std::string& modified) {
  CascadiaPeptideDetails result;
  std::vector<std::string> mods;
  bool terminal_separator = false;
  for (std::size_t i = 0; i < modified.size();) {
    char c = modified[i];
    if (c >= 'A' && c <= 'Z') {
      if (terminal_separator && !result.sequence.empty()) result.valid = false;
      terminal_separator = false;
      result.sequence += c;
      ++i;
    } else if (c == '[') {
      auto end = modified.find(']', i);
      if (end == std::string::npos) throw std::runtime_error("Unclosed modification token");
      std::string position = result.sequence.empty() ? "N-term" : terminal_separator ? "C-term" : std::to_string(result.sequence.size());
      mods.push_back(modified.substr(i+1, end-i-1) + "@" + position);
      i = end + 1;
    } else if (c == '-') {
      terminal_separator = !result.sequence.empty();
      ++i;
    } else {
      result.valid = false;
      ++i;
    }
  }
  for (const auto& mod : mods) {
    if (!result.modifications.empty()) result.modifications += ", ";
    result.modifications += mod;
  }
  if (result.modifications.empty()) result.modifications = "none";
  return result;
}

inline std::string cascadia_number(double value) {
  std::ostringstream out;
  out << std::setprecision(10) << value;
  return out.str();
}

inline std::string cascadia_header_value(const std::string& value) {
  std::ostringstream out;
  for (unsigned char c : value) {
    if (std::isspace(c) || c == '%') out << '%' << std::uppercase << std::hex << std::setw(2) << std::setfill('0') << static_cast<unsigned int>(c) << std::dec;
    else out << c;
  }
  return out.str();
}

inline std::string cascadia_fasta_header(const std::vector<std::string>& ssl_fields, const std::string& modifications) {
  static const char* names[] = {"file", "scan", "charge", "sequence", "score-type", "score", "retention-time", "start-time", "end-time"};
  if (ssl_fields.size() != 9) throw std::runtime_error("Expected nine SSL fields");
  std::string header;
  for (std::size_t i = 0; i < ssl_fields.size(); ++i) {
    if (i == 3) continue;
    if (!header.empty()) header += '_';
    header += std::string(names[i]) + '=' + cascadia_header_value(ssl_fields[i]);
  }
  return header + " [modification: " + modifications + "]";
}
