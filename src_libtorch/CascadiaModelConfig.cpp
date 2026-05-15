#include "CascadiaModelConfig.h"

#include <algorithm>
#include <cctype>
#include <fstream>
#include <sstream>
#include <stdexcept>
#include <string>
#include <unordered_map>


namespace {

/*
* @fn trim
* @brief Removes leading and trailing whitespace from a string.
* @signature std::string trim(std::string value);
* @param value: string to trim.
* @throws None.
* @return Trimmed string.
*/
std::string trim(std::string value) {
  const auto is_space = [](unsigned char ch) { return std::isspace(ch) != 0; };
  value.erase(value.begin(),
              std::find_if_not(value.begin(), value.end(), is_space));
  value.erase(std::find_if_not(value.rbegin(), value.rend(), is_space).base(),
              value.end());
  return value;
}

/*
* @fn strip_comment
* @brief Removes TOML comment text while preserving hash characters inside quoted strings.
* @signature std::string strip_comment(const std::string& line);
* @param line: raw TOML line.
* @throws None.
* @return Line without trailing comment text.
*/
std::string strip_comment(const std::string& line) {
  bool in_string = false;
  char quote = '\0';
  for (std::size_t i = 0; i < line.size(); ++i) {
    const char ch = line[i];
    if ((ch == '"' || ch == '\'') && (i == 0 || line[i - 1] != '\\')) {
      if (!in_string) {
        in_string = true;
        quote = ch;
      } else if (quote == ch) {
        in_string = false;
      }
    } else if (ch == '#' && !in_string) {
      return line.substr(0, i);
    }
  }
  return line;
}

/*
* @fn unquote
* @brief Removes matching single or double quotes from a TOML string value.
* @signature std::string unquote(std::string value);
* @param value: raw TOML value.
* @throws None.
* @return Unquoted string value.
*/
std::string unquote(std::string value) {
  value = trim(std::move(value));
  if (value.size() >= 2 &&
      ((value.front() == '"' && value.back() == '"') ||
       (value.front() == '\'' && value.back() == '\''))) {
    return value.substr(1, value.size() - 2);
  }
  return value;
}

/*
* @fn lower
* @brief Converts a string to lowercase for case-insensitive TOML keys.
* @signature std::string lower(std::string value);
* @param value: input string.
* @throws None.
* @return Lowercase string.
*/
std::string lower(std::string value) {
  std::transform(value.begin(), value.end(), value.begin(), [](unsigned char ch) {
    return static_cast<char>(std::tolower(ch));
  });
  return value;
}

/*
* @fn parse_int
* @brief Parses an integer TOML value or returns a fallback.
* @signature int parse_int(const std::unordered_map<std::string, std::string>& values, const std::string& key, int fallback);
* @param values: parsed TOML key-value map.
* @param key: key to read.
* @param fallback: value returned when the key is absent.
* @throws std::invalid_argument or std::out_of_range when the value cannot be parsed by std::stoi.
* @return Parsed integer or fallback.
*/
int parse_int(const std::unordered_map<std::string, std::string>& values,
              const std::string& key, int fallback) {
  const auto it = values.find(key);
  if (it == values.end()) {
    return fallback;
  }
  return std::stoi(it->second);
}

/*
* @fn parse_double
* @brief Parses a floating-point TOML value or returns a fallback.
* @signature double parse_double(const std::unordered_map<std::string, std::string>& values, const std::string& key, double fallback);
* @param values: parsed TOML key-value map.
* @param key: key to read.
* @param fallback: value returned when the key is absent.
* @throws std::invalid_argument or std::out_of_range when the value cannot be parsed by std::stod.
* @return Parsed double or fallback.
*/
double parse_double(const std::unordered_map<std::string, std::string>& values,
                    const std::string& key, double fallback) {
  const auto it = values.find(key);
  if (it == values.end()) {
    return fallback;
  }
  return std::stod(it->second);
}

/*
* @fn parse_string
* @brief Reads and unquotes a string TOML value or returns a fallback.
* @signature std::string parse_string(const std::unordered_map<std::string, std::string>& values, const std::string& key, const std::string& fallback);
* @param values: parsed TOML key-value map.
* @param key: key to read.
* @param fallback: value returned when the key is absent.
* @throws None.
* @return Parsed string or fallback.
*/
std::string parse_string(const std::unordered_map<std::string, std::string>& values,
                         const std::string& key,
                         const std::string& fallback) {
  const auto it = values.find(key);
  if (it == values.end()) {
    return fallback;
  }
  return unquote(it->second);
}

/*
* @fn resolve_path
* @brief Resolves a TOML path relative to the config file location.
* @signature std::filesystem::path resolve_path(const std::filesystem::path& config_path, const std::string& raw_path);
* @param config_path: path to the TOML file.
* @param raw_path: raw path string from the TOML value.
* @throws None.
* @return Absolute paths unchanged, relative paths rooted at config_path.parent_path(), or an empty path.
*/
std::filesystem::path resolve_path(const std::filesystem::path& config_path,
                                   const std::string& raw_path) {
  if (raw_path.empty()) {
    return {};
  }

  std::filesystem::path path(raw_path);
  if (path.is_absolute()) {
    return path;
  }
  return config_path.parent_path() / path;
}

}  // namespace

/*
* @fn from_toml
* @brief Loads Cascadia model and sequence inference settings from a TOML file.
* @signature CascadiaModelConfig CascadiaModelConfig::from_toml(const std::filesystem::path& config_path);
* @param config_path: path to the TOML configuration file.
* @throws std::runtime_error when the config file cannot be opened or contains invalid assignments.
* @return Parsed CascadiaModelConfig.
*/
CascadiaModelConfig CascadiaModelConfig::from_toml(
    const std::filesystem::path& config_path) {
  std::ifstream input(config_path);
  if (!input) {
    throw std::runtime_error("Cannot open config file: " + config_path.string());
  }

  std::unordered_map<std::string, std::string> values;
  std::string section;
  std::string line;
  int line_number = 0;

  while (std::getline(input, line)) {
    ++line_number;
    line = trim(strip_comment(line));
    if (line.empty()) {
      continue;
    }

    if (line.front() == '[' && line.back() == ']') {
      section = lower(trim(line.substr(1, line.size() - 2)));
      continue;
    }

    const auto equals = line.find('=');
    if (equals == std::string::npos) {
      throw std::runtime_error("Invalid TOML assignment at " +
                               config_path.string() + ":" +
                               std::to_string(line_number));
    }

    const auto key = lower(trim(line.substr(0, equals)));
    const auto value = trim(line.substr(equals + 1));
    values[section.empty() ? key : section + "." + key] = value;
  }

  CascadiaModelConfig config;
  config.model_path = resolve_path(
      config_path,
      parse_string(values, "model.path",
                   parse_string(values, "model_path",
                                config.model_path.string())));
  config.spectrum_path = resolve_path(
      config_path,
      parse_string(values, "sequence.spectrum_path",
                   config.spectrum_path.string()));
  config.device = parse_string(values, "runtime.device",
                               parse_string(values, "model.device",
                                            config.device));
  config.tokenizer = parse_string(values, "model.tokenizer", config.tokenizer);
  config.modifications_path = resolve_path(
      config_path,
      parse_string(values, "model.modifications_path",
                   config.modifications_path.string()));

  config.d_model = parse_int(values, "model.d_model", config.d_model);
  config.n_layers = parse_int(values, "model.n_layers", config.n_layers);
  config.n_head = parse_int(values, "model.n_head", config.n_head);
  config.dim_feedforward =
      parse_int(values, "model.dim_feedforward", config.dim_feedforward);
  config.dropout = parse_double(values, "model.dropout", config.dropout);
  config.rt_width = parse_double(values, "model.rt_width", config.rt_width);
  config.max_charge = parse_int(values, "model.max_charge", config.max_charge);

  config.batch_size = parse_int(values, "sequence.batch_size",
                                config.batch_size);
  config.augmentation_width = parse_int(values, "sequence.augmentation_width",
                                        config.augmentation_width);
  config.candidate_max_charge = parse_int(values, "sequence.candidate_max_charge",
                                          config.candidate_max_charge);
  config.scan_width = parse_int(values, "sequence.scan_width",
                                config.scan_width);
  config.top_n_peaks = parse_int(values, "sequence.top_n_peaks",
                                 config.top_n_peaks);
  config.max_sequence_length = parse_int(values, "sequence.max_sequence_length",
                                         config.max_sequence_length);
  config.score_threshold = parse_double(values, "sequence.score_threshold",
                                        config.score_threshold);

  return config;
}

/*
* @fn summary
* @brief Builds a human-readable summary of the configured model and sequence inference settings.
* @signature std::string CascadiaModelConfig::summary() const;
* @throws None.
* @return Configuration summary text.
*/
std::string CascadiaModelConfig::summary() const {
  std::ostringstream out;
  out << "Cascadia config\n";
  out << "Model path: " << model_path.string() << '\n';
  if (!spectrum_path.empty()) {
    out << "Spectrum path: " << spectrum_path.string() << '\n';
  }
  out << "Device: " << device << '\n';
  out << "Tokenizer: " << tokenizer << '\n';
  if (!modifications_path.empty()) {
    out << "Modifications path: " << modifications_path.string() << '\n';
  }
  out << "Architecture: d_model=" << d_model << ", n_layers=" << n_layers
      << ", n_head=" << n_head
      << ", dim_feedforward=" << dim_feedforward
      << ", dropout=" << dropout << ", rt_width=" << rt_width
      << ", max_charge=" << max_charge << '\n';
  out << "Sequence inference: batch_size=" << batch_size
      << ", augmentation_width=" << augmentation_width
      << ", candidate_max_charge=" << candidate_max_charge
      << ", scan_width=" << scan_width << ", top_n_peaks=" << top_n_peaks
      << ", max_sequence_length=" << max_sequence_length
      << ", score_threshold=" << score_threshold << '\n';
  return out.str();
}
