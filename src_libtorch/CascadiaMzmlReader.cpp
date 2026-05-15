#include "CascadiaMzmlReader.h"

#include <zlib.h>

#include <algorithm>
#include <array>
#include <cctype>
#include <cmath>
#include <cstring>
#include <fstream>
#include <iterator>
#include <limits>
#include <numeric>
#include <stdexcept>
#include <string>
#include <unordered_map>

namespace {

constexpr double kProtonMass = 1.007276;

/*
* @fn read_file
* @brief Reads an entire mzML file into memory.
* @signature std::string read_file(const std::filesystem::path& path);
* @param path: mzML file path.
* @throws std::runtime_error when the file cannot be opened.
* @return File contents as a string.
*/
std::string read_file(const std::filesystem::path& path) {
  std::ifstream input(path, std::ios::binary);
  if (!input) {
    throw std::runtime_error("Cannot open mzML file: " + path.string());
  }
  return {std::istreambuf_iterator<char>(input),
          std::istreambuf_iterator<char>()};
}

/*
* @fn attr_value
* @brief Extracts an XML attribute value from a tag string.
* @signature std::string attr_value(const std::string& tag, const std::string& attr);
* @param tag: XML tag text.
* @param attr: attribute name to read.
* @throws None.
* @return Attribute value, or an empty string when missing.
*/
std::string attr_value(const std::string& tag, const std::string& attr) {
  const auto needle = attr + "=\"";
  const auto begin = tag.find(needle);
  if (begin == std::string::npos) {
    return {};
  }
  const auto value_begin = begin + needle.size();
  const auto value_end = tag.find('"', value_begin);
  if (value_end == std::string::npos) {
    return {};
  }
  return tag.substr(value_begin, value_end - value_begin);
}

/*
* @fn parse_double_or
* @brief Parses a double value or returns a fallback for empty input.
* @signature double parse_double_or(const std::string& value, double fallback = 0.0);
* @param value: string value to parse.
* @param fallback: value returned when input is empty.
* @throws std::invalid_argument or std::out_of_range when std::stod cannot parse a non-empty value.
* @return Parsed double or fallback.
*/
double parse_double_or(const std::string& value, double fallback = 0.0) {
  if (value.empty()) {
    return fallback;
  }
  return std::stod(value);
}

/*
* @fn parse_int_or
* @brief Parses an integer value or returns a fallback for empty input.
* @signature int parse_int_or(const std::string& value, int fallback = 0);
* @param value: string value to parse.
* @param fallback: value returned when input is empty.
* @throws std::invalid_argument or std::out_of_range when std::stoi cannot parse a non-empty value.
* @return Parsed integer or fallback.
*/
int parse_int_or(const std::string& value, int fallback = 0) {
  if (value.empty()) {
    return fallback;
  }
  return std::stoi(value);
}

/*
* @fn find_cv_value
* @brief Finds a controlled-vocabulary parameter value by name in an XML block.
* @signature std::string find_cv_value(const std::string& xml, const std::string& name, std::size_t offset = 0);
* @param xml: XML text to search.
* @param name: cvParam name attribute to match.
* @param offset: search offset in the XML text.
* @throws None.
* @return Matched value attribute, or an empty string when missing.
*/
std::string find_cv_value(const std::string& xml, const std::string& name,
                          std::size_t offset = 0) {
  const auto pos = xml.find("name=\"" + name + "\"", offset);
  if (pos == std::string::npos) {
    return {};
  }
  const auto tag_begin = xml.rfind('<', pos);
  const auto tag_end = xml.find('>', pos);
  if (tag_begin == std::string::npos || tag_end == std::string::npos) {
    return {};
  }
  return attr_value(xml.substr(tag_begin, tag_end - tag_begin + 1), "value");
}

/*
* @fn has_name
* @brief Tests whether an XML block contains a name attribute with the requested value.
* @signature bool has_name(const std::string& xml, const std::string& name);
* @param xml: XML text to search.
* @param name: name attribute value to match.
* @throws None.
* @return True when the name is present.
*/
bool has_name(const std::string& xml, const std::string& name) {
  return xml.find("name=\"" + name + "\"") != std::string::npos;
}

/*
* @fn between
* @brief Extracts text between an opening and closing XML marker.
* @signature std::string between(const std::string& xml, const std::string& open, const std::string& close, std::size_t offset = 0);
* @param xml: XML text to search.
* @param open: opening marker.
* @param close: closing marker.
* @param offset: search offset in the XML text.
* @throws None.
* @return Text between markers, or an empty string when missing.
*/
std::string between(const std::string& xml, const std::string& open,
                    const std::string& close, std::size_t offset = 0) {
  const auto begin = xml.find(open, offset);
  if (begin == std::string::npos) {
    return {};
  }
  const auto content_begin = xml.find('>', begin);
  if (content_begin == std::string::npos) {
    return {};
  }
  const auto end = xml.find(close, content_begin);
  if (end == std::string::npos) {
    return {};
  }
  return xml.substr(content_begin + 1, end - content_begin - 1);
}

/*
* @fn base64_decode
* @brief Decodes base64 mzML binary payload text.
* @signature std::vector<unsigned char> base64_decode(const std::string& input);
* @param input: base64-encoded text.
* @throws None.
* @return Decoded bytes.
*/
std::vector<unsigned char> base64_decode(const std::string& input) {
  static constexpr std::array<int, 256> table = [] {
    std::array<int, 256> out{};
    out.fill(-1);
    const std::string chars =
        "ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz0123456789+/";
    for (int i = 0; i < static_cast<int>(chars.size()); ++i) {
      out[static_cast<unsigned char>(chars[i])] = i;
    }
    return out;
  }();

  std::vector<unsigned char> bytes;
  int value = 0;
  int bits = -8;
  for (const unsigned char ch : input) {
    if (ch == '=') {
      break;
    }
    if (std::isspace(ch)) {
      continue;
    }
    const int decoded = table[ch];
    if (decoded < 0) {
      continue;
    }
    value = (value << 6) + decoded;
    bits += 6;
    if (bits >= 0) {
      bytes.push_back(static_cast<unsigned char>((value >> bits) & 0xff));
      bits -= 8;
    }
  }
  return bytes;
}

/*
* @fn inflate_zlib
* @brief Decompresses zlib-compressed mzML binary payload bytes.
* @signature std::vector<unsigned char> inflate_zlib(const std::vector<unsigned char>& data, std::size_t expected_size);
* @param data: compressed bytes.
* @param expected_size: expected decompressed byte count used for allocation.
* @throws std::runtime_error when zlib initialization or decompression fails.
* @return Decompressed bytes.
*/
std::vector<unsigned char> inflate_zlib(const std::vector<unsigned char>& data,
                                        std::size_t expected_size) {
  z_stream stream{};
  stream.next_in = const_cast<Bytef*>(data.data());
  stream.avail_in = static_cast<uInt>(data.size());

  if (inflateInit(&stream) != Z_OK) {
    throw std::runtime_error("zlib inflateInit failed.");
  }

  std::vector<unsigned char> output;
  output.reserve(expected_size);
  std::array<unsigned char, 16384> buffer{};
  int status = Z_OK;
  while (status == Z_OK) {
    stream.next_out = buffer.data();
    stream.avail_out = static_cast<uInt>(buffer.size());
    status = inflate(&stream, Z_NO_FLUSH);
    const auto produced = buffer.size() - stream.avail_out;
    output.insert(output.end(), buffer.begin(), buffer.begin() + produced);
  }

  inflateEnd(&stream);
  if (status != Z_STREAM_END) {
    throw std::runtime_error("zlib inflate failed while decoding mzML binary.");
  }
  return output;
}

/*
* @fn read_little_endian
* @brief Reads a little-endian scalar value from mzML binary bytes.
* @signature template <typename T> T read_little_endian(const unsigned char* data);
* @param data: pointer to little-endian bytes.
* @throws None.
* @return Decoded scalar value.
*/
template <typename T>
T read_little_endian(const unsigned char* data) {
  T value{};
  std::memcpy(&value, data, sizeof(T));
  return value;
}

/*
* @fn decode_binary_array
* @brief Decodes an mzML binaryDataArray into double values.
* @signature std::vector<double> decode_binary_array(const std::string& binary_xml, std::size_t default_array_length);
* @param binary_xml: binaryDataArray XML block.
* @param default_array_length: mzML default array length used to reserve decompression output.
* @throws std::runtime_error when the binary precision is unsupported or zlib decompression fails.
* @return Decoded m/z or intensity values, or an empty vector for other array types.
*/
std::vector<double> decode_binary_array(const std::string& binary_xml,
                                        std::size_t default_array_length) {
  const bool is_mz = has_name(binary_xml, "m/z array");
  const bool is_intensity = has_name(binary_xml, "intensity array");
  if (!is_mz && !is_intensity) {
    return {};
  }

  const bool is_float64 = has_name(binary_xml, "64-bit float");
  const bool is_float32 = has_name(binary_xml, "32-bit float");
  const bool is_zlib = has_name(binary_xml, "zlib compression");
  const auto binary_text = between(binary_xml, "<binary>", "</binary>");
  if (binary_text.empty()) {
    return {};
  }

  auto bytes = base64_decode(binary_text);
  const std::size_t bytes_per_value = is_float32 ? sizeof(float) : sizeof(double);
  if (is_zlib) {
    bytes = inflate_zlib(bytes, default_array_length * bytes_per_value);
  }

  if (!is_float64 && !is_float32) {
    throw std::runtime_error("Unsupported mzML binary precision.");
  }

  const std::size_t count = bytes.size() / bytes_per_value;
  std::vector<double> values;
  values.reserve(count);
  for (std::size_t i = 0; i < count; ++i) {
    const auto* ptr = bytes.data() + i * bytes_per_value;
    values.push_back(is_float32 ? read_little_endian<float>(ptr)
                                : read_little_endian<double>(ptr));
  }
  return values;
}

/*
* @fn top_peaks_sorted_by_mz
* @brief Selects the most intense peaks, normalizes intensities, and returns them sorted by m/z.
* @signature std::vector<std::pair<double, double>> top_peaks_sorted_by_mz(const CascadiaMzmlReader::Spectrum& spectrum, std::size_t top_n, bool double_sqrt);
* @param spectrum: spectrum containing m/z and intensity arrays.
* @param top_n: maximum number of peaks to keep.
* @param double_sqrt: whether to apply a second square-root intensity transform.
* @throws None.
* @return Vector of normalized m/z-intensity peak pairs sorted by m/z.
*/
std::vector<std::pair<double, double>> top_peaks_sorted_by_mz(
    const CascadiaMzmlReader::Spectrum& spectrum,
    std::size_t top_n,
    bool double_sqrt) {
  std::vector<std::size_t> indices(spectrum.intensity.size());
  std::iota(indices.begin(), indices.end(), 0);
  const auto keep = std::min(top_n, indices.size());
  std::partial_sort(indices.begin(), indices.begin() + keep, indices.end(),
                    [&](std::size_t lhs, std::size_t rhs) {
                      return spectrum.intensity[lhs] > spectrum.intensity[rhs];
                    });
  indices.resize(keep);
  std::sort(indices.begin(), indices.end(), [&](std::size_t lhs,
                                                std::size_t rhs) {
    return spectrum.mz[lhs] < spectrum.mz[rhs];
  });

  std::vector<std::pair<double, double>> peaks;
  peaks.reserve(indices.size());
  double max_intensity = 0.0;
  for (const auto index : indices) {
    double intensity = std::sqrt(std::max(0.0, spectrum.intensity[index]));
    if (double_sqrt) {
      intensity = std::sqrt(intensity);
    }
    max_intensity = std::max(max_intensity, intensity);
    peaks.emplace_back(spectrum.mz[index], intensity);
  }
  if (max_intensity > 0.0) {
    for (auto& peak : peaks) {
      peak.second /= max_intensity;
    }
  }
  return peaks;
}

struct Center {
  double mz = 0.0;
  double rt = 0.0;
  double lower_offset = 0.0;
  double upper_offset = 0.0;
};

}  // namespace

/*
* @fn read
* @brief Reads spectra and metadata from an mzML file.
* @signature std::vector<CascadiaMzmlReader::Spectrum> CascadiaMzmlReader::read(const std::filesystem::path& mzml_path) const;
* @param mzml_path: path to the mzML input file.
* @throws std::runtime_error when the mzML file cannot be opened or contains unsupported binary precision.
* @return Vector of parsed spectra with m/z and intensity arrays.
*/
std::vector<CascadiaMzmlReader::Spectrum> CascadiaMzmlReader::read(
    const std::filesystem::path& mzml_path) const {
  const auto xml = read_file(mzml_path);
  std::vector<Spectrum> spectra;

  std::size_t cursor = 0;
  while (true) {
    const auto begin = xml.find("<spectrum ", cursor);
    if (begin == std::string::npos) {
      break;
    }
    const auto tag_end = xml.find('>', begin);
    const auto end = xml.find("</spectrum>", tag_end);
    if (tag_end == std::string::npos || end == std::string::npos) {
      break;
    }

    const auto open_tag = xml.substr(begin, tag_end - begin + 1);
    const auto block = xml.substr(begin, end + 11 - begin);
    cursor = end + 11;

    Spectrum spectrum;
    spectrum.id = attr_value(open_tag, "id");
    spectrum.ms_level = parse_int_or(find_cv_value(block, "ms level"), 0);
    spectrum.retention_time_seconds =
        60.0 * parse_double_or(find_cv_value(block, "scan start time"), 0.0);

    if (spectrum.ms_level == 2) {
      spectrum.isolation_target_mz =
          parse_double_or(find_cv_value(block, "isolation window target m/z"));
      spectrum.isolation_lower_offset =
          parse_double_or(find_cv_value(block, "isolation window lower offset"));
      spectrum.isolation_upper_offset =
          parse_double_or(find_cv_value(block, "isolation window upper offset"));
      spectrum.selected_ion_mz =
          parse_double_or(find_cv_value(block, "selected ion m/z"),
                          spectrum.isolation_target_mz);
      spectrum.precursor_charge =
          parse_int_or(find_cv_value(block, "charge state"), 0);
    }

    const auto default_array_length =
        static_cast<std::size_t>(parse_int_or(attr_value(open_tag,
                                                         "defaultArrayLength")));
    std::size_t array_cursor = 0;
    while (true) {
      const auto array_begin = block.find("<binaryDataArray", array_cursor);
      if (array_begin == std::string::npos) {
        break;
      }
      const auto array_end = block.find("</binaryDataArray>", array_begin);
      if (array_end == std::string::npos) {
        break;
      }
      const auto array_block =
          block.substr(array_begin, array_end + 18 - array_begin);
      array_cursor = array_end + 18;
      auto values = decode_binary_array(array_block, default_array_length);
      if (values.empty()) {
        continue;
      }
      if (has_name(array_block, "m/z array")) {
        spectrum.mz = std::move(values);
      } else if (has_name(array_block, "intensity array")) {
        spectrum.intensity = std::move(values);
      }
    }

    if (!spectrum.mz.empty() && spectrum.mz.size() == spectrum.intensity.size()) {
      spectra.push_back(std::move(spectrum));
    }
  }

  return spectra;
}

/*
* @fn build_augmented_spectra
* @brief Builds charge-candidate augmented spectra from parsed MS1/MS2 spectra.
* @signature std::vector<CascadiaMzmlReader::AugmentedSpectrum> CascadiaMzmlReader::build_augmented_spectra(const std::vector<Spectrum>& spectra, const Options& options) const;
* @param spectra: parsed mzML spectra.
* @param options: peak, scan-width, and charge-candidate settings.
* @throws None.
* @return Vector of augmented spectra ready for tensor conversion.
*/
std::vector<CascadiaMzmlReader::AugmentedSpectrum>
CascadiaMzmlReader::build_augmented_spectra(
    const std::vector<Spectrum>& spectra,
    const Options& options) const {
  std::vector<Center> centers;
  double previous_ms1_rt = 0.0;
  double cycle_time = 0.0;
  for (const auto& spectrum : spectra) {
    if (spectrum.ms_level == 1) {
      if (previous_ms1_rt > 0.0) {
        cycle_time = spectrum.retention_time_seconds - previous_ms1_rt;
      }
      previous_ms1_rt = spectrum.retention_time_seconds;
    } else if (spectrum.ms_level == 2 && spectrum.isolation_target_mz > 0.0) {
      centers.push_back({spectrum.isolation_target_mz,
                         spectrum.retention_time_seconds,
                         spectrum.isolation_lower_offset,
                         spectrum.isolation_upper_offset});
    }
  }

  if (cycle_time <= 0.0) {
    cycle_time = 2.0;
  }
  const double time_width = (options.scan_width + 1) * cycle_time;

  std::vector<AugmentedSpectrum> augmented;
  for (const auto& center : centers) {
    std::vector<std::pair<double, std::vector<std::pair<double, double>>>> ms2;
    std::vector<std::pair<double, std::vector<std::pair<double, double>>>> ms1;

    for (const auto& spectrum : spectra) {
      const double rt_delta = spectrum.retention_time_seconds - center.rt;
      if (std::abs(rt_delta) >= time_width) {
        continue;
      }

      if (spectrum.ms_level == 2) {
        const bool in_window =
            center.mz > spectrum.isolation_target_mz -
                            spectrum.isolation_lower_offset &&
            center.mz < spectrum.isolation_target_mz +
                            spectrum.isolation_upper_offset;
        if (in_window) {
          ms2.emplace_back(rt_delta,
                           top_peaks_sorted_by_mz(spectrum, options.top_n,
                                                  false));
        }
      } else if (spectrum.ms_level == 1) {
        ms1.emplace_back(rt_delta,
                         top_peaks_sorted_by_mz(spectrum, options.top_n, true));
      }
    }

    if (ms1.empty() || ms2.empty()) {
      continue;
    }

    const auto nearest = [](const auto& lhs, const auto& rhs) {
      return std::abs(lhs.first) < std::abs(rhs.first);
    };
    std::sort(ms1.begin(), ms1.end(), nearest);
    std::sort(ms2.begin(), ms2.end(), nearest);
    ms1.resize(std::min<std::size_t>(options.scan_width, ms1.size()));
    ms2.resize(std::min<std::size_t>(options.scan_width, ms2.size()));

    const double ms1_window = std::max(center.lower_offset, center.upper_offset);
    for (int charge = 1; charge <= options.max_charge; ++charge) {
      AugmentedSpectrum item;
      item.precursor_mz = center.mz;
      item.charge = charge;
      item.retention_time_seconds = center.rt;

      for (const auto& [rt_delta, peaks] : ms2) {
        for (const auto& [mz, intensity] : peaks) {
          item.peaks.push_back({static_cast<float>(mz),
                                static_cast<float>(intensity),
                                static_cast<float>(rt_delta), 2.0F});
        }
      }
      for (const auto& [rt_delta, peaks] : ms1) {
        for (const auto& [mz, intensity] : peaks) {
          if (std::abs(mz - center.mz) < ms1_window + 1.0) {
            item.peaks.push_back({static_cast<float>(mz),
                                  static_cast<float>(intensity),
                                  static_cast<float>(rt_delta), 1.0F});
          }
        }
      }

      if (!item.peaks.empty()) {
        augmented.push_back(std::move(item));
      }
    }
  }

  return augmented;
}

/*
* @fn to_tensors
* @brief Converts augmented spectra into LibTorch tensors and matching candidate metadata.
* @signature CascadiaMzmlReader::TensorBatch CascadiaMzmlReader::to_tensors(const std::vector<AugmentedSpectrum>& spectra) const;
* @param spectra: augmented spectra to batch.
* @throws None.
* @return TensorBatch containing spectra, precursor tensors, retention times, precursor m/z values, and charges.
*/
CascadiaMzmlReader::TensorBatch CascadiaMzmlReader::to_tensors(
    const std::vector<AugmentedSpectrum>& spectra) const {
  TensorBatch batch;
  const auto batch_size = static_cast<int64_t>(spectra.size());
  std::size_t max_peaks = 0;
  for (const auto& spectrum : spectra) {
    max_peaks = std::max(max_peaks, spectrum.peaks.size());
  }

  batch.spectra =
      torch::zeros({batch_size, static_cast<int64_t>(max_peaks), 4},
                   torch::TensorOptions().dtype(torch::kFloat32));
  batch.precursors = torch::zeros({batch_size, 2},
                                  torch::TensorOptions().dtype(torch::kFloat32));

  auto spectra_acc = batch.spectra.accessor<float, 3>();
  auto precursor_acc = batch.precursors.accessor<float, 2>();
  for (int64_t i = 0; i < batch_size; ++i) {
    const auto& spectrum = spectra[static_cast<std::size_t>(i)];
    for (std::size_t j = 0; j < spectrum.peaks.size(); ++j) {
      for (int k = 0; k < 4; ++k) {
        spectra_acc[i][static_cast<int64_t>(j)][k] = spectrum.peaks[j][k];
      }
    }
    precursor_acc[i][0] =
        static_cast<float>((spectrum.precursor_mz - kProtonMass) *
                           spectrum.charge);
    precursor_acc[i][1] = static_cast<float>(spectrum.charge);
    batch.retention_times.push_back(spectrum.retention_time_seconds);
    batch.precursor_mz.push_back(spectrum.precursor_mz);
    batch.charges.push_back(spectrum.charge);
  }

  return batch;
}
