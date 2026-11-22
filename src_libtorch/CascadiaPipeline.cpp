#include "CascadiaPipeline.h"
#include "CascadiaOutput.h"
#include "CascadiaRunStats.h"
#include "CascadiaMzmlReader.h"
#include "CascadiaSequenceForward.h"
#include <torch/torch.h>
#include <algorithm>
#include <cerrno>
#include <cmath>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <regex>
#include <stdexcept>
#include <vector>
#ifndef _WIN32
#include <fcntl.h>
#include <sys/wait.h>
#include <unistd.h>
#endif

namespace {
namespace fs = std::filesystem;
void help() {
  std::cout << "Usage: megapx_sequence --sequence INPUT.mzML --model BUNDLE.pt --out DIRECTORY [options]\n"
    "  --search --megapx /path/to/megapx -- [MegaPX multi_indexing arguments]\n"
    "  --batch-size 32 --threads 4 --device auto|cpu|cuda --score-threshold 0.8\n"
    "  --max-length 32 --max-charge 4 --scan-width 1 --top-peaks 150\n"
    "  --output-mode python|token  default: python (original scoring/filtering)\n"
    "  --jit-optimize on|off       runtime JIT profiling/optimization (default: on)\n"
    "  --limit N              process first N candidates (smoke tests)\n"
    "  --dump-tensors FILE    save first batch for Python parity checking\n"
    "  --search-mode MODE     default: multi_indexing\n"
    "Output: candidates.tsv, predictions.ssl, peptides.fasta, model.json, memory.txt; with search: search.tsv, search.log\n"
    "Export Lightning .ckpt once using export_cascadia_torchscript.py. Inference runs in C++.\n";
}

void run_search(const fs::path& binary, const std::string& mode,
                const std::vector<std::string>& options, const fs::path& out, CascadiaRunStats& stats) {
#ifndef _WIN32
  std::vector<std::string> command{binary.string(), mode};
  command.insert(command.end(), options.begin(), options.end());
  command.insert(command.end(), {"-f", (out / "peptides.fasta").string(), "-F", (out / "search.tsv").string()});
  std::vector<char*> args;
  for (auto& arg : command) args.push_back(arg.data());
  args.push_back(nullptr);
  int log = ::open((out / "search.log").c_str(), O_WRONLY | O_CREAT | O_EXCL, 0644);
  if (log < 0) throw std::runtime_error("Cannot create search.log");
  stats.search_started();
  const auto pid = fork();
  if (pid == 0) {
    if (chdir(out.c_str()) != 0 || dup2(log, STDOUT_FILENO) < 0 || dup2(log, STDERR_FILENO) < 0) _exit(126);
    close(log);
    execv(binary.c_str(), args.data());
    _exit(127);
  }
  close(log);
  if (pid < 0) throw std::runtime_error("Cannot launch MegaPX search");
  int status;
  rusage usage{};
  while (wait4(pid, &status, 0, &usage) < 0) if (errno != EINTR) throw std::runtime_error("Cannot wait for MegaPX");
  stats.search_finished(usage);
  if (!WIFEXITED(status) || WEXITSTATUS(status) != 0)
    throw std::runtime_error("MegaPX search failed; sequencing outputs retained. See " + (out / "search.log").string());
  if (!fs::exists(out / "search.tsv")) throw std::runtime_error("MegaPX did not write search.tsv; see search.log");
#else
  throw std::runtime_error("MegaPX search subprocess currently requires Linux or macOS");
#endif
}
}

int cascadia_cli(int argc, const char* const* argv) {
  CascadiaRunStats stats;
  try {
    CascadiaModelConfig config;
    config.max_sequence_length = 32;
    config.device = "auto";
    CascadiaMzmlReader::Options options;
    fs::path out, binary, dump;
    bool search = false, has_model = false;
    int threads = 4;
    std::size_t limit = 0;
    std::string mode = "multi_indexing";
    std::vector<std::string> search_args;
    if (argc == 1) { help(); return 0; }
    for (int i = 1; i < argc; ++i) {
      std::string arg = argv[i];
      auto value = [&]() -> std::string { if (++i >= argc) throw std::runtime_error("Missing value for " + arg); return argv[i]; };
      auto positive = [&]() { auto text = value(); std::size_t pos; int n = std::stoi(text, &pos); if (pos != text.size() || n < 1) throw std::runtime_error("Expected positive integer for " + arg); return n; };
      if (arg == "--help" || arg == "-h") { help(); return 0; }
      else if (arg == "--sequence") config.spectrum_path = fs::absolute(value());
      else if (arg == "--model") { config.model_path = fs::absolute(value()); has_model = true; }
      else if (arg == "--out" || arg == "-o") out = fs::absolute(value());
      else if (arg == "--jit-optimize") {
        auto setting = value();
        if (setting != "on" && setting != "off") throw std::runtime_error("--jit-optimize must be on or off");
        config.jit_optimize = setting == "on";
      }
      else if (arg == "--output-mode") config.output_mode = value();
      else if (arg == "--device") config.device = value();
      else if (arg == "--batch-size") config.batch_size = positive();
      else if (arg == "--threads") threads = positive();
      else if (arg == "--max-length") config.max_sequence_length = positive();
      else if (arg == "--max-charge") options.max_charge = positive();
      else if (arg == "--scan-width") options.scan_width = positive();
      else if (arg == "--top-peaks") options.top_n = positive();
      else if (arg == "--limit") limit = positive();
      else if (arg == "--dump-tensors") dump = fs::absolute(value());
      else if (arg == "--score-threshold") { auto text = value(); std::size_t pos; config.score_threshold = std::stod(text, &pos); if (pos != text.size()) throw std::runtime_error("Invalid score threshold"); }
      else if (arg == "--search") search = true;
      else if (arg == "--megapx") binary = fs::absolute(value());
      else if (arg == "--search-mode") mode = value();
      else if (arg == "--") { for (++i; i < argc; ++i) search_args.emplace_back(argv[i]); break; }
      else throw std::runtime_error("Unknown sequencing argument: " + arg);
    }
    if (config.spectrum_path.empty() || out.empty() || !has_model) throw std::runtime_error("--sequence, --model and --out are required");
    if (config.output_mode != "python" && config.output_mode != "token") throw std::runtime_error("--output-mode must be python or token");
    if (!std::isfinite(config.score_threshold) || config.score_threshold < 0 || config.score_threshold > 1)
      throw std::runtime_error("--score-threshold must be between 0 and 1");
    if (config.model_path.extension() == ".ckpt") throw std::runtime_error("Lightning .ckpt is not executable by LibTorch. Export it once with export_cascadia_torchscript.py --checkpoint MODEL.ckpt --output MODEL.pt");
    if (!fs::is_regular_file(config.spectrum_path) || !fs::is_regular_file(config.model_path)) throw std::runtime_error("Input spectrum or model file does not exist");
    if (search) {
      if (binary.empty()) throw std::runtime_error("Supply --megapx /path/to/megapx for optional multi_indexing search");
      if (!fs::is_regular_file(binary)) throw std::runtime_error("Supply --megapx /absolute/path/to/megapx");
      if (mode != "multi_indexing") throw std::runtime_error("Integrated search currently supports multi_indexing; other MegaPX modes can use peptides.fasta directly");
      bool reference = false, matrix = false;
      for (std::size_t i = 0; i < search_args.size(); ++i) {
        auto& arg = search_args[i];
        if (arg == "-f" || arg == "--queryFileName" || arg == "-F" || arg == "--filterFileName" || arg == "-o" || arg == "--outputDir" || arg.rfind("--queryFileName=", 0) == 0 || arg.rfind("--filterFileName=", 0) == 0)
          throw std::runtime_error("Query and output paths are managed by --sequence; omit -f, -F and -o from search arguments");
        if (arg == "-i" || arg == "--InputFasta" || arg == "-m" || arg == "--Matrix" || arg == "-b" || arg == "--blackListFile") {
          bool is_ref = arg == "-i" || arg == "--InputFasta";
          bool is_matrix = arg == "-m" || arg == "--Matrix";
          if (++i == search_args.size()) throw std::runtime_error("Missing search path");
          search_args[i] = fs::absolute(search_args[i]).string();
          if (!fs::is_regular_file(search_args[i])) throw std::runtime_error("Search input not found: " + search_args[i]);
          reference |= is_ref; matrix |= is_matrix;
        }
      }
      if (!reference || !matrix) throw std::runtime_error("Search needs reference -i and substitution matrix -m after --");
    } else if (!search_args.empty()) throw std::runtime_error("Search arguments require --search");
    for (const auto& name : {"candidates.tsv", "predictions.ssl", "peptides.fasta", "model.json", "search.tsv", "search.log", "memory.txt"})
      if (fs::exists(out / name)) throw std::runtime_error("Output already exists: " + (out / name).string() + "; choose a new --out directory");
    stats.start(out);
    torch::set_num_threads(threads);
    CascadiaSequenceForward forward(config);
    std::cout << "Loading Cascadia bundle " << config.model_path << '\n' << std::flush;
    forward.load();
    const auto& tokenizer = forward.tokenizer();
    CascadiaMzmlReader reader;
    auto raw = reader.read(config.spectrum_path);
    auto augmented = reader.build_augmented_spectra(raw, options);
    double cycle_time = 0.0, last_ms1 = 0.0, final_isolation_width = 0.0;
    for (const auto& scan : raw) if (scan.ms_level == 2) final_isolation_width = scan.isolation_lower_offset + scan.isolation_upper_offset;
    for (const auto& scan : raw) if (scan.ms_level == 1) { if (last_ms1 > 0) cycle_time = scan.retention_time_seconds - last_ms1; last_ms1 = scan.retention_time_seconds; }
    raw.clear(); raw.shrink_to_fit();
    if (augmented.empty()) throw std::runtime_error("No augmented candidates: input needs DIA MS2 isolation windows and neighboring MS1 spectra");
    if (limit && augmented.size() > limit) augmented.resize(limit);
    fs::create_directories(out);
    std::ofstream candidates(out / "candidates.tsv"), ssl(out / "predictions.ssl"), fasta(out / "peptides.fasta"), metadata(out / "model.json");
    if (!candidates || !ssl || !fasta || !metadata) throw std::runtime_error("Cannot create sequencing outputs");
    metadata << forward.metadata() << '\n';
    candidates << "candidate\tcharge\tprecursor_mz\tretention_time_seconds\tsequence\tunmodified_sequence\tmodifications\tconfidence\tpredicted_mz\tstatus\n";
    ssl << "file\tscan\tcharge\tsequence\tscore-type\tscore\tretention-time\tstart-time\tend-time\n";
    std::size_t accepted = 0;
    const bool python_mode = config.output_mode == "python";
    const bool write_minutes = !python_mode || std::any_of(augmented.begin(), augmented.end(), [](const auto& s) { return static_cast<float>(s.retention_time_seconds) > 500.0F; });
    const double retention_scale = write_minutes ? 60.0 : 1.0;
    const double rt_width = 2 * cycle_time / retention_scale;
    std::cout << "Sequencing " << augmented.size() << " candidates; encoder cached once per batch\n" << std::flush;
    for (std::size_t begin = 0; begin < augmented.size(); begin += config.batch_size) {
      auto end = std::min(augmented.size(), begin + config.batch_size);
      auto batch = reader.to_tensors({augmented.begin() + begin, augmented.begin() + end});
      if (begin == 0 && !dump.empty()) { torch::serialize::OutputArchive archive; archive.write("spectra", batch.spectra); archive.write("precursors", batch.precursors); archive.save_to(dump.string()); }
      auto result = forward.greedy_decode(batch.spectra, batch.precursors, tokenizer, config.max_sequence_length);
      for (std::size_t i = 0; i < result.sequences.size(); ++i) {
        const auto& item = augmented[begin + i];
        const auto& seq = result.sequences[i];
        auto details = cascadia_peptide_details(seq);
        const auto& plain = details.sequence;
        auto legacy_unmodified = std::regex_replace(seq, std::regex("\\[.*?\\]"), "");
        bool terminal_ok = python_mode ? (seq.find('-') == std::string::npos || legacy_unmodified.rfind('-') == 0) : details.valid;
        double confidence = python_mode ? static_cast<double>(std::exp(result.peptide_log_scores[i].item<float>())) : std::exp(result.peptide_log_scores[i].item<double>());
        double predicted_mz = tokenizer.peptide_mz(result.token_ids[i], item.charge);
        bool stopped = (result.token_ids[i] == tokenizer.stop_token_id()).any().item<bool>();
        double reference_mz = item.precursor_mz;
        if (python_mode) {
          // Match the inverse of the legacy model-input mass transform, using
          // the same float32 precursor values that Python's writer receives.
          const double mass = batch.precursors[i][0].item<float>();
          const float quotient = static_cast<float>(mass) / static_cast<float>(item.charge);
          reference_mz = ((static_cast<double>(quotient) + 1.007276) / item.charge + 1.007276);
        }
        const double window = python_mode ? final_isolation_width : item.isolation_width;
        std::string status = (python_mode ? legacy_unmodified.empty() : plain.empty()) ? "empty" : !terminal_ok ? "invalid_terminal_modification" : (!python_mode && !stopped) ? "length_limit" : !(confidence > config.score_threshold) ? "low_confidence" : !(std::abs(predicted_mz - reference_mz) < 2 * window) ? "precursor_mass" : "accepted";
        candidates << std::setprecision(10) << begin+i << '\t' << item.charge << '\t' << item.precursor_mz << '\t' << item.retention_time_seconds << '\t' << seq << '\t' << plain << '\t' << details.modifications << '\t' << confidence << '\t' << predicted_mz << '\t' << status << '\n';
        if (status == "accepted") {
          const double retention_time = python_mode ? static_cast<double>(static_cast<float>(item.retention_time_seconds) / static_cast<float>(retention_scale)) : item.retention_time_seconds / retention_scale;
          std::vector<std::string> fields{config.spectrum_path.string(), std::to_string(begin+i), std::to_string(item.charge), tokenizer.ssl_sequence(result.token_ids[i]), "UNKNOWN", cascadia_number(confidence), cascadia_number(retention_time), cascadia_number(retention_time-rt_width), cascadia_number(retention_time+rt_width)};
          for (std::size_t field = 0; field < fields.size(); ++field) { if (field) ssl << '\t'; ssl << fields[field]; }
          ssl << '\n';
          fasta << '>' << cascadia_fasta_header(fields, details.modifications) << '\n' << plain << '\n';
          ++accepted;
        }
      }
      std::cout << "Decoded " << end << '/' << augmented.size() << "; accepted " << accepted << '\n' << std::flush;
    }
    candidates.close(); ssl.close(); fasta.close(); metadata.close();
    if (!candidates || !ssl || !fasta || !metadata) throw std::runtime_error("Failed to finish writing sequencing output files");
    stats.sequence_finished();
    if (search) {
      if (!accepted) throw std::runtime_error("No accepted peptides to search; inspect candidates.tsv");
      run_search(binary, mode, search_args, out, stats);
    }
    std::cout << "Results written to " << out << '\n';
    stats.finish(true);
    return 0;
  } catch (const std::exception& error) {
    std::cerr << "Sequencing error: " << error.what() << '\n';
    stats.finish(false);
    return 1;
  }
}
