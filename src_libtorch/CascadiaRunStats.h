#pragma once
#include <algorithm>
#include <chrono>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <sys/resource.h>

// ru_maxrss is bytes on macOS and KiB on Linux. Peaks are per process,
// so parent and search-child peaks are reported separately, never added.
inline double cascadia_peak_mib(const rusage& usage) {
#ifdef __APPLE__
  return static_cast<double>(usage.ru_maxrss) / (1024.0 * 1024.0);
#else
  return static_cast<double>(usage.ru_maxrss) / 1024.0;
#endif
}

inline double cascadia_cpu_seconds(const rusage& usage) {
  return usage.ru_utime.tv_sec + usage.ru_stime.tv_sec +
         (usage.ru_utime.tv_usec + usage.ru_stime.tv_usec) / 1e6;
}

class CascadiaRunStats {
  using Clock = std::chrono::steady_clock;
 public:
  void start(const std::filesystem::path& directory) {
    directory_ = directory;
    getrusage(RUSAGE_SELF, &initial_);
    started_ = Clock::now();
    active_ = true;
  }
  void sequence_finished() {
    sequence_seconds_ = seconds_since(started_);
    sequence_complete_ = true;
  }
  void search_started() {
    search_started_ = Clock::now();
    search_ran_ = true;
  }
  void search_finished(const rusage& usage) {
    search_seconds_ = seconds_since(search_started_);
    child_ = usage;
  }
  void finish(bool success) {
    if (!active_) return;
    active_ = false;
    const double elapsed = seconds_since(started_);
    rusage current{};
    getrusage(RUSAGE_SELF, &current);
    const double self_cpu = cascadia_cpu_seconds(current) - cascadia_cpu_seconds(initial_);
    std::ostringstream report;
    report << std::fixed << std::setprecision(3)
           << "\n**************** MegaPX Sequencing Usage Report ****************\n"
           << "Status                 : " << (success ? "complete" : "failed") << '\n'
           << "Real time (total)      : " << elapsed << " sec\n"
           << "Real time (sequencing) : " << (sequence_complete_ ? sequence_seconds_ : elapsed) << " sec\n"
           << "CPU time (total)       : " << self_cpu + cascadia_cpu_seconds(child_) << " sec\n"
           << "CPU time (sequencer)   : " << self_cpu << " sec\n"
           << "Peak RSS (sequencer)   : " << cascadia_peak_mib(current) << " MiB\n";
    if (search_ran_) {
      report << "Real time (search)     : " << search_seconds_ << " sec\n"
             << "CPU time (search)      : " << cascadia_cpu_seconds(child_) << " sec\n"
             << "Peak RSS (search)      : " << cascadia_peak_mib(child_) << " MiB\n";
    }
    report << "***************************************************************\n";
    std::cout << report.str() << std::flush;
    // Validation failures do not create result directories just for metrics.
    if (std::filesystem::is_directory(directory_)) {
      std::ofstream file(directory_ / "memory.txt");
      file << report.str();
      file.close();
      if (!file) std::cerr << "Unable to save usage report to " << directory_ / "memory.txt" << '\n';
    }
  }
 private:
  static double seconds_since(Clock::time_point start) {
    return std::chrono::duration<double>(Clock::now() - start).count();
  }
  bool active_ = false, sequence_complete_ = false, search_ran_ = false;
  Clock::time_point started_, search_started_;
  std::filesystem::path directory_;
  rusage initial_{}, child_{};
  double sequence_seconds_ = 0, search_seconds_ = 0;
};
