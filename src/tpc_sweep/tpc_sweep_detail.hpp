#pragma once

#include <condition_variable>
#include <fstream>
#include <mutex>
#include <ostream>
#include <string>
#include <vector>

#include "newcode/tpc_sweep_runner.hpp"

namespace tpc_sweep {
namespace detail {

struct SweepScenario {
  std::string name;
  int bitgen_seed = 0;
  int channel_seed = 0;
  float ebn0_db = newcode::DEFAULT_EBN0_DB;
};

struct ScenarioOutput {
  std::size_t idx{};
  std::string name;
  int bitgen_seed = 0;
  int channel_seed = 0;
  float ebn0_db = newcode::DEFAULT_EBN0_DB;
  newcode::tpc::TpcPipelineResult result;
};

class DualOut {
 public:
  DualOut(std::ostream& console, const std::string& filepath, bool mirror_console = true);

  template <typename T>
  DualOut& operator<<(const T& value) {
    if (console_) {
      *console_ << value;
    }
    if (file_) {
      file_ << value;
    }
    return *this;
  }

  DualOut& operator<<(std::ostream& (*pf)(std::ostream&));

 private:
  std::ostream* console_;
  std::ofstream file_;
};

class Semaphore {
 public:
  explicit Semaphore(std::size_t count);

  void acquire();
  void release();

 private:
  std::mutex mutex_;
  std::condition_variable cv_;
  std::size_t count_;
};

void ensure_csv_header(const std::string& csv_path);
std::string join_vec(const std::vector<float>& values, char sep, int precision);
void write_csv_row(std::ostream& csv,
                   const std::string& timestamp,
                   const std::string& run_id,
                   const SweepScenario& scenario,
                   const newcode::tpc::TpcPipelineResult& result,
                   const SweepParameterConfig& config);

std::vector<int> generate_random_seeds(int count);
std::vector<float> build_ebn0_values(const SweepParameterConfig& config);

std::vector<SweepScenario> build_scenarios(const SweepParameterConfig& config,
                                           const std::vector<float>& ebn0_candidates,
                                           const std::vector<int>& bitgen_seeds,
                                           const std::vector<int>& channel_seeds);

std::vector<ScenarioOutput> run_scenarios_parallel(
    const std::vector<SweepScenario>& scenarios,
    const SweepParameterConfig& config,
    unsigned max_workers_hint = 0,
    const std::string& stage_tag = std::string{},
    DualOut* log = nullptr);

unsigned resolve_worker_count(const SweepParameterConfig& config);

}  // namespace detail
}  // namespace tpc_sweep
