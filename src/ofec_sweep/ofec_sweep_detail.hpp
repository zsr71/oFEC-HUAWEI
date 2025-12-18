#pragma once

#include <condition_variable>
#include <fstream>
#include <mutex>
#include <ostream>
#include <string>
#include <vector>

#include "newcode/ofec_sweep_runner.hpp"

namespace ofec_sweep {
namespace detail {

struct SweepScenario {
  std::string name;
  std::vector<float> alpha_list;
  std::vector<float> beta_list;
  float alpha_start = 0.0f;
  float alpha_step = 0.0f;
  float beta_start = 0.0f;
  float beta_step = 0.0f;
  float alpha_low = 0.0f;
  float alpha_high = 0.0f;
  float gamma_alpha = 1.0f;
  float beta_low = 0.0f;
  float beta_high = 0.0f;
  float gamma_beta = 1.0f;
  int chase_L = 0;
  int chase_n_test = 0;
  int bitgen_seed = 0;
  int channel_seed = 0;
  float ebn0_db = newcode::DEFAULT_EBN0_DB;
};

struct ScenarioOutput {
  std::size_t idx{};
  std::string name;
  std::vector<float> alpha_list;
  std::vector<float> beta_list;
  float alpha_start = 0.0f;
  float alpha_step = 0.0f;
  float beta_start = 0.0f;
  float beta_step = 0.0f;
  int chase_L = 0;
  int chase_n_test = 0;
  int bitgen_seed = 0;
  int channel_seed = 0;
  float ebn0_db = newcode::DEFAULT_EBN0_DB;
  newcode::PipelineResult result;
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
void ensure_csv_header_v2(const std::string& csv_path);

std::vector<float> generate_sequence(float start, float step, std::size_t length);
float infer_step(const std::vector<float>& values);
std::vector<int> generate_random_seeds(int count);

newcode::PipelineConfig make_pipeline_config(const SweepParameterConfig& config);

std::vector<SweepScenario> build_scenarios(const SweepParameterConfig& config,
                                           const std::vector<float>& ebn0_candidates,
                                           const std::vector<int>& bitgen_seeds,
                                           const std::vector<int>& channel_seeds);

std::vector<float> build_ebn0_values(const SweepParameterConfig& config);

std::string join_vec(const std::vector<float>& values, char sep, int precision);
std::string join_vec(const std::vector<double>& values, char sep, int precision);
double mean(const std::vector<double>& values);
enum class CsvFormat {
  Basic,
  Extended
};

void write_csv_row(std::ostream& csv,
                   const std::string& timestamp,
                   const std::string& run_id,
                   const std::string& stage_tag,
                   std::size_t num_bits,
                   const SweepScenario& scenario,
                   const newcode::PipelineResult& result,
                   CsvFormat format);

std::vector<ScenarioOutput> run_scenarios_parallel(
    const std::vector<SweepScenario>& scenarios,
    const SweepParameterConfig& config,
    unsigned max_workers_hint = 0,
    const std::string& stage_tag = std::string{},
    DualOut* log = nullptr);

unsigned resolve_worker_count(const SweepParameterConfig& config);

}  // namespace detail
}  // namespace ofec_sweep
