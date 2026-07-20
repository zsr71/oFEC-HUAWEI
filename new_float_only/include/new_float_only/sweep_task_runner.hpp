#pragma once

#include <cstddef>
#include <string>
#include <vector>

#include "new_float_only/params.hpp"
#include "new_float_only/rx/ber/ber.hpp"

namespace new_float_only {

struct SweepPattern {
  std::size_t pattern_index = 0;
  std::string pattern_label;
  std::vector<float> alpha_list;
  std::vector<float> beta_list;
};

struct SweepSeedPair {
  std::size_t seed_index = 0;
  int bitgen_seed = 0;
  int channel_seed = 0;
};

struct SweepTask {
  std::size_t task_index = 0;
  std::size_t pattern_index = 0;
  std::size_t seed_index = 0;
  std::string pattern_label;
  std::vector<float> alpha_list;
  std::vector<float> beta_list;
  int bitgen_seed = 0;
  int channel_seed = 0;
};

struct SweepTaskRunnerConfig {
  std::string label = "sweep";
  float ebn0_db = 3.24f;
  unsigned bits_per_symbol = 2;
  bool generate_random_bits = true;
  bool normalize_extrinsic = true;
  Params base_decoder{};

  unsigned max_workers_override = 0;
  bool quiet_pipeline = true;
  bool quiet_logs = false;
};

struct SweepTaskResult {
  std::size_t task_index = 0;
  std::size_t pattern_index = 0;
  std::size_t seed_index = 0;
  std::string pattern_label;
  std::vector<float> alpha_list;
  std::vector<float> beta_list;
  int bitgen_seed = 0;
  int channel_seed = 0;
  BerStats pre_fec;
  BerStats post_fec;
};

struct SweepPatternSummary {
  std::size_t pattern_index = 0;
  std::string pattern_label;
  std::vector<float> alpha_list;
  std::vector<float> beta_list;
  std::size_t trials_completed = 0;
  std::size_t pre_frame_error_trials = 0;
  std::size_t post_frame_error_trials = 0;
  BerStats pre_fec;
  BerStats post_fec;
};

Params normalize_sweep_decoder_config(const DecoderConfig& decoder);
unsigned resolve_sweep_worker_count(const SweepTaskRunnerConfig& config);
std::vector<SweepSeedPair> build_sweep_seed_schedule(std::size_t trial_count,
                                                     int bitgen_seed_base,
                                                     int channel_seed_base,
                                                     const std::vector<int>& bitgen_seeds = {},
                                                     const std::vector<int>& channel_seeds = {});
std::vector<SweepTask> build_sweep_tasks(const std::vector<SweepPattern>& patterns,
                                         const std::vector<SweepSeedPair>& seed_schedule);
std::vector<SweepTaskResult> run_sweep_tasks(const SweepTaskRunnerConfig& config,
                                             const std::vector<SweepTask>& tasks);
std::vector<SweepPatternSummary> aggregate_sweep_task_results(
    const std::vector<SweepPattern>& patterns,
    const std::vector<SweepTaskResult>& task_results);

}  // namespace new_float_only
