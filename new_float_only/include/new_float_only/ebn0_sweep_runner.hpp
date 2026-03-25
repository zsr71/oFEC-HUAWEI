#pragma once

#include <cstddef>
#include <string>
#include <vector>

#include "new_float_only/params.hpp"
#include "new_float_only/rx/ber/ber.hpp"

namespace new_float_only {

struct Ebn0SweepTrialResult {
  std::size_t task_index = 0;
  std::size_t ebn0_index = 0;
  std::size_t trial_index = 0;
  float ebn0_db = 0.0f;
  int bitgen_seed = 0;
  int channel_seed = 0;
  BerStats pre_fec;
  BerStats post_fec;
};

struct Ebn0SweepPointResult {
  std::size_t ebn0_index = 0;
  float ebn0_db = 0.0f;
  std::size_t trials_requested = 0;
  std::size_t trials_completed = 0;
  BerStats pre_fec;
  BerStats post_fec;
  std::size_t pre_frame_error_trials = 0;
  std::size_t post_frame_error_trials = 0;
};

struct Ebn0SweepConfig {
  std::string label = "ofec_ebn0_sweep_float";
  float ebn0_start = 3.24f;
  float ebn0_end = 3.24f;
  int ebn0_points = 1;
  unsigned bits_per_symbol = 2;
  bool generate_random_bits = true;
  bool normalize_extrinsic = true;
  DecoderConfig decoder{};

  std::size_t trial_count = 0;
  int bitgen_seed_base = 56456;
  int channel_seed_base = 57112;
  std::vector<int> bitgen_seeds;
  std::vector<int> channel_seeds;

  unsigned max_workers_override = 0;
  bool quiet_pipeline = true;
  bool quiet_logs = false;
  bool write_summary_csv = true;
  bool write_trial_csv = true;
};

struct Ebn0SweepResult {
  std::size_t points_requested = 0;
  std::size_t points_completed = 0;
  std::vector<float> ebn0_values;
  std::vector<Ebn0SweepPointResult> point_results;
  std::vector<Ebn0SweepTrialResult> trial_results;
  std::string summary_csv_path;
  std::string trial_csv_path;
};

Ebn0SweepResult run_ebn0_sweep(const Ebn0SweepConfig& config);

}  // namespace new_float_only
