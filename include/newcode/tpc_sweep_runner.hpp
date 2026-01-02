#pragma once

#include <string>
#include <vector>

#include "newcode/params.hpp"
#include "newcode/pipeline_runner.hpp"
#include "newcode/tpc_pipeline.hpp"

namespace tpc_sweep {

struct SweepParameterConfig {
  newcode::Params base_params{};
  unsigned bits_per_symbol = 1;
  int max_iters = 4;
  int num_blocks = 1;
  std::vector<float> alpha_schedule;
  std::vector<float> beta_schedule;

  int bitgen_seed_count = 0;
  int channel_seed_count = 0;
  std::vector<int> bitgen_seed_candidates;
  std::vector<int> channel_seed_candidates;

  float ebn0_start = newcode::DEFAULT_EBN0_DB;
  float ebn0_end = newcode::DEFAULT_EBN0_DB;
  int ebn0_points = 1;
  std::vector<float> ebn0_candidates;

  bool generate_random_bits = true;
  bool quiet_pipeline = true;
  bool quiet_logs = false;
};

int run_sweep(const SweepParameterConfig& config);

}  // namespace tpc_sweep
