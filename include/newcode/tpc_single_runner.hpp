#pragma once

#include <string>
#include <vector>

#include "newcode/params.hpp"
#include "newcode/tpc_pipeline.hpp"

namespace tpc_single {

struct Config {
  std::string label;
  float ebn0_db;
  int max_iters;
  int num_blocks;
  unsigned bits_per_symbol;
  int bitgen_seed;
  int channel_seed;
  bool generate_random_bits = true;
  std::size_t llr_bits = 16;
  std::vector<float> alpha_schedule;
  std::vector<float> beta_schedule;
};

int run_tpc_single(const Config& config);

namespace detail {

newcode::Params build_params(const Config& cfg);

newcode::tpc::TpcPipelineConfig build_pipeline_config(const Config& cfg);

void log_run_overview(const Config& cfg, const newcode::Params& params);

void log_pipeline_results(const newcode::tpc::TpcPipelineResult& result);

}  // namespace detail

}  // namespace tpc_single
