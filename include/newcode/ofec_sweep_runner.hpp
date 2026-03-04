#pragma once

#include <string>
#include <vector>

#include "newcode/params.hpp"
#include "newcode/pipeline_runner.hpp"

namespace ofec_sweep {

// 显式 α/β 参数模式
struct ExplicitAlphaBetaPattern {
  std::string label;
  std::vector<float> alpha_list;
  std::vector<float> beta_list;
};

// 配置结构体，包含所有可扫参数及基础参数
struct SweepParameterConfig {
  newcode::Params base_params{};
  std::string interleaver_name = "identity";
  std::string decoder_name = "plain";
  bool normalize_extrinsic = true;
  unsigned bits_per_symbol = 2;

  std::vector<float> alpha_start_candidates;
  std::vector<float> alpha_step_candidates;
  std::vector<float> beta_start_candidates;
  std::vector<float> beta_step_candidates;
  std::vector<int> siso_active_list;
  int mux_group_g = 1;
  std::vector<int> chase_l_candidates;

  int bitgen_seed_count = 0;
  int channel_seed_count = 0;
  std::vector<int> bitgen_seed_candidates;
  std::vector<int> channel_seed_candidates;

  float ebn0_start = newcode::DEFAULT_EBN0_DB;
  float ebn0_end = newcode::DEFAULT_EBN0_DB;
  int ebn0_points = 1;
  std::vector<float> ebn0_candidates;

  std::vector<ExplicitAlphaBetaPattern> explicit_patterns;

  unsigned max_workers_override = 0;
  bool generate_random_bits = true;
  bool normalize_known_prefix_tail = true;
  float quant_clip_ratio = 0.0f;
  bool quiet_pipeline = true;
  bool quiet_logs = false;
};

int run_sweep(const SweepParameterConfig& config);

}  // namespace ofec_sweep
