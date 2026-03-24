#pragma once

#include <limits>
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
  std::vector<float> early_stop_action_sign_beta_list;
};

// 配置结构体，包含所有可扫参数及基础参数
struct SweepParameterConfig {
  newcode::Params base_params{};
  std::string interleaver_name = "identity";
  std::string decoder_name = "chase_baseline";
  std::vector<std::string> decoder_name_candidates;
  bool enable_early_stop = true;
  int early_stop_condition_mode = 1;
  int early_stop_action_mode = 1;
  bool early_stop_cond_v1_require_bch = true;
  bool early_stop_cond_v1_require_overall = true;
  float early_stop_v2_llr_abs_threshold = 0.5f;
  int early_stop_v2_max_unreliable_bits = 8;
  bool early_stop_cond_v2_include_overall = true;
  float early_stop_action_sign_beta_fill =
      std::numeric_limits<float>::quiet_NaN();
  float early_stop_action_residual_divisor = 1.0f;
  float early_stop_action_hard_llr_mag = 1.0f;
  int chase_n_test = 64;
  int chase_topk_keep = 8;
  int chase_group_minima_bits = 3;
  bool normalize_extrinsic = true;
  unsigned bits_per_symbol = 2;

  std::vector<float> alpha_start_candidates;
  std::vector<float> alpha_step_candidates;
  std::vector<float> beta_start_candidates;
  std::vector<float> beta_step_candidates;
  std::vector<float> early_stop_action_beta_start_candidates;
  std::vector<float> early_stop_action_beta_step_candidates;
  std::vector<float> early_stop_action_hard_llr_mag_candidates;
  std::vector<int> early_stop_condition_candidates;
  std::vector<int> early_stop_action_candidates;
  std::vector<bool> early_stop_cond_v1_require_bch_candidates;
  std::vector<bool> early_stop_cond_v1_require_overall_candidates;
  std::vector<float> early_stop_v2_llr_abs_threshold_candidates;
  std::vector<int> early_stop_v2_max_unreliable_bits_candidates;
  std::vector<int> siso_active_list;
  int mux_group_g = 1;
  bool mux_enable_reconfig = false;
  int mux_bypass_scheme = 0;
  std::vector<newcode::mux::MuxEdge> mux_extra_bypass_edges;
  std::vector<int> chase_l_candidates;
  std::vector<int> chase_n_test_candidates;
  std::vector<int> chase_topk_keep_candidates;
  std::vector<int> chase_group_minima_bits_candidates;

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
