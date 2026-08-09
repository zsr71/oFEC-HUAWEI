#pragma once

#include <filesystem>
#include <fstream>
#include <iostream>
#include <limits>
#include <optional>
#include <string>
#include <vector>

#include "newcode/params.hpp"
#include "newcode/pipeline_runner.hpp"
#include "newcode/io/dualwriter.hpp"
namespace ofec_single {

struct Config {
  std::string label;
  float ebn0_db;
  int chaseL_override;
  int chase_n_test_override = -1;
  int chase_topk_keep = 8;
  int chase_group_minima_bits = 3;
  bool normalize_extrinsic;
  unsigned bits_per_symbol;
  int bitgen_seed;
  int channel_seed;
  bool enable_early_stop = true;
  std::vector<int> early_stop_enable_list;
  int early_stop_condition_mode = 1;
  std::vector<int> early_stop_condition_mode_list;
  int early_stop_action_mode = 1;
  std::vector<int> early_stop_action_mode_list;
  int early_stop_bind_group_size = 1;
  std::vector<int> early_stop_bind_group_size_list;
  bool early_stop_cond_v1_require_bch = true;
  bool early_stop_cond_v1_require_overall = true;
  float early_stop_v2_llr_abs_threshold = 0.5f;
  int early_stop_v2_max_unreliable_bits = 8;
  bool early_stop_cond_v2_include_overall = true;
  float alpha_fill = 1.0f;
  float beta_fill = 0.35f;
  std::vector<float> alpha_explicit;
  std::vector<float> beta_explicit;
  float early_stop_action_sign_beta_fill =
      std::numeric_limits<float>::quiet_NaN();
  std::vector<float> early_stop_action_sign_beta_explicit;
  float early_stop_action_residual_divisor = 1.0f;
  float early_stop_action_hard_llr_mag = 1.0f;
  std::vector<int> siso_active_list;
  std::vector<int> hiho_active_list;
  int mux_group_g = 1;
  int mux_scheduling_mode = 0;
  int mux_early_stop_priority_rule = 0;
  bool mux_enable_reconfig = false;
  std::vector<newcode::mux::MuxEdge> mux_extra_bypass_edges;
  bool hybrid_enable = false;
  std::vector<int> hybrid_enable_list;
  float hybrid_hard_llr_mag = 99.0f;
  std::vector<float> hybrid_hard_llr_mag_list;
  newcode::HybridClassifierMode hybrid_classifier_mode =
      newcode::HybridClassifierMode::LegacyHardDecode;
  newcode::HybridSisoBackfillMode hybrid_siso_backfill_mode =
      newcode::HybridSisoBackfillMode::Disabled;
  bool hybrid_normalize_soft_only = false;
  bool level56_shared_enable = false;
  bool level56_temporal_lookahead_enable = false;
  int level56_temporal_group_load_threshold = 16;
  int level56_shared_hiso_active = 32;
  int level56_shared_siso_active = 32;
  newcode::Level56PriorityMode level56_priority_mode =
      newcode::Level56PriorityMode::Level5First;
  newcode::Level56ScheduleMode level56_schedule_mode =
      newcode::Level56ScheduleMode::GlobalPriority;
  newcode::Level56EarlyStopGroupUpdateMode level56_early_stop_group_update_mode =
      newcode::Level56EarlyStopGroupUpdateMode::AllGroups;
  bool level56_single_level_select_enable = false;
  bool level56_unselected_early_stop_action_enable = false;
  bool dump_level56_schedule_stats = false;
  std::string level56_schedule_rounds_output_path;
  std::string level56_schedule_codes_output_path;
  std::string interleaver_name;
  std::string decoder_name;
  bool generate_random_bits = true;
  bool normalize_known_prefix_tail = true;
  float quant_clip_ratio = 0.0f;
  std::size_t llr_bits = 16;
  bool dump_quantized_llr = false;
  std::string quantized_llr_output_path;
  bool dump_work_llr = false;
  std::string work_llr_output_path;
  bool dump_tile_early_stop_samples = false;
  std::string tile_early_stop_samples_output_path;
  bool dump_tile_early_stop_group_bind_debug_samples = false;
  std::string tile_early_stop_group_bind_debug_samples_output_path;
  newcode::Params::DebugTraceConfig debug_trace;
};



int run_ofec_single(const Config& config);

namespace detail {



std::optional<newcode::Params> build_params(const Config& cfg,
                                            io::DualWriter& log);

newcode::PipelineConfig build_pipeline_config(const Config& cfg);

void log_run_overview(const Config& cfg,
                      const newcode::Params& params,
                      io::DualWriter& log);

void log_pipeline_results(const newcode::PipelineResult& result,
                          io::DualWriter& log);

} 

}  
