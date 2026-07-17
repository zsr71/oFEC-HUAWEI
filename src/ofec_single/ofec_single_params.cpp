#include "newcode/ofec_single_runner.hpp"
#include "newcode/io/dualwriter.hpp"
#include "newcode/ofec/mux/mux_config_validate.hpp"
#include "newcode/ofec/mux/mux_group_config_validate.hpp"
#include <algorithm>
#include <cmath>

namespace ofec_single {
namespace detail {

std::optional<newcode::Params> build_params(const Config& cfg,
                                            io::DualWriter& log) {
  newcode::Params params;
  params.BITGEN_SEED = cfg.bitgen_seed;
  params.CHANNEL_SEED = cfg.channel_seed;
  params.BITGEN_RANDOM_BITS = cfg.generate_random_bits;
  params.NORMALIZE_KNOWN_PREFIX_TAIL = cfg.normalize_known_prefix_tail;
  params.ENABLE_EARLY_STOP = cfg.enable_early_stop;
  params.EARLY_STOP_ENABLE_LIST = cfg.early_stop_enable_list;
  params.EARLY_STOP_CONDITION_MODE = cfg.early_stop_condition_mode;
  params.EARLY_STOP_CONDITION_MODE_LIST = cfg.early_stop_condition_mode_list;
  params.EARLY_STOP_ACTION_MODE = cfg.early_stop_action_mode;
  params.EARLY_STOP_ACTION_MODE_LIST = cfg.early_stop_action_mode_list;
  params.EARLY_STOP_BIND_GROUP_SIZE = cfg.early_stop_bind_group_size;
  params.EARLY_STOP_BIND_GROUP_SIZE_LIST = cfg.early_stop_bind_group_size_list;
  params.EARLY_STOP_COND_V1_REQUIRE_BCH = cfg.early_stop_cond_v1_require_bch;
  params.EARLY_STOP_COND_V1_REQUIRE_OVERALL = cfg.early_stop_cond_v1_require_overall;
  params.EARLY_STOP_V2_LLR_ABS_THRESHOLD = cfg.early_stop_v2_llr_abs_threshold;
  params.EARLY_STOP_V2_MAX_UNRELIABLE_BITS = cfg.early_stop_v2_max_unreliable_bits;
  params.EARLY_STOP_COND_V2_INCLUDE_OVERALL = cfg.early_stop_cond_v2_include_overall;
  params.EARLY_STOP_ACTION_RESIDUAL_DIVISOR =
      cfg.early_stop_action_residual_divisor;
  params.EARLY_STOP_ACTION_HARD_LLR_MAG =
      cfg.early_stop_action_hard_llr_mag;
  params.CHASE_TOPK_KEEP = cfg.chase_topk_keep;
  params.CHASE_GROUP_MINIMA_BITS = cfg.chase_group_minima_bits;
  params.debug_trace = cfg.debug_trace;
  params.LLR_BITS = cfg.llr_bits;
  params.DUMP_WORK_LLR = cfg.dump_work_llr;
  params.WORK_LLR_OUTPUT_PATH = cfg.work_llr_output_path;
  if (cfg.llr_bits < 2 || cfg.llr_bits > 16) {
    log << "[ERROR] LLR_BITS 必须在 [2,16]，16 表示浮点，其余使用 qfloat::qfloat<N>\n";
    return std::nullopt;
  }
  if (cfg.early_stop_condition_mode != 1 && cfg.early_stop_condition_mode != 2) {
    log << "[ERROR] early_stop_condition_mode 必须是 1 或 2\n";
    return std::nullopt;
  }
  if (cfg.early_stop_action_mode != 1 &&
      cfg.early_stop_action_mode != 2 &&
      cfg.early_stop_action_mode != 3 &&
      cfg.early_stop_action_mode != 4 &&
      cfg.early_stop_action_mode != 5 &&
      cfg.early_stop_action_mode != 6 &&
      cfg.early_stop_action_mode != 7 &&
      cfg.early_stop_action_mode != 8) {
    log << "[ERROR] early_stop_action_mode 目前必须是 1、2、3、4、5、6、7 或 8\n";
    return std::nullopt;
  }
  if (cfg.early_stop_v2_llr_abs_threshold < 0.0f) {
    log << "[ERROR] early_stop_v2_llr_abs_threshold 必须 >= 0\n";
    return std::nullopt;
  }
  if (cfg.early_stop_bind_group_size < 1) {
    log << "[ERROR] early_stop_bind_group_size 必须 >= 1\n";
    return std::nullopt;
  }
  for (int mode : cfg.early_stop_condition_mode_list) {
    if (mode != 1 && mode != 2) {
      log << "[ERROR] early_stop_condition_mode_list 的元素必须是 1 或 2\n";
      return std::nullopt;
    }
  }
  for (int mode : cfg.early_stop_action_mode_list) {
    if (mode != 1 && mode != 2 && mode != 3 &&
        mode != 4 && mode != 5 && mode != 6 &&
        mode != 7 && mode != 8) {
      log << "[ERROR] early_stop_action_mode_list 的元素必须是 1、2、3、4、5、6、7 或 8\n";
      return std::nullopt;
    }
  }
  for (int bind_group_size : cfg.early_stop_bind_group_size_list) {
    if (bind_group_size < 1) {
      log << "[ERROR] early_stop_bind_group_size_list 的元素必须 >= 1\n";
      return std::nullopt;
    }
  }
  if (cfg.early_stop_v2_max_unreliable_bits < 0 ||
      cfg.early_stop_v2_max_unreliable_bits > static_cast<int>(newcode::Params::BCH_N)) {
    log << "[ERROR] early_stop_v2_max_unreliable_bits 必须在 [0, BCH_N] 范围内\n";
    return std::nullopt;
  }
  if (!(cfg.early_stop_action_residual_divisor > 0.0f)) {
    log << "[ERROR] early_stop_action_residual_divisor 必须 > 0\n";
    return std::nullopt;
  }
  if (!std::isfinite(cfg.early_stop_action_hard_llr_mag)) {
    log << "[ERROR] early_stop_action_hard_llr_mag 必须是有限数\n";
    return std::nullopt;
  }
  if (cfg.chase_topk_keep < 1) {
    log << "[ERROR] chase_topk_keep 必须 >= 1\n";
    return std::nullopt;
  }
  if (cfg.chase_n_test_override == 0 || cfg.chase_n_test_override < -1) {
    log << "[ERROR] chase_n_test_override 必须是 -1 或 >= 1\n";
    return std::nullopt;
  }
  if (cfg.chase_group_minima_bits < 0) {
    log << "[ERROR] chase_group_minima_bits 必须 >= 0\n";
    return std::nullopt;
  }
  params.LLR_CLIP_RATIO = std::clamp(cfg.quant_clip_ratio, 0.0f, 1.0f);
  if (cfg.chaseL_override >= 0) {
    params.CHASE_L = cfg.chaseL_override;
  }
  if (cfg.chase_n_test_override >= 0) {
    params.CHASE_NTEST = cfg.chase_n_test_override;
  } else if (cfg.chaseL_override >= 0) {
    params.CHASE_NTEST = 1 << params.CHASE_L;
  }

  const std::size_t tiles = params.TILES_PER_WIN;

  if (!params.EARLY_STOP_ENABLE_LIST.empty() &&
      params.EARLY_STOP_ENABLE_LIST.size() != tiles) {
    log << "[ERROR] early_stop_enable_list 长度必须等于 TILES_PER_WIN\n";
    return std::nullopt;
  }
  if (!params.EARLY_STOP_CONDITION_MODE_LIST.empty() &&
      params.EARLY_STOP_CONDITION_MODE_LIST.size() != tiles) {
    log << "[ERROR] early_stop_condition_mode_list 长度必须等于 TILES_PER_WIN\n";
    return std::nullopt;
  }
  if (!params.EARLY_STOP_ACTION_MODE_LIST.empty() &&
      params.EARLY_STOP_ACTION_MODE_LIST.size() != tiles) {
    log << "[ERROR] early_stop_action_mode_list 长度必须等于 TILES_PER_WIN\n";
    return std::nullopt;
  }
  if (!params.EARLY_STOP_BIND_GROUP_SIZE_LIST.empty() &&
      params.EARLY_STOP_BIND_GROUP_SIZE_LIST.size() != tiles) {
    log << "[ERROR] early_stop_bind_group_size_list 长度必须等于 TILES_PER_WIN\n";
    return std::nullopt;
  }

  if (!cfg.alpha_explicit.empty()) {
    if (cfg.alpha_explicit.size() != tiles) {
      log << "[ERROR] kAlpha_explicit 长度必须等于 TILES_PER_WIN\n";
      return std::nullopt;
    }
    params.ALPHA_LIST = cfg.alpha_explicit;
  } else {
    params.ALPHA_LIST.assign(tiles, cfg.alpha_fill);
  }
  if (!params.ALPHA_LIST.empty()) {
    params.ALPHA = params.ALPHA_LIST.front();
  }

  if (!cfg.beta_explicit.empty()) {
    if (cfg.beta_explicit.size() != tiles) {
      log << "[ERROR] kBeta_explicit 长度必须等于 TILES_PER_WIN\n";
      return std::nullopt;
    }
    params.beta_list = cfg.beta_explicit;
  } else {
    params.beta_list.assign(tiles, cfg.beta_fill);
  }
  if (!params.beta_list.empty()) {
    params.beta = params.beta_list.front();
  }

  if (!cfg.early_stop_action_sign_beta_explicit.empty()) {
    if (cfg.early_stop_action_sign_beta_explicit.size() != tiles) {
      log << "[ERROR] early_stop_action_sign_beta_explicit 长度必须等于 TILES_PER_WIN\n";
      return std::nullopt;
    }
    params.EARLY_STOP_ACTION_SIGN_BETA_LIST =
        cfg.early_stop_action_sign_beta_explicit;
  } else if (std::isfinite(cfg.early_stop_action_sign_beta_fill)) {
    params.EARLY_STOP_ACTION_SIGN_BETA_LIST.assign(
        tiles, cfg.early_stop_action_sign_beta_fill);
  } else {
    // 兼容旧行为：如果顶层没有单独配置 early-stop beta，则沿用 Chase 的 beta_list。
    params.EARLY_STOP_ACTION_SIGN_BETA_LIST = params.beta_list;
  }
  if (!params.EARLY_STOP_ACTION_SIGN_BETA_LIST.empty()) {
    params.EARLY_STOP_ACTION_SIGN_BETA =
        params.EARLY_STOP_ACTION_SIGN_BETA_LIST.front();
  } else {
    params.EARLY_STOP_ACTION_SIGN_BETA = params.beta;
  }

  if (!cfg.siso_active_list.empty()) {
  params.SISO_ACTIVE_LIST = cfg.siso_active_list;
  }
  if (!cfg.hiho_active_list.empty()) {
    params.HIHO_ACTIVE_LIST = cfg.hiho_active_list;
  }
  params.MUX_GROUP_G = cfg.mux_group_g;
  params.MUX_GROUP_G_LIST = cfg.mux_group_g_list;
  params.MUX_SCHEDULING_MODE = cfg.mux_scheduling_mode;
  params.MUX_EARLY_STOP_PRIORITY_RULE = cfg.mux_early_stop_priority_rule;
  params.MUX_ENABLE_RECONFIG = cfg.mux_enable_reconfig;
  params.MUX_EXTRA_BYPASS_EDGES = cfg.mux_extra_bypass_edges;
  params.HYBRID_ENABLE = cfg.hybrid_enable;
  params.HYBRID_ENABLE_LIST = cfg.hybrid_enable_list;
  params.HYBRID_HARD_LLR_MAG = cfg.hybrid_hard_llr_mag;
  params.HYBRID_HARD_LLR_MAG_LIST = cfg.hybrid_hard_llr_mag_list;
  params.HYBRID_CLASSIFIER_MODE = cfg.hybrid_classifier_mode;
  params.HYBRID_SISO_BACKFILL_MODE = cfg.hybrid_siso_backfill_mode;
  params.HYBRID_USE_FAST_CLASSIFIER =
      cfg.hybrid_classifier_mode != newcode::HybridClassifierMode::LegacyHardDecode;
  params.HYBRID_NORMALIZE_SOFT_ONLY = cfg.hybrid_normalize_soft_only;
  if (!params.HYBRID_ENABLE_LIST.empty() &&
      params.HYBRID_ENABLE_LIST.size() != tiles) {
    log << "[ERROR] hybrid_enable_list 长度必须等于 TILES_PER_WIN\n";
    return std::nullopt;
  }
  if (!std::isfinite(params.HYBRID_HARD_LLR_MAG)) {
    log << "[ERROR] hybrid_hard_llr_mag 必须是有限数\n";
    return std::nullopt;
  }
  if (!params.HYBRID_HARD_LLR_MAG_LIST.empty() &&
      params.HYBRID_HARD_LLR_MAG_LIST.size() != tiles) {
    log << "[ERROR] hybrid_hard_llr_mag_list 长度必须等于 TILES_PER_WIN\n";
    return std::nullopt;
  }
  for (float hard_mag : params.HYBRID_HARD_LLR_MAG_LIST) {
    if (!std::isfinite(hard_mag)) {
      log << "[ERROR] hybrid_hard_llr_mag_list 的元素必须是有限数\n";
      return std::nullopt;
    }
  }
  const auto mux_ok =
      newcode::mux::validate_siso_active_list(params.SISO_ACTIVE_LIST, tiles);
  if (!mux_ok.ok) {
    log << "[ERROR] " << mux_ok.error << "\n";
    return std::nullopt;
  }
  const auto hiho_ok =
      newcode::mux::validate_hiho_active_list(params.HIHO_ACTIVE_LIST, tiles);
  if (!hiho_ok.ok) {
    log << "[ERROR] " << hiho_ok.error << "\n";
    return std::nullopt;
  }
  const std::size_t rows_to_decode =
      static_cast<std::size_t>(params.CHASE_SBR) *
      newcode::Params::BITS_PER_SUBBLOCK_DIM;
  if (!params.MUX_GROUP_G_LIST.empty() &&
      params.MUX_GROUP_G_LIST.size() != tiles) {
    log << "[ERROR] mux_group_g_list 长度必须等于 TILES_PER_WIN\n";
    return std::nullopt;
  }
  for (std::size_t tile_idx = 0; tile_idx < tiles; ++tile_idx) {
    const int group_g = newcode::mux::pick_mux_group_g_for_tile(
        params.MUX_GROUP_G_LIST, tile_idx, params.MUX_GROUP_G);
    const auto group_ok =
        newcode::mux::validate_mux_group_g(group_g, rows_to_decode);
    if (!group_ok.ok) {
      log << "[ERROR] tile " << tile_idx << ": " << group_ok.error << "\n";
      return std::nullopt;
    }
  }
  if (params.MUX_SCHEDULING_MODE != 0 && params.MUX_SCHEDULING_MODE != 1) {
    log << "[ERROR] mux_scheduling_mode 目前必须是 0 或 1\n";
    return std::nullopt;
  }
  if (params.MUX_EARLY_STOP_PRIORITY_RULE != 0 &&
      params.MUX_EARLY_STOP_PRIORITY_RULE != 1) {
    log << "[ERROR] mux_early_stop_priority_rule 目前必须是 0 或 1\n";
    return std::nullopt;
  }
  if (params.MUX_ENABLE_RECONFIG) {
    for (std::size_t tile_idx = 0; tile_idx < params.SISO_ACTIVE_LIST.size(); ++tile_idx) {
      const int group_g = newcode::mux::pick_mux_group_g_for_tile(
          params.MUX_GROUP_G_LIST, tile_idx, params.MUX_GROUP_G);
      const auto reconfig_ok = newcode::mux::validate_mux_reconfig_runtime(
          group_g, params.SISO_ACTIVE_LIST[tile_idx], rows_to_decode);
      if (!reconfig_ok.ok) {
        log << "[ERROR] tile " << tile_idx << ": " << reconfig_ok.error << "\n";
        return std::nullopt;
      }
    }
  }

  return params;
}

newcode::PipelineConfig build_pipeline_config(const Config& cfg) {
  newcode::PipelineConfig pipeline_cfg;
  pipeline_cfg.decoder_name = cfg.decoder_name;
  pipeline_cfg.interleaver_name = cfg.interleaver_name;
  pipeline_cfg.normalize_extrinsic = cfg.normalize_extrinsic;
  pipeline_cfg.bits_per_symbol = cfg.bits_per_symbol;
  pipeline_cfg.dump_quantized_llr = cfg.dump_quantized_llr;
  pipeline_cfg.quantized_llr_output_path = cfg.quantized_llr_output_path;
  pipeline_cfg.dump_work_llr = cfg.dump_work_llr;
  pipeline_cfg.work_llr_output_path = cfg.work_llr_output_path;
  if (pipeline_cfg.dump_quantized_llr && pipeline_cfg.quantized_llr_output_path.empty()) {
    pipeline_cfg.quantized_llr_output_path = "data/llr/" + cfg.label + "_quantized_llr.txt";
  }
  if (pipeline_cfg.dump_work_llr && pipeline_cfg.work_llr_output_path.empty()) {
    pipeline_cfg.work_llr_output_path = "data/llr/" + cfg.label + "_work_llr.txt";
  }
  return pipeline_cfg;
}

}  // namespace detail
}  // namespace ofec_single
