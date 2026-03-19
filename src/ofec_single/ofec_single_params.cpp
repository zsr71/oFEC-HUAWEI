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
  params.EARLY_STOP_CONDITION_MODE = cfg.early_stop_condition_mode;
  params.EARLY_STOP_ACTION_MODE = cfg.early_stop_action_mode;
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
      cfg.early_stop_action_mode != 3) {
    log << "[ERROR] early_stop_action_mode 目前必须是 1、2 或 3\n";
    return std::nullopt;
  }
  if (cfg.early_stop_v2_llr_abs_threshold < 0.0f) {
    log << "[ERROR] early_stop_v2_llr_abs_threshold 必须 >= 0\n";
    return std::nullopt;
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
  params.MUX_GROUP_G = cfg.mux_group_g;
  params.MUX_ENABLE_RECONFIG = cfg.mux_enable_reconfig;
  params.MUX_EXTRA_BYPASS_EDGES = cfg.mux_extra_bypass_edges;
  const auto mux_ok =
      newcode::mux::validate_siso_active_list(params.SISO_ACTIVE_LIST, tiles);
  if (!mux_ok.ok) {
    log << "[ERROR] " << mux_ok.error << "\n";
    return std::nullopt;
  }
  const std::size_t rows_to_decode =
      static_cast<std::size_t>(params.CHASE_SBR) *
      newcode::Params::BITS_PER_SUBBLOCK_DIM;
  const auto group_ok =
      newcode::mux::validate_mux_group_g(params.MUX_GROUP_G, rows_to_decode);
  if (!group_ok.ok) {
    log << "[ERROR] " << group_ok.error << "\n";
    return std::nullopt;
  }
  if (params.MUX_ENABLE_RECONFIG) {
    for (std::size_t tile_idx = 0; tile_idx < params.SISO_ACTIVE_LIST.size(); ++tile_idx) {
      const auto reconfig_ok = newcode::mux::validate_mux_reconfig_runtime(
          params.MUX_GROUP_G, params.SISO_ACTIVE_LIST[tile_idx], rows_to_decode);
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
