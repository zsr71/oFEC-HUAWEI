#include "newcode/ofec_single_runner.hpp"

namespace ofec_single {
namespace detail {

std::optional<newcode::Params> build_params(const Config& cfg,
                                            DualWriter& log) {
  newcode::Params params;
  params.BITGEN_SEED = cfg.bitgen_seed;
  params.CHANNEL_SEED = cfg.channel_seed;
  params.BITGEN_RANDOM_BITS = cfg.generate_random_bits;
  params.NORMALIZE_KNOWN_PREFIX_TAIL = cfg.normalize_known_prefix_tail;
  params.debug_trace = cfg.debug_trace;
  if (cfg.chaseL_override >= 0) {
    params.CHASE_L = cfg.chaseL_override;
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

  return params;
}

newcode::PipelineConfig build_pipeline_config(const Config& cfg) {
  newcode::PipelineConfig pipeline_cfg;
  pipeline_cfg.decoder_name = cfg.decoder_name;
  pipeline_cfg.interleaver_name = cfg.interleaver_name;
  pipeline_cfg.normalize_extrinsic = cfg.normalize_extrinsic;
  pipeline_cfg.bits_per_symbol = cfg.bits_per_symbol;
  return pipeline_cfg;
}

}  // namespace detail
}  // namespace ofec_single
