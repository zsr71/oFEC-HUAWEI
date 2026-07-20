#include "newcode/tpc_single_runner.hpp"

#include <iostream>

namespace tpc_single {

int run_tpc_single(const Config& config) {
  newcode::Params params = detail::build_params(config);
  detail::log_run_overview(config, params);
  newcode::tpc::TpcPipelineConfig pipeline_cfg =
      detail::build_pipeline_config(config);
  const auto result =
      newcode::tpc::run_tpc_pipeline(params, pipeline_cfg, config.label,
                                     config.ebn0_db);
  detail::log_pipeline_results(result);
  return 0;
}

namespace detail {

newcode::Params build_params(const Config& cfg) {
  newcode::Params params;
  params.BITGEN_SEED = cfg.bitgen_seed;
  params.CHANNEL_SEED = cfg.channel_seed;
  params.BITGEN_RANDOM_BITS = cfg.generate_random_bits;
  params.LLR_BITS = cfg.llr_bits;
  return params;
}

newcode::tpc::TpcPipelineConfig build_pipeline_config(const Config& cfg) {
  newcode::tpc::TpcPipelineConfig pipeline_cfg;
  pipeline_cfg.bits_per_symbol = cfg.bits_per_symbol;
  pipeline_cfg.max_iters = cfg.max_iters;
  pipeline_cfg.num_blocks = cfg.num_blocks;
  pipeline_cfg.alpha_schedule = cfg.alpha_schedule;
  pipeline_cfg.beta_schedule = cfg.beta_schedule;
  pipeline_cfg.quiet = true;
  return pipeline_cfg;
}

void log_run_overview(const Config& cfg, const newcode::Params& params) {
  std::cout << "[INFO] run_tpc_pipeline(label=" << cfg.label
            << ", Eb/N0=" << cfg.ebn0_db << " dB)\n";
  std::cout << "[INFO] RNG seeds (bitgen/channel) = " << params.BITGEN_SEED
            << "/" << params.CHANNEL_SEED << "\n";
  std::cout << "[INFO] LLR bits = " << params.LLR_BITS
            << " (" << (params.LLR_BITS == 16 ? "float" : "qfloat") << ")\n";
  std::cout << "[INFO] TPC iters = " << cfg.max_iters
            << " | blocks=" << cfg.num_blocks << "\n";
  std::cout << "[INFO] TPC schedules: alpha=" << cfg.alpha_schedule.size()
            << ", beta=" << cfg.beta_schedule.size() << "\n";
}

void log_pipeline_results(const newcode::tpc::TpcPipelineResult& result) {
  std::cout << "[RESULT] Pre-FEC BER=" << result.pre_fec.ber
            << " (errs=" << result.pre_fec.errors << "/"
            << result.pre_fec.total << ")\n";
  std::cout << "[RESULT] Post-FEC BER=" << result.post_fec.ber
            << " (errs=" << result.post_fec.errors << "/"
            << result.post_fec.total << ")\n";
  std::cout << "[INFO] TPC iters=" << result.iters
            << ", avg_iters=" << result.avg_iters
            << ", early_stop=" << (result.early_stop ? "true" : "false")
            << ", early_stop_blocks=" << result.early_stop_blocks
            << "/" << result.blocks << "\n";
  if (!result.early_stop_start_pct.empty()) {
    std::cout << "[INFO] early_stop_start_pct: ";
    for (std::size_t i = 0; i < result.early_stop_start_pct.size(); ++i) {
      std::cout << result.early_stop_start_pct[i];
      if (i + 1 < result.early_stop_start_pct.size()) {
        std::cout << ", ";
      }
    }
    std::cout << "\n";
  }
  if (!result.per_block_ber_iter3.empty()) {
    std::cout << "[INFO] ber_after_iter3: ";
    for (std::size_t i = 0; i < result.per_block_ber_iter3.size(); ++i) {
      std::cout << result.per_block_ber_iter3[i];
      if (i + 1 < result.per_block_ber_iter3.size()) {
        std::cout << ", ";
      }
    }
    std::cout << "\n";
  }
  if (!result.per_block_early_stop_iter4.empty()) {
    std::cout << "[INFO] early_stop_at_iter4_start: ";
    for (std::size_t i = 0; i < result.per_block_early_stop_iter4.size(); ++i) {
      std::cout << result.per_block_early_stop_iter4[i];
      if (i + 1 < result.per_block_early_stop_iter4.size()) {
        std::cout << ", ";
      }
    }
    std::cout << "\n";
  }
  if (!result.per_block_ber_iter4.empty()) {
    std::cout << "[INFO] ber_after_iter4: ";
    for (std::size_t i = 0; i < result.per_block_ber_iter4.size(); ++i) {
      std::cout << result.per_block_ber_iter4[i];
      if (i + 1 < result.per_block_ber_iter4.size()) {
        std::cout << ", ";
      }
    }
    std::cout << "\n";
  }
}

}  // namespace detail
}  // namespace tpc_single
