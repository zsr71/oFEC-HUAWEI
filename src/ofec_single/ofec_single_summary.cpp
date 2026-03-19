#include "newcode/ofec_single_runner.hpp"
#include "newcode/io/dualwriter.hpp"
#include <algorithm>
#include <iomanip>
#include <sstream>

namespace ofec_single {
namespace {

constexpr std::size_t kDefaultPositionsToPrint = 10;

std::string format_positions(const std::vector<std::size_t>& positions,
                             std::size_t max_count = kDefaultPositionsToPrint) {
  std::ostringstream oss;
  const std::size_t count =
      (max_count == 0) ? positions.size()
                       : std::min(max_count, positions.size());
  oss << "[";
  for (std::size_t i = 0; i < count; ++i) {
    if (i) {
      oss << ", ";
    }
    oss << positions[i];
  }
  if (count < positions.size()) {
    if (count > 0) {
      oss << ", ";
    }
    oss << "... (+"
        << (positions.size() - count)
        << " more)";
  }
  oss << "]";
  return oss.str();
}

}  // namespace

namespace detail {

void log_run_overview(const Config& cfg,
                      const newcode::Params& params,
                      io::DualWriter& log) {
  log << "[INFO] run_pipeline(label=" << cfg.label
      << ", decoder=" << cfg.decoder_name
      << ", Eb/N0=" << cfg.ebn0_db
      << " dB, CHASE_L=" << params.CHASE_L
      << ", CHASE_NTEST=" << params.CHASE_NTEST
      << ", CHASE_TOPK_KEEP=" << params.CHASE_TOPK_KEEP
      << ", CHASE_GROUP_MINIMA_BITS=" << params.CHASE_GROUP_MINIMA_BITS
      << ")\n";
  log << "[INFO] RNG seeds (bitgen/channel) = "
      << params.BITGEN_SEED << "/"
      << params.CHANNEL_SEED << "\n";
  log << "[INFO] LLR bits = " << params.LLR_BITS
      << " (" << (params.LLR_BITS == 16 ? "float" : "qfloat") << ")\n";
  log << "[INFO] Dump quantized LLR = "
      << (cfg.dump_quantized_llr ? "ON" : "OFF");
  if (cfg.dump_quantized_llr) {
    log << " -> " << cfg.quantized_llr_output_path;
  }
  log << "\n";
}

void log_pipeline_results(const newcode::PipelineResult& result,
                          io::DualWriter& log) {
  log << "[RESULT] Pre-FEC BER=" << result.pre_fec.ber
      << " (errs=" << result.pre_fec.errors
      << "/" << result.pre_fec.total << ")";
  if (result.has_pre_fec_quantized_hard) {
    log << " | Pre-FEC BER (quantized hard)="
        << result.pre_fec_quantized_hard.ber
        << " (errs=" << result.pre_fec_quantized_hard.errors
        << "/" << result.pre_fec_quantized_hard.total << ")";
  }
  log << " | Post-FEC BER=" << result.post_fec.ber
      << " (errs=" << result.post_fec.errors
      << "/" << result.post_fec.total << ")\n";

  log << "[DETAIL] Pre-FEC error positions: "
      << format_positions(result.pre_fec_error_positions) << "\n";
  if (result.has_pre_fec_quantized_hard) {
    log << "[DETAIL] Pre-FEC error positions (quantized hard): "
        << format_positions(result.pre_fec_quantized_hard_error_positions) << "\n";
  }
  log << "[DETAIL] Post-FEC error positions: "
      << format_positions(result.post_fec_error_positions) << "\n";

  if (!result.tile_early_stop_pct.empty()) {
    std::ostringstream oss;
    oss.setf(std::ios::fixed);
    oss << std::setprecision(1);
    for (std::size_t i = 0; i < result.tile_early_stop_pct.size(); ++i) {
      oss << result.tile_early_stop_pct[i];
      if (i + 1 < result.tile_early_stop_pct.size()) {
        oss << ", ";
      }
    }
    log << "[RESULT] EarlyStop hit rates (%): "
        << oss.str() << "\n";
  }
  if (!result.tile_row_early_stop_pct.empty()) {
    std::ostringstream oss;
    oss.setf(std::ios::fixed);
    oss << std::setprecision(1);
    for (std::size_t i = 0; i < result.tile_row_early_stop_pct.size(); ++i) {
      oss << result.tile_row_early_stop_pct[i];
      if (i + 1 < result.tile_row_early_stop_pct.size()) {
        oss << ", ";
      }
    }
    log << "[RESULT] EarlyStop (per-row) hit rates (%): "
        << oss.str() << "\n";
  }

  if (!result.dequantized_llr_path.empty()) {
    log << "[INFO] Dequantized LLR saved to " << result.dequantized_llr_path << "\n";
  }
  if (!result.float_llr_path.empty()) {
    log << "[INFO] Float (pre-quant) LLR saved to " << result.float_llr_path << "\n";
  }
  if (!result.quantized_codes_path.empty()) {
    log << "[INFO] Quantized codes saved to " << result.quantized_codes_path << "\n";
  }
}

}  // namespace detail
}  // namespace ofec_single
