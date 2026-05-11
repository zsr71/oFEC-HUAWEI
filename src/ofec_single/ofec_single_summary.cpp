#include "newcode/ofec_single_runner.hpp"
#include "newcode/ofec/hybrid/hybrid_classifier.hpp"
#include "newcode/io/dualwriter.hpp"
#include <algorithm>
#include <iomanip>
#include <sstream>

namespace ofec_single {
namespace {

constexpr std::size_t kDefaultPositionsToPrint = 10;

const char* hybrid_classifier_mode_name(newcode::HybridClassifierMode mode) {
  switch (mode) {
    case newcode::HybridClassifierMode::LegacyHardDecode:
      return "legacy_hard_decode";
    case newcode::HybridClassifierMode::RepoFastClassifier:
      return "repo_fast_classifier";
    case newcode::HybridClassifierMode::FriendS1S3Classifier:
      return "friend_s1s3_classifier";
  }
  return "unknown";
}

std::string format_compact_int_list(const std::vector<int>& values) {
  std::ostringstream oss;
  oss << "[";
  for (std::size_t i = 0; i < values.size(); ++i) {
    if (i) {
      oss << ", ";
    }
    oss << values[i];
  }
  oss << "]";
  return oss.str();
}

std::string format_compact_size_t_list(const std::vector<std::size_t>& values) {
  std::ostringstream oss;
  oss << "[";
  for (std::size_t i = 0; i < values.size(); ++i) {
    if (i) {
      oss << ", ";
    }
    oss << values[i];
  }
  oss << "]";
  return oss.str();
}

std::string format_compact_double_list(const std::vector<double>& values,
                                       int precision = 1) {
  std::ostringstream oss;
  oss.setf(std::ios::fixed);
  oss << std::setprecision(precision);
  oss << "[";
  for (std::size_t i = 0; i < values.size(); ++i) {
    if (i) {
      oss << ", ";
    }
    oss << values[i];
  }
  oss << "]";
  return oss.str();
}

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
  log << "[INFO] Early-stop = " << (params.ENABLE_EARLY_STOP ? "ON" : "OFF")
      << ", per-tile enable list = "
      << format_compact_int_list(params.EARLY_STOP_ENABLE_LIST)
      << ", per-tile condition list = "
      << format_compact_int_list(params.EARLY_STOP_CONDITION_MODE_LIST)
      << ", per-tile action list = "
      << format_compact_int_list(params.EARLY_STOP_ACTION_MODE_LIST)
      << ", per-tile bind list = "
      << format_compact_int_list(params.EARLY_STOP_BIND_GROUP_SIZE_LIST)
      << ", bind_group_size = " << params.EARLY_STOP_BIND_GROUP_SIZE << "\n";
  log << "[INFO] MUX group/scheduling/rule = "
      << params.MUX_GROUP_G << "/"
      << params.MUX_SCHEDULING_MODE << "/"
      << params.MUX_EARLY_STOP_PRIORITY_RULE
      << ", reconfig = " << (params.MUX_ENABLE_RECONFIG ? "ON" : "OFF") << "\n";
  log << "[INFO] Hybrid prepass = "
      << (params.HYBRID_ENABLE ? "ON" : "OFF")
      << ", per-tile enable list = "
      << format_compact_int_list(params.HYBRID_ENABLE_LIST)
      << ", classifier_mode = "
      << hybrid_classifier_mode_name(
             newcode::detail::effective_hybrid_classifier_mode(params))
      << ", normalize_soft_only = "
      << (params.HYBRID_NORMALIZE_SOFT_ONLY ? "ON" : "OFF") << "\n";
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
    log << "[RESULT] EarlyStop (per-row) hit rates (%): "
        << format_compact_double_list(result.tile_row_early_stop_pct) << "\n";
  }
  if (!result.tile_hard_finish_count.empty()) {
    log << "[RESULT] Hybrid hard-finish counts: "
        << format_compact_size_t_list(result.tile_hard_finish_count) << "\n";
  }
  if (!result.tile_hard_finish_pct.empty()) {
    log << "[RESULT] Hybrid hard-finish rates (% of rows): "
        << format_compact_double_list(result.tile_hard_finish_pct) << "\n";
  }

  if (!result.tile_unscheduled_count.empty()) {
    log << "[RESULT] Unscheduled counts: "
        << format_compact_size_t_list(result.tile_unscheduled_count) << "\n";
  }
  if (!result.tile_unscheduled_pct.empty()) {
    log << "[RESULT] Unscheduled rates (% of rows): "
        << format_compact_double_list(result.tile_unscheduled_pct) << "\n";
  }
  if (!result.tile_unscheduled_among_need_pct.empty()) {
    log << "[RESULT] Unscheduled rates (% of NeedSiso rows): "
        << format_compact_double_list(result.tile_unscheduled_among_need_pct) << "\n";
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
