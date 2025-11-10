#include "newcode/ofec_single_runner.hpp"

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
                      DualWriter& log) {
  log << "[INFO] run_pipeline(label=" << cfg.label
      << ", Eb/N0=" << cfg.ebn0_db
      << " dB, CHASE_L=" << params.CHASE_L << ")\n";
  log << "[INFO] RNG seeds (bitgen/channel) = "
      << params.BITGEN_SEED << "/"
      << params.CHANNEL_SEED << "\n";
}

void log_pipeline_results(const newcode::PipelineResult& result,
                          DualWriter& log) {
  log << "[RESULT] Pre-FEC BER=" << result.pre_fec.ber
      << " (errs=" << result.pre_fec.errors
      << "/" << result.pre_fec.total << ")";
  log << " | Post-FEC BER=" << result.post_fec.ber
      << " (errs=" << result.post_fec.errors
      << "/" << result.post_fec.total << ")\n";

  log << "[DETAIL] Pre-FEC error positions: "
      << format_positions(result.pre_fec_error_positions) << "\n";
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
}

}  // namespace detail
}  // namespace ofec_single
