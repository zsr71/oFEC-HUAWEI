#include "new_float_only/single_runner.hpp"

#include <algorithm>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <optional>
#include <sstream>
#include <stdexcept>

#include "new_float_only/io/dualwriter.hpp"
#include "new_float_only/io/prepare_log_file.hpp"
#include "new_float_only/pipeline_runner.hpp"

namespace new_float_only {
namespace {

std::string format_positions(const std::vector<std::size_t>& positions,
                             std::size_t max_count = 10) {
  std::ostringstream oss;
  const std::size_t count = std::min(max_count, positions.size());
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
    oss << "... (+" << (positions.size() - count) << " more)";
  }
  oss << "]";
  return oss.str();
}

/**
 * 规范化解码配置。
 * 关键语句：按 CHASE_L 自动回填 CHASE_NTEST，并检查 alpha/beta 列表长度。
 */
Params normalize_decoder_config(const SingleRunConfig& config) {
  Params params = config.decoder;
  params.CHASE_NTEST = std::max(params.CHASE_NTEST, 1 << params.CHASE_L);

  if (!params.valid()) {
    throw std::invalid_argument("DecoderConfig is invalid");
  }
  if (!params.ALPHA_LIST.empty() && params.ALPHA_LIST.size() != params.TILES_PER_WIN) {
    throw std::invalid_argument("ALPHA_LIST size must equal TILES_PER_WIN");
  }
  if (!params.beta_list.empty() && params.beta_list.size() != params.TILES_PER_WIN) {
    throw std::invalid_argument("beta_list size must equal TILES_PER_WIN");
  }
  if (!params.ALPHA_LIST.empty()) {
    params.ALPHA = params.ALPHA_LIST.front();
  } else {
    params.ALPHA_LIST.assign(params.TILES_PER_WIN, params.ALPHA);
  }
  if (!params.beta_list.empty()) {
    params.beta = params.beta_list.front();
  } else {
    params.beta_list.assign(params.TILES_PER_WIN, params.beta);
  }
  return params;
}

PipelineConfig build_pipeline_config(const SingleRunConfig& config) {
  PipelineConfig pipeline;
  pipeline.normalize_extrinsic = config.normalize_extrinsic;
  pipeline.bits_per_symbol = config.bits_per_symbol;
  pipeline.bitgen_seed = config.bitgen_seed;
  pipeline.channel_seed = config.channel_seed;
  pipeline.generate_random_bits = config.generate_random_bits;
  pipeline.quiet = false;
  return pipeline;
}

void log_run_overview(const SingleRunConfig& config,
                      const Params& params,
                      io::DualWriter& log) {
  log << "[INFO] run_single(label=" << config.label
      << ", Eb/N0=" << config.ebn0_db
      << " dB, CHASE_L=" << params.CHASE_L << ")\n";
  log << "[INFO] RNG seeds (bitgen/channel) = "
      << config.bitgen_seed << "/"
      << config.channel_seed << "\n";
  log << "[INFO] bits_per_symbol=" << config.bits_per_symbol
      << ", normalize_extrinsic=" << (config.normalize_extrinsic ? "ON" : "OFF")
      << "\n";
  log << "[INFO] pipeline=direct no-permutation path\n";
}

void log_pipeline_result(const PipelineResult& result, io::DualWriter& log) {
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
  if (!result.work_llr_path.empty()) {
    log << "[INFO] work_llr saved to " << result.work_llr_path << "\n";
  }
}

}  // namespace

SingleRunResult run_single(const SingleRunConfig& config) {
  const Params params = normalize_decoder_config(config);

  const std::filesystem::path data_dir = "data";
  std::string log_path;
  std::ofstream log_file = io::prepare_log_file(data_dir, log_path);
  io::DualWriter log(log_file);

  log_run_overview(config, params, log);
  const PipelineResult pipeline_result =
      run_pipeline(params, build_pipeline_config(config), config.label, config.ebn0_db);
  log_pipeline_result(pipeline_result, log);
  log << "[INFO] log saved at " << log_path << "\n";

  SingleRunResult result;
  result.ebn0_db = pipeline_result.ebn0_db;
  result.pre_fec = pipeline_result.pre_fec;
  result.post_fec = pipeline_result.post_fec;
  result.pre_fec_error_positions = pipeline_result.pre_fec_error_positions;
  result.post_fec_error_positions = pipeline_result.post_fec_error_positions;
  result.log_path = log_path;
  result.work_llr_path = pipeline_result.work_llr_path;
  return result;
}

}  // namespace new_float_only
