#pragma once

#include <limits>
#include <string>
#include <vector>

#include "newcode/ber.hpp"
#include "newcode/decoder_api.hpp"
#include "newcode/params.hpp"

namespace newcode {

struct PipelineConfig {
  std::string decoder_name = "ebchPF";
  std::string interleaver_name = "ofec";
  bool normalize_extrinsic = true;
  unsigned bits_per_symbol = 2;
  bool quiet = false;
};

struct PipelineResult {
  float ebn0_db = std::numeric_limits<float>::quiet_NaN();
  BerStats pre_fec;
  BerStats post_fec;
  std::vector<std::size_t> pre_fec_error_positions;
  std::vector<std::size_t> post_fec_error_positions;
  std::vector<double> tile_early_stop_pct;
};

inline constexpr float DEFAULT_EBN0_DB = 3.24f;

PipelineResult run_pipeline(const Params& params,
                            const std::string& label,
                            float ebn0_dB = DEFAULT_EBN0_DB);

PipelineResult run_pipeline(const Params& params,
                            const PipelineConfig& config,
                            const std::string& label,
                            float ebn0_dB = DEFAULT_EBN0_DB);

} // namespace newcode
