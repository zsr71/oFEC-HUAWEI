#pragma once

#include <limits>
#include <string>
#include <vector>

#include "newcode/ber.hpp"
#include "newcode/params.hpp"

namespace newcode {

struct PipelineResult {
  float ebn0_db = std::numeric_limits<float>::quiet_NaN();
  BerStats pre_fec;
  BerStats post_fec;
  std::vector<double> tile_early_stop_pct;
};

inline constexpr float DEFAULT_EBN0_DB = 3.24f;

PipelineResult run_pipeline(const Params& params,
                            const std::string& label,
                            float ebn0_dB = DEFAULT_EBN0_DB);

} // namespace newcode
