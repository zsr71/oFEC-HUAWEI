#pragma once

#include <string>
#include <vector>

#include "newcode/ber.hpp"
#include "newcode/params.hpp"

namespace newcode {

struct PipelineResult {
  BerStats pre_fec;
  BerStats post_fec;
  std::vector<double> tile_early_stop_pct;
};

PipelineResult run_pipeline(const Params& params, const std::string& label);

} // namespace newcode
