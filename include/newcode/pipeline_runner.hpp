#pragma once

#include <string>

#include "newcode/ber.hpp"
#include "newcode/params.hpp"

namespace newcode {

struct PipelineResult {
  BerStats pre_fec;
  BerStats post_fec;
};

PipelineResult run_pipeline(const Params& params, const std::string& label);

} // namespace newcode

