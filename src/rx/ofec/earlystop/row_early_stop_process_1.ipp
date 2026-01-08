#pragma once

#include <cstddef>

#include "newcode/common/qfloat/llr_utils.hpp"
#include "newcode/params.hpp"

namespace newcode {

template <typename LLR>
void row_early_stop_process_1(const LLR* lin256,
                              float* y2_256,
                              const newcode::Params& p)
{
  const float mag = p.beta;
  for (size_t j = 0; j < newcode::Params::BCH_N; ++j) {
    const float v = qfloat::llr_to_float(lin256[j]);
    y2_256[j] = (v < 0.0f) ? -mag : mag;
  }
}

} // namespace newcode
