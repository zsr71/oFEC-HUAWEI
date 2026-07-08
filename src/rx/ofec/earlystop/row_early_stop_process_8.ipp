#pragma once

#include <cmath>
#include <cstddef>

#include "newcode/common/qfloat/llr_utils.hpp"
#include "newcode/params.hpp"

namespace newcode {

template <typename LLR>
void row_early_stop_process_8(const LLR* lin256,
                              const LLR* lch256,
                              float* y2_256,
                              const newcode::Params& p)
{
  const float mag = std::fabs(p.EARLY_STOP_ACTION_SIGN_BETA);
  const float alpha = p.ALPHA;
  const float inv_alpha =
      (std::fabs(alpha) > 1e-8f) ? (1.0f / alpha) : 1.0f;
  for (size_t j = 0; j < newcode::Params::BCH_N; ++j) {
    const float v = qfloat::llr_to_float(lin256[j]);
    const float target = (v < 0.0f) ? -mag : mag;
    const float channel = qfloat::llr_to_float(lch256[j]);
    y2_256[j] = (target - channel) * inv_alpha;
  }
}

} // namespace newcode
