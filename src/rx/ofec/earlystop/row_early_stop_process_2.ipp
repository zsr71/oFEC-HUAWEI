#pragma once

#include <cstddef>

#include "newcode/common/qfloat/llr_utils.hpp"
#include "newcode/params.hpp"

namespace newcode {

template <typename LLR>
void row_early_stop_process_2(const LLR* lin256,
                              const LLR* lch256,
                              float* y2_256,
                              const newcode::Params& p)
{
  const float divisor = p.EARLY_STOP_ACTION_RESIDUAL_DIVISOR;
  for (size_t j = 0; j < newcode::Params::BCH_N; ++j) {
    y2_256[j] =
        (qfloat::llr_to_float(lin256[j]) - qfloat::llr_to_float(lch256[j])) / divisor;
  }
}

} // namespace newcode
