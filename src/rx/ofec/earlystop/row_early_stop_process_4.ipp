#pragma once

#include <cmath>
#include <cstddef>

#include "newcode/common/qfloat/llr_utils.hpp"
#include "newcode/params.hpp"

namespace newcode {

template <typename LLR>
void row_early_stop_process_4(const LLR* lin256,
                              const LLR* lch256,
                              float* y2_256,
                              const newcode::Params& p)
{
  const float mag = p.EARLY_STOP_ACTION_SIGN_BETA;
  const float alpha = p.ALPHA;
  const float inv_alpha =
      (std::fabs(alpha) > 1e-8f) ? (1.0f / alpha) : 1.0f;
  for (size_t j = 0; j < newcode::Params::BCH_N; ++j) {
    const float v = qfloat::llr_to_float(lin256[j]);
    const float raw =
        (qfloat::llr_to_float(lin256[j]) - qfloat::llr_to_float(lch256[j])) +
        ((v < 0.0f) ? -mag : mag);
    // 预先除以当前 tile 的 alpha，使公共后处理乘回 alpha 后，
    // 最终写回值等效于“动作1输出不经过 alpha 缩放”。
    y2_256[j] = raw * inv_alpha;
  }
}

} // namespace newcode
