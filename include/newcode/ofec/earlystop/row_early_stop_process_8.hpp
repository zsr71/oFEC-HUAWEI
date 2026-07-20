#pragma once

#include "newcode/params.hpp"

namespace newcode {

/**
 * early-stop 动作 8：
 * - 把 ±EARLY_STOP_ACTION_SIGN_BETA 视作目标 posterior；
 * - 先减去 channel LLR 得到 extrinsic；
 * - 再预除 alpha，抵消公共后处理的 alpha 缩放。
 */
template <typename LLR>
void row_early_stop_process_8(const LLR* lin256,
                              const LLR* lch256,
                              float* y2_256,
                              const newcode::Params& p);

}  // namespace newcode

#include "ofec/earlystop/row_early_stop_process_8.ipp"
