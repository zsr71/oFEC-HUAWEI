#pragma once

#include "newcode/params.hpp"

namespace newcode {

/**
 * early-stop 动作 6：
 * - 不做 BCH hard-decode；
 * - 直接根据当前 lin256 的符号输出 ±EARLY_STOP_ACTION_SIGN_BETA；
 * - beta 使用 early-stop 专用的 sign beta，而不是普通 Chase beta。
 */
template <typename LLR>
void row_early_stop_process_6(const LLR* lin256,
                              const LLR* lch256,
                              float* y2_256,
                              const newcode::Params& p);

}  // namespace newcode

#include "ofec/earlystop/row_early_stop_process_6.ipp"
