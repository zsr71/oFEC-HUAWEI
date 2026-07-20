#pragma once

#include "newcode/params.hpp"

namespace newcode {

/**
 * early-stop 动作 7：
 * - 不做 BCH hard-decode；
 * - 直接根据当前 lin256 的符号输出 ±EARLY_STOP_ACTION_SIGN_BETA；
 * - 在动作内预除 alpha，抵消公共后处理的 alpha 缩放。
 */
template <typename LLR>
void row_early_stop_process_7(const LLR* lin256,
                              const LLR* lch256,
                              float* y2_256,
                              const newcode::Params& p);

}  // namespace newcode

#include "ofec/earlystop/row_early_stop_process_7.ipp"
