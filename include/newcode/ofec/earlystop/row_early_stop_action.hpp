#pragma once

#include "newcode/params.hpp"

namespace newcode {

/**
 * 根据顶层配置选择早停命中后的行级输出动作。
 * 当前支持：
 * 1. residual + sign * beta
 * 2. residual only / scaled
 */
template <typename LLR>
void apply_row_early_stop_action(const LLR* lin256,
                                 const LLR* lch256,
                                 float* y2_256,
                                 const newcode::Params& p);

}  // namespace newcode

#include "ofec/earlystop/row_early_stop_action.ipp"
