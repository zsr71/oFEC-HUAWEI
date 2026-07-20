#pragma once

#include "newcode/params.hpp"

namespace newcode {

/**
 * 根据顶层配置选择早停命中后的行级输出动作。
 * 当前支持：
 * 1. residual + sign * beta
 * 2. residual only / scaled
 * 3. BCH 硬解成功后直接输出符号化 LLR；失败则本轮不产出输出
 * 4. residual + sign * beta，但在动作内预除 alpha 抵消公共缩放
 * 5. residual 预除 alpha 后再加 sign * beta
 * 6. 仅按 sign(lin) 输出 ±early-stop beta，不做 BCH 硬解
 * 7. 仅按 sign(lin) 输出 ±early-stop beta，并预除 alpha
 * 8. 将 ±early-stop beta 视作目标 posterior，减 channel 后预除 alpha
 */
template <typename LLR>
bool apply_row_early_stop_action(const LLR* lin256,
                                 const LLR* lch256,
                                 float* y2_256,
                                 const newcode::Params& p);

}  // namespace newcode

#include "ofec/earlystop/row_early_stop_action.ipp"
