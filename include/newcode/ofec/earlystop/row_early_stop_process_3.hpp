#pragma once

#include "newcode/params.hpp"

namespace newcode {

/**
 * early-stop 动作 3：
 * 1. 对当前 256 位码字做硬判并送入 BCH(255,239) 硬解码；
 * 2. 解码成功后重建 overall parity；
 * 3. 直接输出仅由硬解结果决定的 ±|hard_llr_mag| 外信息。
 *
 * 返回值：
 * - true：硬解成功，y2_256 已被完整填充
 * - false：硬解失败，调用方应视为本轮不产出输出
 */
template <typename LLR>
bool row_early_stop_process_3(const LLR* lin256,
                              const LLR* lch256,
                              float* y2_256,
                              const newcode::Params& p);

}  // namespace newcode

#include "ofec/earlystop/row_early_stop_process_3.ipp"
