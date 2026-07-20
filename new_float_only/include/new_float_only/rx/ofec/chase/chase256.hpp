#pragma once

#include "new_float_only/params.hpp"

namespace chase {

/**
 * Plain Chase(256) 软输入软输出译码。
 * Lin256 为“信道 + 历史”的当前输入，Lch256 为原始信道 LLR，
 * Y2_256 输出的是 extrinsic，而不是最终的 posterior。
 */
void chase_decode_256_plain(const float* lin256,
                            const float* lch256,
                            float* y2_256,
                            const new_float_only::Params& config);

/**
 * 两参版本为兼容包装：把同一份数组同时当作 lin 与 lch 输入。
 */
void chase_decode_256_plain(const float* y256,
                            float* y2_256,
                            const new_float_only::Params& config);

}  // namespace chase
