#pragma once

#include <array>

#include "new_float_only/params.hpp"

namespace new_float_only {

/**
 * 硬判决 oFEC 行译码入口。
 * 输入为 256 维行级 LLR，输出为 256 维 extrinsic。
 * 返回值表示 BCH 硬译码是否成功。
 */
bool perform_hard_decode(const std::array<float, 256>& lin256,
                         const std::array<float, 256>& lch256,
                         std::array<float, 256>& y2,
                         const Params& config);

}  // namespace new_float_only
