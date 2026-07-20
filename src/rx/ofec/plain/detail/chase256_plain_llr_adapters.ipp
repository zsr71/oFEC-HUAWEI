#pragma once

#include <cmath>
#include <limits>
#include <type_traits>

namespace newcode {
namespace detail {

// ----- LLR adapters -----
template<typename LLR>
inline float llr_to_float(LLR x)
{
    // 输入:
    // - x: 任意 LLR 表示，可能是 float/int/qfloat。
    // 输出:
    // - 对应的 float 值。
    // 用途:
    // - Plain Chase 内部统一在 float 域计算，先把各种 LLR 类型转成 float。
    return static_cast<float>(x);
}

template<typename LLR, typename std::enable_if<std::is_floating_point<LLR>::value, int>::type = 0>
inline LLR llr_from_float(float x)
{
    // 输入:
    // - x: float 域结果。
    // 输出:
    // - 同数值语义的浮点 LLR。
    // 用途:
    // - 浮点路径无需额外量化，直接做类型转换返回。
    return static_cast<LLR>(x);
}

template<typename LLR, typename std::enable_if<std::is_integral<LLR>::value && std::is_signed<LLR>::value, int>::type = 0>
inline LLR llr_from_float(float x)
{
    // 输入:
    // - x: float 域结果。
    // 输出:
    // - 饱和并四舍五入后的整数 LLR。
    // 用途:
    // - 给整数量化路径做安全写回，避免越界。
    const float lo = (float)std::numeric_limits<LLR>::min();
    const float hi = (float)std::numeric_limits<LLR>::max();
    if (x < lo) x = lo;
    if (x > hi) x = hi;
    return (LLR)std::lrintf(x);
}

template<typename LLR, typename std::enable_if<!std::is_arithmetic<LLR>::value, int>::type = 0>
inline LLR llr_from_float(float x)
{
    // 输入:
    // - x: float 域结果。
    // 输出:
    // - 自定义 LLR 类型对象。
    // 用途:
    // - 对 qfloat 之类的自定义量化类型，复用其 from_float() 完成量化。
    return LLR::from_float(x);
}

} // namespace detail
} // namespace newcode
