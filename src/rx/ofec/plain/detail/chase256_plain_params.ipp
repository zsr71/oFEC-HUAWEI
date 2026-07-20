#pragma once

namespace newcode {
namespace detail {

// ----- select coefficients from newcode::Params -----
static inline void pick_cp(const newcode::Params& p, float& beta, float& alpha)
{
    // 输入:
    // - p: 解码参数集合。
    // - beta/alpha: 输出参数引用。
    // 输出:
    // - beta 取自 p.beta，alpha 取自 p.ALPHA。
    // 用途:
    // - 从统一参数对象中提取 Plain Chase/Pyndiah 需要的核心系数。
    //   这里 beta 复用为 fallback 外信息幅度。
    beta = p.beta;
    alpha = p.ALPHA;
}

} // namespace detail
} // namespace newcode
