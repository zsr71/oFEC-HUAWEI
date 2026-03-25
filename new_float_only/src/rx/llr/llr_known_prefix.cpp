#include "new_float_only/llr_known_prefix.hpp"

#include <algorithm>
#include <limits>
#include <cmath>   // std::fabs

namespace new_float_only {

// 对接收端 LLR 矩阵应用“已知前缀为 0”的约束，并可选地对其余区域做归一化。
//
// 这里的处理分两步：
// 1. 把最前面的 known_rows 行强制写成“可靠的比特 0”
// 2. 如果开启了 NORMALIZE_KNOWN_PREFIX_TAIL，则只对 known_rows 之后的尾部区域
//    做一次整体幅度归一化，避免前缀强制写值后影响真实数据区的幅度尺度
void apply_known_zero_prefix(matrix::Matrix<float>& llr_mat, const Params& p)
{
    const size_t R = llr_mat.rows();
    const size_t C = llr_mat.cols();
    if (R == 0 || C == 0) return;

    // 已知前缀行数由 tile 高度决定；如果矩阵本身更短，就按矩阵总行数截断。
    size_t known_rows = std::min(p.tile_height_rows(), R);

    // 第一步：已知前缀行强制设为比特 0。
    // 这里直接写一个固定的正 LLR，用来告诉后续解码器“这些位置应当被当作 0”。
    if (known_rows > 0)
    {
        constexpr float bit0_llr = 2.0f;
        for (size_t r = 0; r < known_rows; ++r)
            for (size_t c = 0; c < C; ++c)
                llr_mat[r][c] = bit0_llr;
    }

    // 第二步：只对 known_rows 之外的尾部区域做整体幅度归一化。
    // 这里统计尾部所有元素的绝对值均值，再把尾部整体缩放到平均绝对值约为 1。
    // 这样做的目的是让真正参与解码的数据区保持稳定量纲，而不是让前缀强制写入的值
    // 参与均值统计。
    if (known_rows < R && p.NORMALIZE_KNOWN_PREFIX_TAIL)
    {
        double acc_abs = 0.0;
        const size_t n_elems = (R - known_rows) * C;

        for (size_t r = known_rows; r < R; ++r)
            for (size_t c = 0; c < C; ++c)
                acc_abs += std::fabs(static_cast<double>(llr_mat[r][c]));

        if (n_elems > 0)
        {
            const double llrMean = acc_abs / static_cast<double>(n_elems);
            if (llrMean > 0.0)
            {
                const float invMean = static_cast<float>(1.0 / llrMean);
                // 只对尾部切片做缩放，不修改前面已知前缀区域。
                auto tail = llr_mat.slice(known_rows, R, 0, C);
                scale(tail, invMean);
            }
            // 若 llrMean==0 则不缩放，保持原值（避免除零）
        }
    }
}

} // namespace new_float_only
