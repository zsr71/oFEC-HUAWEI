#include "newcode/llr_known_prefix.hpp"

#include <algorithm>
#include <limits>
#include <cmath>   // std::fabs

namespace newcode {

void apply_known_zero_prefix(matrix::Matrix<float>& llr_mat, const Params& p)
{
    const size_t R = llr_mat.rows();
    const size_t C = llr_mat.cols();
    if (R == 0 || C == 0) return;

    // 已知前缀行数：取“前缀 tile 数 × 单 tile 高度”与总行数的较小者
    size_t known_rows = std::min(p.known_prefix_rows(), R);

    // 1) 已知前缀行强制为比特0（大正 LLR）
    if (known_rows > 0)
    {
        constexpr float bit0_llr = 2.0f;
        for (size_t r = 0; r < known_rows; ++r)
            for (size_t c = 0; c < C; ++c)
                llr_mat[r][c] = bit0_llr;
    }

    // 2) 对 “known_rows 之外的区域” 做整体均值归一化：
    //    llrMean = mean(abs(llrMat(:)))；仅统计/缩放 [known_rows..R-1, 0..C-1]
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
                auto tail = llr_mat.slice(known_rows, R, 0, C);
                scale(tail, invMean);
            }
            // 若 llrMean==0 则不缩放，保持原值（避免除零）
        }
    }
}

} // namespace newcode
