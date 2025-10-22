#include "newcode/llr_known_prefix.hpp"

#include <algorithm>
#include <limits>

namespace newcode {

void apply_known_zero_prefix(Matrix<float>& llr_mat, const Params& p)
{
    if (llr_mat.rows() == 0 || llr_mat.cols() == 0) {
        return;
    }

    size_t known_rows = p.win_height_rows();
    known_rows = std::min(known_rows, llr_mat.rows());

    if (known_rows == 0) {
        return;
    }

    const float bit0_llr = 100.0f; // 一个很大的正值，表示强烈相信是比特 0
    for (size_t r = 0; r < known_rows; ++r) {
        for (size_t c = 0; c < llr_mat.cols(); ++c) {
            llr_mat[r][c] = bit0_llr;
        }
    }
}

} // namespace newcode

