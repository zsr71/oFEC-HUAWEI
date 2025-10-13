#include "newcode/llr_known_prefix.hpp"

#include <algorithm>
#include <limits>

namespace newcode {

void apply_known_zero_prefix(Matrix<float>& llr_mat, const Params& p)
{
    if (llr_mat.rows() == 0 || llr_mat.cols() == 0) {
        return;
    }

    const size_t bits_per_row = Params::BITS_PER_SUBBLOCK_DIM;
    const size_t known_subrows = (p.NUM_GUARD_SUBROWS * 2u) + Params::INFO_SUBROWS_PER_CODE;
    size_t known_rows = known_subrows * bits_per_row;
    known_rows = std::min(known_rows, llr_mat.rows());

    if (known_rows == 0) {
        return;
    }

    const float bit0_llr = std::numeric_limits<float>::infinity();
    for (size_t r = 0; r < known_rows; ++r) {
        for (size_t c = 0; c < llr_mat.cols(); ++c) {
            llr_mat[r][c] = bit0_llr;
        }
    }
}

} // namespace newcode

