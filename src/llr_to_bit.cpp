#include "newcode/llr_to_bit.hpp"
#include "newcode/matrix.hpp"   // 如果你的 Matrix 头不在这个路径，请改成实际路径
#include <cstddef>
#include <cstdint>

namespace newcode {

static inline uint8_t hard_decide(float L) noexcept { return L >= 0.f ? 0u : 1u; }

Matrix<uint8_t> llr_to_bit(const Matrix<float>& llr_mat) {
    Matrix<uint8_t> bit_mat(llr_mat.rows(), llr_mat.cols());
    for (size_t r = 0; r < llr_mat.rows(); ++r) {
        for (size_t c = 0; c < llr_mat.cols(); ++c) {
            bit_mat[r][c] = hard_decide(llr_mat[r][c]);
        }
    }
    return bit_mat;
}

} // namespace newcode
