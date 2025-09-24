#include "newcode/ofec_llr_matrix.hpp"

namespace newcode {

Matrix<float> llr_to_matrix_row_major(const std::vector<float>& llr,
                                      size_t rows,
                                      size_t cols)
{
    const size_t need = rows * cols;
    if (llr.size() != need) {
        throw std::invalid_argument(
            "[llr_to_matrix_row_major] size mismatch: llr.size()="
            + std::to_string(llr.size()) + ", expected=" + std::to_string(need));
    }

    Matrix<float> M(rows, cols);
    size_t idx = 0;
    for (size_t r = 0; r < rows; ++r) {
        for (size_t c = 0; c < cols; ++c) {
            M[r][c] = llr[idx++];
        }
    }
    return M;
}

void fill_llr_matrix_row_major(const std::vector<float>& llr,
                               Matrix<float>& M)
{
    const size_t rows = M.rows();
    const size_t cols = M.cols();
    const size_t need = rows * cols;
    if (llr.size() != need) {
        throw std::invalid_argument(
            "[fill_llr_matrix_row_major] size mismatch: llr.size()="
            + std::to_string(llr.size()) + ", expected=" + std::to_string(need));
    }

    size_t idx = 0;
    for (size_t r = 0; r < rows; ++r) {
        for (size_t c = 0; c < cols; ++c) {
            M[r][c] = llr[idx++];
        }
    }
}

} // namespace newcode
