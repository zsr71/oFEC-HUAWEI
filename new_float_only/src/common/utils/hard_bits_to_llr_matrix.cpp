#include "new_float_only/common/matrix/hard_bits_to_llr_matrix.hpp"

namespace matrix {

Matrix<float> hard_bits_to_llr_matrix(const Matrix<uint8_t>& bits_mat, float A)
{
  Matrix<float> m(bits_mat.rows(), bits_mat.cols());
  for (size_t r = 0; r < bits_mat.rows(); ++r)
    for (size_t c = 0; c < bits_mat.cols(); ++c)
      m[r][c] = (bits_mat[r][c] ? -A : +A);
  return m;
}

} // namespace matrix
