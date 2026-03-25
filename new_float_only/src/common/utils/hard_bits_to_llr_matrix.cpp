#include "new_float_only/common/matrix/hard_bits_to_llr_matrix.hpp"

namespace matrix {

/**
 * @brief 把硬判比特矩阵转换成同尺寸的浮点 LLR 矩阵。
 *
 * 映射规则很直接：
 * - 比特 `0` 映射成 `+A`
 * - 比特 `1` 映射成 `-A`
 *
 * 这个函数常用于从发送端码字矩阵构造参考 LLR，也可用于把硬判结果转成
 * 后续模块可复用的 LLR 形式。输出矩阵的尺寸与输入矩阵完全一致。
 *
 * @param bits_mat 输入的硬判比特矩阵。
 * @param A        映射后的 LLR 幅度绝对值。
 * @return Matrix<float> 与输入同尺寸的 LLR 矩阵。
 */
Matrix<float> hard_bits_to_llr_matrix(const Matrix<uint8_t>& bits_mat, float A)
{
  // 先分配一个与输入比特矩阵同尺寸的输出矩阵。
  Matrix<float> m(bits_mat.rows(), bits_mat.cols());

  // 逐元素做 0 -> +A、1 -> -A 的符号映射。
  for (size_t r = 0; r < bits_mat.rows(); ++r)
    for (size_t c = 0; c < bits_mat.cols(); ++c)
      m[r][c] = (bits_mat[r][c] ? -A : +A);
  return m;
}

} // namespace matrix
