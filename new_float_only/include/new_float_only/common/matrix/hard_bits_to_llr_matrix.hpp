#pragma once

#include <cstdint>

#include "new_float_only/common/matrix/matrix.hpp"

namespace matrix {

Matrix<float> hard_bits_to_llr_matrix(const Matrix<uint8_t>& bits_mat, float A = 50.0f);

} // namespace matrix