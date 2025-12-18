#pragma once

#include <cstdint>

#include "newcode/matrix.hpp"

namespace newcode {

Matrix<float> hard_bits_to_llr_matrix(const Matrix<uint8_t>& bits_mat, float A = 50.0f);

} // namespace newcode
