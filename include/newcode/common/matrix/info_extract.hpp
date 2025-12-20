#pragma once
#include <vector>
#include <cstdint>
#include <stdexcept>

#include "newcode/common/matrix/matrix.hpp"
#include "newcode/params.hpp"

namespace matrix {


std::vector<uint8_t> rx_info_from_bit_llr(const matrix::Matrix<float>& bit_llr_mat, const newcode::Params& p);

} // namespace matrix
