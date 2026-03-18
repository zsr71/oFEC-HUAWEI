#pragma once
#include <vector>
#include <cstdint>
#include <stdexcept>

#include "new_float_only/common/matrix/matrix.hpp"
#include "new_float_only/params.hpp"

namespace matrix {


std::vector<uint8_t> rx_info_from_bit_llr(const matrix::Matrix<float>& bit_llr_mat, const new_float_only::Params& p);

} // namespace matrix
