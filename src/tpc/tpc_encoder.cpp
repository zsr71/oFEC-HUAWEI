#include "tpc_internal.hpp"

#include <array>
#include <stdexcept>

#include "newcode/common/bch/bch_255_239.hpp"

namespace newcode {
namespace tpc {

matrix::Matrix<uint8_t> tpc_encode(const std::vector<uint8_t>& info_bits,
                                   const Params& /*params*/) {
  const std::size_t info_bits_needed = kTpcInfoDim * kTpcInfoDim;
  if (info_bits.size() != info_bits_needed) {
    throw std::invalid_argument("tpc_encode: info_bits size mismatch.");
  }

  // 239x239 信息块按行填入矩阵
  matrix::Matrix<uint8_t> info_mat(kTpcInfoDim, kTpcInfoDim);
  std::size_t idx = 0;
  for (std::size_t r = 0; r < kTpcInfoDim; ++r) {
    for (std::size_t c = 0; c < kTpcInfoDim; ++c) {
      info_mat[r][c] = static_cast<uint8_t>(info_bits[idx++] & 1u);
    }
  }

  // 行编码：每行 239 -> 256
  matrix::Matrix<uint8_t> row_encoded(kTpcInfoDim, kTpcCodeDim);
  for (std::size_t r = 0; r < kTpcInfoDim; ++r) {
    std::vector<uint8_t> row_bits;
    row_bits.reserve(kTpcInfoDim);
    for (std::size_t c = 0; c < kTpcInfoDim; ++c) {
      row_bits.push_back(info_mat[r][c]);
    }
    const auto cw = bch::bch_255_239_encode(row_bits);
    for (std::size_t c = 0; c < kTpcCodeDim; ++c) {
      row_encoded[r][c] = cw[c];
    }
  }

  // 列编码：对 239 行结果逐列编码成 256x256
  matrix::Matrix<uint8_t> code_mat(kTpcCodeDim, kTpcCodeDim);
  for (std::size_t c = 0; c < kTpcCodeDim; ++c) {
    std::vector<uint8_t> col_bits;
    col_bits.reserve(kTpcInfoDim);
    for (std::size_t r = 0; r < kTpcInfoDim; ++r) {
      col_bits.push_back(row_encoded[r][c]);
    }
    const auto cw = bch::bch_255_239_encode(col_bits);
    for (std::size_t r = 0; r < kTpcCodeDim; ++r) {
      code_mat[r][c] = cw[r];
    }
  }

  return code_mat;
}

}  // namespace tpc
}  // namespace newcode
