#pragma once

#include <vector>
#include <cstdint>

#include "newcode/matrix.hpp"
#include "newcode/params.hpp"
#include "newcode/qfloat.hpp"

namespace newcode {

template <typename LLR>
struct DecoderCoreResult {
  Matrix<LLR> lout;
  std::vector<bool> produced_rows;
};

template <typename LLR>
DecoderCoreResult<LLR> Decoder_Core(const Matrix<LLR>& lin_matrix,
                                    const Matrix<LLR>& lch_matrix,
                                    bool use_hard_decode,
                                    const Params& p);

extern template DecoderCoreResult<float> Decoder_Core<float>(const Matrix<float>&,
                                                             const Matrix<float>&,
                                                             bool,
                                                             const Params&);
extern template DecoderCoreResult<int8_t> Decoder_Core<int8_t>(const Matrix<int8_t>&,
                                                               const Matrix<int8_t>&,
                                                               bool,
                                                               const Params&);
extern template DecoderCoreResult<qfloat<4>> Decoder_Core<qfloat<4>>(const Matrix<qfloat<4>>&,
                                                                     const Matrix<qfloat<4>>&,
                                                                     bool,
                                                                     const Params&);
extern template DecoderCoreResult<qfloat<5>> Decoder_Core<qfloat<5>>(const Matrix<qfloat<5>>&,
                                                                     const Matrix<qfloat<5>>&,
                                                                     bool,
                                                                     const Params&);

} // namespace newcode
