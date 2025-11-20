#pragma once

#include <vector>
#include <cstdint>

#include "newcode/matrix.hpp"
#include "newcode/params.hpp"
#include "newcode/qfloat.hpp"

namespace newcode {

template <typename LLR>
struct DecoderCoreResult {
  Matrix<float> lout;
  std::vector<bool> produced_rows;
};

template <typename LLR>
DecoderCoreResult<LLR> Decoder_Core_plain(const Matrix<LLR>& lin_matrix,
                                          const Matrix<LLR>& lch_matrix,
                                          bool use_hard_decode,
                                          const Params& p);

template <typename LLR>
DecoderCoreResult<LLR> Decoder_Core_ebchPF(const Matrix<LLR>& lin_matrix,
                                           const Matrix<LLR>& lch_matrix,
                                           bool use_hard_decode,
                                           const Params& p);

extern template DecoderCoreResult<float> Decoder_Core_plain<float>(const Matrix<float>&,
                                                                   const Matrix<float>&,
                                                                   bool,
                                                                   const Params&);
extern template DecoderCoreResult<int8_t> Decoder_Core_plain<int8_t>(const Matrix<int8_t>&,
                                                                     const Matrix<int8_t>&,
                                                                     bool,
                                                                     const Params&);
extern template DecoderCoreResult<float> Decoder_Core_ebchPF<float>(const Matrix<float>&,
                                                                    const Matrix<float>&,
                                                                    bool,
                                                                    const Params&);
extern template DecoderCoreResult<int8_t> Decoder_Core_ebchPF<int8_t>(const Matrix<int8_t>&,
                                                                      const Matrix<int8_t>&,
                                                                      bool,
                                                                      const Params&);

#define DECLARE_DECODER_CORE_QFLOAT(N) \
extern template DecoderCoreResult<qfloat<N>> Decoder_Core_plain<qfloat<N>>( \
    const Matrix<qfloat<N>>&, const Matrix<qfloat<N>>&, bool, const Params&); \
extern template DecoderCoreResult<qfloat<N>> Decoder_Core_ebchPF<qfloat<N>>( \
    const Matrix<qfloat<N>>&, const Matrix<qfloat<N>>&, bool, const Params&);

DECLARE_DECODER_CORE_QFLOAT(2)
DECLARE_DECODER_CORE_QFLOAT(3)
DECLARE_DECODER_CORE_QFLOAT(4)
DECLARE_DECODER_CORE_QFLOAT(5)
DECLARE_DECODER_CORE_QFLOAT(6)
DECLARE_DECODER_CORE_QFLOAT(7)
DECLARE_DECODER_CORE_QFLOAT(8)
DECLARE_DECODER_CORE_QFLOAT(9)
DECLARE_DECODER_CORE_QFLOAT(10)
DECLARE_DECODER_CORE_QFLOAT(11)
DECLARE_DECODER_CORE_QFLOAT(12)
DECLARE_DECODER_CORE_QFLOAT(13)
DECLARE_DECODER_CORE_QFLOAT(14)
DECLARE_DECODER_CORE_QFLOAT(15)

#undef DECLARE_DECODER_CORE_QFLOAT

} // namespace newcode
