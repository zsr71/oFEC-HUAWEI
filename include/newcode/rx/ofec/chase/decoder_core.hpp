#pragma once

#include <vector>
#include <cstdint>

#include "newcode/common/matrix/matrix.hpp"
#include "newcode/params.hpp"
#include "newcode/common/qfloat/qfloat.hpp"

namespace chase {

template <typename LLR>
struct DecoderCoreResult {
  matrix::Matrix<float> lout;
  std::vector<bool> produced_rows;
};

template <typename LLR>
DecoderCoreResult<LLR> Decoder_Core_plain(const matrix::Matrix<LLR>& lin_matrix,
                                          const matrix::Matrix<LLR>& lch_matrix,
                                          bool use_hard_decode,
                                          const newcode::Params& p,
                                          const std::vector<bool>* early_stop_row_flags,
                                          const std::vector<uint8_t>* mux_state);

template <typename LLR>
DecoderCoreResult<LLR> Decoder_Core_ebchPF(const matrix::Matrix<LLR>& lin_matrix,
                                           const matrix::Matrix<LLR>& lch_matrix,
                                           bool use_hard_decode,
                                           const newcode::Params& p,
                                           const std::vector<bool>* early_stop_row_flags,
                                           const std::vector<uint8_t>* mux_state);

extern template DecoderCoreResult<float> Decoder_Core_plain<float>(const matrix::Matrix<float>&,
                                                                   const matrix::Matrix<float>&,
                                                                   bool,
                                                                   const newcode::Params&,
                                                                   const std::vector<bool>*,
                                                                   const std::vector<uint8_t>*);
extern template DecoderCoreResult<int8_t> Decoder_Core_plain<int8_t>(const matrix::Matrix<int8_t>&,
                                                                     const matrix::Matrix<int8_t>&,
                                                                     bool,
                                                                     const newcode::Params&,
                                                                     const std::vector<bool>*,
                                                                     const std::vector<uint8_t>*);
extern template DecoderCoreResult<float> Decoder_Core_ebchPF<float>(const matrix::Matrix<float>&,
                                                                    const matrix::Matrix<float>&,
                                                                    bool,
                                                                    const newcode::Params&,
                                                                    const std::vector<bool>*,
                                                                    const std::vector<uint8_t>*);
extern template DecoderCoreResult<int8_t> Decoder_Core_ebchPF<int8_t>(const matrix::Matrix<int8_t>&,
                                                                      const matrix::Matrix<int8_t>&,
                                                                      bool,
                                                                      const newcode::Params&,
                                                                      const std::vector<bool>*,
                                                                      const std::vector<uint8_t>*);

#define DECLARE_DECODER_CORE_QFLOAT(N) \
extern template DecoderCoreResult<qfloat::qfloat<N>> Decoder_Core_plain<qfloat::qfloat<N>>( \
    const matrix::Matrix<qfloat::qfloat<N>>& lin_matrix, const matrix::Matrix<qfloat::qfloat<N>>& lch_matrix, bool, const newcode::Params&, \
    const std::vector<bool>*, const std::vector<uint8_t>*); \
extern template DecoderCoreResult<qfloat::qfloat<N>> Decoder_Core_ebchPF<qfloat::qfloat<N>>( \
    const matrix::Matrix<qfloat::qfloat<N>>& lin_matrix, const matrix::Matrix<qfloat::qfloat<N>>& lch_matrix, bool, const newcode::Params&, \
    const std::vector<bool>*, const std::vector<uint8_t>*);

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
