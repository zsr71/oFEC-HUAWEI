#include "newcode/decoder_core.hpp"

#include "newcode/chase256.hpp"
#include "newcode/ofec_decoder_hard.hpp"

#include <array>
#include <stdexcept>

namespace newcode {
namespace {

template<typename LLR>
using ChaseFn = void (*)(const LLR*, const LLR*, float*, const Params&);

template<typename LLR>
DecoderCoreResult<LLR> Decoder_Core_impl(const Matrix<LLR>& lin_matrix,
                                         const Matrix<LLR>& lch_matrix,
                                         bool use_hard_decode,
                                         const Params& p,
                                         ChaseFn<LLR> chase_fn)
{
  const size_t rows = lin_matrix.rows();
  const size_t cols = lin_matrix.cols();
  const size_t expected_cols =
      2 * Params::NUM_SUBBLOCK_COLS * Params::BITS_PER_SUBBLOCK_DIM;
  if (cols != expected_cols) {
    throw std::invalid_argument("Decoder_Core: unexpected column count");
  }

  DecoderCoreResult<LLR> result{
      Matrix<float>(rows, cols),
      std::vector<bool>(rows, false)};

  for (size_t row = 0; row < rows; ++row) {
    std::array<LLR, Params::BCH_N> LinVec{};
    std::array<LLR, Params::BCH_N> LchVec{};

    for (size_t col = 0; col < cols; ++col) {
      LinVec[col] = lin_matrix[row][col];
      LchVec[col] = lch_matrix[row][col];
    }

    std::array<float, Params::BCH_N> Y2{};
    bool produced = false;

    if (use_hard_decode) {
      produced = perform_hard_decode<LLR>(LinVec, LchVec, Y2, p);
    } else {
      chase_fn(LinVec.data(), LchVec.data(), Y2.data(), p);
      produced = true;
    }

    if (produced) {
      result.produced_rows[row] = true;
      for (size_t col = 0; col < cols; ++col) {
        result.lout[row][col] = Y2[col];
      }
    }
  }

  return result;
}

} // namespace

template<typename LLR>
DecoderCoreResult<LLR> Decoder_Core_plain(const Matrix<LLR>& lin_matrix,
                                          const Matrix<LLR>& lch_matrix,
                                          bool use_hard_decode,
                                          const Params& p)
{
  return Decoder_Core_impl(lin_matrix, lch_matrix, use_hard_decode, p,
                           &chase_decode_256_plain<LLR>);
}

template<typename LLR>
DecoderCoreResult<LLR> Decoder_Core_ebchPF(const Matrix<LLR>& lin_matrix,
                                           const Matrix<LLR>& lch_matrix,
                                           bool use_hard_decode,
                                           const Params& p)
{
  return Decoder_Core_impl(lin_matrix, lch_matrix, use_hard_decode, p,
                           &chase_decode_256_ebchPF<LLR>);
}

template DecoderCoreResult<float> Decoder_Core_plain<float>(const Matrix<float>&,
                                                            const Matrix<float>&,
                                                            bool,
                                                            const Params&);
template DecoderCoreResult<int8_t> Decoder_Core_plain<int8_t>(const Matrix<int8_t>&,
                                                              const Matrix<int8_t>&,
                                                              bool,
                                                              const Params&);

template DecoderCoreResult<float> Decoder_Core_ebchPF<float>(const Matrix<float>&,
                                                             const Matrix<float>&,
                                                             bool,
                                                             const Params&);
template DecoderCoreResult<int8_t> Decoder_Core_ebchPF<int8_t>(const Matrix<int8_t>&,
                                                               const Matrix<int8_t>&,
                                                               bool,
                                                               const Params&);

#define INSTANTIATE_DECODER_CORE_QFLOAT(N) \
template DecoderCoreResult<qfloat<N>> Decoder_Core_plain<qfloat<N>>( \
    const Matrix<qfloat<N>>&, const Matrix<qfloat<N>>&, bool, const Params&); \
template DecoderCoreResult<qfloat<N>> Decoder_Core_ebchPF<qfloat<N>>( \
    const Matrix<qfloat<N>>&, const Matrix<qfloat<N>>&, bool, const Params&);

INSTANTIATE_DECODER_CORE_QFLOAT(2)
INSTANTIATE_DECODER_CORE_QFLOAT(3)
INSTANTIATE_DECODER_CORE_QFLOAT(4)
INSTANTIATE_DECODER_CORE_QFLOAT(5)
INSTANTIATE_DECODER_CORE_QFLOAT(6)
INSTANTIATE_DECODER_CORE_QFLOAT(7)
INSTANTIATE_DECODER_CORE_QFLOAT(8)
INSTANTIATE_DECODER_CORE_QFLOAT(9)
INSTANTIATE_DECODER_CORE_QFLOAT(10)
INSTANTIATE_DECODER_CORE_QFLOAT(11)
INSTANTIATE_DECODER_CORE_QFLOAT(12)
INSTANTIATE_DECODER_CORE_QFLOAT(13)
INSTANTIATE_DECODER_CORE_QFLOAT(14)
INSTANTIATE_DECODER_CORE_QFLOAT(15)

#undef INSTANTIATE_DECODER_CORE_QFLOAT

} // namespace newcode
