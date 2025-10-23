#include "newcode/decoder_core.hpp"

#include "newcode/chase256.hpp"
#include "newcode/ofec_decoder_hard.hpp"

#include <array>
#include <stdexcept>

namespace newcode {

template<typename LLR>
void chase_decode_256(const LLR* Lin256, const LLR* Lch256, LLR* Y2_256, const Params& p);

template <typename LLR>
DecoderCoreResult<LLR> Decoder_Core(const Matrix<LLR>& lin_matrix,
                                    const Matrix<LLR>& lch_matrix,
                                    bool use_hard_decode,
                                    const Params& p)
{
  const size_t rows = lin_matrix.rows();
  const size_t cols = lin_matrix.cols();
  const size_t expected_cols =
      2 * Params::NUM_SUBBLOCK_COLS * Params::BITS_PER_SUBBLOCK_DIM;
  if (cols != expected_cols) {
    throw std::invalid_argument("Decoder_Core: unexpected column count");
  }

  DecoderCoreResult<LLR> result{
      Matrix<LLR>(rows, cols),
      std::vector<bool>(rows, false)};

  for (size_t row = 0; row < rows; ++row) {
    std::array<LLR, Params::BCH_N> LinVec{};
    std::array<LLR, Params::BCH_N> LchVec{};

    for (size_t col = 0; col < cols; ++col) {
      LinVec[col] = lin_matrix[row][col];
      LchVec[col] = lch_matrix[row][col];
    }

    std::array<LLR, Params::BCH_N> Y2{};
    bool produced = false;

    if (use_hard_decode) {
      produced = perform_hard_decode<LLR>(LinVec, LchVec, Y2, p);
    } else {
      chase_decode_256<LLR>(LinVec.data(), LchVec.data(), Y2.data(), p);
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

template DecoderCoreResult<float> Decoder_Core<float>(const Matrix<float>&,
                                                      const Matrix<float>&,
                                                      bool,
                                                      const Params&);
template DecoderCoreResult<int8_t> Decoder_Core<int8_t>(const Matrix<int8_t>&,
                                                        const Matrix<int8_t>&,
                                                        bool,
                                                        const Params&);
template DecoderCoreResult<qfloat<4>> Decoder_Core<qfloat<4>>(const Matrix<qfloat<4>>&,
                                                              const Matrix<qfloat<4>>&,
                                                              bool,
                                                              const Params&);
template DecoderCoreResult<qfloat<5>> Decoder_Core<qfloat<5>>(const Matrix<qfloat<5>>&,
                                                              const Matrix<qfloat<5>>&,
                                                              bool,
                                                              const Params&);

} // namespace newcode

