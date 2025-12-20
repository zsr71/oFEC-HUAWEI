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
DecoderCoreResult<LLR> Decoder_Core_impl(const matrix::Matrix<LLR>& lin_matrix,
                                         const matrix::Matrix<LLR>& lch_matrix,
                                         bool use_hard_decode,
                                         const Params& p,
                                         ChaseFn<LLR> chase_fn,
                                         const std::vector<bool>* early_stop_row_flags)
{
  const size_t rows = lin_matrix.rows();
  const size_t cols = lin_matrix.cols();
  const size_t expected_cols =
      2 * Params::NUM_SUBBLOCK_COLS * Params::BITS_PER_SUBBLOCK_DIM;
  if (cols != expected_cols) {
    throw std::invalid_argument("Decoder_Core: unexpected column count");
  }

  DecoderCoreResult<LLR> result{
      matrix::Matrix<float>(rows, cols),
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

    Params row_params = p;
    const auto& trace_cfg = p.debug_trace;
    row_params.debug_trace.active_chase_entries.clear();
    row_params.debug_trace.chase_expected_bits = trace_cfg.chase_expected_bits;
    if (trace_cfg.chase_expected_bits &&
        static_cast<size_t>(row) < trace_cfg.chase_expected_bits->size()) {
      row_params.debug_trace.chase_expected_bits_row =
          &(*trace_cfg.chase_expected_bits)[row];
    } else {
      row_params.debug_trace.chase_expected_bits_row = nullptr;
    }
    for (const auto& entry : trace_cfg.active_chase_entries) {
      if (entry.row_index == static_cast<int>(row)) {
        row_params.debug_trace.active_chase_entries.push_back(entry);
      }
    }
    if (!row_params.debug_trace.active_chase_entries.empty()) {
      row_params.debug_trace.chase_decoder_row = static_cast<int>(row);
      row_params.debug_trace.chase_decoder_col =
          row_params.debug_trace.active_chase_entries.front().k;
    } else if (trace_cfg.chase_decoder_row >= 0 &&
               static_cast<int>(row) == trace_cfg.chase_decoder_row) {
      row_params.debug_trace.chase_decoder_row = static_cast<int>(row);
      row_params.debug_trace.chase_decoder_col = trace_cfg.chase_decoder_col;
    } else {
      row_params.debug_trace.chase_decoder_row = -1;
      row_params.debug_trace.chase_decoder_col = -1;
    }

    if (use_hard_decode) {
      produced = perform_hard_decode<LLR>(LinVec, LchVec, Y2, p);
    } else {
      if (early_stop_row_flags &&
          row < early_stop_row_flags->size() &&
          (*early_stop_row_flags)[row]) {
        //produced = false;
        //continue;
      }
      chase_fn(LinVec.data(), LchVec.data(), Y2.data(), row_params);
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
DecoderCoreResult<LLR> Decoder_Core_plain(const matrix::Matrix<LLR>& lin_matrix,
                                          const matrix::Matrix<LLR>& lch_matrix,
                                          bool use_hard_decode,
                                          const Params& p,
                                          const std::vector<bool>* early_stop_row_flags)
{
  return Decoder_Core_impl(lin_matrix, lch_matrix, use_hard_decode, p,
                           &chase_decode_256_plain<LLR>,
                           early_stop_row_flags);
}

template<typename LLR>
DecoderCoreResult<LLR> Decoder_Core_ebchPF(const matrix::Matrix<LLR>& lin_matrix,
                                           const matrix::Matrix<LLR>& lch_matrix,
                                           bool use_hard_decode,
                                           const Params& p,
                                           const std::vector<bool>* early_stop_row_flags)
{
  return Decoder_Core_impl(lin_matrix, lch_matrix, use_hard_decode, p,
                           &chase_decode_256_ebchPF<LLR>,
                           early_stop_row_flags);
}

template DecoderCoreResult<float> Decoder_Core_plain<float>(const matrix::Matrix<float>&,
                                                            const matrix::Matrix<float>&,
                                                            bool,
                                                            const Params&,
                                                            const std::vector<bool>*);
template DecoderCoreResult<int8_t> Decoder_Core_plain<int8_t>(const matrix::Matrix<int8_t>&,
                                                              const matrix::Matrix<int8_t>&,
                                                              bool,
                                                              const Params&,
                                                              const std::vector<bool>*);

template DecoderCoreResult<float> Decoder_Core_ebchPF<float>(const matrix::Matrix<float>&,
                                                             const matrix::Matrix<float>&,
                                                             bool,
                                                             const Params&,
                                                             const std::vector<bool>*);
template DecoderCoreResult<int8_t> Decoder_Core_ebchPF<int8_t>(const matrix::Matrix<int8_t>&,
                                                               const matrix::Matrix<int8_t>&,
                                                               bool,
                                                               const Params&,
                                                               const std::vector<bool>*);

#define INSTANTIATE_DECODER_CORE_QFLOAT(N) \
template DecoderCoreResult<qfloat::qfloat<N>> Decoder_Core_plain<qfloat::qfloat<N>>( \
    const matrix::Matrix<qfloat::qfloat<N>>& lin_matrix, const matrix::Matrix<qfloat::qfloat<N>>& lch_matrix, bool, const Params&, \
    const std::vector<bool>*); \
template DecoderCoreResult<qfloat::qfloat<N>> Decoder_Core_ebchPF<qfloat::qfloat<N>>( \
    const matrix::Matrix<qfloat::qfloat<N>>& lin_matrix, const matrix::Matrix<qfloat::qfloat<N>>& lch_matrix, bool, const Params&, \
    const std::vector<bool>*);

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
