#include "newcode/shared_bch_hard_decode64.hpp"

#include <array>
#include <stdexcept>

#include "newcode/ofec_decoder_hard.hpp"

namespace newcode {
namespace {

template <typename CoreLLR>
constexpr std::size_t kCodeBits = Params::BCH_N;

template <typename CoreLLR>
void validate_shared_batch_plan(const SharedBatchPlan<CoreLLR>& batch) {
  if (batch.lin_matrix.rows() == 0 ||
      batch.lin_matrix.cols() != kCodeBits<CoreLLR>) {
    throw std::invalid_argument(
        "shared_bch_hard_decode_core: lin_matrix must have shape R x 256");
  }
  if (batch.lch_matrix.rows() != batch.lin_matrix.rows() ||
      batch.lch_matrix.cols() != batch.lin_matrix.cols()) {
    throw std::invalid_argument(
        "shared_bch_hard_decode_core: lch_matrix shape must match lin_matrix");
  }
  for (const auto& slice : batch.slices) {
    if (slice.row_offset + slice.row_count > batch.lin_matrix.rows()) {
      throw std::invalid_argument(
          "shared_bch_hard_decode_core: slice range exceeds batch row count");
    }
  }
}

template <typename CoreLLR>
void load_row(const matrix::Matrix<CoreLLR>& matrix_in,
              std::size_t row,
              std::array<CoreLLR, Params::BCH_N>* out) {
  for (std::size_t col = 0; col < Params::BCH_N; ++col) {
    (*out)[col] = matrix_in[row][col];
  }
}

}  // namespace

template <typename CoreLLR>
chase::DecoderCoreResult<CoreLLR> shared_bch_hard_decode_core(
    const SharedBatchPlan<CoreLLR>& batch,
    SharedDecoderStats* stats) {
  validate_shared_batch_plan(batch);

  chase::DecoderCoreResult<CoreLLR> result;
  result.lout =
      matrix::Matrix<float>(batch.lin_matrix.rows(), batch.lin_matrix.cols());
  result.produced_rows.assign(batch.lin_matrix.rows(), false);

  for (std::size_t row = 0; row < batch.lin_matrix.rows(); ++row) {
    std::array<CoreLLR, Params::BCH_N> lin_vec{};
    std::array<CoreLLR, Params::BCH_N> lch_vec{};
    std::array<float, Params::BCH_N> y2{};
    load_row(batch.lin_matrix, row, &lin_vec);
    load_row(batch.lch_matrix, row, &lch_vec);

    const bool produced =
        perform_hard_decode<CoreLLR>(lin_vec, lch_vec, y2, batch.params_for_core);
    result.produced_rows[row] = produced;
    if (!produced) {
      continue;
    }
    for (std::size_t col = 0; col < Params::BCH_N; ++col) {
      result.lout[row][col] = y2[col];
    }
  }

  if (stats) {
    stats->total_rows = result.produced_rows.size();
    stats->produced_rows = 0;
    for (bool produced : result.produced_rows) {
      if (produced) {
        ++stats->produced_rows;
      }
    }
    stats->failed_rows = stats->total_rows - stats->produced_rows;
    stats->slices.clear();
    stats->slices.reserve(batch.slices.size());
    for (std::size_t slice_index = 0; slice_index < batch.slices.size();
         ++slice_index) {
      const auto& slice = batch.slices[slice_index];
      SharedDecoderSliceStats slice_stats;
      slice_stats.slice_index = slice_index;
      slice_stats.stream_id = slice.stream_id;
      slice_stats.row_offset = slice.row_offset;
      slice_stats.row_count = slice.row_count;
      for (std::size_t row = 0; row < slice.row_count; ++row) {
        if (result.produced_rows[slice.row_offset + row]) {
          ++slice_stats.produced_rows;
        }
      }
      slice_stats.failed_rows = slice_stats.row_count - slice_stats.produced_rows;
      stats->slices.push_back(slice_stats);
    }
  }

  return result;
}

template <typename CoreLLR>
chase::DecoderCoreResult<CoreLLR> shared_bch_hard_decode64(
    const matrix::Matrix<CoreLLR>& lin_matrix,
    const matrix::Matrix<CoreLLR>& lch_matrix,
    const Params& p) {
  if (lin_matrix.rows() != 64 || lin_matrix.cols() != kCodeBits<CoreLLR>) {
    throw std::invalid_argument(
        "shared_bch_hard_decode64: lin_matrix must be 64 x 256");
  }
  if (lch_matrix.rows() != 64 || lch_matrix.cols() != kCodeBits<CoreLLR>) {
    throw std::invalid_argument(
        "shared_bch_hard_decode64: lch_matrix must be 64 x 256");
  }

  SharedBatchPlan<CoreLLR> batch;
  batch.lin_matrix = lin_matrix;
  batch.lch_matrix = lch_matrix;
  batch.params_for_core = p;
  batch.slices.push_back(SharedBatchSlice{
      .row_offset = 0,
      .row_count = lin_matrix.rows(),
      .stream_id = 0,
  });
  return shared_bch_hard_decode_core(batch, nullptr);
}

template chase::DecoderCoreResult<float> shared_bch_hard_decode_core<float>(
    const SharedBatchPlan<float>& batch,
    SharedDecoderStats* stats);

template chase::DecoderCoreResult<int8_t> shared_bch_hard_decode_core<int8_t>(
    const SharedBatchPlan<int8_t>& batch,
    SharedDecoderStats* stats);

template chase::DecoderCoreResult<float> shared_bch_hard_decode64<float>(
    const matrix::Matrix<float>& lin_matrix,
    const matrix::Matrix<float>& lch_matrix,
    const Params& p);

template chase::DecoderCoreResult<int8_t> shared_bch_hard_decode64<int8_t>(
    const matrix::Matrix<int8_t>& lin_matrix,
    const matrix::Matrix<int8_t>& lch_matrix,
    const Params& p);

}  // namespace newcode
