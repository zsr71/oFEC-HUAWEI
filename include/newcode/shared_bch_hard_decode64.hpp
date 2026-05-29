#pragma once

#include <cstddef>
#include <vector>

#include "newcode/common/matrix/matrix.hpp"
#include "newcode/params.hpp"
#include "newcode/rx/ofec/chase/decoder_core.hpp"

namespace newcode {

struct SharedBatchSlice {
  std::size_t row_offset = 0;
  std::size_t row_count = 0;
  std::size_t stream_id = 0;
};

template <typename CoreLLR>
struct SharedBatchPlan {
  matrix::Matrix<CoreLLR> lin_matrix;
  matrix::Matrix<CoreLLR> lch_matrix;
  Params params_for_core;
  std::vector<SharedBatchSlice> slices;
};

struct SharedDecoderSliceStats {
  std::size_t slice_index = 0;
  std::size_t stream_id = 0;
  std::size_t row_offset = 0;
  std::size_t row_count = 0;
  std::size_t produced_rows = 0;
  std::size_t failed_rows = 0;
};

struct SharedDecoderStats {
  std::size_t total_rows = 0;
  std::size_t produced_rows = 0;
  std::size_t failed_rows = 0;
  std::vector<SharedDecoderSliceStats> slices;
};

template <typename CoreLLR>
chase::DecoderCoreResult<CoreLLR> shared_bch_hard_decode_core(
    const SharedBatchPlan<CoreLLR>& batch,
    SharedDecoderStats* stats = nullptr);

template <typename CoreLLR>
chase::DecoderCoreResult<CoreLLR> shared_bch_hard_decode64(
    const matrix::Matrix<CoreLLR>& lin_matrix,
    const matrix::Matrix<CoreLLR>& lch_matrix,
    const Params& p);

extern template chase::DecoderCoreResult<float> shared_bch_hard_decode_core<float>(
    const SharedBatchPlan<float>& batch,
    SharedDecoderStats* stats);

extern template chase::DecoderCoreResult<int8_t> shared_bch_hard_decode_core<int8_t>(
    const SharedBatchPlan<int8_t>& batch,
    SharedDecoderStats* stats);

extern template chase::DecoderCoreResult<float> shared_bch_hard_decode64<float>(
    const matrix::Matrix<float>& lin_matrix,
    const matrix::Matrix<float>& lch_matrix,
    const Params& p);

extern template chase::DecoderCoreResult<int8_t> shared_bch_hard_decode64<int8_t>(
    const matrix::Matrix<int8_t>& lin_matrix,
    const matrix::Matrix<int8_t>& lch_matrix,
    const Params& p);

}  // namespace newcode
