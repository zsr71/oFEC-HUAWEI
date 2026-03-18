#include "new_float_only/ofec_decoder_hard.hpp"
#include "new_float_only/params.hpp"
#include "new_float_only/rx/ofec/chase/chase256.hpp"
#include "new_float_only/rx/ofec/chase/decoder_core.hpp"

#include <array>
#include <stdexcept>

namespace chase {
namespace {

/**
 * 行级核心实现。
 * 每一行要么进入硬判决 BCH，要么进入 plain Chase；
 * 这里不再有 early-stop fallback，也不再接收 MUX 状态。
 */
DecoderCoreResult Decoder_Core_impl(const matrix::Matrix<float>& lin_matrix,
                                    const matrix::Matrix<float>& lch_matrix,
                                    bool use_hard_decode,
                                    const new_float_only::Params& config) {
  const std::size_t rows = lin_matrix.rows();
  const std::size_t cols = lin_matrix.cols();
  const std::size_t expected_cols = 2 * new_float_only::Params::NUM_SUBBLOCK_COLS *
                                    new_float_only::Params::BITS_PER_SUBBLOCK_DIM;
  if (cols != expected_cols) {
    throw std::invalid_argument("Decoder_Core_plain: unexpected column count");
  }

  DecoderCoreResult result{matrix::Matrix<float>(rows, cols), std::vector<bool>(rows, false)};

  for (std::size_t row = 0; row < rows; ++row) {
    std::array<float, new_float_only::Params::BCH_N> lin_vec{};
    std::array<float, new_float_only::Params::BCH_N> lch_vec{};
    for (std::size_t col = 0; col < cols; ++col) {
      lin_vec[col] = lin_matrix[row][col];
      lch_vec[col] = lch_matrix[row][col];
    }

    std::array<float, new_float_only::Params::BCH_N> y2{};
    bool produced = false;

    // 调试信息按“当前行”裁剪，只保留与本行有关的 Chase 跟踪项。
    new_float_only::Params row_params = config;
    row_params.debug_trace.active_chase_entries.clear();
    row_params.debug_trace.chase_expected_bits = config.debug_trace.chase_expected_bits;
    if (config.debug_trace.chase_expected_bits &&
        row < config.debug_trace.chase_expected_bits->size()) {
      row_params.debug_trace.chase_expected_bits_row =
          &(*config.debug_trace.chase_expected_bits)[row];
    } else {
      row_params.debug_trace.chase_expected_bits_row = nullptr;
    }
    for (const auto& entry : config.debug_trace.active_chase_entries) {
      if (entry.row_index == static_cast<int>(row)) {
        row_params.debug_trace.active_chase_entries.push_back(entry);
      }
    }
    if (!row_params.debug_trace.active_chase_entries.empty()) {
      row_params.debug_trace.chase_decoder_row = static_cast<int>(row);
      row_params.debug_trace.chase_decoder_col =
          row_params.debug_trace.active_chase_entries.front().k;
    } else {
      row_params.debug_trace.chase_decoder_row = -1;
      row_params.debug_trace.chase_decoder_col = -1;
    }

    if (use_hard_decode) {
      produced = new_float_only::perform_hard_decode(lin_vec, lch_vec, y2, row_params);
    } else {
      chase_decode_256_plain(lin_vec.data(), lch_vec.data(), y2.data(), row_params);
      produced = true;
    }

    if (produced) {
      result.produced_rows[row] = true;
      for (std::size_t col = 0; col < cols; ++col) {
        result.lout[row][col] = y2[col];
      }
    }
  }

  return result;
}

}  // namespace

DecoderCoreResult Decoder_Core_plain(const matrix::Matrix<float>& lin_matrix,
                                     const matrix::Matrix<float>& lch_matrix,
                                     bool use_hard_decode,
                                     const new_float_only::Params& config) {
  return Decoder_Core_impl(lin_matrix, lch_matrix, use_hard_decode, config);
}

} // namespace chase
