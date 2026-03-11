#pragma once

#include "ofec_window_impl.ipp"

#include "newcode/decoder_api.hpp"
#include "newcode/common/qfloat/llr_utils.hpp"
#include "newcode/ofec/mux/mux_config_validate.hpp"
#include "newcode/ofec/mux/mux_group_config_validate.hpp"
#include "newcode/quantized_llr_dump.hpp"

#include <string>
#include <utility>
#include <vector>

namespace newcode {
namespace detail {

template <typename LLR>
matrix::Matrix<LLR> ofec_decode_llr_impl(const matrix::Matrix<LLR>& llr_mat, const newcode::Params& p,
                                 std::vector<TileEarlyStopCounter>* tile_stats,
                                 bool normalize_extrinsic,
                                 const matrix::Matrix<float>* tx_llr_ref,
                                 CoreFn<typename LinMatrixAdapter<LLR>::core_type> core_fn)
{
  const size_t N = newcode::Params::NUM_SUBBLOCK_COLS * newcode::Params::BITS_PER_SUBBLOCK_DIM;

  const size_t RROWS = llr_mat.rows();
  const size_t CCOLS = llr_mat.cols();
  if (CCOLS != N)
      throw std::invalid_argument("ofec_decode_llr: llr_mat cols != N.");

  assert(p.valid());
  const auto mux_ok =
      newcode::mux::validate_siso_active_list(p.SISO_ACTIVE_LIST, p.TILES_PER_WIN);
  if (!mux_ok.ok) {
    throw std::invalid_argument("ofec_decode_llr: " + mux_ok.error);
  }
  const std::size_t rows_to_decode =
      static_cast<std::size_t>(p.CHASE_SBR) *
      newcode::Params::BITS_PER_SUBBLOCK_DIM;
  const auto group_ok =
      newcode::mux::validate_mux_group_g(p.MUX_GROUP_G, rows_to_decode);
  if (!group_ok.ok) {
    throw std::invalid_argument("ofec_decode_llr: " + group_ok.error);
  }
  if (p.MUX_ENABLE_RECONFIG) {
    for (std::size_t tile_idx = 0; tile_idx < p.SISO_ACTIVE_LIST.size(); ++tile_idx) {
      const auto reconfig_ok = newcode::mux::validate_mux_reconfig_runtime(
          p.MUX_GROUP_G, p.SISO_ACTIVE_LIST[tile_idx], rows_to_decode);
      if (!reconfig_ok.ok) {
        throw std::invalid_argument(
            "ofec_decode_llr: tile " + std::to_string(tile_idx) + ": " +
            reconfig_ok.error);
      }
    }
  }

  const size_t TILE_HEIGHT_ROWS = p.tile_height_rows();
  const size_t TILE_STRIDE_ROWS = p.tile_stride_rows();
  const size_t WIN_HEIGHT_ROWS  = p.win_height_rows();
  const size_t POP_PUSH_ROWS    = p.pop_push_rows();
  const size_t TILES_PER_WIN    = p.TILES_PER_WIN;


  matrix::Matrix<LLR> channel_llr = llr_mat;
  matrix::Matrix<LLR> work_llr(RROWS, N);
  for (size_t r = 0; r < RROWS; ++r)
      for (size_t c = 0; c < N; ++c)
          work_llr[r][c] = qfloat::llr_from_float<LLR>(0.0f);
  matrix::Matrix<float> last_tile_history_llr(RROWS, N);

  if (RROWS < WIN_HEIGHT_ROWS) {
    if (tile_stats) {
      tile_stats->assign(p.TILES_PER_WIN, TileEarlyStopCounter{});
    }
    return channel_llr;
  }

  size_t win_start     = p.initial_win_start_rows();
  const size_t last_ws = RROWS - WIN_HEIGHT_ROWS;

  std::vector<TileEarlyStopCounter> local_tile_stats;
  std::vector<TileEarlyStopCounter>* stats_ptr = nullptr;
  if (tile_stats) {
    local_tile_stats.assign(p.TILES_PER_WIN, TileEarlyStopCounter{});
    stats_ptr = &local_tile_stats;
  }

  while (win_start <= last_ws) {
    const size_t win_end = win_start + WIN_HEIGHT_ROWS - 1;
    process_window_impl<LLR>(work_llr, channel_llr,
                             win_start, win_end, p,
                             TILE_HEIGHT_ROWS, TILE_STRIDE_ROWS, TILES_PER_WIN,
                             stats_ptr,
                             normalize_extrinsic,
                             tx_llr_ref,
                             core_fn,
                             &last_tile_history_llr);

    win_start += POP_PUSH_ROWS;
  }

  if (tile_stats) {
    *tile_stats = std::move(local_tile_stats);
  }

  matrix::Matrix<LLR> out(RROWS, N);
  for (size_t r = 0; r < RROWS; ++r) {
    for (size_t c = 0; c < N; ++c) {
      const float sum = qfloat::llr_to_float(channel_llr[r][c]) + last_tile_history_llr[r][c];
      out[r][c] = qfloat::llr_from_float<LLR>(sum);
    }
  }

  // 可选：保存窗口累积后的 work_llr（解码前）
  if (p.DUMP_WORK_LLR) {
    matrix::Matrix<float> work_float(RROWS, N);
    for (size_t r = 0; r < RROWS; ++r)
      for (size_t c = 0; c < N; ++c)
        work_float[r][c] = qfloat::llr_to_float(work_llr[r][c]);

    DecodeRequest dump_req{
        .label = "work_llr",
        .channel_llr = work_float,
        .tx_llr_ref = nullptr,
        .params = p,
        .format = LlrFormat::Float,
        .quant_bits = 16,
        .quant_clip = 0.0f,
        .normalize_extrinsic = true,
        .quiet = false,
        .dump_quantized_llr = true,
        .quantized_llr_output_path = p.WORK_LLR_OUTPUT_PATH.empty()
            ? std::string("data/llr/quantized_llr_work.txt")
            : p.WORK_LLR_OUTPUT_PATH,
        .dump_float_llr = false,
        .float_llr_output_path = {},
        .dump_quantized_codes = false,
        .quantized_codes_output_path = {},
        .dump_work_llr = true,
        .work_llr_output_path = {}
    };
    (void)dump_quantized_llr(work_float, dump_req, true,
                             dump_req.quantized_llr_output_path, "_work");
  }

  return out;
}

} // namespace detail
} // namespace newcode
