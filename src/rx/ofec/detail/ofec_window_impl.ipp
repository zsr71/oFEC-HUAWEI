#pragma once

#include "ofec_tile_impl.ipp"

#include "newcode/ofec_decoder.hpp"
#include "newcode/common/qfloat/llr_utils.hpp"

#include <vector>

namespace newcode {
namespace detail {

template <typename LLR>
void process_window_impl(matrix::Matrix<LLR>& work_llr,
                         const matrix::Matrix<LLR>& channel_llr,
                         size_t win_start, size_t win_end, const newcode::Params& p,
                         size_t tile_height_rows, size_t tile_stride_rows, size_t TILES_PER_WIN,
                         std::vector<TileEarlyStopCounter>* tile_stats,
                         bool normalize_extrinsic,
                         const matrix::Matrix<float>* tx_llr_ref,
                         CoreFn<typename LinMatrixAdapter<LLR>::core_type> core_fn,
                         matrix::Matrix<float>* last_tile_history_accum)
{
  (void)win_start;
  const auto& trace_cfg = p.debug_trace;
  const bool trace_has_coords = (trace_cfg.row >= 0 && trace_cfg.col >= 0);
  const bool trace_mismatch =
      trace_cfg.enable && trace_cfg.log_mismatch && trace_has_coords;
  const long trace_row = trace_has_coords ? trace_cfg.row : -1;
  const long trace_col = trace_has_coords ? trace_cfg.col : -1;
  const size_t N = newcode::Params::NUM_SUBBLOCK_COLS * newcode::Params::BITS_PER_SUBBLOCK_DIM;
  static size_t chase_invocation_counter = 0;

  auto pick_float = [](const std::vector<float>& tbl, size_t idx, float fallback) -> float {
      return (idx < tbl.size()) ? tbl[idx] : fallback;
  };
  auto pick_int = [](const std::vector<int>& tbl, size_t idx, int fallback) -> int {
      return (idx < tbl.size()) ? tbl[idx] : fallback;
  };
  std::vector<bool> hard_tile_mask(TILES_PER_WIN);
  int last_soft_tile_idx = -1;
  for (size_t t = 0; t < TILES_PER_WIN; ++t) {
      const bool is_hard = pick_int(p.HARD_TILE_LIST, t, p.HARD_DECODE_DEFAULT ? 1 : 0) != 0;
      hard_tile_mask[t] = is_hard;
      if (!is_hard) {
          last_soft_tile_idx = static_cast<int>(t);
      }
  }

  for (size_t t = 0; t < TILES_PER_WIN; ++t)
  {
        const size_t tile_bottom_row = win_end  - t * tile_stride_rows;
        const size_t tile_top_row    = tile_bottom_row + 1 - tile_height_rows;

        const size_t tile_height_rows_actual = tile_bottom_row - tile_top_row + 1;
        matrix::Matrix<LLR> tile_in(tile_height_rows_actual, N);
        matrix::Matrix<LLR> ch_tile(tile_height_rows_actual, N);

        const bool use_hard = hard_tile_mask[t];
        const bool use_history_input =
            use_hard && last_tile_history_accum &&
            last_soft_tile_idx >= 0 &&
            static_cast<int>(t) > last_soft_tile_idx;

        for (size_t r = 0; r < tile_height_rows_actual; ++r) {
            for (size_t c = 0; c < N; ++c) {
                const size_t global_row = tile_top_row + r;
                if (use_history_input) {
                    const float hist = (*last_tile_history_accum)[global_row][c];
                    tile_in[r][c] = qfloat::llr_from_float<LLR>(hist);
                } else {
                    tile_in[r][c]  = work_llr[global_row][c];
                }
                ch_tile[r][c]  = channel_llr[global_row][c];
            }
        }

        newcode::Params tile_params = p;
        tile_params.beta = pick_float(p.beta_list, t, p.beta);
        tile_params.ALPHA = pick_float(p.ALPHA_LIST, t, p.ALPHA);
        tile_params.debug_trace.chase_tile_index = static_cast<int>(t);
        tile_params.debug_trace.chase_invocation =
            static_cast<int>(++chase_invocation_counter);

        const bool capture_history =
            last_tile_history_accum &&
            (use_hard || static_cast<int>(t) == last_soft_tile_idx);

    TileProcessResult<LLR> tile_result = process_tile_impl<LLR>(tile_in, ch_tile, tile_params,
                                                                /*tile_top_row_global=*/tile_top_row,
                                                                /*use_hard_decode=*/use_hard,
                                                                /*normalize_extrinsic=*/normalize_extrinsic,
                                                                tx_llr_ref,
                                                                core_fn,
                                                                last_tile_history_accum,
                                                                capture_history);

    if (tile_stats && t < tile_stats->size()) {
      auto& counter = (*tile_stats)[t];
      counter.total += 1;
      if (tile_result.early_stop_triggered) {
        counter.triggered += 1;
      }
      counter.row_total += tile_result.rows_total;
      counter.row_triggered += tile_result.rows_early_stop;
    }

    for (size_t r = 0; r < tile_height_rows_actual; ++r) {
      const size_t global_row = tile_top_row + r;
      for (size_t c = 0; c < work_llr.cols(); ++c) {
        const auto incoming = tile_result.tile_out[r][c];
        if (trace_mismatch &&
            static_cast<long>(global_row) == trace_row &&
            static_cast<long>(c) == trace_col) {
          const float existing_val = qfloat::llr_to_float(work_llr[global_row][c]);
          const float incoming_val = qfloat::llr_to_float(incoming);
          const float channel_val  = qfloat::llr_to_float(channel_llr[global_row][c]);
          if (existing_val != incoming_val) {
            std::cout << "Mismatch at work_llr[" << trace_row << "][" << trace_col
                      << "]: tile index " << t
                      << " incoming=" << incoming_val << '\n'
                      << " channel =" << channel_val << '\n';
          }
        }
        work_llr[global_row][c] = incoming;
      }
    }
  }
}

} // namespace detail
} // namespace newcode
