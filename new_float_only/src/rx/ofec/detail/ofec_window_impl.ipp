#pragma once

#include "ofec_tile_impl.ipp"

#include "new_float_only/ofec_decoder.hpp"

#include <vector>

namespace new_float_only {
namespace detail {

/**
 * 滑窗内部实现。
 * 关键语句：
 * 1. tile_bottom_row/tile_top_row 负责把窗口索引换成全局矩阵坐标；
 * 2. use_history_input 只在 hard tile 且位于最后一个 soft tile 之上时生效；
 * 3. 每个 Tile 的输出会直接覆盖 work_llr 对应范围。
 */
void process_window_impl(matrix::Matrix<float>& work_llr,
                         const matrix::Matrix<float>& channel_llr,
                         size_t win_start, size_t win_end, const new_float_only::Params& p,
                         size_t tile_height_rows, size_t tile_stride_rows, size_t TILES_PER_WIN,
                         bool normalize_extrinsic,
                         const matrix::Matrix<float>* tx_llr_ref,
                         CoreFn core_fn,
                         matrix::Matrix<float>* last_tile_history_accum)
{
  (void)win_start;
  const auto& trace_cfg = p.debug_trace;
  const bool trace_has_coords = (trace_cfg.row >= 0 && trace_cfg.col >= 0);
  const bool trace_mismatch =
      trace_cfg.enable && trace_cfg.log_mismatch && trace_has_coords;
  const long trace_row = trace_has_coords ? trace_cfg.row : -1;
  const long trace_col = trace_has_coords ? trace_cfg.col : -1;
  const size_t N = new_float_only::Params::NUM_SUBBLOCK_COLS * new_float_only::Params::BITS_PER_SUBBLOCK_DIM;
  // Chase 调试编号只用于 trace 文件命名；这里改成 thread_local，
  // 避免并行 sweep 时多个线程同时更新同一个静态计数器。
  thread_local size_t chase_invocation_counter = 0;

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
        // 当前 tile 在全局矩阵中的上下边界。window 从底部往上逐个取 tile，
        // 相邻 tile 之间按 tile_stride_rows 错开。
        const size_t tile_bottom_row = win_end  - t * tile_stride_rows;
        const size_t tile_top_row    = tile_bottom_row + 1 - tile_height_rows;

        // 当前实现下 tile 高度通常是固定值，这里仍按实际边界重新算一次，
        // 便于后续构造与当前 tile 对齐的局部输入矩阵。
        const size_t tile_height_rows_actual = tile_bottom_row - tile_top_row + 1;
        matrix::Matrix<float> tile_in(tile_height_rows_actual, N);
        matrix::Matrix<float> ch_tile(tile_height_rows_actual, N);

        const bool use_hard = hard_tile_mask[t];
        // 对 hard tile，如果它位于最后一个 soft tile 之后，就不再从当前 work_llr
        // 取历史，而是改用 last_tile_history_accum 中累积的最后 soft tile history。
        // 这样做是为了让后续 hard tile 直接吃到“最后一层软译码写回后的历史结果”。
        const bool use_history_input =
            use_hard && last_tile_history_accum &&
            last_soft_tile_idx >= 0 &&
            static_cast<int>(t) > last_soft_tile_idx;

        for (size_t r = 0; r < tile_height_rows_actual; ++r) {
            for (size_t c = 0; c < N; ++c) {
                const size_t global_row = tile_top_row + r;
                if (use_history_input) {
                    // hard tile 的 history 输入改取最后 soft tile 的累计历史。
                    const float hist = (*last_tile_history_accum)[global_row][c];
                    tile_in[r][c] = hist;
                } else {
                    // 普通情况下，tile 输入 history 直接来自当前工作矩阵。
                    tile_in[r][c]  = work_llr[global_row][c];
                }
                // 信道项始终直接从全局 channel_llr 中按相同行列切片。
                ch_tile[r][c]  = channel_llr[global_row][c];
            }
        }

        new_float_only::Params tile_params = p;
        tile_params.beta = pick_float(p.beta_list, t, p.beta);
        tile_params.ALPHA = pick_float(p.ALPHA_LIST, t, p.ALPHA);
        tile_params.debug_trace.chase_tile_index = static_cast<int>(t);
        tile_params.debug_trace.chase_invocation =
            static_cast<int>(++chase_invocation_counter);

        const bool capture_history =
            last_tile_history_accum &&
            (use_hard || static_cast<int>(t) == last_soft_tile_idx);

    TileProcessResult tile_result = process_tile_impl(tile_in, ch_tile, tile_params,
                                                      /*tile_top_row_global=*/tile_top_row,
                                                      /*use_hard_decode=*/use_hard,
                                                      /*normalize_extrinsic=*/normalize_extrinsic,
                                                      tx_llr_ref,
                                                      core_fn,
                                                      last_tile_history_accum,
                                                      capture_history);

    for (size_t r = 0; r < tile_height_rows_actual; ++r) {
      const size_t global_row = tile_top_row + r;
      for (size_t c = 0; c < work_llr.cols(); ++c) {
        const auto incoming = tile_result.tile_out[r][c];
        if (trace_mismatch &&
            static_cast<long>(global_row) == trace_row &&
            static_cast<long>(c) == trace_col) {
          const float existing_val = work_llr[global_row][c];
          const float incoming_val = incoming;
          const float channel_val  = channel_llr[global_row][c];
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
} // namespace new_float_only
