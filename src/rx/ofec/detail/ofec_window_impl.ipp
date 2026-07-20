#pragma once

#include "ofec_tile_impl.ipp"

#include "newcode/ofec_decoder.hpp"
#include "newcode/common/qfloat/llr_utils.hpp"
#include "newcode/ofec/mux/mux_siso_budget.hpp"

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
  // 输入:
  // - work_llr: 当前全局工作矩阵，保存已经累积的外信息/历史信息，会被原地更新。
  // - channel_llr: 原始信道矩阵，只读。
  // - [win_start, win_end]: 当前处理窗口的全局行范围。
  // - p: 参数集合。
  // - tile_height_rows/tile_stride_rows/TILES_PER_WIN: 当前窗口的 tile 划分方式。
  // - tile_stats: 可选 early-stop 统计输出。
  // - normalize_extrinsic: 是否归一化每个 tile 的外信息。
  // - tx_llr_ref: 可选参考矩阵，用于调试。
  // - core_fn: Chase core 回调。
  // - last_tile_history_accum: 记录“最后一个有效 tile”的历史值，供窗口外层合成输出。
  // 输出:
  // - 无返回值；通过原地更新 work_llr / tile_stats / last_tile_history_accum 生效。
  // 用途:
  // - 在一个滑动窗口内部，按 tile 顺序切片、解码、写回，并把结果覆盖回全局工作矩阵。
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
      // 先根据配置判断每个 tile 走软译码还是硬判回退。
      const bool is_hard = pick_int(p.HARD_TILE_LIST, t, p.HARD_DECODE_DEFAULT ? 1 : 0) != 0;
      hard_tile_mask[t] = is_hard;
      if (!is_hard) {
          last_soft_tile_idx = static_cast<int>(t);
      }
  }

  for (size_t t = 0; t < TILES_PER_WIN; ++t)
  {
        // 当前 tile 在窗口中的全局行范围。
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
                    // 某些硬判 tile 直接吃“最后 soft tile 留下的历史值”。
                    const float hist = (*last_tile_history_accum)[global_row][c];
                    tile_in[r][c] = qfloat::llr_from_float<LLR>(hist);
                } else {
                    // 常规路径从 work_llr 读取当前先验/外信息。
                    tile_in[r][c]  = work_llr[global_row][c];
                }
                ch_tile[r][c]  = channel_llr[global_row][c];
            }
        }

        newcode::Params tile_params = p;
        tile_params.beta = pick_float(p.beta_list, t, p.beta);
        tile_params.EARLY_STOP_ACTION_SIGN_BETA =
            pick_float(p.EARLY_STOP_ACTION_SIGN_BETA_LIST,
                       t,
                       p.EARLY_STOP_ACTION_SIGN_BETA);
        tile_params.ENABLE_EARLY_STOP =
            pick_int(p.EARLY_STOP_ENABLE_LIST,
                     t,
                     p.ENABLE_EARLY_STOP ? 1 : 0) != 0;
        tile_params.EARLY_STOP_CONDITION_MODE =
            pick_int(p.EARLY_STOP_CONDITION_MODE_LIST,
                     t,
                     p.EARLY_STOP_CONDITION_MODE);
        tile_params.EARLY_STOP_ACTION_MODE =
            pick_int(p.EARLY_STOP_ACTION_MODE_LIST,
                     t,
                     p.EARLY_STOP_ACTION_MODE);
        tile_params.EARLY_STOP_BIND_GROUP_SIZE =
            pick_int(p.EARLY_STOP_BIND_GROUP_SIZE_LIST,
                     t,
                     p.EARLY_STOP_BIND_GROUP_SIZE);
        tile_params.HYBRID_ENABLE =
            pick_int(p.HYBRID_ENABLE_LIST,
                     t,
                     p.HYBRID_ENABLE ? 1 : 0) != 0;
        tile_params.HYBRID_HARD_LLR_MAG =
            pick_float(p.HYBRID_HARD_LLR_MAG_LIST,
                       t,
                       p.HYBRID_HARD_LLR_MAG);
        tile_params.ALPHA = pick_float(p.ALPHA_LIST, t, p.ALPHA);
        tile_params.debug_trace.chase_tile_index = static_cast<int>(t);
        tile_params.debug_trace.chase_invocation =
            static_cast<int>(++chase_invocation_counter);
        const int siso_active_for_tile =
            newcode::mux::pick_siso_active_for_tile(p.SISO_ACTIVE_LIST, t);
        const int hiho_active_for_tile =
            newcode::mux::pick_hiho_active_for_tile(p.HIHO_ACTIVE_LIST, t);

        const bool capture_history =
            last_tile_history_accum &&
            (use_hard || static_cast<int>(t) == last_soft_tile_idx);

    TileProcessResult<LLR> tile_result = process_tile_impl<LLR>(tile_in, ch_tile, tile_params,
                                                                /*tile_top_row_global=*/tile_top_row,
                                                                siso_active_for_tile,
                                                                hiho_active_for_tile,
                                                                /*use_hard_decode=*/use_hard,
                                                                /*normalize_extrinsic=*/normalize_extrinsic,
                                                                tx_llr_ref,
                                                                core_fn,
                                                                last_tile_history_accum,
                                                                capture_history);

    if (tile_stats && t < tile_stats->size()) {
      // 将本 tile 的 early-stop 统计累加到窗口级统计数组。
      auto& counter = (*tile_stats)[t];
      counter.total += 1;
      if (tile_result.early_stop_triggered) {
        counter.triggered += 1;
      }
      counter.row_total += tile_result.rows_total;
      counter.row_triggered += tile_result.rows_early_stop;
      counter.row_hard_finish += tile_result.rows_hard_finish;
      counter.row_need_siso_before_mux += tile_result.rows_need_siso_before_mux;
      counter.row_unscheduled += tile_result.rows_unscheduled;
      counter.samples.push_back(TileEarlyStopSample{
          .invocation = counter.total,
          .tile_index = t,
          .rows_total = tile_result.rows_total,
          .rows_passed = tile_result.rows_early_stop,
          .rows_hard_finish = tile_result.rows_hard_finish,
          .rows_need_siso_before_mux = tile_result.rows_need_siso_before_mux,
          .rows_unscheduled = tile_result.rows_unscheduled,
      });
      if (tile_result.has_group_bind_debug_sample) {
        auto group_bind_sample = tile_result.group_bind_debug_sample;
        group_bind_sample.invocation = counter.total;
        group_bind_sample.tile_index = t;
        counter.group_bind_debug_samples.push_back(std::move(group_bind_sample));
      }
      counter.hybrid_class_counts.push_back(tile_result.hybrid_class_count);
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
        // tile 解码写回后的结果覆盖到全局工作矩阵，供后续 tile/窗口继续使用。
        work_llr[global_row][c] = incoming;
      }
    }
  }
}

} // namespace detail
} // namespace newcode
