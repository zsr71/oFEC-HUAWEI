#pragma once

#include "ofec_tile_impl.ipp"
#include "ofec_level56_shared.ipp"

#include "newcode/ofec_decoder.hpp"
#include "newcode/common/qfloat/llr_utils.hpp"
#include "newcode/ofec/mux/mux_config_validate.hpp"
#include "newcode/ofec/mux/mux_siso_budget.hpp"

#include <iterator>
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
                         matrix::Matrix<float>* last_tile_history_accum,
                         std::size_t* level56_shared_invocation,
                         Level56TemporalState<LLR>* level56_temporal_state = nullptr)
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
  validate_level56_shared_config(p);
  const std::size_t mux_tile_count =
      p.LEVEL56_SHARED_ENABLE ? kLevel5TileIndex : TILES_PER_WIN;
  const auto mux_ok = p.LEVEL56_SHARED_ENABLE
      ? newcode::mux::validate_siso_active_prefix(p.SISO_ACTIVE_LIST,
                                                  mux_tile_count)
      : newcode::mux::validate_siso_active_list(p.SISO_ACTIVE_LIST,
                                                mux_tile_count);
  if (!mux_ok.ok) {
    throw std::invalid_argument("process_window_impl: " + mux_ok.error);
  }
  const auto hiho_ok = p.LEVEL56_SHARED_ENABLE
      ? newcode::mux::validate_hiho_active_prefix(p.HIHO_ACTIVE_LIST,
                                                  mux_tile_count)
      : newcode::mux::validate_hiho_active_list(p.HIHO_ACTIVE_LIST,
                                                mux_tile_count);
  if (!hiho_ok.ok) {
    throw std::invalid_argument("process_window_impl: " + hiho_ok.error);
  }
  const auto& trace_cfg = p.debug_trace;
  const bool trace_has_coords = (trace_cfg.row >= 0 && trace_cfg.col >= 0);
  const bool trace_mismatch =
      trace_cfg.enable && trace_cfg.log_mismatch && trace_has_coords;
  const long trace_row = trace_has_coords ? trace_cfg.row : -1;
  const long trace_col = trace_has_coords ? trace_cfg.col : -1;
  const size_t N = newcode::Params::NUM_SUBBLOCK_COLS * newcode::Params::BITS_PER_SUBBLOCK_DIM;
  const size_t rows_to_decode =
      static_cast<size_t>(p.CHASE_SBR) *
      newcode::Params::BITS_PER_SUBBLOCK_DIM;
  static size_t chase_invocation_counter = 0;

  auto pick_float = [](const std::vector<float>& tbl, size_t idx, float fallback) -> float {
      return (idx < tbl.size()) ? tbl[idx] : fallback;
  };
  auto pick_int = [](const std::vector<int>& tbl, size_t idx, int fallback) -> int {
      return (idx < tbl.size()) ? tbl[idx] : fallback;
  };
  auto make_tile_params = [&](size_t t) {
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
      return tile_params;
  };

  auto copy_tile_from_global = [&](size_t top_row) {
    matrix::Matrix<LLR> tile(tile_height_rows, N);
    for (size_t r = 0; r < tile_height_rows; ++r) {
      for (size_t c = 0; c < N; ++c) {
        tile[r][c] = work_llr[top_row + r][c];
      }
    }
    return tile;
  };

  auto make_temporal_batch = [&](size_t batch_win_start)
      -> std::optional<Level56TemporalBatch<LLR>> {
    const size_t batch_win_end =
        batch_win_start + (win_end - win_start);
    if (batch_win_end >= channel_llr.rows()) {
      return std::nullopt;
    }
    const size_t bottom5 = batch_win_end - kLevel5TileIndex * tile_stride_rows;
    const size_t top5 = bottom5 + 1 - tile_height_rows;
    const size_t bottom6 = batch_win_end - kLevel6TileIndex * tile_stride_rows;
    const size_t top6 = bottom6 + 1 - tile_height_rows;
    Level56TemporalBatch<LLR> batch;
    batch.tile_top5 = top5;
    batch.tile_top6 = top6;
    batch.tile_in5 = copy_tile_from_global(top5);
    batch.tile_in6 = copy_tile_from_global(top6);
    batch.ch_tile5 = matrix::Matrix<LLR>(tile_height_rows, N);
    batch.ch_tile6 = matrix::Matrix<LLR>(tile_height_rows, N);
    for (size_t r = 0; r < tile_height_rows; ++r) {
      for (size_t c = 0; c < N; ++c) {
        batch.ch_tile5[r][c] = channel_llr[top5 + r][c];
        batch.ch_tile6[r][c] = channel_llr[top6 + r][c];
      }
    }
    const auto params5 = make_tile_params(kLevel5TileIndex);
    const auto params6 = make_tile_params(kLevel6TileIndex);
    auto prep5 = prepare_tile_inputs(batch.tile_in5, batch.ch_tile5, params5,
                                     top5, params5.CHASE_SBR, rows_to_decode,
                                     tx_llr_ref);
    auto prep6 = prepare_tile_inputs(batch.tile_in6, batch.ch_tile6, params6,
                                     top6, params6.CHASE_SBR, rows_to_decode,
                                     tx_llr_ref);
    batch.early5 = detect_level56_early_stop(prep5, params5);
    batch.early6 = detect_level56_early_stop(prep6, params6);
    append_level56_early_stop_only_entries(
        prep5, batch.early5.effective, 5, &batch.entries);
    append_level56_early_stop_only_entries(
        prep6, batch.early6.effective, 6, &batch.entries);
    return batch;
  };

  std::optional<Level56TemporalBatch<LLR>> lookahead_batch;
  auto update_tile_stats = [&](size_t t,
                               const TileProcessResult<LLR>& tile_result) {
    if (!tile_stats || t >= tile_stats->size()) {
      return;
    }
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
  };
  auto write_tile_to_work = [&](size_t t,
                                size_t tile_top_row,
                                const matrix::Matrix<LLR>& tile_out) {
    for (size_t r = 0; r < tile_out.rows(); ++r) {
      const size_t global_row = tile_top_row + r;
      for (size_t c = 0; c < work_llr.cols(); ++c) {
        const auto incoming = tile_out[r][c];
        if (trace_mismatch &&
            static_cast<long>(global_row) == trace_row &&
            static_cast<long>(c) == trace_col) {
          const float existing_val = qfloat::llr_to_float(work_llr[global_row][c]);
          const float incoming_val = qfloat::llr_to_float(incoming);
          const float channel_val = qfloat::llr_to_float(channel_llr[global_row][c]);
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
        if (p.LEVEL56_SHARED_ENABLE && t == kLevel5TileIndex) {
          const size_t bottom5 = win_end - kLevel5TileIndex * tile_stride_rows;
          const size_t top5 = bottom5 + 1 - tile_height_rows;
          const size_t bottom6 = win_end - kLevel6TileIndex * tile_stride_rows;
          const size_t top6 = bottom6 + 1 - tile_height_rows;
          matrix::Matrix<LLR> tile_in5(tile_height_rows, N);
          matrix::Matrix<LLR> tile_in6(tile_height_rows, N);
          matrix::Matrix<LLR> ch_tile5(tile_height_rows, N);
          matrix::Matrix<LLR> ch_tile6(tile_height_rows, N);
          for (size_t r = 0; r < tile_height_rows; ++r) {
            for (size_t c = 0; c < N; ++c) {
              tile_in5[r][c] = work_llr[top5 + r][c];
              tile_in6[r][c] = work_llr[top6 + r][c];
              ch_tile5[r][c] = channel_llr[top5 + r][c];
              ch_tile6[r][c] = channel_llr[top6 + r][c];
            }
          }
          auto params5 = make_tile_params(kLevel5TileIndex);
          auto params6 = make_tile_params(kLevel6TileIndex);
          const std::size_t invocation = level56_shared_invocation
              ? (*level56_shared_invocation)++
              : 0u;
          if (level56_temporal_state && p.LEVEL56_TEMPORAL_LOOKAHEAD_ENABLE) {
            // Level 1-4 have prepared the relevant rows by this point. Build
            // the new future batch before any current Level 5/6 writeback.
            lookahead_batch =
                make_temporal_batch(win_start + p.pop_push_rows());
            // t=1 is the batch being decoded now. Rebuild it from the latest
            // work_llr so prior-window Level5/6 writeback is visible. The
            // lookahead snapshot is only a prediction input for this window.
            auto current_batch = make_temporal_batch(win_start);
            if (!current_batch) {
              throw std::logic_error(
                  "LEVEL56 temporal current batch is unavailable");
            }
            const bool have_previous =
                level56_temporal_state->previous.has_value();
            const std::size_t x = have_previous
                                      ? level56_temporal_pending_group_count(
                                            level56_temporal_state->previous->entries)
                                      : 0;
            const std::size_t k1 =
                level56_temporal_candidate_group_count(current_batch->entries);
            const std::size_t k2 = lookahead_batch
                                      ? level56_temporal_candidate_group_count(
                                            lookahead_batch->entries)
                                      : 0;
            const auto temporal_branch = select_level56_temporal_branch(
                have_previous, lookahead_batch.has_value(), x, k1, k2,
                static_cast<std::size_t>(
                    p.LEVEL56_TEMPORAL_GROUP_LOAD_THRESHOLD));
            const bool supplement_history =
                temporal_branch == Level56TemporalBranch::SupplementHistory;
            const std::size_t t0_pending_before =
                have_previous
                    ? level56_temporal_pending_code_count(
                          level56_temporal_state->previous->entries)
                    : 0;
            const std::size_t t1_pending_before =
                level56_temporal_pending_code_count(current_batch->entries);
            std::vector<bool> t0_was_pending;
            if (supplement_history) {
              t0_was_pending.reserve(
                  level56_temporal_state->previous->entries.size());
              for (const auto& entry :
                   level56_temporal_state->previous->entries) {
                t0_was_pending.push_back(is_level56_pending_decode(entry));
              }
            }
            Level56SharedResult<LLR> current_result;
            std::optional<Level56ScheduleSample> history_schedule_sample;
            std::size_t used_t0 = 0;

            auto run_batch = [&](Level56TemporalBatch<LLR>* batch,
                                 std::vector<Level56DispatchEntry>* prior,
                                 std::size_t budget,
                                 std::size_t slot_offset,
                                 bool preserve) {
              const auto& input5 = batch->tile_in5;
              const auto& input6 = batch->tile_in6;
              const auto& channel5 = batch->ch_tile5;
              const auto& channel6 = batch->ch_tile6;
              auto out5 = copy_tile_from_global(batch->tile_top5);
              auto out6 = copy_tile_from_global(batch->tile_top6);
              auto result = process_level56_shared(
                  input5, channel5, input6, channel6,
                  make_tile_params(kLevel5TileIndex),
                  make_tile_params(kLevel6TileIndex),
                  batch->tile_top5, batch->tile_top6, invocation,
                  normalize_extrinsic, tx_llr_ref, core_fn,
                  last_tile_history_accum, prior, budget, slot_offset, preserve,
                  &batch->early5, &batch->early6, &out5, &out6);
              batch->entries = result.dispatch;
              for (size_t r = 0; r < result.level5.tile_out.rows(); ++r) {
                for (size_t c = 0; c < N; ++c) {
                  work_llr[batch->tile_top5 + r][c] =
                      result.level5.tile_out[r][c];
                  work_llr[batch->tile_top6 + r][c] =
                      result.level6.tile_out[r][c];
                }
              }
              return result;
            };

            if (supplement_history) {
              auto previous_result = run_batch(
                  &*level56_temporal_state->previous,
                  &level56_temporal_state->previous->entries,
                  kLevel56MaxGroupEntries, 0, true);
              if (previous_result.has_schedule_sample) {
                history_schedule_sample =
                    std::move(previous_result.schedule_sample);
              }
              // The scheduler sample is the authoritative count after planning.
              // A later pass uses the number of assigned slots in the cached batch.
              for (const auto& entry : level56_temporal_state->previous->entries) {
                if (entry.assigned_entry_slot >= 0) {
                  used_t0 = std::max(
                      used_t0,
                      static_cast<std::size_t>(entry.assigned_entry_slot + 1));
                }
              }
              current_result = run_batch(
                  &*current_batch, nullptr,
                  kLevel56MaxGroupEntries - used_t0, used_t0, false);
            } else {
              current_result = run_batch(
                  &*current_batch, nullptr,
                  kLevel56MaxGroupEntries, 0, false);
            }
            const std::size_t t0_pending_after =
                have_previous
                    ? level56_temporal_pending_code_count(
                          level56_temporal_state->previous->entries)
                    : 0;
            const std::size_t t1_pending_after =
                level56_temporal_pending_code_count(current_result.dispatch);
            std::size_t t0_new_produced = 0;
            if (supplement_history) {
              const auto& history_entries =
                  level56_temporal_state->previous->entries;
              for (std::size_t index = 0; index < history_entries.size();
                   ++index) {
                t0_new_produced +=
                    t0_was_pending[index] && history_entries[index].produced;
              }
            }
            if (current_result.has_schedule_sample) {
              auto& sample = current_result.schedule_sample;
              const auto t1_group_entry_counts = sample.group_entry_counts;
              sample.temporal_lookahead_enabled = true;
              sample.temporal_has_history = have_previous;
              sample.temporal_has_future = lookahead_batch.has_value();
              sample.temporal_x = x;
              sample.temporal_k1 = k1;
              sample.temporal_k2 = k2;
              sample.temporal_branch = temporal_branch;
              sample.temporal_t0_group_entries = supplement_history ? used_t0 : 0;
              sample.temporal_t1_group_entries =
                  sample.total_group_entries - sample.temporal_t0_group_entries;
              sample.temporal_t0_pending_before = t0_pending_before;
              sample.temporal_t0_pending_after = t0_pending_after;
              sample.temporal_t1_pending_before = t1_pending_before;
              sample.temporal_t1_pending_after = t1_pending_after;
              sample.temporal_t0_new_produced = t0_new_produced;
              sample.temporal_t1_group_entry_counts = t1_group_entry_counts;

              for (auto& round : sample.rounds) {
                round.time_index = 1;
              }
              if (history_schedule_sample) {
                sample.temporal_t0_group_entry_counts =
                    history_schedule_sample->group_entry_counts;
                for (auto& round : history_schedule_sample->rounds) {
                  round.time_index = 0;
                }
                sample.rounds.insert(
                    sample.rounds.begin(),
                    std::make_move_iterator(
                        history_schedule_sample->rounds.begin()),
                    std::make_move_iterator(
                        history_schedule_sample->rounds.end()));
                for (std::size_t group = 0;
                     group < sample.group_entry_counts.size(); ++group) {
                  sample.group_entry_counts[group] +=
                      history_schedule_sample->group_entry_counts[group];
                }
                sample.planned_hiso_count +=
                    history_schedule_sample->planned_hiso_count;
                sample.planned_siso_count +=
                    history_schedule_sample->planned_siso_count;
              }

              std::vector<Level56ScheduleCodeSample> temporal_codes;
              temporal_codes.reserve(
                  (have_previous ? kLevel56GroupedCodeCount : 0u) +
                  kLevel56GroupedCodeCount +
                  (lookahead_batch ? kLevel56GroupedCodeCount : 0u));
              const auto append_temporal_codes =
                  [&](const std::vector<Level56DispatchEntry>& entries,
                      std::size_t time_index,
                      Level56TemporalInfoType info_type) {
                    for (const auto& entry : entries) {
                      temporal_codes.push_back(make_level56_schedule_code_sample(
                          entry, time_index, info_type, true));
                    }
                  };
              if (have_previous) {
                append_temporal_codes(
                    level56_temporal_state->previous->entries, 0,
                    Level56TemporalInfoType::DecodeInfo);
              }
              append_temporal_codes(
                  current_result.dispatch, 1,
                  Level56TemporalInfoType::EarlyStopInfo);
              if (lookahead_batch) {
                append_temporal_codes(
                    lookahead_batch->entries, 2,
                    Level56TemporalInfoType::EarlyStopInfo);
              }
              sample.codes = std::move(temporal_codes);
            }
            update_tile_stats(kLevel5TileIndex, current_result.level5);
            update_tile_stats(kLevel6TileIndex, current_result.level6);
            if (tile_stats && current_result.has_schedule_sample) {
              (*tile_stats)[kLevel5TileIndex].level56_schedule_samples.push_back(
                  std::move(current_result.schedule_sample));
            }
            level56_temporal_state->previous = std::move(*current_batch);
            level56_temporal_state->current.reset();
            ++t;
            continue;
          }
          auto shared = process_level56_shared(
              tile_in5, ch_tile5, tile_in6, ch_tile6,
              params5, params6, top5, top6, invocation,
              normalize_extrinsic, tx_llr_ref, core_fn,
              last_tile_history_accum);
          update_tile_stats(kLevel5TileIndex, shared.level5);
          update_tile_stats(kLevel6TileIndex, shared.level6);
          if (tile_stats && shared.has_schedule_sample) {
            (*tile_stats)[kLevel5TileIndex].level56_schedule_samples.push_back(
                std::move(shared.schedule_sample));
          }
          write_tile_to_work(kLevel5TileIndex, top5, shared.level5.tile_out);
          write_tile_to_work(kLevel6TileIndex, top6, shared.level6.tile_out);
          ++t;
          continue;
        }

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

        newcode::Params tile_params = make_tile_params(t);
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

    update_tile_stats(t, tile_result);
    write_tile_to_work(t, tile_top_row, tile_result.tile_out);
  }
}

} // namespace detail
} // namespace newcode
