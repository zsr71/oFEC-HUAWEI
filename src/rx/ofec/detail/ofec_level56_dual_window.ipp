#pragma once

#include <algorithm>
#include <iterator>
#include <limits>
#include <optional>
#include <stdexcept>
#include <utility>
#include <vector>

namespace newcode {
namespace detail {

struct Level56DualWindowSlot {
  bool valid = false;
  std::size_t batch_id = 0;
  std::size_t global_tile_top = 0;
};

struct Level56DualWindowStage {
  bool classified = false;
  bool terminal = false;
  Level56DualWindowCompletionReason reason =
      Level56DualWindowCompletionReason::None;
  Level56EarlyStopEvaluation early_stop;
  std::vector<Level56DispatchEntry> entries;
  // 每轮调度结束后的 entry-slot 快照；下一轮仍会重新参加 mux，历史只用于
  // 状态核查和观测，不参与资源锁定。
  std::vector<std::vector<int>> mux_history;
};

struct Level56DualWindowBatch {
  std::size_t batch_id = 0;
  std::size_t arrival_time = 0;
  std::size_t global_tile_top = 0;
  newcode::Params params5;
  newcode::Params params6;
  Level56DualWindowStage level5;
  Level56DualWindowStage level6;
};

struct Level56DualWindowState {
  bool initialized = false;
  std::size_t service_time = 0;
  std::size_t next_batch_id = 0;
  std::size_t s5 = 0;
  std::size_t s6 = 0;
  std::size_t recovery_debt5 = 0;
  std::size_t recovery_debt6 = 0;
  std::vector<Level56DualWindowSlot> slots;
  std::vector<Level56DualWindowBatch> batches;
};

inline void initialize_level56_dual_window_state(
    Level56DualWindowState* state,
    const newcode::Params& p) {
  if (!state || state->initialized) {
    return;
  }
  const std::size_t rows = p.LEVEL56_BUFFER_ROWS + 44u;
  if ((rows % 2u) != 0u || rows < 44u) {
    throw std::invalid_argument(
        "LEVEL56 dual-window SRAM length must be an even value >= 44");
  }
  state->slots.assign(rows / 2u, Level56DualWindowSlot{});
  state->s6 = p.LEVEL56_BUFFER_ROWS;
  state->s5 = p.LEVEL56_BUFFER_ROWS + 22u;
  state->initialized = true;
}

inline std::size_t level56_dual_decode_slot(std::size_t window_start) {
  return (window_start + 20u) / 2u;
}

inline bool level56_dual_windows_legal(
    const Level56DualWindowState& state,
    std::size_t s5,
    std::size_t s6) {
  const std::size_t rows = state.slots.size() * 2u;
  return (s5 % 2u) == 0u && (s6 % 2u) == 0u &&
         s6 + 22u <= s5 && s5 + 22u <= rows;
}

inline Level56DualWindowBatch* level56_dual_batch(
    Level56DualWindowState* state,
    std::size_t batch_id) {
  if (!state || batch_id >= state->batches.size() ||
      state->batches[batch_id].batch_id != batch_id) {
    throw std::logic_error("LEVEL56 dual-window batch mapping is invalid");
  }
  return &state->batches[batch_id];
}

inline const Level56DualWindowBatch* level56_dual_batch(
    const Level56DualWindowState& state,
    std::size_t batch_id) {
  if (batch_id >= state.batches.size() ||
      state.batches[batch_id].batch_id != batch_id) {
    throw std::logic_error("LEVEL56 dual-window batch mapping is invalid");
  }
  return &state.batches[batch_id];
}

inline std::optional<std::size_t> level56_dual_visible_batch_id(
    const Level56DualWindowState& state,
    std::size_t level) {
  const std::size_t start = level == 5u ? state.s5 : state.s6;
  const std::size_t slot_index = level56_dual_decode_slot(start);
  if (slot_index >= state.slots.size() || !state.slots[slot_index].valid) {
    return std::nullopt;
  }
  const std::size_t batch_id = state.slots[slot_index].batch_id;
  const auto* batch = level56_dual_batch(state, batch_id);
  if (level == 6u && !batch->level5.terminal) {
    return std::nullopt;
  }
  return batch_id;
}

inline Level56DualWindowStage* level56_dual_visible_stage(
    Level56DualWindowState* state,
    std::size_t level) {
  const auto id = level56_dual_visible_batch_id(*state, level);
  if (!id) {
    return nullptr;
  }
  auto* batch = level56_dual_batch(state, *id);
  return level == 5u ? &batch->level5 : &batch->level6;
}

inline bool level56_dual_code_complete(
    const Level56DispatchEntry& entry) {
  return entry.already_decoded && entry.writeback_complete &&
         !entry.forced_evicted;
}

inline bool level56_dual_stage_complete(
    const Level56DualWindowStage& stage) {
  return stage.classified && !stage.entries.empty() &&
         std::all_of(stage.entries.begin(), stage.entries.end(),
                     level56_dual_code_complete);
}

inline bool level56_dual_stage_full_early_stop(
    const Level56DualWindowStage& stage) {
  return stage.classified && !stage.entries.empty() &&
         std::all_of(stage.entries.begin(), stage.entries.end(),
                     [](const Level56DispatchEntry& entry) {
                       return entry.early_stop_hit;
                     });
}

inline void level56_dual_fill_step_position(
    const Level56DualWindowState& state,
    Level56DualWindowStepSample* step) {
  step->s5 = state.s5;
  step->s6 = state.s6;
  if (const auto h5 = level56_dual_visible_batch_id(state, 5u)) {
    step->h5_valid = true;
    step->h5_batch_id = *h5;
  }
  if (const auto h6 = level56_dual_visible_batch_id(state, 6u)) {
    step->h6_valid = true;
    step->h6_batch_id = *h6;
  }
  step->recovery_debt5 = state.recovery_debt5;
  step->recovery_debt6 = state.recovery_debt6;
}

inline void level56_dual_record_step(
    const Level56DualWindowState& state,
    Level56DualWindowStepType type,
    std::size_t step_index,
    std::vector<Level56DualWindowStepSample>* steps,
    bool level5_active = false,
    bool level6_active = false,
    std::size_t entries_before = 0,
    std::size_t entries_used = 0,
    Level56DualWindowCompletionReason reason5 =
        Level56DualWindowCompletionReason::None,
    Level56DualWindowCompletionReason reason6 =
        Level56DualWindowCompletionReason::None) {
  if (!steps) {
    return;
  }
  Level56DualWindowStepSample step;
  step.service_time = state.service_time;
  step.step_index = step_index;
  step.type = type;
  step.level5_active = level5_active;
  step.level6_active = level6_active;
  step.entries_before = entries_before;
  step.entries_used = entries_used;
  step.entries_after = entries_before + entries_used;
  step.level5_reason = reason5;
  step.level6_reason = reason6;
  level56_dual_fill_step_position(state, &step);
  steps->push_back(std::move(step));
}

template <typename ForcedFn = Level56NoopRetirement>
inline void level56_dual_force_stage(
    Level56DualWindowState* state,
    std::size_t level,
    Level56DualWindowTimeSample* time_sample,
    std::size_t* step_index,
    std::vector<Level56DualWindowStepSample>* steps,
    ForcedFn finalize_forced = {}) {
  auto* stage = level56_dual_visible_stage(state, level);
  if (!stage || stage->terminal) {
    return;
  }
  const auto batch_id = level56_dual_visible_batch_id(*state, level);
  const auto* batch = level56_dual_batch(state, *batch_id);
  const auto& params = level == 6u ? batch->params6 : batch->params5;
  for (auto& entry : stage->entries) {
    if (entry.already_decoded && entry.writeback_complete) {
      continue;
    }
    entry.forced_evicted = true;
    entry.pending = false;
    entry.remaining_for_schedule = false;
    entry.needs_execution = false;
    entry.final_action = Level56FinalAction::Unscheduled;
    entry.decode_status = Level56DecodeStatus::ForcedEvicted;
    finalize_forced(entry, params);
  }
  stage->terminal = true;
  stage->reason = Level56DualWindowCompletionReason::ForcedEvicted;
  if (level == 5u) {
    state->recovery_debt5 += 2u;
    ++time_sample->forced_evicted5;
  } else {
    state->recovery_debt6 += 2u;
    ++time_sample->forced_evicted6;
  }
  level56_dual_record_step(
      *state, Level56DualWindowStepType::ForcedEviction, (*step_index)++,
      steps, level == 5u, level == 6u, 0u, 0u,
      level == 5u ? Level56DualWindowCompletionReason::ForcedEvicted
                  : Level56DualWindowCompletionReason::None,
      level == 6u ? Level56DualWindowCompletionReason::ForcedEvicted
                  : Level56DualWindowCompletionReason::None);
}

inline Level56DualWindowSlot level56_dual_make_arrival(
    Level56DualWindowState* state,
    std::size_t global_tile_top,
    const newcode::Params& params5,
    const newcode::Params& params6) {
  Level56DualWindowBatch batch;
  batch.batch_id = state->next_batch_id++;
  batch.arrival_time = state->service_time;
  batch.global_tile_top = global_tile_top;
  batch.params5 = params5;
  batch.params6 = params6;
  state->batches.push_back(std::move(batch));
  return Level56DualWindowSlot{
      .valid = true,
      .batch_id = state->batches.back().batch_id,
      .global_tile_top = global_tile_top,
  };
}

template <typename ForcedFn = Level56NoopRetirement>
inline void level56_dual_advance_boundary(
    Level56DualWindowState* state,
    const Level56DualWindowSlot& incoming,
    Level56DualWindowTimeSample* time_sample,
    std::size_t* step_index,
    std::vector<Level56DualWindowStepSample>* steps,
    ForcedFn finalize_forced = {}) {
  if (state->service_time == 0u) {
    state->slots.back() = incoming;
    level56_dual_record_step(
        *state, Level56DualWindowStepType::BoundaryAdvance, (*step_index)++,
        steps);
    return;
  }

  const auto old_h5 = level56_dual_visible_batch_id(*state, 5u);
  const auto old_h6 = level56_dual_visible_batch_id(*state, 6u);
  bool follow5 = old_h5 &&
                 !level56_dual_batch(*state, *old_h5)->level5.terminal;
  bool follow6 = old_h6 &&
                 !level56_dual_batch(*state, *old_h6)->level6.terminal;

  if (follow6 && state->s6 < 2u) {
    level56_dual_force_stage(
        state, 6u, time_sample, step_index, steps, finalize_forced);
    follow6 = false;
  }

  std::size_t next_s5 = follow5 ? state->s5 - 2u : state->s5;
  std::size_t next_s6 = follow6 ? state->s6 - 2u : state->s6;
  if (next_s6 + 22u > next_s5) {
    if (!follow5) {
      throw std::logic_error(
          "LEVEL56 dual-window boundary cannot preserve non-overlap");
    }
    level56_dual_force_stage(
        state, 5u, time_sample, step_index, steps, finalize_forced);
    follow5 = false;
    next_s5 = state->s5;
  }
  if (!level56_dual_windows_legal(*state, next_s5, next_s6)) {
    throw std::logic_error(
        "LEVEL56 dual-window boundary produced an illegal window pair");
  }

  if (!state->slots.empty() && state->slots.front().valid) {
    time_sample->output_valid = true;
    time_sample->output_batch_id = state->slots.front().batch_id;
  }
  state->slots.erase(state->slots.begin());
  state->slots.push_back(incoming);
  state->s5 = next_s5;
  state->s6 = next_s6;

  if (follow5) {
    const auto h5 = level56_dual_visible_batch_id(*state, 5u);
    if (!h5 || h5 != old_h5) {
      throw std::logic_error(
          "LEVEL56 Level 5 necessary following lost its physical batch");
    }
  } else if (state->recovery_debt5 >= 2u &&
             level56_dual_visible_batch_id(*state, 5u)) {
    state->recovery_debt5 -= 2u;
  }
  if (follow6) {
    const auto h6 = level56_dual_visible_batch_id(*state, 6u);
    if (!h6 || h6 != old_h6) {
      throw std::logic_error(
          "LEVEL56 Level 6 necessary following lost its physical batch");
    }
  } else if (state->recovery_debt6 >= 2u &&
             level56_dual_visible_batch_id(*state, 6u)) {
    state->recovery_debt6 -= 2u;
  }
  level56_dual_record_step(
      *state, Level56DualWindowStepType::BoundaryAdvance, (*step_index)++,
      steps);
}

inline bool level56_dual_recovery_destination_valid(
    const Level56DualWindowState& state,
    std::size_t level,
    std::size_t start) {
  const std::size_t slot_index = level56_dual_decode_slot(start);
  if (slot_index >= state.slots.size() || !state.slots[slot_index].valid) {
    return false;
  }
  const auto* batch =
      level56_dual_batch(state, state.slots[slot_index].batch_id);
  return level != 6u || batch->level5.terminal;
}

inline bool level56_dual_apply_recovery(
    Level56DualWindowState* state,
    std::size_t* step_index,
    std::vector<Level56DualWindowStepSample>* steps) {
  bool moved = false;
  const auto move = [&](bool move5, bool move6) {
    const std::size_t candidate5 = state->s5 + (move5 ? 2u : 0u);
    const std::size_t candidate6 = state->s6 + (move6 ? 2u : 0u);
    if (!level56_dual_windows_legal(*state, candidate5, candidate6)) {
      return false;
    }
    if (move5 && !level56_dual_recovery_destination_valid(
                     *state, 5u, candidate5)) {
      return false;
    }
    if (move6 && !level56_dual_recovery_destination_valid(
                     *state, 6u, candidate6)) {
      return false;
    }
    state->s5 = candidate5;
    state->s6 = candidate6;
    if (move5) state->recovery_debt5 -= 2u;
    if (move6) state->recovery_debt6 -= 2u;
    level56_dual_record_step(
        *state, Level56DualWindowStepType::Recovery, (*step_index)++, steps,
        move5, move6);
    return true;
  };

  while (state->recovery_debt5 >= 2u &&
         state->recovery_debt6 >= 2u && move(true, true)) {
    moved = true;
  }
  while (state->recovery_debt5 >= 2u && move(true, false)) {
    moved = true;
  }
  while (state->recovery_debt6 >= 2u && move(false, true)) {
    moved = true;
  }
  return moved;
}

template <typename LLR>
matrix::Matrix<LLR> level56_dual_copy_tile(
    const matrix::Matrix<LLR>& source,
    std::size_t top,
    std::size_t rows) {
  matrix::Matrix<LLR> result(rows, source.cols());
  for (std::size_t row = 0; row < rows; ++row) {
    for (std::size_t col = 0; col < source.cols(); ++col) {
      result[row][col] = source[top + row][col];
    }
  }
  return result;
}

template <typename LLR>
void level56_dual_write_tile(
    matrix::Matrix<LLR>* destination,
    std::size_t top,
    const matrix::Matrix<LLR>& tile) {
  for (std::size_t row = 0; row < tile.rows(); ++row) {
    for (std::size_t col = 0; col < tile.cols(); ++col) {
      (*destination)[top + row][col] = tile[row][col];
    }
  }
}

template <typename LLR>
void classify_level56_dual_stage(
    Level56DualWindowBatch* batch,
    std::size_t level,
    const matrix::Matrix<LLR>& work_llr,
    const matrix::Matrix<LLR>& channel_llr,
    const matrix::Matrix<float>* tx_llr_ref) {
  auto* stage = level == 5u ? &batch->level5 : &batch->level6;
  if (stage->classified || stage->terminal) {
    return;
  }
  const auto& params = level == 5u ? batch->params5 : batch->params6;
  const std::size_t tile_rows = params.tile_height_rows();
  auto tile = level56_dual_copy_tile(
      work_llr, batch->global_tile_top, tile_rows);
  auto channel = level56_dual_copy_tile(
      channel_llr, batch->global_tile_top, tile_rows);
  const std::size_t rows_to_decode =
      static_cast<std::size_t>(params.CHASE_SBR) *
      newcode::Params::BITS_PER_SUBBLOCK_DIM;
  auto prep = prepare_tile_inputs(
      tile, channel, params, batch->global_tile_top, params.CHASE_SBR,
      rows_to_decode, tx_llr_ref);
  stage->early_stop = detect_level56_early_stop(prep, params);
  append_level56_entries(
      prep, stage->early_stop.effective, level, &stage->entries);
  for (std::size_t index = 0; index < stage->entries.size(); ++index) {
    stage->entries[index].shared_row = index;
  }
  stage->classified = true;
}

template <typename LLR>
std::size_t run_level56_dual_shared_round(
    matrix::Matrix<LLR>& work_llr,
    const matrix::Matrix<LLR>& channel_llr,
    bool level5_active,
    bool level6_active,
    std::size_t max_entries,
    std::size_t entry_offset,
    bool normalize_extrinsic,
    const matrix::Matrix<float>* tx_llr_ref,
    CoreFn<typename LinMatrixAdapter<LLR>::core_type> core_fn,
    matrix::Matrix<float>* last_tile_history_accum,
    std::size_t* level56_shared_invocation,
    Level56DualWindowState* state,
    std::vector<TileEarlyStopCounter>* tile_stats) {
  const auto h5 = level56_dual_visible_batch_id(*state, 5u);
  const auto h6 = level56_dual_visible_batch_id(*state, 6u);
  if ((level5_active && !h5) || (level6_active && !h6) ||
      (!level5_active && !level6_active)) {
    throw std::logic_error("LEVEL56 dual-window round has no active object");
  }
  auto* batch5 = h5 ? level56_dual_batch(state, *h5) : nullptr;
  auto* batch6 = h6 ? level56_dual_batch(state, *h6) : nullptr;
  auto* reference = batch5 ? batch5 : batch6;
  const std::size_t top5 = batch5 ? batch5->global_tile_top
                                  : reference->global_tile_top;
  const std::size_t top6 = batch6 ? batch6->global_tile_top
                                  : reference->global_tile_top;
  const auto& params5 = batch5 ? batch5->params5 : reference->params5;
  const auto& params6 = batch6 ? batch6->params6 : reference->params6;
  const std::size_t tile_rows = params5.tile_height_rows();
  auto input5 = level56_dual_copy_tile(work_llr, top5, tile_rows);
  auto input6 = level56_dual_copy_tile(work_llr, top6, tile_rows);
  auto channel5 = level56_dual_copy_tile(channel_llr, top5, tile_rows);
  auto channel6 = level56_dual_copy_tile(channel_llr, top6, tile_rows);
  auto output5 = input5;
  auto output6 = input6;

  std::vector<Level56DispatchEntry> combined;
  combined.reserve(kLevel56GroupedCodeCount);
  const auto append_stage = [&](const Level56DualWindowStage* stage,
                                std::size_t level) {
    if (stage) {
      for (auto entry : stage->entries) {
        entry.shared_row = combined.size();
        entry.source_level = level;
        combined.push_back(std::move(entry));
      }
      return;
    }
    for (std::size_t row = 0; row < kLevel56GroupedCodeCount / 2u; ++row) {
      Level56DispatchEntry entry;
      entry.shared_row = combined.size();
      entry.source_level = level;
      entry.source_local_row = row;
      entry.already_decoded = true;
      entry.needs_execution = false;
      entry.decode_status = Level56DecodeStatus::CompletedNoWriteback;
      combined.push_back(std::move(entry));
    }
  };
  append_stage(level5_active ? &batch5->level5 : nullptr, 5u);
  append_stage(level6_active ? &batch6->level6 : nullptr, 6u);

  const std::size_t invocation =
      level56_shared_invocation ? (*level56_shared_invocation)++ : 0u;
  auto result = process_level56_shared(
      input5, channel5, input6, channel6, params5, params6, top5, top6,
      invocation, normalize_extrinsic, tx_llr_ref, core_fn,
      last_tile_history_accum, &combined, max_entries, entry_offset, true,
      level5_active ? &batch5->level5.early_stop : nullptr,
      level6_active ? &batch6->level6.early_stop : nullptr,
      &output5, &output6, true,
      level5_active ? batch5->batch_id : batch6->batch_id,
      state->service_time, level5_active, level6_active);

  if (level5_active) {
    batch5->level5.entries.assign(
        result.dispatch.begin(), result.dispatch.begin() + 32);
    std::vector<int> slots;
    slots.reserve(batch5->level5.entries.size());
    for (std::size_t index = 0; index < batch5->level5.entries.size(); ++index) {
      batch5->level5.entries[index].shared_row = index;
      slots.push_back(batch5->level5.entries[index].assigned_entry_slot);
    }
    batch5->level5.mux_history.push_back(std::move(slots));
    level56_dual_write_tile(&work_llr, top5, result.level5.tile_out);
  }
  if (level6_active) {
    batch6->level6.entries.assign(
        result.dispatch.begin() + 32, result.dispatch.end());
    std::vector<int> slots;
    slots.reserve(batch6->level6.entries.size());
    for (std::size_t index = 0; index < batch6->level6.entries.size(); ++index) {
      batch6->level6.entries[index].shared_row = index;
      slots.push_back(batch6->level6.entries[index].assigned_entry_slot);
    }
    batch6->level6.mux_history.push_back(std::move(slots));
    level56_dual_write_tile(&work_llr, top6, result.level6.tile_out);
  }

  const auto update_stats = [&](std::size_t tile_index,
                                const TileProcessResult<LLR>& tile_result) {
    if (!tile_stats || tile_index >= tile_stats->size()) return;
    auto& counter = (*tile_stats)[tile_index];
    ++counter.total;
    counter.triggered += tile_result.early_stop_triggered;
    counter.row_total += tile_result.rows_total;
    counter.row_triggered += tile_result.rows_early_stop;
    counter.row_hard_finish += tile_result.rows_hard_finish;
    counter.row_need_siso_before_mux += tile_result.rows_need_siso_before_mux;
    counter.row_unscheduled += tile_result.rows_unscheduled;
    counter.hybrid_class_counts.push_back(tile_result.hybrid_class_count);
  };
  if (level5_active) update_stats(kLevel5TileIndex, result.level5);
  if (level6_active) update_stats(kLevel6TileIndex, result.level6);
  if (tile_stats && result.has_schedule_sample) {
    (*tile_stats)[kLevel5TileIndex].level56_schedule_samples.push_back(
        std::move(result.schedule_sample));
  }

  std::size_t used = 0;
  for (const auto& entry : result.dispatch) {
    if (entry.assigned_entry_slot >= static_cast<int>(entry_offset)) {
      used = std::max(
          used, static_cast<std::size_t>(entry.assigned_entry_slot) + 1u -
                    entry_offset);
    }
  }
  if (used > max_entries) {
    throw std::logic_error("LEVEL56 dual-window round exceeded entry budget");
  }
  return used;
}

template <typename LLR>
void service_level56_dual_window_time(
    matrix::Matrix<LLR>& work_llr,
    const matrix::Matrix<LLR>& channel_llr,
    const newcode::Params& p,
    const std::optional<std::size_t>& arriving_global_tile_top,
    const newcode::Params& params5,
    const newcode::Params& params6,
    std::vector<TileEarlyStopCounter>* tile_stats,
    bool normalize_extrinsic,
    const matrix::Matrix<float>* tx_llr_ref,
    CoreFn<typename LinMatrixAdapter<LLR>::core_type> core_fn,
    matrix::Matrix<float>* last_tile_history_accum,
    std::size_t* level56_shared_invocation,
    Level56DualWindowState* state,
    bool drain) {
  initialize_level56_dual_window_state(state, p);
  Level56DualWindowTimeSample time_sample;
  time_sample.service_time = state->service_time;
  time_sample.drain = drain;
  time_sample.s5_before = state->s5;
  time_sample.s6_before = state->s6;
  std::vector<Level56DualWindowStepSample> steps;
  std::size_t step_index = 0;

  Level56DualWindowSlot incoming;
  if (arriving_global_tile_top) {
    incoming = level56_dual_make_arrival(
        state, *arriving_global_tile_top, params5, params6);
    time_sample.input_valid = true;
    time_sample.input_batch_id = incoming.batch_id;
  }
  level56_dual_advance_boundary(
      state, incoming, &time_sample, &step_index, &steps,
      [&](const Level56DispatchEntry& entry, const newcode::Params& params) {
        finalize_level56_forced_prior(
            work_llr, last_tile_history_accum, entry, params);
      });

  std::size_t entries_used = 0;
  std::size_t guard = 0;
  while (++guard <= state->slots.size() * 4u + 16u) {
    bool classified_any = false;
    const auto h5 = level56_dual_visible_batch_id(*state, 5u);
    const auto h6 = level56_dual_visible_batch_id(*state, 6u);
    if (h5) {
      auto* batch = level56_dual_batch(state, *h5);
      if (!batch->level5.classified && !batch->level5.terminal) {
        classify_level56_dual_stage(
            batch, 5u, work_llr, channel_llr, tx_llr_ref);
        classified_any = true;
      }
    }
    if (h6) {
      auto* batch = level56_dual_batch(state, *h6);
      if (!batch->level6.classified && !batch->level6.terminal) {
        classify_level56_dual_stage(
            batch, 6u, work_llr, channel_llr, tx_llr_ref);
        classified_any = true;
      }
    }
    if (classified_any) {
      level56_dual_record_step(
          *state, Level56DualWindowStepType::Classify, step_index++, &steps,
          h5.has_value(), h6.has_value(), entries_used);
    }

    auto* stage5 = level56_dual_visible_stage(state, 5u);
    auto* stage6 = level56_dual_visible_stage(state, 6u);
    const bool full5 = stage5 && !stage5->terminal &&
                       level56_dual_stage_full_early_stop(*stage5);
    const bool full6 = stage6 && !stage6->terminal &&
                       level56_dual_stage_full_early_stop(*stage6);
    if (full5 || full6) {
      (void)run_level56_dual_shared_round(
          work_llr, channel_llr, full5, full6, 0u, entries_used,
          normalize_extrinsic, tx_llr_ref, core_fn, last_tile_history_accum,
          level56_shared_invocation, state, tile_stats);
      Level56DualWindowCompletionReason reason5 =
          Level56DualWindowCompletionReason::None;
      Level56DualWindowCompletionReason reason6 =
          Level56DualWindowCompletionReason::None;
      if (full5) {
        stage5->terminal = true;
        stage5->reason = Level56DualWindowCompletionReason::FullEarlyStop;
        state->recovery_debt5 += 2u;
        ++time_sample.c5;
        ++time_sample.full_early_stop5;
        reason5 = Level56DualWindowCompletionReason::FullEarlyStop;
      }
      if (full6) {
        stage6->terminal = true;
        stage6->reason = Level56DualWindowCompletionReason::FullEarlyStop;
        state->recovery_debt6 += 2u;
        ++time_sample.c6;
        ++time_sample.full_early_stop6;
        reason6 = Level56DualWindowCompletionReason::FullEarlyStop;
      }
      level56_dual_record_step(
          *state, Level56DualWindowStepType::FullEarlyStop, step_index++,
          &steps, full5, full6, entries_used, 0u, reason5, reason6);
      (void)level56_dual_apply_recovery(state, &step_index, &steps);
      continue;
    }

    if (entries_used >= kLevel56MaxGroupEntries) {
      break;
    }
    stage5 = level56_dual_visible_stage(state, 5u);
    stage6 = level56_dual_visible_stage(state, 6u);
    const bool active5 = stage5 && stage5->classified && !stage5->terminal;
    const bool active6 = stage6 && stage6->classified && !stage6->terminal;
    if (!active5 && !active6) {
      break;
    }

    const std::size_t remaining = kLevel56MaxGroupEntries - entries_used;
    const std::size_t used = run_level56_dual_shared_round(
        work_llr, channel_llr, active5, active6, remaining, entries_used,
        normalize_extrinsic, tx_llr_ref, core_fn, last_tile_history_accum,
        level56_shared_invocation, state, tile_stats);
    ++time_sample.ordinary_rounds;
    entries_used += used;
    Level56DualWindowCompletionReason reason5 =
        Level56DualWindowCompletionReason::None;
    Level56DualWindowCompletionReason reason6 =
        Level56DualWindowCompletionReason::None;
    stage5 = level56_dual_visible_stage(state, 5u);
    stage6 = level56_dual_visible_stage(state, 6u);
    if (active5 && stage5 && !stage5->terminal &&
        level56_dual_stage_complete(*stage5)) {
      stage5->terminal = true;
      stage5->reason = Level56DualWindowCompletionReason::Normal;
      state->recovery_debt5 += 2u;
      ++time_sample.c5;
      reason5 = Level56DualWindowCompletionReason::Normal;
    }
    if (active6 && stage6 && !stage6->terminal &&
        level56_dual_stage_complete(*stage6)) {
      stage6->terminal = true;
      stage6->reason = Level56DualWindowCompletionReason::Normal;
      state->recovery_debt6 += 2u;
      ++time_sample.c6;
      reason6 = Level56DualWindowCompletionReason::Normal;
    }
    level56_dual_record_step(
        *state, Level56DualWindowStepType::OrdinarySchedule, step_index++,
        &steps, active5, active6, entries_used - used, used, reason5, reason6);
    const bool completed = reason5 != Level56DualWindowCompletionReason::None ||
                           reason6 != Level56DualWindowCompletionReason::None;
    (void)level56_dual_apply_recovery(state, &step_index, &steps);
    if (used == 0u && !completed) {
      break;
    }
  }
  if (guard > state->slots.size() * 4u + 16u) {
    throw std::logic_error("LEVEL56 dual-window microstep loop did not converge");
  }

  time_sample.ordinary_entries_used = entries_used;
  time_sample.s5_after = state->s5;
  time_sample.s6_after = state->s6;
  if (const auto h5 = level56_dual_visible_batch_id(*state, 5u)) {
    time_sample.h5_valid = true;
    time_sample.h5_batch_id = *h5;
  }
  if (const auto h6 = level56_dual_visible_batch_id(*state, 6u)) {
    time_sample.h6_valid = true;
    time_sample.h6_batch_id = *h6;
  }
  time_sample.recovery_debt5 = state->recovery_debt5;
  time_sample.recovery_debt6 = state->recovery_debt6;
  if (tile_stats && tile_stats->size() > kLevel5TileIndex) {
    auto& counter = (*tile_stats)[kLevel5TileIndex];
    counter.level56_dual_window_time_samples.push_back(std::move(time_sample));
    counter.level56_dual_window_step_samples.insert(
        counter.level56_dual_window_step_samples.end(),
        std::make_move_iterator(steps.begin()),
        std::make_move_iterator(steps.end()));
  }
  ++state->service_time;
}

inline bool level56_dual_has_real_slots(
    const Level56DualWindowState& state) {
  return std::any_of(state.slots.begin(), state.slots.end(),
                     [](const Level56DualWindowSlot& slot) {
                       return slot.valid;
                     });
}

inline bool level56_dual_all_real_stages_terminal(
    const Level56DualWindowState& state) {
  return std::all_of(
      state.batches.begin(), state.batches.end(),
      [](const Level56DualWindowBatch& batch) {
        return batch.level5.terminal && batch.level6.terminal;
      });
}

template <typename LLR>
void drain_level56_dual_window_at_frame_end(
    matrix::Matrix<LLR>& work_llr,
    const matrix::Matrix<LLR>& channel_llr,
    const newcode::Params& p,
    std::vector<TileEarlyStopCounter>* tile_stats,
    bool normalize_extrinsic,
    const matrix::Matrix<float>* tx_llr_ref,
    CoreFn<typename LinMatrixAdapter<LLR>::core_type> core_fn,
    matrix::Matrix<float>* last_tile_history_accum,
    std::size_t* level56_shared_invocation,
    Level56DualWindowState* state) {
  if (!state || !state->initialized) {
    return;
  }
  newcode::Params params5 = p;
  newcode::Params params6 = p;
  std::size_t guard = 0;
  const std::size_t max_times = state->slots.size() + 2u;
  while (level56_dual_has_real_slots(*state)) {
    if (++guard > max_times) {
      throw std::logic_error("LEVEL56 dual-window drain did not terminate");
    }
    service_level56_dual_window_time(
        work_llr, channel_llr, p, std::nullopt, params5, params6, tile_stats,
        normalize_extrinsic, tx_llr_ref, core_fn, last_tile_history_accum,
        level56_shared_invocation, state, true);
  }
  if (!level56_dual_all_real_stages_terminal(*state)) {
    throw std::logic_error(
        "LEVEL56 dual-window drain emptied SRAM with unfinished stages");
  }
}

}  // namespace detail
}  // namespace newcode
