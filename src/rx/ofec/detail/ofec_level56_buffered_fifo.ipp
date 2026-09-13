#pragma once

#include <algorithm>
#include <deque>
#include <limits>
#include <stdexcept>
#include <vector>

namespace newcode {
namespace detail {

template <typename LLR>
struct Level56BufferedBatch {
  std::size_t batch_id = 0;
  std::size_t arrival_time = 0;
  std::size_t tile_top5 = 0;
  std::size_t tile_top6 = 0;
  std::size_t window_start_at_arrival = 0;
  newcode::Params params5;
  newcode::Params params6;
  Level56EarlyStopEvaluation early5;
  Level56EarlyStopEvaluation early6;
  std::vector<Level56DispatchEntry> entries;
  bool classified = false;
};

template <typename LLR>
struct Level56BufferedFifoState {
  std::deque<Level56BufferedBatch<LLR>> fifo;
  std::size_t next_batch_id = 0;
  std::size_t service_time = 0;
  std::size_t window_start = 0;
  bool initialized = false;
};

// Split-buffer experiment: a persistent FIFO item belongs to exactly one
// level and therefore contains exactly 32 codes.  The 64-code object required
// by the shared Group4 scheduler is assembled only for one scheduling round
// and is split back into these two independent items immediately afterwards.
template <typename LLR>
struct Level56SplitBufferedBatch {
  std::size_t batch_id = 0;
  std::size_t arrival_time = 0;
  std::size_t source_level = 0;
  std::size_t tile_top = 0;
  std::size_t window_start_at_arrival = 0;
  newcode::Params params;
  Level56EarlyStopEvaluation early;
  std::vector<Level56DispatchEntry> entries;
  std::size_t ordinary_service_count = 0;
  bool classified = false;
};

template <typename LLR>
struct Level56SplitBufferedFifoState {
  std::deque<Level56SplitBufferedBatch<LLR>> fifo5;
  std::deque<Level56SplitBufferedBatch<LLR>> fifo6;
  std::size_t next_batch_id5 = 0;
  std::size_t next_batch_id6 = 0;
  std::size_t service_time = 0;
  // These are S5(t) and S6(t): the remaining distance to each level's
  // eviction boundary, in sub-block rows.
  std::size_t window_start5 = 0;
  std::size_t window_start6 = 0;
  bool initialized = false;
};

struct Level56SplitBufferedServiceOutcome {
  struct LevelOutcome {
    bool had_head_before = false;
    std::size_t head_batch_id_before = 0;
    std::size_t completed_batches = 0;
    std::size_t full_early_stop_batches = 0;
    std::size_t forced_evicted_batches = 0;
    std::vector<Level56BufferedRetirementSample> retirements;
    std::vector<std::size_t> forced_evicted_global_rows;
  };

  struct OrdinaryRound {
    std::size_t entry_budget = 0;
    std::size_t entry_slot_offset = 0;
    std::size_t entries_used = 0;
    bool served_level5 = false;
    bool served_level6 = false;
    std::size_t batch_id5 = std::numeric_limits<std::size_t>::max();
    std::size_t batch_id6 = std::numeric_limits<std::size_t>::max();
    bool completed_level5 = false;
    bool completed_level6 = false;
  };

  LevelOutcome level5;
  LevelOutcome level6;
  std::size_t ordinary_entries_used = 0;
  std::vector<OrdinaryRound> ordinary_rounds;
};

constexpr std::size_t kLevel56CodesPerLevel =
    kLevel56GroupedCodeCount / 2u;
constexpr std::size_t kLevel56SplitFixedDelaySubblockRows = 2u;

struct Level56BufferedServiceOutcome {
  struct OrdinaryService {
    std::size_t batch_id = 0;
    std::size_t entry_budget = 0;
    std::size_t entry_slot_offset = 0;
    std::size_t entries_used = 0;
    bool completed = false;
  };

  bool had_head_before = false;
  std::size_t head_batch_id_before = 0;
  bool ordinary_service_used = false;
  std::size_t ordinary_entries_used = 0;
  std::vector<OrdinaryService> ordinary_services;
  std::size_t completed_batches = 0;
  std::size_t full_early_stop_batches = 0;
  std::size_t forced_evicted_batches = 0;
  std::vector<Level56BufferedRetirementSample> retirements;
  std::vector<std::size_t> forced_evicted_global_rows;
};

inline bool level56_buffered_code_complete(
    const Level56DispatchEntry& entry) {
  return entry.already_decoded || entry.forced_evicted;
}

inline std::size_t level56_buffered_pending_count(
    const std::vector<Level56DispatchEntry>& entries) {
  return static_cast<std::size_t>(std::count_if(
      entries.begin(), entries.end(), [](const auto& entry) {
        return entry.pending && !level56_buffered_code_complete(entry);
      }));
}

inline bool level56_buffered_batch_complete(
    const std::vector<Level56DispatchEntry>& entries) {
  return entries.size() == kLevel56GroupedCodeCount &&
         std::all_of(entries.begin(), entries.end(), [](const auto& entry) {
           return entry.already_decoded && entry.writeback_complete &&
                  !entry.forced_evicted;
         });
}

inline bool level56_buffered_batch_all_early_stop(
    const std::vector<Level56DispatchEntry>& entries) {
  return entries.size() == kLevel56GroupedCodeCount &&
         std::all_of(entries.begin(), entries.end(), [](const auto& entry) {
           return entry.early_stop_hit;
         });
}

inline bool level56_split_buffered_batch_complete(
    const std::vector<Level56DispatchEntry>& entries) {
  return entries.size() == kLevel56CodesPerLevel &&
         std::all_of(entries.begin(), entries.end(), [](const auto& entry) {
           return entry.already_decoded && entry.writeback_complete &&
                  !entry.forced_evicted;
         });
}

inline bool level56_split_buffered_batch_all_early_stop(
    const std::vector<Level56DispatchEntry>& entries) {
  return entries.size() == kLevel56CodesPerLevel &&
         std::all_of(entries.begin(), entries.end(), [](const auto& entry) {
           return entry.early_stop_hit;
         });
}

inline std::size_t level56_buffered_next_window_start(
    std::size_t current,
    std::size_t completed_batches,
    std::size_t buffer_rows) {
  current = std::min(current, buffer_rows);
  const std::size_t distance_to_upper_clip = buffer_rows - current + 2u;
  const std::size_t completions_to_upper_clip =
      (distance_to_upper_clip + 1u) / 2u;
  if (completed_batches >= completions_to_upper_clip) {
    return buffer_rows;
  }
  const std::size_t before_clip = current + completed_batches * 2u;
  return before_clip > 2u ? before_clip - 2u : 0u;
}

inline void mark_level56_full_early_stop_complete(
    std::vector<Level56DispatchEntry>* entries) {
  for (auto& entry : *entries) {
    if (!entry.early_stop_hit) {
      throw std::logic_error(
          "LEVEL56 full EarlyStop retirement contains a non-EarlyStop code");
    }
    entry.planned_hiso = false;
    entry.planned_siso = false;
    entry.remaining_for_schedule = false;
    entry.final_action = Level56FinalAction::EarlyStopAction;
    entry.assigned_entry_slot = -1;
    entry.assigned_core = -1;
    entry.produced = true;
    entry.needs_execution = false;
    entry.already_decoded = true;
    entry.pending = false;
    entry.forced_evicted = false;
    entry.writeback_complete = true;
    entry.decode_status = Level56DecodeStatus::Produced;
  }
}

inline void mark_level56_batch_forced_evicted(
    std::vector<Level56DispatchEntry>* entries) {
  for (auto& entry : *entries) {
    if (entry.already_decoded) {
      continue;
    }
    entry.planned_hiso = false;
    entry.planned_siso = false;
    entry.remaining_for_schedule = false;
    entry.final_action = Level56FinalAction::Unscheduled;
    entry.assigned_entry_slot = -1;
    entry.assigned_core = -1;
    entry.produced = false;
    entry.needs_execution = false;
    entry.pending = false;
    entry.forced_evicted = true;
    entry.writeback_complete = false;
    entry.decode_status = Level56DecodeStatus::ForcedEvicted;
  }
}

template <typename LLR>
void initialize_level56_buffered_fifo_state(
    Level56BufferedFifoState<LLR>* state,
    const newcode::Params& p) {
  if (!state->initialized) {
    state->window_start = p.LEVEL56_BUFFER_ROWS;
    state->initialized = true;
  }
}

template <typename LLR>
void initialize_level56_split_buffered_fifo_state(
    Level56SplitBufferedFifoState<LLR>* state,
    const newcode::Params& p) {
  if (!state->initialized) {
    state->window_start5 =
        p.level5_buffer_rows() + kLevel56SplitFixedDelaySubblockRows;
    state->window_start6 = p.level6_buffer_rows();
    state->initialized = true;
  }
}

inline Level56DispatchEntry make_level56_split_bypass_entry(
    std::size_t source_level, std::size_t source_local_row) {
  Level56DispatchEntry entry;
  entry.shared_row = source_level == 5u
                         ? source_local_row
                         : kLevel56CodesPerLevel + source_local_row;
  entry.source_level = source_level;
  entry.source_local_row = source_local_row;
  entry.eligibility = Level56Eligibility::None;
  entry.final_action = Level56FinalAction::Unscheduled;
  entry.needs_execution = false;
  entry.already_decoded = true;
  entry.pending = false;
  entry.writeback_complete = true;
  entry.decode_status = Level56DecodeStatus::Produced;
  return entry;
}

inline std::vector<Level56DispatchEntry> make_level56_split_dispatch(
    const std::vector<Level56DispatchEntry>* entries5,
    const std::vector<Level56DispatchEntry>* entries6) {
  if ((entries5 && entries5->size() != kLevel56CodesPerLevel) ||
      (entries6 && entries6->size() != kLevel56CodesPerLevel)) {
    throw std::invalid_argument(
        "LEVEL56 split FIFO requires 32 cached entries per present level");
  }
  std::vector<Level56DispatchEntry> dispatch;
  dispatch.reserve(kLevel56GroupedCodeCount);
  for (std::size_t row = 0; row < kLevel56CodesPerLevel; ++row) {
    auto entry = entries5 ? (*entries5)[row]
                          : make_level56_split_bypass_entry(5u, row);
    entry.shared_row = row;
    entry.source_level = 5u;
    entry.source_local_row = row;
    dispatch.push_back(std::move(entry));
  }
  for (std::size_t row = 0; row < kLevel56CodesPerLevel; ++row) {
    auto entry = entries6 ? (*entries6)[row]
                          : make_level56_split_bypass_entry(6u, row);
    entry.shared_row = kLevel56CodesPerLevel + row;
    entry.source_level = 6u;
    entry.source_local_row = row;
    dispatch.push_back(std::move(entry));
  }
  return dispatch;
}

inline void split_level56_dispatch(
    const std::vector<Level56DispatchEntry>& dispatch,
    std::vector<Level56DispatchEntry>* entries5,
    std::vector<Level56DispatchEntry>* entries6) {
  if (dispatch.size() != kLevel56GroupedCodeCount) {
    throw std::invalid_argument(
        "LEVEL56 split FIFO scheduler returned a non-64-code dispatch");
  }
  if (entries5) {
    entries5->assign(dispatch.begin(),
                     dispatch.begin() +
                         static_cast<std::ptrdiff_t>(kLevel56CodesPerLevel));
    for (std::size_t row = 0; row < entries5->size(); ++row) {
      (*entries5)[row].shared_row = row;
      (*entries5)[row].source_level = 5u;
      (*entries5)[row].source_local_row = row;
    }
  }
  if (entries6) {
    entries6->assign(
        dispatch.begin() +
            static_cast<std::ptrdiff_t>(kLevel56CodesPerLevel),
        dispatch.end());
    for (std::size_t row = 0; row < entries6->size(); ++row) {
      (*entries6)[row].shared_row = kLevel56CodesPerLevel + row;
      (*entries6)[row].source_level = 6u;
      (*entries6)[row].source_local_row = row;
    }
  }
}

template <typename LLR, typename ClassifyFn, typename OrdinaryServiceFn,
          typename FullEarlyStopFn>
Level56BufferedServiceOutcome service_level56_buffered_fifo_time(
    Level56BufferedFifoState<LLR>* state,
    bool force_original_head_at_boundary,
    std::size_t ordinary_entry_capacity,
    ClassifyFn classify,
    OrdinaryServiceFn ordinary_service,
    FullEarlyStopFn full_early_stop) {
  Level56BufferedServiceOutcome outcome;
  if (state->fifo.empty()) {
    return outcome;
  }

  outcome.had_head_before = true;
  outcome.head_batch_id_before = state->fifo.front().batch_id;
  std::size_t remaining_entries = ordinary_entry_capacity;
  const auto retire_head = [&](Level56BufferedBatchRetirement reason) {
    const auto& head = state->fifo.front();
    outcome.retirements.push_back(Level56BufferedRetirementSample{
        .batch_id = head.batch_id,
        .arrival_time = head.arrival_time,
        .reason = reason,
    });
    state->fifo.pop_front();
  };

  while (!state->fifo.empty()) {
    auto& head = state->fifo.front();
    const bool was_classified = head.classified;
    classify(&head);

    if (level56_buffered_batch_all_early_stop(head.entries)) {
      full_early_stop(&head);
      if (!level56_buffered_batch_complete(head.entries)) {
        throw std::logic_error(
            "LEVEL56 full EarlyStop action did not complete every code");
      }
      ++outcome.completed_batches;
      ++outcome.full_early_stop_batches;
      retire_head(Level56BufferedBatchRetirement::FullEarlyStop);
      continue;
    }

    if (outcome.ordinary_services.size() >= 2u || remaining_entries == 0u) {
      // A newly exposed non-FullEarlyStop head may be inspected only to decide
      // whether the zero-cost fast path can continue. It has not received an
      // ordinary service, so do not cache this early classification across the
      // next window's preceding tile updates.
      if (!was_classified) {
        head.classified = false;
        head.entries.clear();
      }
      break;
    }
    outcome.ordinary_service_used = true;
    const std::size_t slot_offset =
        ordinary_entry_capacity - remaining_entries;
    const std::size_t entries_used =
        ordinary_service(&head, remaining_entries, slot_offset);
    if (entries_used == 0u || entries_used > remaining_entries) {
      throw std::logic_error(
          "LEVEL56 ordinary service returned an invalid entry count");
    }
    remaining_entries -= entries_used;
    outcome.ordinary_entries_used += entries_used;
    const bool completed = level56_buffered_batch_complete(head.entries);
    outcome.ordinary_services.push_back(
        Level56BufferedServiceOutcome::OrdinaryService{
            .batch_id = head.batch_id,
            .entry_budget = remaining_entries + entries_used,
            .entry_slot_offset = slot_offset,
            .entries_used = entries_used,
            .completed = completed,
        });
    if (!completed) {
      if (entries_used < remaining_entries + entries_used) {
        throw std::logic_error(
            "LEVEL56 unfinished ordinary batch did not consume its budget");
      }
      break;
    }
    ++outcome.completed_batches;
    retire_head(Level56BufferedBatchRetirement::Normal);
    // The completed head has already written back. The next loop iteration
    // classifies the new head from the updated shared-memory image. After two
    // ordinary batches, only a run of full-EarlyStop heads may still retire.
  }

  if (force_original_head_at_boundary && outcome.had_head_before &&
      !state->fifo.empty() &&
      state->fifo.front().batch_id == outcome.head_batch_id_before) {
    for (const auto& entry : state->fifo.front().entries) {
      if (!entry.already_decoded) {
        outcome.forced_evicted_global_rows.push_back(entry.source_global_row);
      }
    }
    mark_level56_batch_forced_evicted(&state->fifo.front().entries);
    ++outcome.forced_evicted_batches;
    retire_head(Level56BufferedBatchRetirement::ForcedEvicted);
  }
  return outcome;
}

template <typename LLR, typename ClassifyFn, typename OrdinaryServiceFn,
          typename FullEarlyStopFn>
Level56SplitBufferedServiceOutcome service_level56_split_buffered_fifo_time(
    Level56SplitBufferedFifoState<LLR>* state,
    bool force_level5_original_head_at_boundary,
    bool force_level6_original_head_at_boundary,
    std::size_t ordinary_entry_capacity,
    ClassifyFn classify,
    OrdinaryServiceFn ordinary_service,
    FullEarlyStopFn full_early_stop) {
  Level56SplitBufferedServiceOutcome outcome;
  if (!state->fifo5.empty()) {
    outcome.level5.had_head_before = true;
    outcome.level5.head_batch_id_before = state->fifo5.front().batch_id;
  }
  if (!state->fifo6.empty()) {
    outcome.level6.had_head_before = true;
    outcome.level6.head_batch_id_before = state->fifo6.front().batch_id;
  }

  const auto retire_head = [&](auto* fifo, auto* level_outcome,
                               Level56BufferedBatchRetirement reason) {
    const auto& head = fifo->front();
    level_outcome->retirements.push_back(Level56BufferedRetirementSample{
        .batch_id = head.batch_id,
        .arrival_time = head.arrival_time,
        .ordinary_service_count = head.ordinary_service_count,
        .reason = reason,
    });
    fifo->pop_front();
  };
  const auto clean_full_early_stop = [&](auto* fifo, auto* level_outcome) {
    while (!fifo->empty()) {
      auto& head = fifo->front();
      const bool was_classified = head.classified;
      classify(&head);
      if (!level56_split_buffered_batch_all_early_stop(head.entries)) {
        // Classifying a newly exposed ordinary head is permitted only for the
        // zero-cost check.  If no ordinary round subsequently consumes it,
        // the classification must not survive intervening SRAM updates.
        if (!was_classified) {
          head.classified = false;
          head.entries.clear();
        }
        break;
      }
      full_early_stop(&head);
      if (!level56_split_buffered_batch_complete(head.entries)) {
        throw std::logic_error(
            "LEVEL56 split FullEarlyStop did not complete all 32 codes");
      }
      ++level_outcome->completed_batches;
      ++level_outcome->full_early_stop_batches;
      retire_head(fifo, level_outcome,
                  Level56BufferedBatchRetirement::FullEarlyStop);
    }
  };

  std::size_t remaining_entries = ordinary_entry_capacity;
  for (std::size_t round = 0; round < 2u; ++round) {
    // The agreed order is stable: clean Level 5, then Level 6, before every
    // ordinary round.  This ordering affects logs only; neither path consumes
    // a Group4 entry.
    clean_full_early_stop(&state->fifo5, &outcome.level5);
    clean_full_early_stop(&state->fifo6, &outcome.level6);
    if (remaining_entries == 0u ||
        (state->fifo5.empty() && state->fifo6.empty())) {
      break;
    }

    auto* head5 = state->fifo5.empty() ? nullptr : &state->fifo5.front();
    auto* head6 = state->fifo6.empty() ? nullptr : &state->fifo6.front();
    if (head5) classify(head5);
    if (head6) classify(head6);
    const std::size_t slot_offset = ordinary_entry_capacity - remaining_entries;
    const std::size_t batch_id5 = head5
        ? head5->batch_id : std::numeric_limits<std::size_t>::max();
    const std::size_t batch_id6 = head6
        ? head6->batch_id : std::numeric_limits<std::size_t>::max();
    const std::size_t entries_used = ordinary_service(
        head5, head6, remaining_entries, slot_offset);
    if (entries_used == 0u || entries_used > remaining_entries) {
      throw std::logic_error(
          "LEVEL56 split ordinary round returned an invalid entry count");
    }
    remaining_entries -= entries_used;
    if (head5) ++head5->ordinary_service_count;
    if (head6) ++head6->ordinary_service_count;
    outcome.ordinary_entries_used += entries_used;

    const bool completed5 =
        head5 && level56_split_buffered_batch_complete(head5->entries);
    const bool completed6 =
        head6 && level56_split_buffered_batch_complete(head6->entries);
    outcome.ordinary_rounds.push_back(
        Level56SplitBufferedServiceOutcome::OrdinaryRound{
            .entry_budget = remaining_entries + entries_used,
            .entry_slot_offset = slot_offset,
            .entries_used = entries_used,
            .served_level5 = head5 != nullptr,
            .served_level6 = head6 != nullptr,
            .batch_id5 = batch_id5,
            .batch_id6 = batch_id6,
            .completed_level5 = completed5,
            .completed_level6 = completed6,
        });
    if (completed5) {
      ++outcome.level5.completed_batches;
      retire_head(&state->fifo5, &outcome.level5,
                  Level56BufferedBatchRetirement::Normal);
    }
    if (completed6) {
      ++outcome.level6.completed_batches;
      retire_head(&state->fifo6, &outcome.level6,
                  Level56BufferedBatchRetirement::Normal);
    }
    if (!completed5 && !completed6) {
      if (entries_used < remaining_entries + entries_used) {
        throw std::logic_error(
            "LEVEL56 split unfinished round did not consume its budget");
      }
      break;
    }
  }

  // A zero-cost run may still retire newly exposed FullEarlyStop heads after
  // the second ordinary round, exactly as in the previous two-batch scheme.
  clean_full_early_stop(&state->fifo5, &outcome.level5);
  clean_full_early_stop(&state->fifo6, &outcome.level6);

  const auto force_original_head = [&](auto* fifo, auto* level_outcome,
                                       bool force) {
    if (!force || !level_outcome->had_head_before || fifo->empty() ||
        fifo->front().batch_id != level_outcome->head_batch_id_before) {
      return;
    }
    // The original head must have been classified before it can reach this
    // point (ordinary service or the fast-path check).
    for (const auto& entry : fifo->front().entries) {
      if (!entry.already_decoded) {
        level_outcome->forced_evicted_global_rows.push_back(
            entry.source_global_row);
      }
    }
    mark_level56_batch_forced_evicted(&fifo->front().entries);
    ++level_outcome->forced_evicted_batches;
    retire_head(fifo, level_outcome,
                Level56BufferedBatchRetirement::ForcedEvicted);
  };
  force_original_head(&state->fifo5, &outcome.level5,
                      force_level5_original_head_at_boundary);
  force_original_head(&state->fifo6, &outcome.level6,
                      force_level6_original_head_at_boundary);
  return outcome;
}

template <typename LLR>
void classify_level56_buffered_batch(
    Level56BufferedBatch<LLR>* batch,
    const matrix::Matrix<LLR>& tile_in5,
    const matrix::Matrix<LLR>& ch_tile5,
    const matrix::Matrix<LLR>& tile_in6,
    const matrix::Matrix<LLR>& ch_tile6,
    std::size_t rows_to_decode,
    const matrix::Matrix<float>* tx_llr_ref) {
  if (batch->classified) {
    return;
  }
  auto prep5 = prepare_tile_inputs(
      tile_in5, ch_tile5, batch->params5, batch->tile_top5,
      batch->params5.CHASE_SBR, rows_to_decode, tx_llr_ref);
  auto prep6 = prepare_tile_inputs(
      tile_in6, ch_tile6, batch->params6, batch->tile_top6,
      batch->params6.CHASE_SBR, rows_to_decode, tx_llr_ref);
  batch->early5 = detect_level56_early_stop(prep5, batch->params5);
  batch->early6 = detect_level56_early_stop(prep6, batch->params6);
  batch->entries.clear();
  batch->entries.reserve(kLevel56GroupedCodeCount);
  append_level56_entries(prep5, batch->early5.effective, 5,
                         &batch->entries);
  append_level56_entries(prep6, batch->early6.effective, 6,
                         &batch->entries);
  if (batch->entries.size() != kLevel56GroupedCodeCount) {
    throw std::logic_error(
        "LEVEL56 buffered batch classification did not produce 64 codes");
  }
  for (auto& entry : batch->entries) {
    entry.needs_execution = false;
    entry.pending = false;
    entry.already_decoded = false;
    entry.forced_evicted = false;
    entry.writeback_complete = false;
  }
  batch->classified = true;
}

template <typename LLR>
void classify_level56_split_buffered_batch(
    Level56SplitBufferedBatch<LLR>* batch,
    const matrix::Matrix<LLR>& tile_in,
    const matrix::Matrix<LLR>& ch_tile,
    std::size_t rows_to_decode,
    const matrix::Matrix<float>* tx_llr_ref) {
  if (batch->classified) {
    return;
  }
  if (batch->source_level != 5u && batch->source_level != 6u) {
    throw std::invalid_argument(
        "LEVEL56 split FIFO batch source level must be 5 or 6");
  }
  auto prep = prepare_tile_inputs(
      tile_in, ch_tile, batch->params, batch->tile_top,
      batch->params.CHASE_SBR, rows_to_decode, tx_llr_ref);
  batch->early = detect_level56_early_stop(prep, batch->params);
  batch->entries.clear();
  batch->entries.reserve(kLevel56CodesPerLevel);
  append_level56_entries(prep, batch->early.effective, batch->source_level,
                         &batch->entries);
  if (batch->entries.size() != kLevel56CodesPerLevel) {
    throw std::logic_error(
        "LEVEL56 split FIFO classification did not produce 32 codes");
  }
  for (std::size_t row = 0; row < batch->entries.size(); ++row) {
    auto& entry = batch->entries[row];
    entry.shared_row = batch->source_level == 5u
                           ? row : kLevel56CodesPerLevel + row;
    entry.needs_execution = false;
    entry.pending = false;
    entry.already_decoded = false;
    entry.forced_evicted = false;
    entry.writeback_complete = false;
  }
  batch->classified = true;
}

}  // namespace detail
}  // namespace newcode
