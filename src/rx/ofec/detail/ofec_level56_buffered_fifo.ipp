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
  bool has_last_used_cycle = false;
  std::size_t last_used_cycle = 0;
  std::size_t release_cycle = 0;
};

struct Level56BufferedLatencyProtection {
  std::size_t source_batch_id = 0;
  std::size_t last_used_cycle = 0;
  std::size_t release_cycle = 0;
};

template <typename LLR>
struct Level56BufferedFifoState {
  std::deque<Level56BufferedBatch<LLR>> fifo;
  std::size_t next_batch_id = 0;
  std::size_t service_time = 0;
  std::size_t window_start = 0;
  std::vector<Level56BufferedLatencyProtection> latency_protections;
  bool initialized = false;
};

struct Level56BufferedServiceOutcome {
  struct OrdinaryService {
    std::size_t batch_id = 0;
    std::size_t entry_budget = 0;
    std::size_t entry_slot_offset = 0;
    std::size_t entries_used = 0;
    std::size_t first_used_cycle = 0;
    std::size_t last_used_cycle = 0;
    std::size_t release_cycle = 0;
    bool completed = false;
  };

  bool had_head_before = false;
  std::size_t head_batch_id_before = 0;
  bool ordinary_service_used = false;
  std::size_t interval_begin_cycle = 0;
  std::size_t interval_end_cycle = 0;
  std::size_t ordinary_entries_used = 0;
  std::size_t latency_blocked_cycles = 0;
  std::size_t effective_release_cycle = 0;
  std::size_t latency_candidate_batch_id = 0;
  bool had_latency_candidate = false;
  std::vector<std::size_t> blocking_source_batches;
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

inline std::size_t level56_buffered_effective_release_cycle(
    const std::vector<Level56BufferedLatencyProtection>& protections,
    std::size_t candidate_batch_id,
    std::vector<std::size_t>* blocking_sources) {
  std::size_t effective_release_cycle = 0;
  if (blocking_sources) {
    blocking_sources->clear();
  }
  for (const auto& protection : protections) {
    if (candidate_batch_id < protection.source_batch_id + 3u) {
      continue;
    }
    if (blocking_sources) {
      blocking_sources->push_back(protection.source_batch_id);
    }
    effective_release_cycle =
        std::max(effective_release_cycle, protection.release_cycle);
  }
  return effective_release_cycle;
}

inline void update_level56_buffered_latency_protection(
    std::vector<Level56BufferedLatencyProtection>* protections,
    std::size_t source_batch_id,
    std::size_t last_used_cycle,
    std::size_t release_cycle) {
  const auto found = std::find_if(
      protections->begin(), protections->end(),
      [source_batch_id](const auto& protection) {
        return protection.source_batch_id == source_batch_id;
      });
  if (found == protections->end()) {
    protections->push_back(Level56BufferedLatencyProtection{
        .source_batch_id = source_batch_id,
        .last_used_cycle = last_used_cycle,
        .release_cycle = release_cycle,
    });
    return;
  }
  found->last_used_cycle = last_used_cycle;
  found->release_cycle = release_cycle;
}

template <typename LLR, typename ClassifyFn, typename OrdinaryServiceFn,
          typename FullEarlyStopFn>
Level56BufferedServiceOutcome service_level56_buffered_fifo_time(
    Level56BufferedFifoState<LLR>* state,
    bool force_original_head_at_boundary,
    std::size_t ordinary_entry_capacity,
    std::size_t siso_decoder_latency,
    ClassifyFn classify,
    OrdinaryServiceFn ordinary_service,
    FullEarlyStopFn full_early_stop) {
  Level56BufferedServiceOutcome outcome;
  outcome.interval_begin_cycle =
      state->service_time * kLevel56MaxGroupEntries;
  outcome.interval_end_cycle =
      outcome.interval_begin_cycle + ordinary_entry_capacity;
  std::erase_if(state->latency_protections, [&](const auto& protection) {
    return protection.release_cycle <= outcome.interval_begin_cycle;
  });
  if (state->fifo.empty()) {
    return outcome;
  }

  outcome.had_head_before = true;
  outcome.head_batch_id_before = state->fifo.front().batch_id;
  std::size_t cursor = outcome.interval_begin_cycle;
  std::size_t service_index = 0;
  const auto retire_at = [&](std::size_t index,
                             Level56BufferedBatchRetirement reason) {
    const auto& batch = state->fifo[index];
    outcome.retirements.push_back(Level56BufferedRetirementSample{
        .batch_id = batch.batch_id,
        .arrival_time = batch.arrival_time,
        .reason = reason,
    });
    state->fifo.erase(state->fifo.begin() +
                      static_cast<std::ptrdiff_t>(index));
  };

  while (service_index < state->fifo.size()) {
    if (cursor >= outcome.interval_end_cycle) {
      break;
    }
    auto& batch = state->fifo[service_index];
    const bool was_classified = batch.classified;
    classify(&batch);

    if (level56_buffered_batch_all_early_stop(batch.entries)) {
      full_early_stop(&batch);
      if (!level56_buffered_batch_complete(batch.entries)) {
        throw std::logic_error(
            "LEVEL56 full EarlyStop action did not complete every code");
      }
      ++outcome.completed_batches;
      ++outcome.full_early_stop_batches;
      retire_at(service_index, Level56BufferedBatchRetirement::FullEarlyStop);
      continue;
    }

    std::vector<std::size_t> blocking_sources;
    const std::size_t effective_release_cycle =
        level56_buffered_effective_release_cycle(
            state->latency_protections, batch.batch_id, &blocking_sources);
    if (effective_release_cycle > cursor) {
      outcome.had_latency_candidate = true;
      outcome.latency_candidate_batch_id = batch.batch_id;
      outcome.effective_release_cycle = effective_release_cycle;
      outcome.blocking_source_batches = blocking_sources;
      const std::size_t next_cursor =
          std::min(effective_release_cycle, outcome.interval_end_cycle);
      outcome.latency_blocked_cycles += next_cursor - cursor;
      cursor = next_cursor;
    }
    if (cursor >= outcome.interval_end_cycle) {
      // Classification before a blocked ordinary service is speculative: the
      // shared SRAM may change before the next interval reaches this batch.
      if (!was_classified) {
        batch.classified = false;
        batch.entries.clear();
      }
      break;
    }

    const std::size_t entry_budget = outcome.interval_end_cycle - cursor;
    const std::size_t slot_offset = cursor - outcome.interval_begin_cycle;
    const std::size_t entries_used =
        ordinary_service(&batch, entry_budget, slot_offset);
    if (entries_used > entry_budget) {
      throw std::logic_error(
          "LEVEL56 ordinary service exceeded its remaining entry budget");
    }
    if (entries_used == 0u) {
      // A non-fast batch that makes no ordinary progress must stop this
      // interval; otherwise an unchanged batch could cause a busy loop.
      if (!was_classified) {
        batch.classified = false;
        batch.entries.clear();
      }
      break;
    }
    outcome.ordinary_service_used = true;
    const std::size_t first_used_cycle = cursor;
    const std::size_t last_used_cycle = cursor + entries_used - 1u;
    cursor += entries_used;
    batch.has_last_used_cycle = true;
    batch.last_used_cycle = last_used_cycle;
    batch.release_cycle = last_used_cycle + siso_decoder_latency;
    update_level56_buffered_latency_protection(
        &state->latency_protections, batch.batch_id,
        batch.last_used_cycle, batch.release_cycle);
    outcome.ordinary_entries_used += entries_used;

    const bool completed = level56_buffered_batch_complete(batch.entries);
    outcome.ordinary_services.push_back(
        Level56BufferedServiceOutcome::OrdinaryService{
            .batch_id = batch.batch_id,
            .entry_budget = entry_budget,
            .entry_slot_offset = slot_offset,
            .entries_used = entries_used,
            .first_used_cycle = first_used_cycle,
            .last_used_cycle = last_used_cycle,
            .release_cycle = batch.release_cycle,
            .completed = completed,
        });
    if (completed) {
      ++outcome.completed_batches;
      retire_at(service_index, Level56BufferedBatchRetirement::Normal);
    } else {
      // 本时刻临时向后服务；下一时刻仍会从 FIFO 最早未完成 batch 开始。
      ++service_index;
    }
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
    retire_at(0, Level56BufferedBatchRetirement::ForcedEvicted);
  }
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

}  // namespace detail
}  // namespace newcode
