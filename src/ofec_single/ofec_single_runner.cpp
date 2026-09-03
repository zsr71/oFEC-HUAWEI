#include "newcode/ofec_single_runner.hpp"
#include "newcode/io/dualwriter.hpp"
#include "newcode/io/prepare_log_file.hpp"

#include <array>
#include <chrono>
#include <ctime>
#include <iomanip>
#include <map>
#include <sstream>

namespace {

std::string make_run_id() {
  const auto now = std::chrono::system_clock::now();
  const std::time_t now_time = std::chrono::system_clock::to_time_t(now);
  std::tm tm{};
#if defined(_WIN32)
  localtime_s(&tm, &now_time);
#else
  localtime_r(&now_time, &tm);
#endif
  std::ostringstream oss;
  oss << std::put_time(&tm, "%Y%m%d_%H%M%S");
  return oss.str();
}

template <typename Range>
std::string join_numbers(const Range& values, char separator = '|') {
  std::ostringstream oss;
  bool first = true;
  for (const auto value : values) {
    if (!first) {
      oss << separator;
    }
    first = false;
    oss << value;
  }
  return oss.str();
}

template <typename Range>
std::string join_one_based(const Range& values, char separator = '|') {
  std::ostringstream oss;
  bool first = true;
  for (const auto value : values) {
    if (!first) {
      oss << separator;
    }
    first = false;
    oss << (value + 1u);
  }
  return oss.str();
}

const char* level56_branch_name(newcode::Level56ScheduleBranch branch) {
  switch (branch) {
    case newcode::Level56ScheduleBranch::K0:
      return "K0";
    case newcode::Level56ScheduleBranch::KLessThan8:
      return "K_LT_8";
    case newcode::Level56ScheduleBranch::KEqual8:
      return "K_EQ_8";
    case newcode::Level56ScheduleBranch::KGreaterThan8:
      return "K_GT_8";
  }
  return "UNKNOWN";
}

const char* level56_temporal_branch_name(
    newcode::Level56TemporalBranch branch) {
  switch (branch) {
    case newcode::Level56TemporalBranch::Disabled: return "Disabled";
    case newcode::Level56TemporalBranch::NoHistory: return "NoHistory";
    case newcode::Level56TemporalBranch::SupplementHistory:
      return "SupplementHistory";
    case newcode::Level56TemporalBranch::CurrentFirst: return "CurrentFirst";
    case newcode::Level56TemporalBranch::NoFuture: return "NoFuture";
  }
  return "Unknown";
}

const char* level56_info_type_name(newcode::Level56TemporalInfoType type) {
  switch (type) {
    case newcode::Level56TemporalInfoType::DecodeInfo: return "DecodeInfo";
    case newcode::Level56TemporalInfoType::EarlyStopInfo:
      return "EarlyStopInfo";
  }
  return "Unknown";
}

const char* level56_decode_status_name(newcode::Level56DecodeStatus status) {
  switch (status) {
    case newcode::Level56DecodeStatus::NotDecoded: return "NotDecoded";
    case newcode::Level56DecodeStatus::Produced: return "Produced";
    case newcode::Level56DecodeStatus::ActionFailed: return "ActionFailed";
    case newcode::Level56DecodeStatus::CompletedNoWriteback:
      return "CompletedNoWriteback";
    case newcode::Level56DecodeStatus::ForcedEvicted: return "ForcedEvicted";
  }
  return "Unknown";
}

const char* level56_buffered_retirement_name(
    newcode::Level56BufferedBatchRetirement reason) {
  switch (reason) {
    case newcode::Level56BufferedBatchRetirement::Normal: return "Normal";
    case newcode::Level56BufferedBatchRetirement::FullEarlyStop:
      return "FullEarlyStop";
    case newcode::Level56BufferedBatchRetirement::ForcedEvicted:
      return "ForcedEvicted";
  }
  return "Unknown";
}

const char* level56_hybrid_class_name(uint8_t value) {
  switch (value) {
    case 0: return "None";
    case 1: return "BchHardDecoded";
    case 2: return "Clean";
    case 3: return "ParityOnly";
    case 4: return "OneMain";
    case 5: return "OneMainPlusParity";
    case 6: return "TwoMain";
    case 7: return "Suspicious";
    case 8: return "HardFail";
    default: return "Unknown";
  }
}

const char* level56_eligibility_name(uint8_t value) {
  switch (value) {
    case 0: return "None";
    case 1: return "HisoOnly";
    case 2: return "SisoOnly";
    case 3: return "HisoOrSiso";
    default: return "Unknown";
  }
}

const char* level56_action_name(uint8_t value) {
  switch (value) {
    case 0: return "EarlyStopAction";
    case 1: return "HisoDecode";
    case 2: return "SisoDecode";
    case 3: return "Unscheduled";
    default: return "Unknown";
  }
}

void dump_tile_early_stop_samples_csv(
    const std::vector<newcode::TileEarlyStopSample>& samples,
    const std::filesystem::path& output_path) {
  std::filesystem::path parent = output_path.parent_path();
  if (!parent.empty()) {
    std::filesystem::create_directories(parent);
  }
  std::ofstream out(output_path);
  out << "invocation,tile_index,rows_total,rows_passed,"
         "rows_hard_finish,rows_need_siso_before_mux,rows_unscheduled\n";
  for (const auto& sample : samples) {
    out << sample.invocation << ','
        << sample.tile_index << ','
        << sample.rows_total << ','
        << sample.rows_passed << ','
        << sample.rows_hard_finish << ','
        << sample.rows_need_siso_before_mux << ','
        << sample.rows_unscheduled << '\n';
  }
}

void dump_tile_early_stop_group_bind_debug_samples_csv(
    const std::vector<newcode::TileEarlyStopGroupBindDebugSample>& samples,
    const std::filesystem::path& output_path,
    const std::string& run_id,
    const std::string& label) {
  std::filesystem::path parent = output_path.parent_path();
  if (!parent.empty()) {
    std::filesystem::create_directories(parent);
  }
  std::ofstream out(output_path);
  out << "run_id,label,invocation,tile_index,condition_mode,bind_group_size,"
         "raw_early_stop_flags,bound_early_stop_flags\n";
  for (const auto& sample : samples) {
    out << run_id << ','
        << label << ','
        << sample.invocation << ','
        << sample.tile_index << ','
        << sample.condition_mode << ','
        << sample.bind_group_size << ','
        << sample.raw_early_stop_flags << ','
        << sample.bound_early_stop_flags << '\n';
  }
}

void dump_level56_schedule_rounds_csv(
    const std::vector<newcode::Level56ScheduleSample>& samples,
    const std::filesystem::path& output_path,
    const std::string& run_id,
    const std::string& label) {
  const auto parent = output_path.parent_path();
  if (!parent.empty()) {
    std::filesystem::create_directories(parent);
  }
  std::ofstream out(output_path);
  out << "run_id,label,invocation,buffered_fifo_enabled,buffered_t,"
         "buffered_batch,temporal_enabled,has_history,has_future,"
         "X,K1,K2,temporal_branch,t0_group_entries,t1_group_entries,"
         "t0_pending_before,t0_pending_after,t0_pending_reduced,"
         "t1_pending_before,t1_pending_after,t1_pending_reduced,"
         "pending_reduced_total,t0_new_produced,"
         "K,branch,time_index,round_index,used_entries_before,"
         "used_entries_after,initial_counts,remaining_before,selected_groups_1based,"
         "remaining_after,group_entry_counts,t0_group_entry_counts,"
         "t1_group_entry_counts,entry_slot_offset,group_entries_used,"
         "total_group_entries,planned_hiso,"
         "planned_siso,idle_hiso_capacity,idle_siso_capacity\n";
  for (const auto& sample : samples) {
    const auto write_row = [&](long round_index,
                               std::size_t time_index,
                               std::size_t used_before,
                               std::size_t used_after,
                               const auto& initial_counts,
                               std::size_t initial_nonzero_groups,
                               newcode::Level56ScheduleBranch branch,
                               const auto& remaining_before,
                               const auto& selected_groups,
                               const auto& remaining_after,
                               const auto& group_entry_counts) {
      out << run_id << ',' << label << ',' << sample.invocation << ','
          << sample.buffered_fifo_enabled << ',';
      if (sample.buffered_fifo_enabled) {
        out << sample.buffered_service_time << ",B"
            << sample.buffered_batch_id;
      } else {
        out << ',';
      }
      out << ',' << sample.temporal_lookahead_enabled << ','
          << sample.temporal_has_history << ','
          << sample.temporal_has_future << ','
          << sample.temporal_x << ',' << sample.temporal_k1 << ','
          << sample.temporal_k2 << ','
          << level56_temporal_branch_name(sample.temporal_branch) << ','
          << sample.temporal_t0_group_entries << ','
          << sample.temporal_t1_group_entries << ','
          << sample.temporal_t0_pending_before << ','
          << sample.temporal_t0_pending_after << ','
          << (sample.temporal_t0_pending_before -
              sample.temporal_t0_pending_after) << ','
          << sample.temporal_t1_pending_before << ','
          << sample.temporal_t1_pending_after << ','
          << (sample.temporal_t1_pending_before -
              sample.temporal_t1_pending_after) << ','
          << (sample.temporal_t0_pending_before -
              sample.temporal_t0_pending_after +
              sample.temporal_t1_pending_before -
              sample.temporal_t1_pending_after) << ','
          << sample.temporal_t0_new_produced << ','
          << initial_nonzero_groups << ','
          << level56_branch_name(branch) << ','
          << time_index << ',' << round_index << ','
          << used_before << ',' << used_after << ','
          << '"' << join_numbers(initial_counts) << "\",\""
          << join_numbers(remaining_before) << "\",\""
          << join_one_based(selected_groups) << "\",\""
          << join_numbers(remaining_after) << "\",\""
          << join_numbers(group_entry_counts) << "\",\""
          << join_numbers(sample.temporal_t0_group_entry_counts) << "\",\""
          << join_numbers(sample.temporal_t1_group_entry_counts) << "\","
          << sample.entry_slot_offset << ','
          << sample.group_entries_used << ','
          << sample.total_group_entries << ','
          << sample.planned_hiso_count << ','
          << sample.planned_siso_count << ','
          << (sample.entry_capacity >= sample.planned_hiso_count
                  ? sample.entry_capacity - sample.planned_hiso_count
                  : 0u) << ','
          << (sample.entry_capacity >= sample.planned_siso_count
                  ? sample.entry_capacity - sample.planned_siso_count
                  : 0u) << '\n';
    };
    if (sample.rounds.empty()) {
      const std::vector<std::size_t> no_groups;
      write_row(-1, 1, 0, 0, sample.initial_counts,
                sample.initial_nonzero_groups, sample.branch,
                sample.initial_counts, no_groups,
                sample.initial_counts, sample.group_entry_counts);
      continue;
    }
    for (std::size_t index = 0; index < sample.rounds.size(); ++index) {
      const auto& round = sample.rounds[index];
      const auto& round_group_entry_counts =
          sample.temporal_lookahead_enabled
              ? (round.time_index == 0
                     ? sample.temporal_t0_group_entry_counts
                     : sample.temporal_t1_group_entry_counts)
              : sample.group_entry_counts;
      write_row(static_cast<long>(index + 1u), round.time_index,
                round.used_entries_before, round.used_entries_after,
                round.initial_counts, round.initial_nonzero_groups,
                round.branch,
                round.remaining_before, round.selected_groups,
                round.remaining_after, round_group_entry_counts);
    }
  }
}

void dump_level56_schedule_codes_csv(
    const std::vector<newcode::Level56ScheduleSample>& samples,
    const std::filesystem::path& output_path,
    const std::string& run_id,
    const std::string& label) {
  const auto parent = output_path.parent_path();
  if (!parent.empty()) {
    std::filesystem::create_directories(parent);
  }
  std::ofstream out(output_path);
  out << "run_id,label,invocation,buffered_fifo_enabled,buffered_t,"
         "buffered_batch,shared_index,time_index,time_offset,info_type,"
         "decode_status,code,source_level,source_local_row,source_global_row,group,"
         "position_in_group,early_stop_hit,hybrid_class,resource_eligibility,"
         "planned_hiso,planned_siso,remaining_for_schedule,final_action,"
         "assigned_entry_slot,assigned_core,produced,already_decoded,pending,"
         "forced_evicted,writeback_complete\n";
  for (const auto& sample : samples) {
    for (const auto& code : sample.codes) {
      out << run_id << ',' << label << ',' << sample.invocation << ','
          << sample.buffered_fifo_enabled << ',';
      if (sample.buffered_fifo_enabled) {
        out << sample.buffered_service_time << ",B"
            << sample.buffered_batch_id;
      } else {
        out << ',';
      }
      out << ',' << code.code_index << ',' << code.time_index << ','
          << code.time_offset << ',' << level56_info_type_name(code.info_type)
          << ',' << level56_decode_status_name(code.decode_status) << ','
          << (code.code_index + 1u) << ',' << code.source_level << ','
          << (code.source_local_row + 1u) << ',' << code.source_global_row
          << ',' << (code.group_index + 1u)
          << ',' << (code.position_in_group + 1u) << ','
          << code.early_stop_hit << ','
          << level56_hybrid_class_name(code.hybrid_class) << ','
          << level56_eligibility_name(code.resource_eligibility) << ','
          << code.planned_hiso << ',' << code.planned_siso << ','
          << code.remaining_for_schedule << ','
          << level56_action_name(code.final_action) << ','
          << code.assigned_entry_slot << ',' << code.assigned_core << ','
          << code.produced << ',' << code.already_decoded << ','
          << code.pending << ',' << code.forced_evicted << ','
          << code.writeback_complete << '\n';
    }
  }
}

void dump_level56_equivalence_observation_csv(
    const std::vector<newcode::Level56ScheduleSample>& samples,
    const std::filesystem::path& output_path,
    const std::string& run_id,
    const std::string& label,
    const std::vector<std::size_t>& post_fec_error_positions) {
  const auto parent = output_path.parent_path();
  if (!parent.empty()) {
    std::filesystem::create_directories(parent);
  }
  std::ofstream out(output_path);
  out << "run_id,label,invocation,buffered_fifo_enabled,buffered_t,"
         "buffered_batch,code,source_level,source_local_row,source_global_row,"
         "early_stop_hit,hybrid_class,resource_eligibility,final_action,"
         "planned_hiso,planned_siso,assigned_entry_slot,assigned_core,"
         "produced,already_decoded,pending,forced_evicted,writeback_complete,"
         "channel_input_hash,decoder_input_hash,decoder_output_available,"
         "decoder_output_hash,tile_row_hash_after_writeback,"
         "level6_history_available,level6_history_hash_after_writeback,"
         "hiso_hard_word_observation_available,"
         "hiso_lin_hard_bit_errors_vs_expected,"
         "hiso_corrected_hard_bit_errors_vs_expected,"
         "hiso_corrected_hard_error_positions,"
         "hiso_corrected_error_info_positions,"
         "hiso_error_info_positions_in_final_post_fec_count\n";
  for (const auto& sample : samples) {
    for (const auto& code : sample.codes) {
      std::size_t final_overlap = 0;
      for (const auto position : code.hiso_corrected_error_info_positions) {
        if (std::binary_search(post_fec_error_positions.begin(),
                               post_fec_error_positions.end(), position)) {
          ++final_overlap;
        }
      }
      out << run_id << ',' << label << ',' << sample.invocation << ','
          << sample.buffered_fifo_enabled << ',';
      if (sample.buffered_fifo_enabled) {
        out << sample.buffered_service_time << ",B"
            << sample.buffered_batch_id;
      } else {
        out << ',';
      }
      out << ',' << (code.code_index + 1u) << ',' << code.source_level << ','
          << (code.source_local_row + 1u) << ',' << code.source_global_row
          << ',' << code.early_stop_hit << ','
          << level56_hybrid_class_name(code.hybrid_class) << ','
          << level56_eligibility_name(code.resource_eligibility) << ','
          << level56_action_name(code.final_action) << ','
          << code.planned_hiso << ',' << code.planned_siso << ','
          << code.assigned_entry_slot << ',' << code.assigned_core << ','
          << code.produced << ',' << code.already_decoded << ','
          << code.pending << ',' << code.forced_evicted << ','
          << code.writeback_complete << ',' << code.channel_input_hash << ','
          << code.decoder_input_hash << ',' << code.decoder_output_available
          << ',' << code.decoder_output_hash << ','
          << code.tile_row_hash_after_writeback << ','
          << code.level6_history_available << ','
          << code.level6_history_hash_after_writeback << ','
          << code.hiso_hard_word_observation_available << ','
          << code.hiso_lin_hard_bit_errors_vs_expected << ','
          << code.hiso_corrected_hard_bit_errors_vs_expected << ','
          << join_numbers(code.hiso_corrected_hard_error_positions) << ','
          << join_numbers(code.hiso_corrected_error_info_positions) << ','
          << final_overlap << '\n';
    }
  }
}

void dump_level56_buffered_times_csv(
    const std::vector<newcode::Level56BufferedTimeSample>& samples,
    const std::filesystem::path& output_path,
    const std::string& run_id,
    const std::string& label) {
  const auto parent = output_path.parent_path();
  if (!parent.empty()) {
    std::filesystem::create_directories(parent);
  }
  std::ofstream out(output_path);
  out << "run_id,label,t,arrived_batch,fifo_depth_before,had_head_before,"
         "head_batch_before,pending_before,ordinary_service_used,"
         "ordinary_batch,ordinary_schedule_invocation,ordinary_service_count,"
         "ordinary_entries_used,ordinary_services,C_t,"
         "full_early_stop_batches,forced_evicted_batches,S_t,S_next,"
         "fifo_depth_after,pending_after,retirements,forced_evicted_global_rows\n";
  for (const auto& sample : samples) {
    std::ostringstream retirements;
    std::ostringstream forced_rows;
    for (std::size_t index = 0; index < sample.retirements.size(); ++index) {
      if (index > 0) {
        retirements << '|';
      }
      const auto& retirement = sample.retirements[index];
      retirements << 'B' << retirement.batch_id << ':'
                  << level56_buffered_retirement_name(retirement.reason);
    }
    for (std::size_t index = 0;
         index < sample.forced_evicted_global_rows.size(); ++index) {
      if (index > 0) {
        forced_rows << '|';
      }
      forced_rows << sample.forced_evicted_global_rows[index];
    }
    out << run_id << ',' << label << ',' << sample.service_time << ",B"
        << sample.arrived_batch_id << ',' << sample.fifo_depth_before << ','
        << sample.had_head_before << ',';
    if (sample.had_head_before) {
      out << 'B' << sample.head_batch_id_before;
    }
    out << ',' << sample.pending_before << ','
        << sample.ordinary_service_used << ',';
    if (sample.ordinary_service_used) {
      out << 'B' << sample.ordinary_batch_id << ','
          << sample.ordinary_schedule_invocation;
    } else {
      out << ',';
    }
    std::ostringstream ordinary_services;
    for (std::size_t index = 0; index < sample.ordinary_services.size();
         ++index) {
      if (index > 0) {
        ordinary_services << '|';
      }
      const auto& service = sample.ordinary_services[index];
      ordinary_services << 'B' << service.batch_id
                        << ":I" << service.schedule_invocation
                        << ":budget" << service.entry_budget
                        << ":offset" << service.entry_slot_offset
                        << ":used" << service.entries_used
                        << ":complete" << service.completed;
    }
    out << ',' << sample.ordinary_service_count << ','
        << sample.ordinary_entries_used << ',' << ordinary_services.str()
        << ',' << sample.completed_batches << ','
        << sample.full_early_stop_batches << ','
        << sample.forced_evicted_batches << ',' << sample.window_start_before
        << ',' << sample.window_start_after << ',' << sample.fifo_depth_after
        << ',' << sample.pending_after << ',' << retirements.str() << ','
        << forced_rows.str() << '\n';
  }
}

void log_level56_buffered_summary(
    const std::vector<newcode::Level56BufferedTimeSample>& samples,
    io::DualWriter& log) {
  if (samples.empty()) {
    log << "[INFO] Level56 buffered FIFO samples: 0\n";
    return;
  }
  std::size_t completed = 0;
  std::size_t full_early_stop = 0;
  std::size_t forced_evicted = 0;
  std::size_t max_fifo_depth = 0;
  std::size_t min_window_start = samples.front().window_start_before;
  std::array<std::size_t, 3> completion_buckets{};
  for (const auto& sample : samples) {
    completed += sample.completed_batches;
    full_early_stop += sample.full_early_stop_batches;
    forced_evicted += sample.forced_evicted_batches;
    max_fifo_depth = std::max(max_fifo_depth, sample.fifo_depth_before);
    min_window_start = std::min(
        min_window_start,
        std::min(sample.window_start_before, sample.window_start_after));
    ++completion_buckets[std::min<std::size_t>(sample.completed_batches, 2u)];
  }
  log << "[INFO] Level56 buffered FIFO: t=" << samples.size()
      << ", completed=" << completed
      << ", full_early_stop=" << full_early_stop
      << ", forced_evicted=" << forced_evicted
      << ", max_fifo_depth=" << max_fifo_depth
      << ", min_S_t=" << min_window_start
      << ", C_t(0/1/2+)=" << completion_buckets[0] << '/'
      << completion_buckets[1] << '/' << completion_buckets[2] << "\n";
}

void log_level56_schedule_summary(
    const std::vector<newcode::Level56ScheduleSample>& samples,
    io::DualWriter& log) {
  std::array<std::size_t, 4> branch_counts{};
  std::array<std::size_t, 5> temporal_branch_counts{};
  std::array<std::size_t, 5> temporal_t0_pending_reduced{};
  std::array<std::size_t, 5> temporal_t1_pending_reduced{};
  std::size_t current_first_with_history_pending_windows = 0;
  std::size_t current_first_history_pending_codes = 0;
  std::size_t supplement_history_new_produced = 0;
  std::size_t temporal_window_count = 0;
  std::array<std::size_t, 16> group_entries{};
  std::array<std::array<std::size_t, 4>, 2> level_actions{};
  std::array<std::array<std::size_t, 4>, 9> class_actions{};
  std::size_t total_entries = 0;
  std::size_t total_hiso = 0;
  std::size_t total_siso = 0;
  std::size_t total_entry_capacity = 0;
  std::map<std::size_t, std::size_t> buffered_capacity_by_service_time;
  for (const auto& sample : samples) {
    ++branch_counts[static_cast<std::size_t>(sample.branch)];
    if (sample.temporal_lookahead_enabled) {
      ++temporal_window_count;
      ++temporal_branch_counts[
          static_cast<std::size_t>(sample.temporal_branch)];
      const auto temporal_index =
          static_cast<std::size_t>(sample.temporal_branch);
      temporal_t0_pending_reduced[temporal_index] +=
          sample.temporal_t0_pending_before -
          sample.temporal_t0_pending_after;
      temporal_t1_pending_reduced[temporal_index] +=
          sample.temporal_t1_pending_before -
          sample.temporal_t1_pending_after;
      if (sample.temporal_branch ==
              newcode::Level56TemporalBranch::CurrentFirst &&
          sample.temporal_x > 0) {
        ++current_first_with_history_pending_windows;
        current_first_history_pending_codes +=
            sample.temporal_t0_pending_after;
      }
      if (sample.temporal_branch ==
          newcode::Level56TemporalBranch::SupplementHistory) {
        supplement_history_new_produced +=
            sample.temporal_t0_new_produced;
      }
    }
    // A temporal sample merges history and current scheduling into one
    // observation, so its legacy absolute end position is the combined use.
    // Buffered FIFO samples remain one invocation each and must use the
    // per-invocation count to avoid counting the second batch's offset twice.
    total_entries += sample.temporal_lookahead_enabled
        ? sample.total_group_entries
        : sample.group_entries_used;
    total_hiso += sample.planned_hiso_count;
    total_siso += sample.planned_siso_count;
    if (sample.buffered_fifo_enabled) {
      auto& capacity =
          buffered_capacity_by_service_time[sample.buffered_service_time];
      capacity = std::max(
          capacity, sample.entry_slot_offset + sample.entry_capacity);
    } else {
      total_entry_capacity +=
          sample.entry_slot_offset + sample.entry_capacity;
    }
    for (std::size_t group = 0; group < group_entries.size(); ++group) {
      group_entries[group] += sample.group_entry_counts[group];
    }
    for (const auto& code : sample.codes) {
      if (sample.temporal_lookahead_enabled && code.time_index != 1) {
        continue;
      }
      if ((code.source_level == 5 || code.source_level == 6) &&
          code.final_action < 4) {
        ++level_actions[code.source_level - 5u][code.final_action];
      }
      if (code.hybrid_class < class_actions.size() &&
          code.final_action < class_actions.front().size()) {
        ++class_actions[code.hybrid_class][code.final_action];
      }
    }
  }
  for (const auto& [service_time, capacity] :
       buffered_capacity_by_service_time) {
    (void)service_time;
    total_entry_capacity += capacity;
  }
  log << "[RESULT] Level56 schedule calls/branches(K0,K<8,K=8,K>8) = "
      << samples.size() << "/[" << branch_counts[0] << ", "
      << branch_counts[1] << ", " << branch_counts[2] << ", "
      << branch_counts[3] << "]\n";
  if (temporal_window_count > 0) {
    log << "[RESULT] Level56 temporal windows/branches("
        << "NoHistory,SupplementHistory,CurrentFirst,NoFuture) = "
        << temporal_window_count << "/["
        << temporal_branch_counts[static_cast<std::size_t>(
               newcode::Level56TemporalBranch::NoHistory)]
        << ", "
        << temporal_branch_counts[static_cast<std::size_t>(
               newcode::Level56TemporalBranch::SupplementHistory)]
        << ", "
        << temporal_branch_counts[static_cast<std::size_t>(
               newcode::Level56TemporalBranch::CurrentFirst)]
        << ", "
        << temporal_branch_counts[static_cast<std::size_t>(
               newcode::Level56TemporalBranch::NoFuture)]
        << "]\n";
    log << "[RESULT] Level56 CurrentFirst windows with X>0/history pending "
        << "codes = " << current_first_with_history_pending_windows << '/'
        << current_first_history_pending_codes << "\n";
    log << "[RESULT] Level56 SupplementHistory t0 newly produced codes = "
        << supplement_history_new_produced << "\n";
    log << "[RESULT] Level56 temporal pending reduced by branch "
        << "(NoHistory,SupplementHistory,CurrentFirst,NoFuture), t0/t1/total = ";
    for (const auto branch : {
             newcode::Level56TemporalBranch::NoHistory,
             newcode::Level56TemporalBranch::SupplementHistory,
             newcode::Level56TemporalBranch::CurrentFirst,
             newcode::Level56TemporalBranch::NoFuture}) {
      const auto index = static_cast<std::size_t>(branch);
      log << (branch == newcode::Level56TemporalBranch::NoHistory ? "[" : ", ")
          << temporal_t0_pending_reduced[index] << '/'
          << temporal_t1_pending_reduced[index] << '/'
          << (temporal_t0_pending_reduced[index] +
              temporal_t1_pending_reduced[index]);
    }
    log << "]\n";
  }
  log << "[RESULT] Level56 entry/HISO/SISO/idle-HISO/idle-SISO totals = "
      << total_entries << '/' << total_hiso << '/' << total_siso << '/'
      << (total_entry_capacity >= total_hiso
              ? total_entry_capacity - total_hiso
              : 0u) << '/'
      << (total_entry_capacity >= total_siso
              ? total_entry_capacity - total_siso
              : 0u) << "\n";
  log << "[RESULT] Level56 per-group entry counts = ["
      << join_numbers(group_entries, ',') << "]\n";
  log << "[RESULT] Level5 actions(EarlyStop,HISO,SISO,Unscheduled) = ["
      << join_numbers(level_actions[0], ',') << "]\n";
  log << "[RESULT] Level6 actions(EarlyStop,HISO,SISO,Unscheduled) = ["
      << join_numbers(level_actions[1], ',') << "]\n";
  for (const uint8_t row_class : {uint8_t{3}, uint8_t{4}, uint8_t{5},
                                  uint8_t{6}, uint8_t{8}}) {
    log << "[RESULT] Level56 class "
        << level56_hybrid_class_name(row_class)
        << " actions(HISO,SISO,Unscheduled) = ["
        << class_actions[row_class][1] << ','
        << class_actions[row_class][2] << ','
        << class_actions[row_class][3] << "]\n";
  }
}

}  // namespace

namespace ofec_single {

int run_ofec_single(const Config& config) {
  const std::filesystem::path data_dir = "data";
  const std::string run_id = make_run_id();
  std::string log_path;
  std::ofstream log_file = io::prepare_log_file(data_dir, log_path);
  io::DualWriter log(log_file);

  std::optional<newcode::Params> params = detail::build_params(config, log);
  if (!params.has_value()) {
    return 2;
  }

  detail::log_run_overview(config, *params, log);
  newcode::PipelineConfig pipeline_cfg = detail::build_pipeline_config(config);
  const newcode::PipelineResult result =
      newcode::run_pipeline(*params, pipeline_cfg, config.label, config.ebn0_db);
  if (!config.post_fec_error_positions_output_path.empty()) {
    const std::filesystem::path output_path =
        config.post_fec_error_positions_output_path;
    const auto parent = output_path.parent_path();
    if (!parent.empty()) {
      std::filesystem::create_directories(parent);
    }
    std::ofstream out(output_path);
    for (const auto position : result.post_fec_error_positions) {
      out << position << '\n';
    }
    log << "[INFO] Full Post-FEC error positions saved to "
        << output_path.string() << "\n";
  }
  if (config.dump_tile_early_stop_samples) {
    std::filesystem::path output_path = config.tile_early_stop_samples_output_path;
    if (output_path.empty()) {
      output_path = std::filesystem::path("data/early_stop_hist") /
                    (config.label + "_tile_early_stop_samples.csv");
    }
    dump_tile_early_stop_samples_csv(result.tile_early_stop_samples, output_path);
    log << "[INFO] tile early-stop samples saved to " << output_path.string() << "\n";
  }
  if (config.dump_tile_early_stop_group_bind_debug_samples) {
    std::filesystem::path output_path =
        config.tile_early_stop_group_bind_debug_samples_output_path;
    if (output_path.empty()) {
      output_path = std::filesystem::path("data/early_stop_debug") /
                    (config.label + "_group_bind_debug.csv");
    }
    dump_tile_early_stop_group_bind_debug_samples_csv(
        result.tile_early_stop_group_bind_debug_samples,
        output_path,
        run_id,
        config.label);
    log << "[INFO] tile early-stop group-bind debug samples saved to "
        << output_path.string() << "\n";
  }
  if (config.dump_level56_schedule_stats) {
    std::filesystem::path rounds_path =
        config.level56_schedule_rounds_output_path;
    if (rounds_path.empty()) {
      rounds_path = std::filesystem::path("data/level56_schedule") /
                    (config.label + "_level56_schedule_rounds.csv");
    }
    std::filesystem::path codes_path =
        config.level56_schedule_codes_output_path;
    if (codes_path.empty()) {
      codes_path = std::filesystem::path("data/level56_schedule") /
                   (config.label + "_level56_schedule_codes.csv");
    }
    dump_level56_schedule_rounds_csv(
        result.level56_schedule_samples, rounds_path, run_id, config.label);
    dump_level56_schedule_codes_csv(
        result.level56_schedule_samples, codes_path, run_id, config.label);
    const std::filesystem::path buffered_times_path =
        std::filesystem::path("data/level56_schedule") /
        (config.label + "_level56_buffered_times.csv");
    dump_level56_buffered_times_csv(
        result.level56_buffered_time_samples, buffered_times_path,
        run_id, config.label);
    log_level56_schedule_summary(result.level56_schedule_samples, log);
    log_level56_buffered_summary(result.level56_buffered_time_samples, log);
    log << "[INFO] Level56 schedule rounds saved to "
        << rounds_path.string() << "\n";
    log << "[INFO] Level56 schedule codes saved to "
        << codes_path.string() << "\n";
    log << "[INFO] Level56 buffered times saved to "
        << buffered_times_path.string() << "\n";
  }
  if (config.dump_level56_equivalence_observation) {
    std::filesystem::path output_path =
        config.level56_equivalence_observation_output_path;
    if (output_path.empty()) {
      output_path = std::filesystem::path("data/level56_observation") /
                    (config.label + "_level56_equivalence_observation.csv");
    }
    dump_level56_equivalence_observation_csv(
        result.level56_schedule_samples, output_path, run_id, config.label,
        result.post_fec_error_positions);
    log << "[INFO] Level56 equivalence observation saved to "
        << output_path.string() << "\n";
  }
  detail::log_pipeline_results(result, log);
  log << "[INFO] log saved at " << log_path << "\n";
  return 0;
}

}  // namespace ofec_single
