#pragma once

#include <array>
#include <optional>
#include <sstream>

namespace newcode {
namespace detail {

constexpr std::size_t kLevel5TileIndex = 4;
constexpr std::size_t kLevel6TileIndex = 5;
constexpr std::size_t kLevel56GroupedCodeCount = 64;
constexpr std::size_t kLevel56GroupedGroupCount = 16;
constexpr std::size_t kLevel56CodesPerGroup = 4;
constexpr std::size_t kLevel56MaxGroupEntries = 8;

enum class Level56Eligibility : uint8_t {
  None = 0,
  HisoOnly,
  SisoOnly,
  HisoOrSiso
};

enum class Level56FinalAction : uint8_t {
  EarlyStopAction = 0,
  HisoDecode,
  SisoDecode,
  Unscheduled
};

struct Level56DispatchEntry {
  std::size_t shared_row = 0;
  std::size_t source_level = 0;
  std::size_t source_local_row = 0;
  std::size_t source_global_row = 0;
  bool early_stop_hit = false;
  HybridRowClass hybrid_class = HybridRowClass::None;
  Level56Eligibility eligibility = Level56Eligibility::None;
  bool planned_hiso = false;
  bool planned_siso = false;
  bool remaining_for_schedule = false;
  Level56FinalAction final_action = Level56FinalAction::Unscheduled;
  int assigned_entry_slot = -1;
  int assigned_core = -1;
  // 完整解码状态中的输出标志；调度阶段为 false，执行阶段完成后更新。
  bool produced = false;
  Level56DecodeStatus decode_status = Level56DecodeStatus::NotDecoded;
  // 历史批次补解时，避免已经完成的 action 被重复执行。
  bool needs_execution = true;
};

struct Level56EarlyStopEvaluation {
  TileEarlyStopResult raw;
  TileEarlyStopResult effective;
};

template <typename LLR>
struct Level56SharedResult {
  TileProcessResult<LLR> level5;
  TileProcessResult<LLR> level6;
  std::vector<Level56DispatchEntry> dispatch;
  bool has_schedule_sample = false;
  Level56ScheduleSample schedule_sample;
};

template <typename LLR>
struct Level56TemporalBatch {
  matrix::Matrix<LLR> tile_in5;
  matrix::Matrix<LLR> ch_tile5;
  matrix::Matrix<LLR> tile_in6;
  matrix::Matrix<LLR> ch_tile6;
  std::size_t tile_top5 = 0;
  std::size_t tile_top6 = 0;
  Level56EarlyStopEvaluation early5;
  Level56EarlyStopEvaluation early6;
  std::vector<Level56DispatchEntry> entries;
};

template <typename LLR>
struct Level56TemporalState {
  std::optional<Level56TemporalBatch<LLR>> previous;
  std::optional<Level56TemporalBatch<LLR>> current;
};

inline bool is_level56_pending_decode(const Level56DispatchEntry& entry) {
  return !entry.early_stop_hit &&
         entry.final_action == Level56FinalAction::Unscheduled &&
         !entry.produced;
}

inline std::size_t level56_temporal_pending_code_count(
    const std::vector<Level56DispatchEntry>& entries) {
  return static_cast<std::size_t>(std::count_if(
      entries.begin(), entries.end(), is_level56_pending_decode));
}

inline std::size_t level56_temporal_pending_group_count(
    const std::vector<Level56DispatchEntry>& entries) {
  std::array<bool, kLevel56GroupedGroupCount> groups{};
  for (const auto& entry : entries) {
    if (is_level56_pending_decode(entry)) {
      groups[entry.shared_row / kLevel56CodesPerGroup] = true;
    }
  }
  return static_cast<std::size_t>(
      std::count(groups.begin(), groups.end(), true));
}

inline std::size_t level56_temporal_candidate_group_count(
    const std::vector<Level56DispatchEntry>& entries) {
  std::array<bool, kLevel56GroupedGroupCount> groups{};
  for (const auto& entry : entries) {
    // K1/K2 count ordinary non-EarlyStop candidates that can enter the
    // original Group4 scheduler. Lookahead entries have no classification
    // yet, but their row-level EarlyStop result is sufficient for this count.
    if (!entry.early_stop_hit) {
      groups[entry.shared_row / kLevel56CodesPerGroup] = true;
    }
  }
  return static_cast<std::size_t>(
      std::count(groups.begin(), groups.end(), true));
}

inline Level56TemporalBranch select_level56_temporal_branch(
    bool has_history,
    bool has_future,
    std::size_t x,
    std::size_t k1,
    std::size_t k2,
    std::size_t group_load_threshold = kLevel56GroupedGroupCount) {
  if (!has_history) {
    return Level56TemporalBranch::NoHistory;
  }
  if (!has_future) {
    return Level56TemporalBranch::NoFuture;
  }
  if (x > kLevel56GroupedGroupCount ||
      k1 > kLevel56GroupedGroupCount ||
      k2 > kLevel56GroupedGroupCount) {
    throw std::invalid_argument(
        "LEVEL56 temporal group counts must be within 0..16");
  }
  if (group_load_threshold > 3 * kLevel56GroupedGroupCount) {
    throw std::invalid_argument(
        "LEVEL56 temporal group load threshold must be within 0..48");
  }
  return x > 0 && x + k1 + k2 < group_load_threshold
             ? Level56TemporalBranch::SupplementHistory
             : Level56TemporalBranch::CurrentFirst;
}

inline Level56ScheduleCodeSample make_level56_schedule_code_sample(
    const Level56DispatchEntry& entry,
    std::size_t time_index,
    Level56TemporalInfoType info_type,
    bool use_temporal_shared_index) {
  return Level56ScheduleCodeSample{
      .code_index = use_temporal_shared_index
                        ? time_index * kLevel56GroupedCodeCount +
                              entry.shared_row
                        : entry.shared_row,
      .time_index = time_index,
      .time_offset = static_cast<int>(time_index) - 1,
      .info_type = info_type,
      .decode_status = entry.decode_status,
      .source_level = entry.source_level,
      .source_local_row = entry.source_local_row,
      .group_index = entry.shared_row / kLevel56CodesPerGroup,
      .position_in_group = entry.shared_row % kLevel56CodesPerGroup,
      .early_stop_hit = entry.early_stop_hit,
      .hybrid_class = static_cast<uint8_t>(entry.hybrid_class),
      .resource_eligibility = static_cast<uint8_t>(entry.eligibility),
      .planned_hiso = entry.planned_hiso,
      .planned_siso = entry.planned_siso,
      .remaining_for_schedule = entry.remaining_for_schedule,
      .final_action = static_cast<uint8_t>(entry.final_action),
      .assigned_entry_slot = entry.assigned_entry_slot,
      .assigned_core = entry.assigned_core,
      .produced = entry.produced,
  };
}

inline int pick_level56_int(const std::vector<int>& values,
                            std::size_t index,
                            int fallback) {
  return index < values.size() ? values[index] : fallback;
}

inline float pick_level56_float(const std::vector<float>& values,
                                std::size_t index,
                                float fallback) {
  return index < values.size() ? values[index] : fallback;
}

inline void validate_level56_shared_config(const newcode::Params& p) {
  if (!p.LEVEL56_SHARED_ENABLE) {
    if (p.LEVEL56_TEMPORAL_LOOKAHEAD_ENABLE) {
      throw std::invalid_argument(
          "LEVEL56 temporal lookahead requires shared mode ON");
    }
    return;
  }
  if (p.LEVEL56_TEMPORAL_LOOKAHEAD_ENABLE &&
      p.LEVEL56_SCHEDULE_MODE !=
          newcode::Level56ScheduleMode::Group4LoadSortedMultiround) {
    throw std::invalid_argument(
        "LEVEL56 temporal lookahead requires Group4 multiround scheduling");
  }
  if (p.LEVEL56_TEMPORAL_LOOKAHEAD_ENABLE &&
      p.LEVEL56_EARLY_STOP_GROUP_UPDATE_MODE !=
          newcode::Level56EarlyStopGroupUpdateMode::AllGroups) {
    throw std::invalid_argument(
        "LEVEL56 temporal lookahead requires AllGroups early-stop updates");
  }
  if (p.LEVEL56_TEMPORAL_LOOKAHEAD_ENABLE) {
    for (const std::size_t tile : {kLevel5TileIndex, kLevel6TileIndex}) {
      if (pick_level56_int(p.EARLY_STOP_BIND_GROUP_SIZE_LIST, tile,
                           p.EARLY_STOP_BIND_GROUP_SIZE) != 1) {
        throw std::invalid_argument(
            "LEVEL56 temporal lookahead requires row-level early-stop "
            "decisions for both Level 5 and Level 6");
      }
    }
    if (p.LEVEL56_SINGLE_LEVEL_SELECT_ENABLE) {
      throw std::invalid_argument(
          "LEVEL56 temporal lookahead does not support single-level "
          "selection; AllGroups requires both levels");
    }
  }
  if (p.TILES_PER_WIN != 6 || p.CHASE_SBR != 2 ||
      p.TILE_OVERLAP_BR != 0 || p.MUX_GROUP_G != 1) {
    throw std::invalid_argument(
        "LEVEL56 shared mode requires TILES_PER_WIN=6, CHASE_SBR=2, "
        "TILE_OVERLAP_BR=0 and MUX_GROUP_G=1");
  }
  if (p.LEVEL56_SHARED_HISO_ACTIVE < 0 ||
      p.LEVEL56_SHARED_HISO_ACTIVE > 64 ||
      p.LEVEL56_SHARED_SISO_ACTIVE < 0 ||
      p.LEVEL56_SHARED_SISO_ACTIVE > 64) {
    throw std::invalid_argument(
        "LEVEL56 shared HISO/SISO capacities must be in [0,64]");
  }
  switch (p.LEVEL56_PRIORITY_MODE) {
    case newcode::Level56PriorityMode::Level5First:
    case newcode::Level56PriorityMode::Level6First:
      break;
    default:
      throw std::invalid_argument(
          "LEVEL56 shared priority mode is invalid");
  }
  switch (p.LEVEL56_SCHEDULE_MODE) {
    case newcode::Level56ScheduleMode::GlobalPriority:
    case newcode::Level56ScheduleMode::Group4LoadSortedMultiround:
      break;
    default:
      throw std::invalid_argument(
          "LEVEL56 shared schedule mode is invalid");
  }
  switch (p.LEVEL56_EARLY_STOP_GROUP_UPDATE_MODE) {
    case newcode::Level56EarlyStopGroupUpdateMode::AllGroups:
    case newcode::Level56EarlyStopGroupUpdateMode::EnteredGroupsOnly:
    case newcode::Level56EarlyStopGroupUpdateMode::FillIdleEntries:
      break;
    default:
      throw std::invalid_argument(
          "LEVEL56 early-stop group update mode is invalid");
  }

  if (p.LEVEL56_SCHEDULE_MODE ==
      newcode::Level56ScheduleMode::Group4LoadSortedMultiround) {
    if (p.LEVEL56_SHARED_HISO_ACTIVE != 8 ||
        p.LEVEL56_SHARED_SISO_ACTIVE != 8) {
      throw std::invalid_argument(
          "LEVEL56 grouped multiround mode requires HISO/SISO capacity 8/8");
    }
    if (p.LEVEL56_SINGLE_LEVEL_SELECT_ENABLE) {
      throw std::invalid_argument(
          "LEVEL56 grouped multiround mode requires single-level selection OFF");
    }
    if (p.LEVEL56_PRIORITY_MODE != newcode::Level56PriorityMode::Level5First) {
      throw std::invalid_argument(
          "LEVEL56 grouped multiround mode requires Level5First tie-breaking");
    }
    if (effective_hybrid_classifier_mode(p) !=
        newcode::HybridClassifierMode::FriendS1S3WithS0Classifier) {
      throw std::invalid_argument(
          "LEVEL56 grouped multiround mode requires "
          "FriendS1S3WithS0Classifier");
    }
  }

  const auto same_int = [&](const std::vector<int>& values, int fallback) {
    return pick_level56_int(values, kLevel5TileIndex, fallback) ==
           pick_level56_int(values, kLevel6TileIndex, fallback);
  };
  const auto same_float = [&](const std::vector<float>& values, float fallback) {
    return pick_level56_float(values, kLevel5TileIndex, fallback) ==
           pick_level56_float(values, kLevel6TileIndex, fallback);
  };
  if (!same_int(p.EARLY_STOP_ENABLE_LIST, p.ENABLE_EARLY_STOP ? 1 : 0) ||
      !same_int(p.EARLY_STOP_CONDITION_MODE_LIST,
                p.EARLY_STOP_CONDITION_MODE) ||
      !same_int(p.EARLY_STOP_ACTION_MODE_LIST, p.EARLY_STOP_ACTION_MODE) ||
      !same_int(p.EARLY_STOP_BIND_GROUP_SIZE_LIST,
                p.EARLY_STOP_BIND_GROUP_SIZE) ||
      !same_int(p.HYBRID_ENABLE_LIST, p.HYBRID_ENABLE ? 1 : 0) ||
      !same_float(p.EARLY_STOP_ACTION_SIGN_BETA_LIST,
                  p.EARLY_STOP_ACTION_SIGN_BETA) ||
      !same_float(p.HYBRID_HARD_LLR_MAG_LIST, p.HYBRID_HARD_LLR_MAG)) {
    throw std::invalid_argument(
        "LEVEL56 shared mode requires identical Level 5/6 parameters "
        "except ALPHA_LIST and beta_list");
  }
  if (pick_level56_int(p.HARD_TILE_LIST, kLevel5TileIndex,
                       p.HARD_DECODE_DEFAULT ? 1 : 0) != 0 ||
      pick_level56_int(p.HARD_TILE_LIST, kLevel6TileIndex,
                       p.HARD_DECODE_DEFAULT ? 1 : 0) != 0) {
    throw std::invalid_argument(
        "LEVEL56 shared mode requires Level 5 and Level 6 to be soft tiles");
  }
  if (pick_level56_int(p.HYBRID_ENABLE_LIST, kLevel5TileIndex,
                       p.HYBRID_ENABLE ? 1 : 0) == 0) {
    throw std::invalid_argument(
        "LEVEL56 shared mode requires hybrid classification on Level 5/6");
  }
  if (pick_level56_int(p.EARLY_STOP_ENABLE_LIST, kLevel5TileIndex,
                       p.ENABLE_EARLY_STOP ? 1 : 0) == 0) {
    throw std::invalid_argument(
        "LEVEL56 shared mode requires early-stop on Level 5/6");
  }
  if (effective_hybrid_classifier_mode(p) ==
      newcode::HybridClassifierMode::LegacyHardDecode) {
    throw std::invalid_argument(
        "LEVEL56 shared mode requires a classify-only hybrid classifier");
  }
}

inline int level56_class_priority(HybridRowClass row_class) {
  switch (row_class) {
    case HybridRowClass::ParityOnly:
      return 0;
    case HybridRowClass::OneMain:
      return 1;
    case HybridRowClass::OneMainPlusParity:
      return 2;
    case HybridRowClass::TwoMain:
      return 3;
    case HybridRowClass::Suspicious:
      return 4;
    case HybridRowClass::HardFail:
      return 5;
    default:
      return 6;
  }
}

inline bool level56_hiso_eligible(HybridRowClass row_class) {
  return row_class == HybridRowClass::ParityOnly ||
         row_class == HybridRowClass::OneMain ||
         row_class == HybridRowClass::OneMainPlusParity ||
         row_class == HybridRowClass::TwoMain;
}

inline std::size_t select_level56_decode_level(
    const TileEarlyStopResult& early5,
    const TileEarlyStopResult& early6,
    const newcode::Params& p) {
  if (!p.LEVEL56_SINGLE_LEVEL_SELECT_ENABLE) {
    return 0;
  }
  if (early5.rows_passed < early6.rows_passed) {
    return 5;
  }
  if (early6.rows_passed < early5.rows_passed) {
    return 6;
  }
  if (p.LEVEL56_PRIORITY_MODE == newcode::Level56PriorityMode::Level6First) {
    return 6;
  }
  if (p.LEVEL56_PRIORITY_MODE == newcode::Level56PriorityMode::Level5First) {
    return 5;
  }
  throw std::invalid_argument("LEVEL56 shared priority mode is invalid");
}

inline std::vector<std::size_t> order_level56_candidates(
    const std::vector<Level56DispatchEntry>& entries,
    newcode::Level56PriorityMode mode,
    std::size_t selected_level = 0) {
  std::vector<std::size_t> ordered;
  ordered.reserve(entries.size());

  for (int priority = 0; priority <= 5; ++priority) {
    std::vector<std::size_t> level5;
    std::vector<std::size_t> level6;
    for (std::size_t index = 0; index < entries.size(); ++index) {
      const auto& entry = entries[index];
      if (entry.early_stop_hit ||
          (selected_level != 0 && entry.source_level != selected_level) ||
          level56_class_priority(entry.hybrid_class) != priority) {
        continue;
      }
      (entry.source_level == 5 ? level5 : level6).push_back(index);
    }
    const auto by_local_row = [&](std::size_t lhs, std::size_t rhs) {
      return entries[lhs].source_local_row < entries[rhs].source_local_row;
    };
    std::sort(level5.begin(), level5.end(), by_local_row);
    std::sort(level6.begin(), level6.end(), by_local_row);

    if (mode == newcode::Level56PriorityMode::Level6First) {
      ordered.insert(ordered.end(), level6.begin(), level6.end());
      ordered.insert(ordered.end(), level5.begin(), level5.end());
    } else if (mode == newcode::Level56PriorityMode::Level5First) {
      ordered.insert(ordered.end(), level5.begin(), level5.end());
      ordered.insert(ordered.end(), level6.begin(), level6.end());
    } else {
      throw std::invalid_argument("LEVEL56 shared priority mode is invalid");
    }
  }
  return ordered;
}

inline void schedule_level56_rows_global_priority(
    std::vector<Level56DispatchEntry>* entries,
    const newcode::Params& p,
    std::size_t selected_level = 0) {
  for (auto& entry : *entries) {
    entry.planned_hiso = false;
    entry.planned_siso = false;
    entry.remaining_for_schedule = !entry.early_stop_hit;
    entry.assigned_entry_slot = -1;
    entry.assigned_core = -1;
  }
  const auto ordered = order_level56_candidates(
      *entries, p.LEVEL56_PRIORITY_MODE, selected_level);

  // The unselected level cannot consume shared HISO/SISO capacity. An
  // explicitly enabled early-stop action is preserved because it uses neither.
  if (selected_level != 0) {
    for (auto& entry : *entries) {
      const bool preserve_early_stop_action =
          p.LEVEL56_UNSELECTED_EARLY_STOP_ACTION_ENABLE &&
          entry.early_stop_hit;
      if (entry.source_level != selected_level &&
          !preserve_early_stop_action) {
        entry.final_action = Level56FinalAction::Unscheduled;
      }
    }
  }

  std::size_t soft_total = 0;
  std::size_t flexible_total = 0;
  for (const auto& entry : *entries) {
    if (!entry.early_stop_hit &&
        (selected_level == 0 || entry.source_level == selected_level)) {
      ++soft_total;
      flexible_total += entry.eligibility == Level56Eligibility::HisoOrSiso;
    }
  }
  const std::size_t siso_capacity =
      static_cast<std::size_t>(p.LEVEL56_SHARED_SISO_ACTIVE);
  const std::size_t hiso_capacity =
      static_cast<std::size_t>(p.LEVEL56_SHARED_HISO_ACTIVE);
  const std::size_t reclaim_needed =
      soft_total > siso_capacity ? soft_total - siso_capacity : 0u;
  const std::size_t reclaim_count =
      std::min({reclaim_needed, flexible_total, hiso_capacity});

  std::size_t reclaimed = 0;
  for (std::size_t index : ordered) {
    auto& entry = (*entries)[index];
    if (reclaimed < reclaim_count &&
        entry.eligibility == Level56Eligibility::HisoOrSiso) {
      entry.final_action = Level56FinalAction::HisoDecode;
      ++reclaimed;
    }
  }

  std::size_t siso_scheduled = 0;
  for (std::size_t index : ordered) {
    auto& entry = (*entries)[index];
    if (entry.final_action == Level56FinalAction::HisoDecode) {
      continue;
    }
    if ((entry.eligibility == Level56Eligibility::SisoOnly ||
         entry.eligibility == Level56Eligibility::HisoOrSiso) &&
        siso_scheduled < siso_capacity) {
      entry.final_action = Level56FinalAction::SisoDecode;
      ++siso_scheduled;
    } else {
      entry.final_action = Level56FinalAction::Unscheduled;
    }
  }

  std::size_t final_hiso_count = 0;
  std::size_t final_siso_count = 0;
  for (auto& entry : *entries) {
    if (selected_level != 0 && entry.source_level != selected_level) {
      const auto expected_action =
          p.LEVEL56_UNSELECTED_EARLY_STOP_ACTION_ENABLE &&
                  entry.early_stop_hit
              ? Level56FinalAction::EarlyStopAction
              : Level56FinalAction::Unscheduled;
      if (entry.final_action != expected_action) {
        throw std::logic_error(
            "LEVEL56 single-level mode did not bypass the unselected level");
      }
      continue;
    }
    if (entry.early_stop_hit !=
        (entry.final_action == Level56FinalAction::EarlyStopAction)) {
      throw std::logic_error(
          "LEVEL56 scheduler violated early-stop action conservation");
    }
    if (entry.final_action == Level56FinalAction::HisoDecode) {
      if (entry.eligibility != Level56Eligibility::HisoOnly &&
          entry.eligibility != Level56Eligibility::HisoOrSiso) {
        throw std::logic_error(
            "LEVEL56 scheduler assigned an ineligible row to HISO");
      }
      entry.planned_hiso = true;
      entry.remaining_for_schedule = false;
      ++final_hiso_count;
    } else if (entry.final_action == Level56FinalAction::SisoDecode) {
      if (entry.eligibility != Level56Eligibility::SisoOnly &&
          entry.eligibility != Level56Eligibility::HisoOrSiso) {
        throw std::logic_error(
            "LEVEL56 scheduler assigned an ineligible row to SISO");
      }
      entry.planned_siso = true;
      entry.remaining_for_schedule = false;
      ++final_siso_count;
    }
  }
  if (final_hiso_count > hiso_capacity ||
      final_siso_count > siso_capacity) {
    throw std::logic_error("LEVEL56 scheduler exceeded shared capacity");
  }
}

inline int level56_group_flexible_priority(HybridRowClass row_class) {
  switch (row_class) {
    case HybridRowClass::TwoMain:
      return 0;
    case HybridRowClass::OneMainPlusParity:
      return 1;
    case HybridRowClass::OneMain:
      return 2;
    case HybridRowClass::ParityOnly:
      return 3;
    default:
      return 4;
  }
}

inline std::size_t level56_group_remaining_count(
    const std::vector<Level56DispatchEntry>& entries,
    std::size_t group_index) {
  const std::size_t begin = group_index * kLevel56CodesPerGroup;
  const std::size_t end = begin + kLevel56CodesPerGroup;
  return static_cast<std::size_t>(std::count_if(
      entries.begin() + static_cast<std::ptrdiff_t>(begin),
      entries.begin() + static_cast<std::ptrdiff_t>(end),
      [](const auto& entry) { return entry.remaining_for_schedule; }));
}

inline void validate_level56_grouped_entries(
    const std::vector<Level56DispatchEntry>& entries) {
  if (entries.size() != kLevel56GroupedCodeCount) {
    throw std::invalid_argument(
        "LEVEL56 grouped multiround mode requires exactly 64 shared rows");
  }
  for (std::size_t index = 0; index < entries.size(); ++index) {
    const auto& entry = entries[index];
    const std::size_t expected_level = index < 32 ? 5u : 6u;
    const std::size_t expected_local_row = index % 32u;
    if (entry.shared_row != index || entry.source_level != expected_level ||
        entry.source_local_row != expected_local_row) {
      throw std::invalid_argument(
          "LEVEL56 grouped multiround mode requires code 1-32 = Level 5 "
          "row 1-32 and code 33-64 = Level 6 row 1-32");
    }
  }
}

inline void plan_level56_group_entry(
    std::vector<Level56DispatchEntry>* entries,
    std::size_t group_index,
    int entry_slot) {
  const std::size_t begin = group_index * kLevel56CodesPerGroup;
  const std::size_t end = begin + kLevel56CodesPerGroup;
  std::vector<std::size_t> siso_only;
  std::vector<std::size_t> flexible;

  for (std::size_t index = begin; index < end; ++index) {
    const auto& entry = (*entries)[index];
    if (!entry.remaining_for_schedule) {
      continue;
    }
    if (entry.eligibility == Level56Eligibility::SisoOnly) {
      if (entry.hybrid_class != HybridRowClass::HardFail) {
        throw std::logic_error(
            "LEVEL56 grouped scheduler received non-HardFail SisoOnly row");
      }
      siso_only.push_back(index);
    } else if (entry.eligibility == Level56Eligibility::HisoOrSiso) {
      if (level56_group_flexible_priority(entry.hybrid_class) >= 4) {
        throw std::logic_error(
            "LEVEL56 grouped scheduler received unsupported flexible class");
      }
      flexible.push_back(index);
    } else {
      throw std::logic_error(
          "LEVEL56 grouped scheduler received an ineligible remaining row");
    }
  }

  const auto by_group_position = [](std::size_t lhs, std::size_t rhs) {
    return lhs < rhs;
  };
  std::sort(siso_only.begin(), siso_only.end(), by_group_position);
  std::sort(flexible.begin(), flexible.end(), [&](std::size_t lhs,
                                                  std::size_t rhs) {
    const int lhs_priority =
        level56_group_flexible_priority((*entries)[lhs].hybrid_class);
    const int rhs_priority =
        level56_group_flexible_priority((*entries)[rhs].hybrid_class);
    return lhs_priority != rhs_priority ? lhs_priority < rhs_priority
                                        : lhs < rhs;
  });

  std::optional<std::size_t> siso_index;
  if (!siso_only.empty()) {
    siso_index = siso_only.front();
  } else if (!flexible.empty()) {
    siso_index = flexible.front();
    flexible.erase(flexible.begin());
  }

  const std::optional<std::size_t> hiso_index =
      flexible.empty() ? std::nullopt
                       : std::optional<std::size_t>(flexible.front());

  if (!siso_index && !hiso_index) {
    throw std::logic_error(
        "LEVEL56 grouped scheduler selected a group without a candidate");
  }

  if (siso_index) {
    auto& entry = (*entries)[*siso_index];
    entry.planned_siso = true;
    entry.remaining_for_schedule = false;
    entry.final_action = Level56FinalAction::SisoDecode;
    entry.assigned_entry_slot = entry_slot;
    entry.assigned_core = entry_slot;
  }
  if (hiso_index) {
    auto& entry = (*entries)[*hiso_index];
    entry.planned_hiso = true;
    entry.remaining_for_schedule = false;
    entry.final_action = Level56FinalAction::HisoDecode;
    entry.assigned_entry_slot = entry_slot;
    entry.assigned_core = entry_slot;
  }
}

inline void schedule_level56_rows_group4_load_sorted_multiround(
    std::vector<Level56DispatchEntry>* entries,
    const newcode::Params& p,
    std::size_t selected_level,
    Level56ScheduleSample* sample,
    std::size_t max_group_entries = kLevel56MaxGroupEntries,
    std::size_t entry_slot_offset = 0,
    bool preserve_entry_state = false) {
  if (!p.LEVEL56_SHARED_ENABLE) {
    throw std::invalid_argument(
        "LEVEL56 grouped multiround mode requires shared mode ON");
  }
  if (selected_level != 0 || p.LEVEL56_SINGLE_LEVEL_SELECT_ENABLE) {
    throw std::invalid_argument(
        "LEVEL56 grouped multiround mode does not support single-level selection");
  }
  if (p.LEVEL56_SHARED_HISO_ACTIVE != 8 ||
      p.LEVEL56_SHARED_SISO_ACTIVE != 8) {
    throw std::invalid_argument(
        "LEVEL56 grouped multiround mode requires HISO/SISO capacity 8/8");
  }
  if (p.LEVEL56_PRIORITY_MODE != newcode::Level56PriorityMode::Level5First) {
    throw std::invalid_argument(
        "LEVEL56 grouped multiround mode requires Level5First tie-breaking");
  }
  if (effective_hybrid_classifier_mode(p) !=
      newcode::HybridClassifierMode::FriendS1S3WithS0Classifier) {
    throw std::invalid_argument(
        "LEVEL56 grouped multiround mode requires "
        "FriendS1S3WithS0Classifier");
  }
  validate_level56_grouped_entries(*entries);

  if (sample) {
    *sample = Level56ScheduleSample{};
  }

  for (auto& entry : *entries) {
    entry.planned_hiso = false;
    entry.planned_siso = false;
    if (!preserve_entry_state) {
      entry.remaining_for_schedule = !entry.early_stop_hit;
      entry.final_action = entry.early_stop_hit
                               ? Level56FinalAction::EarlyStopAction
                               : Level56FinalAction::Unscheduled;
      entry.produced = false;
      entry.decode_status = Level56DecodeStatus::NotDecoded;
    }
    entry.assigned_entry_slot = -1;
    entry.assigned_core = -1;
  }

  std::array<std::size_t, kLevel56GroupedGroupCount> initial_counts{};
  std::vector<std::size_t> initial_nonzero_groups;
  for (std::size_t group = 0; group < kLevel56GroupedGroupCount; ++group) {
    initial_counts[group] = level56_group_remaining_count(*entries, group);
    if (initial_counts[group] > 0) {
      initial_nonzero_groups.push_back(group);
    }
  }

  if (sample) {
    sample->initial_counts = initial_counts;
    sample->initial_nonzero_groups = initial_nonzero_groups.size();
  }

  const auto sort_groups_by_count = [](std::vector<std::size_t>* groups,
                                       const auto& counts) {
    std::sort(groups->begin(), groups->end(), [&](std::size_t lhs,
                                                  std::size_t rhs) {
      return counts[lhs] != counts[rhs] ? counts[lhs] > counts[rhs]
                                        : lhs < rhs;
    });
  };

  if (max_group_entries > kLevel56MaxGroupEntries ||
      entry_slot_offset + max_group_entries > kLevel56MaxGroupEntries) {
    throw std::invalid_argument(
        "LEVEL56 grouped scheduler budget must fit within eight entry slots");
  }
  int used_group_entries = 0;
  std::array<bool, kLevel56GroupedGroupCount> group_entered{};
  const auto plan_groups = [&](const std::vector<std::size_t>& groups,
                               std::size_t count) {
    Level56ScheduleRoundSample round;
    if (sample) {
      round.round_index = sample->rounds.size();
      round.used_entries_before =
          entry_slot_offset + static_cast<std::size_t>(used_group_entries);
      for (std::size_t group = 0; group < kLevel56GroupedGroupCount; ++group) {
        round.remaining_before[group] =
            level56_group_remaining_count(*entries, group);
      }
      round.selected_groups.assign(groups.begin(), groups.begin() +
                                                    static_cast<std::ptrdiff_t>(count));
    }
    for (std::size_t i = 0; i < count; ++i) {
      plan_level56_group_entry(
          entries, groups[i],
          static_cast<int>(entry_slot_offset +
                           static_cast<std::size_t>(used_group_entries)));
      group_entered[groups[i]] = true;
      if (sample) {
        ++sample->group_entry_counts[groups[i]];
      }
      ++used_group_entries;
    }
    if (sample) {
      round.used_entries_after =
          entry_slot_offset + static_cast<std::size_t>(used_group_entries);
      for (std::size_t group = 0; group < kLevel56GroupedGroupCount; ++group) {
        round.remaining_after[group] =
            level56_group_remaining_count(*entries, group);
      }
      sample->rounds.push_back(std::move(round));
    }
  };

  const std::size_t initial_nonzero_count = initial_nonzero_groups.size();
  if (initial_nonzero_count == 0) {
    if (sample) {
      sample->branch = Level56ScheduleBranch::K0;
    }
  } else if (initial_nonzero_count > max_group_entries) {
    if (sample) {
      sample->branch = Level56ScheduleBranch::KGreaterThan8;
    }
    if (max_group_entries > 0) {
      sort_groups_by_count(&initial_nonzero_groups, initial_counts);
      plan_groups(initial_nonzero_groups, max_group_entries);
    }
  } else if (initial_nonzero_count == max_group_entries) {
    if (sample) {
      sample->branch = Level56ScheduleBranch::KEqual8;
    }
    plan_groups(initial_nonzero_groups, initial_nonzero_count);
  } else if (initial_nonzero_count > 0) {
    if (sample) {
      sample->branch = Level56ScheduleBranch::KLessThan8;
    }
    plan_groups(initial_nonzero_groups, initial_nonzero_count);

    while (initial_nonzero_count < max_group_entries &&
           used_group_entries < static_cast<int>(max_group_entries)) {
      std::array<std::size_t, kLevel56GroupedGroupCount> remaining_counts{};
      std::vector<std::size_t> remaining_groups;
      for (std::size_t group = 0; group < kLevel56GroupedGroupCount; ++group) {
        remaining_counts[group] = level56_group_remaining_count(*entries, group);
        if (remaining_counts[group] > 0) {
          remaining_groups.push_back(group);
        }
      }
      if (remaining_groups.empty()) {
        break;
      }
      sort_groups_by_count(&remaining_groups, remaining_counts);
      const std::size_t available_entries =
          max_group_entries -
          static_cast<std::size_t>(used_group_entries);
      plan_groups(remaining_groups,
                  std::min(available_entries, remaining_groups.size()));
    }
  }

  if (p.LEVEL56_EARLY_STOP_GROUP_UPDATE_MODE ==
          newcode::Level56EarlyStopGroupUpdateMode::FillIdleEntries &&
      used_group_entries < static_cast<int>(max_group_entries)) {
    Level56ScheduleRoundSample round;
    if (sample) {
      round.round_index = sample->rounds.size();
      round.used_entries_before =
          entry_slot_offset + static_cast<std::size_t>(used_group_entries);
      for (std::size_t group = 0; group < kLevel56GroupedGroupCount; ++group) {
        round.remaining_before[group] =
            level56_group_remaining_count(*entries, group);
      }
    }
    for (std::size_t group = 0;
         group < kLevel56GroupedGroupCount &&
         used_group_entries < static_cast<int>(max_group_entries);
         ++group) {
      if (group_entered[group]) {
        continue;
      }
      group_entered[group] = true;
      ++used_group_entries;
      if (sample) {
        round.selected_groups.push_back(group);
        ++sample->group_entry_counts[group];
      }
    }
    if (sample && !round.selected_groups.empty()) {
      round.used_entries_after =
          entry_slot_offset + static_cast<std::size_t>(used_group_entries);
      for (std::size_t group = 0; group < kLevel56GroupedGroupCount; ++group) {
        round.remaining_after[group] =
            level56_group_remaining_count(*entries, group);
      }
      sample->rounds.push_back(std::move(round));
    }
  }

  if (p.LEVEL56_EARLY_STOP_GROUP_UPDATE_MODE !=
      newcode::Level56EarlyStopGroupUpdateMode::AllGroups) {
    for (std::size_t group = 0; group < kLevel56GroupedGroupCount; ++group) {
      if (group_entered[group]) {
        continue;
      }
      const std::size_t begin = group * kLevel56CodesPerGroup;
      const std::size_t end = begin + kLevel56CodesPerGroup;
      for (std::size_t index = begin; index < end; ++index) {
        auto& entry = (*entries)[index];
        if (!entry.early_stop_hit) {
          continue;
        }
        entry.final_action = Level56FinalAction::Unscheduled;
        entry.remaining_for_schedule = false;
        entry.assigned_entry_slot = -1;
        entry.assigned_core = -1;
      }
    }
  }

  std::size_t planned_hiso = 0;
  std::size_t planned_siso = 0;
  for (const auto& entry : *entries) {
    const std::size_t group = entry.shared_row / kLevel56CodesPerGroup;
    const bool early_stop_action_expected =
        entry.early_stop_hit &&
        (p.LEVEL56_EARLY_STOP_GROUP_UPDATE_MODE ==
             newcode::Level56EarlyStopGroupUpdateMode::AllGroups ||
         group_entered[group]);
    if (early_stop_action_expected !=
        (entry.final_action == Level56FinalAction::EarlyStopAction)) {
      throw std::logic_error(
          "LEVEL56 grouped scheduler violated early-stop action conservation");
    }
    if (entry.planned_hiso && entry.planned_siso) {
      throw std::logic_error(
          "LEVEL56 grouped scheduler reserved one row twice");
    }
    if (entry.planned_hiso || entry.planned_siso) {
      if (entry.assigned_entry_slot < 0 ||
          entry.assigned_entry_slot >=
              static_cast<int>(kLevel56MaxGroupEntries) ||
          entry.assigned_core != entry.assigned_entry_slot ||
          entry.remaining_for_schedule) {
        throw std::logic_error(
            "LEVEL56 grouped scheduler produced an invalid entry-slot route");
      }
    }
    planned_hiso += entry.planned_hiso;
    planned_siso += entry.planned_siso;
  }
  if (used_group_entries > static_cast<int>(max_group_entries) ||
      planned_hiso > max_group_entries ||
      planned_siso > max_group_entries) {
    throw std::logic_error(
        "LEVEL56 grouped scheduler exceeded the eight-entry budget");
  }

  if (sample) {
    sample->total_group_entries =
        entry_slot_offset + static_cast<std::size_t>(used_group_entries);
    sample->planned_hiso_count = planned_hiso;
    sample->planned_siso_count = planned_siso;
    for (auto& round : sample->rounds) {
      round.initial_counts = sample->initial_counts;
      round.initial_nonzero_groups = sample->initial_nonzero_groups;
      round.branch = sample->branch;
    }
    sample->codes.reserve(entries->size());
    for (const auto& entry : *entries) {
      sample->codes.push_back(make_level56_schedule_code_sample(
          entry, 1, Level56TemporalInfoType::DecodeInfo, false));
    }
  }
}

inline void schedule_level56_rows(std::vector<Level56DispatchEntry>* entries,
                                  const newcode::Params& p,
                                  std::size_t selected_level = 0,
                                  Level56ScheduleSample* sample = nullptr,
                                  std::size_t max_group_entries =
                                      kLevel56MaxGroupEntries,
                                  std::size_t entry_slot_offset = 0,
                                  bool preserve_entry_state = false) {
  if (p.LEVEL56_SCHEDULE_MODE ==
      newcode::Level56ScheduleMode::Group4LoadSortedMultiround) {
    schedule_level56_rows_group4_load_sorted_multiround(
        entries, p, selected_level, sample, max_group_entries,
        entry_slot_offset, preserve_entry_state);
    return;
  }
  schedule_level56_rows_global_priority(entries, p, selected_level);
}

inline void route_level56_g1_global_priority(
    std::vector<Level56DispatchEntry>* entries,
    const newcode::Params& p,
    std::size_t selected_level = 0) {
  const auto ordered = order_level56_candidates(
      *entries, p.LEVEL56_PRIORITY_MODE, selected_level);
  int hiso_core = 0;
  int siso_core = 0;
  for (std::size_t index : ordered) {
    auto& entry = (*entries)[index];
    if (entry.final_action == Level56FinalAction::HisoDecode) {
      entry.assigned_core = hiso_core++;
    } else if (entry.final_action == Level56FinalAction::SisoDecode) {
      entry.assigned_core = siso_core++;
    }
  }
  if (hiso_core > p.LEVEL56_SHARED_HISO_ACTIVE ||
      siso_core > p.LEVEL56_SHARED_SISO_ACTIVE) {
    throw std::logic_error("LEVEL56 G=1 routing exceeded shared capacity");
  }
}

inline void route_level56_g1(std::vector<Level56DispatchEntry>* entries,
                             const newcode::Params& p,
                             std::size_t selected_level = 0) {
  if (p.LEVEL56_SCHEDULE_MODE ==
      newcode::Level56ScheduleMode::Group4LoadSortedMultiround) {
    for (const auto& entry : *entries) {
      if ((entry.planned_hiso || entry.planned_siso) &&
          (entry.assigned_entry_slot < 0 ||
           entry.assigned_core != entry.assigned_entry_slot)) {
        throw std::logic_error(
            "LEVEL56 grouped routing lost its entry-slot assignment");
      }
    }
    return;
  }
  route_level56_g1_global_priority(entries, p, selected_level);
}

template <typename LLR>
Level56EarlyStopEvaluation detect_level56_early_stop(
    const TilePrepared<LLR>& prep,
    const newcode::Params& p) {
  Level56EarlyStopEvaluation evaluation;
  if (p.ENABLE_EARLY_STOP) {
    evaluation.raw = detect_tile_early_stop(prep.lin_matrix, p);
    evaluation.effective = evaluation.raw;
    if (p.EARLY_STOP_CONDITION_MODE == 1 &&
        p.EARLY_STOP_BIND_GROUP_SIZE > 1) {
      evaluation.effective = apply_group_bound_early_stop(
          evaluation.raw, p.EARLY_STOP_BIND_GROUP_SIZE);
    }
  } else {
    const std::size_t rows = prep.lin_matrix.rows();
    evaluation.raw.row_passed_flags.assign(rows, false);
    evaluation.raw.row_details.assign(rows, TileEarlyStopRowDetail{});
    evaluation.raw.rows_passed = 0;
    evaluation.raw.rows_total = rows;
    evaluation.raw.all_rows_passed = false;
    evaluation.effective = evaluation.raw;
  }
  return evaluation;
}

template <typename LLR>
void append_level56_entries(const TilePrepared<LLR>& prep,
                            const TileEarlyStopResult& early_stop,
                            std::size_t source_level,
                            std::vector<Level56DispatchEntry>* entries) {
  using CoreLLR = typename TilePrepared<LLR>::CoreLLR;
  for (std::size_t row = 0; row < prep.lin_matrix.rows(); ++row) {
    Level56DispatchEntry entry;
    entry.shared_row = entries->size();
    entry.source_level = source_level;
    entry.source_local_row = row;
    entry.source_global_row = prep.row_global_lookup[row];
    entry.early_stop_hit =
        row < early_stop.row_passed_flags.size() &&
        early_stop.row_passed_flags[row];
    if (entry.early_stop_hit) {
      entry.final_action = Level56FinalAction::EarlyStopAction;
      entries->push_back(entry);
      continue;
    }

    std::array<CoreLLR, newcode::Params::BCH_N> lin{};
    for (std::size_t col = 0; col < lin.size(); ++col) {
      lin[col] = prep.lin_matrix[row][col];
    }
    HybridRowClass row_class = HybridRowClass::HardFail;
    (void)run_selected_hybrid_classifier_classify_only(
        lin, prep.params_for_core, &row_class);
    if (row_class == HybridRowClass::Clean) {
      std::ostringstream oss;
      oss << "LEVEL56 classify-only returned Clean after early-stop: level="
          << source_level << ", source_local_row=" << row;
      throw std::runtime_error(oss.str());
    }
    entry.hybrid_class = row_class;
    entry.eligibility = level56_hiso_eligible(row_class)
                            ? Level56Eligibility::HisoOrSiso
                            : Level56Eligibility::SisoOnly;
    entries->push_back(entry);
  }
}

template <typename LLR>
void append_level56_early_stop_only_entries(
    const TilePrepared<LLR>& prep,
    const TileEarlyStopResult& early_stop,
    std::size_t source_level,
    std::vector<Level56DispatchEntry>* entries) {
  for (std::size_t row = 0; row < prep.lin_matrix.rows(); ++row) {
    Level56DispatchEntry entry;
    entry.shared_row = entries->size();
    entry.source_level = source_level;
    entry.source_local_row = row;
    entry.source_global_row = prep.row_global_lookup[row];
    entry.early_stop_hit =
        row < early_stop.row_passed_flags.size() &&
        early_stop.row_passed_flags[row];
    // Future lookahead exposes only the row-level EarlyStop decision. It does
    // not classify, reserve a resource, or execute EarlyStopAction yet.
    entry.remaining_for_schedule = !entry.early_stop_hit;
    entry.needs_execution = false;
    entries->push_back(entry);
  }
}

template <typename LLR>
void append_level56_bypassed_entries(
    const TilePrepared<LLR>& prep,
    const TileEarlyStopResult& early_stop,
    bool preserve_early_stop_action,
    std::size_t source_level,
    std::vector<Level56DispatchEntry>* entries) {
  for (std::size_t row = 0; row < prep.lin_matrix.rows(); ++row) {
    Level56DispatchEntry entry;
    entry.shared_row = entries->size();
    entry.source_level = source_level;
    entry.source_local_row = row;
    entry.source_global_row = prep.row_global_lookup[row];
    entry.early_stop_hit =
        preserve_early_stop_action &&
        row < early_stop.row_passed_flags.size() &&
        early_stop.row_passed_flags[row];
    entry.hybrid_class = HybridRowClass::None;
    entry.eligibility = Level56Eligibility::None;
    entry.final_action = entry.early_stop_hit
                             ? Level56FinalAction::EarlyStopAction
                             : Level56FinalAction::Unscheduled;
    entries->push_back(entry);
  }
}

template <typename LLR>
chase::DecoderCoreResult<typename TilePrepared<LLR>::CoreLLR>
execute_level56_slice(
    const TilePrepared<LLR>& prep,
    const std::vector<Level56DispatchEntry>& entries,
    std::size_t source_level,
    bool normalize_extrinsic,
    CoreFn<typename TilePrepared<LLR>::CoreLLR> core_fn) {
  using CoreLLR = typename TilePrepared<LLR>::CoreLLR;
  const std::size_t rows = prep.lin_matrix.rows();
  const std::size_t cols = prep.lin_matrix.cols();
  chase::DecoderCoreResult<CoreLLR> result{
      matrix::Matrix<float>(rows, cols), std::vector<bool>(rows, false)};
  std::vector<uint8_t> soft_state(
      rows, static_cast<uint8_t>(newcode::mux::StateTag::Unscheduled));
  std::vector<bool> normalize_rows(rows, false);

  for (const auto& entry : entries) {
    if (!entry.needs_execution || entry.source_level != source_level ||
        entry.final_action == Level56FinalAction::Unscheduled) {
      continue;
    }
    const std::size_t row = entry.source_local_row;
    std::array<CoreLLR, newcode::Params::BCH_N> lin{};
    std::array<CoreLLR, newcode::Params::BCH_N> lch{};
    std::array<float, newcode::Params::BCH_N> y2{};
    load_row_vectors(prep, row, &lin, &lch);

    if (entry.final_action == Level56FinalAction::EarlyStopAction) {
      if (newcode::apply_row_early_stop_action(
              lin.data(), lch.data(), y2.data(), prep.params_for_core)) {
        result.produced_rows[row] = true;
        for (std::size_t col = 0; col < cols; ++col) {
          result.lout[row][col] = y2[col];
        }
      }
    } else if (entry.final_action == Level56FinalAction::HisoDecode) {
      try {
        execute_hybrid_hard_class(entry.hybrid_class, lin,
                                  prep.params_for_core, &y2);
      } catch (const std::exception& ex) {
        std::ostringstream oss;
        oss << "LEVEL56 HISO executor failed: level=" << source_level
            << ", source_local_row=" << row << ": " << ex.what();
        throw std::runtime_error(oss.str());
      }
      result.produced_rows[row] = true;
      for (std::size_t col = 0; col < cols; ++col) {
        result.lout[row][col] = y2[col];
      }
    } else if (entry.final_action == Level56FinalAction::SisoDecode) {
      soft_state[row] = static_cast<uint8_t>(newcode::mux::StateTag::NeedSiso);
      normalize_rows[row] = true;
    }
  }

  if (std::any_of(normalize_rows.begin(), normalize_rows.end(),
                  [](bool selected) { return selected; })) {
    const auto soft_result = run_decoder_core(
        prep.lin_matrix, prep.lch_matrix, false, prep.params_for_core,
        nullptr, &soft_state, core_fn);
    for (std::size_t row = 0; row < rows; ++row) {
      if (!normalize_rows[row] || row >= soft_result.produced_rows.size() ||
          !soft_result.produced_rows[row]) {
        continue;
      }
      result.produced_rows[row] = true;
      for (std::size_t col = 0; col < cols; ++col) {
        result.lout[row][col] = soft_result.lout[row][col];
      }
    }
  }

  postprocess_decoder_result<LLR>(&result, normalize_extrinsic, false,
                                  prep.params_for_core, &normalize_rows);
  log_target_history(prep.params_for_core, result.lout,
                     result.produced_rows);
  return result;
}

template <typename LLR>
TileProcessResult<LLR> build_level56_tile_result(
    const matrix::Matrix<LLR>& tile_in,
    const TileEarlyStopResult& raw_early_stop,
    const TileEarlyStopResult& early_stop,
    const std::vector<Level56DispatchEntry>& entries,
    std::size_t source_level,
    std::size_t shared_invocation,
    const newcode::Params& params,
    matrix::Matrix<LLR> tile_out) {
  std::size_t hard = 0;
  std::size_t soft_candidates = 0;
  std::size_t unscheduled = 0;
  HybridClassCount class_count{};
  class_count.invocation = shared_invocation;
  class_count.tile_index = source_level == 5 ? kLevel5TileIndex
                                             : kLevel6TileIndex;
  for (const auto& entry : entries) {
    if (entry.source_level != source_level || entry.early_stop_hit) {
      continue;
    }
    ++soft_candidates;
    ++class_count.rows_seen_by_hybrid;
    increment_hybrid_class_count(entry.hybrid_class, &class_count);
    hard += entry.final_action == Level56FinalAction::HisoDecode;
    unscheduled += entry.final_action == Level56FinalAction::Unscheduled;
  }
  (void)tile_in;
  auto result = TileProcessResult<LLR>{
      std::move(tile_out), early_stop.all_rows_passed,
      early_stop.rows_passed, early_stop.rows_total, hard,
      soft_candidates, unscheduled};
  result.hybrid_class_count = std::move(class_count);
  if (params.ENABLE_EARLY_STOP &&
      params.EARLY_STOP_CONDITION_MODE == 1 &&
      params.EARLY_STOP_BIND_GROUP_SIZE > 1) {
    result.has_group_bind_debug_sample = true;
    result.group_bind_debug_sample.condition_mode =
        params.EARLY_STOP_CONDITION_MODE;
    result.group_bind_debug_sample.bind_group_size =
        params.EARLY_STOP_BIND_GROUP_SIZE;
    result.group_bind_debug_sample.raw_early_stop_flags =
        early_stop_flags_to_bitstring(raw_early_stop.row_passed_flags);
    result.group_bind_debug_sample.bound_early_stop_flags =
        early_stop_flags_to_bitstring(early_stop.row_passed_flags);
  }
  return result;
}

template <typename LLR>
Level56SharedResult<LLR> process_level56_shared(
    const matrix::Matrix<LLR>& tile_in5,
    const matrix::Matrix<LLR>& ch_tile5,
    const matrix::Matrix<LLR>& tile_in6,
    const matrix::Matrix<LLR>& ch_tile6,
    const newcode::Params& params5,
    const newcode::Params& params6,
    std::size_t tile_top5,
    std::size_t tile_top6,
    std::size_t shared_invocation,
    bool normalize_extrinsic,
    const matrix::Matrix<float>* tx_llr_ref,
    CoreFn<typename LinMatrixAdapter<LLR>::core_type> core_fn,
    matrix::Matrix<float>* last_tile_history_accum,
    std::vector<Level56DispatchEntry>* prior_entries = nullptr,
    std::size_t max_group_entries = kLevel56MaxGroupEntries,
    std::size_t entry_slot_offset = 0,
    bool preserve_entry_state = false,
    const Level56EarlyStopEvaluation* cached_early5 = nullptr,
    const Level56EarlyStopEvaluation* cached_early6 = nullptr,
    matrix::Matrix<LLR>* output_tile5 = nullptr,
    matrix::Matrix<LLR>* output_tile6 = nullptr) {
  const std::size_t rows = static_cast<std::size_t>(params5.CHASE_SBR) *
                           newcode::Params::BITS_PER_SUBBLOCK_DIM;
  auto prep5 = prepare_tile_inputs(tile_in5, ch_tile5, params5, tile_top5,
                                   params5.CHASE_SBR, rows, tx_llr_ref);
  auto prep6 = prepare_tile_inputs(tile_in6, ch_tile6, params6, tile_top6,
                                   params6.CHASE_SBR, rows, tx_llr_ref);
  const auto early5 = cached_early5 ? *cached_early5
                                   : detect_level56_early_stop(prep5, params5);
  const auto early6 = cached_early6 ? *cached_early6
                                   : detect_level56_early_stop(prep6, params6);
  const std::size_t selected_level = select_level56_decode_level(
      early5.effective, early6.effective, params5);

  std::vector<Level56DispatchEntry> entries;
  entries.reserve(rows * 2u);
  const bool preserve_unselected_early_stop_action =
      params5.LEVEL56_UNSELECTED_EARLY_STOP_ACTION_ENABLE;
  if (prior_entries) {
    if (prior_entries->size() != rows * 2u) {
      throw std::invalid_argument(
          "LEVEL56 temporal history requires exactly 64 cached entries");
    }
    entries = *prior_entries;
    for (std::size_t index = 0; index < rows; ++index) {
      entries[index].source_global_row = prep5.row_global_lookup[index];
      entries[rows + index].source_global_row = prep6.row_global_lookup[index];
    }
    auto refresh_pending = [](const TilePrepared<LLR>& prep,
                              std::size_t begin,
                              std::vector<Level56DispatchEntry>* values) {
      for (std::size_t row = 0; row < prep.lin_matrix.rows(); ++row) {
        auto& entry = (*values)[begin + row];
        const bool pending = !entry.early_stop_hit &&
                             entry.final_action == Level56FinalAction::Unscheduled &&
                             !entry.produced;
        entry.planned_hiso = false;
        entry.planned_siso = false;
        entry.remaining_for_schedule = pending;
        entry.needs_execution = false;
        if (!pending) {
          continue;
        }
        std::array<typename TilePrepared<LLR>::CoreLLR,
                   newcode::Params::BCH_N> lin{};
        for (std::size_t col = 0; col < lin.size(); ++col) {
          lin[col] = prep.lin_matrix[row][col];
        }
        HybridRowClass row_class = HybridRowClass::HardFail;
        (void)run_selected_hybrid_classifier_classify_only(
            lin, prep.params_for_core, &row_class);
        if (row_class == HybridRowClass::Clean) {
          throw std::runtime_error(
              "LEVEL56 temporal history classify-only returned Clean");
        }
        entry.hybrid_class = row_class;
        entry.eligibility = level56_hiso_eligible(row_class)
                                ? Level56Eligibility::HisoOrSiso
                                : Level56Eligibility::SisoOnly;
      }
    };
    refresh_pending(prep5, 0, &entries);
    refresh_pending(prep6, rows, &entries);
  } else {
    const bool preserve_unselected_early_stop_action =
        params5.LEVEL56_UNSELECTED_EARLY_STOP_ACTION_ENABLE;
    if (selected_level != 0 && selected_level != 5) {
      append_level56_bypassed_entries(
          prep5, early5.effective, preserve_unselected_early_stop_action,
          5, &entries);
    } else {
      append_level56_entries(prep5, early5.effective, 5, &entries);
    }
    if (selected_level != 0 && selected_level != 6) {
      append_level56_bypassed_entries(
          prep6, early6.effective, preserve_unselected_early_stop_action,
          6, &entries);
    } else {
      append_level56_entries(prep6, early6.effective, 6, &entries);
    }
    for (auto& entry : entries) {
      entry.produced = false;
      entry.needs_execution = true;
    }
  }
  Level56ScheduleSample schedule_sample;
  Level56ScheduleSample* schedule_sample_ptr =
      params5.LEVEL56_SCHEDULE_OBSERVABILITY_ENABLE
          ? &schedule_sample
          : nullptr;
  schedule_level56_rows(
      &entries, params5, selected_level, schedule_sample_ptr,
      max_group_entries, entry_slot_offset, preserve_entry_state);
  route_level56_g1(&entries, params5, selected_level);

  for (auto& entry : entries) {
    if (preserve_entry_state) {
      entry.needs_execution = entry.planned_hiso || entry.planned_siso;
    } else {
      entry.needs_execution = entry.final_action != Level56FinalAction::Unscheduled;
    }
  }

  using SharedCoreLLR = typename TilePrepared<LLR>::CoreLLR;
  chase::DecoderCoreResult<SharedCoreLLR> decoded5{
      matrix::Matrix<float>(prep5.lin_matrix.rows(), prep5.lin_matrix.cols()),
      std::vector<bool>(prep5.lin_matrix.rows(), false)};
  if (selected_level == 0 || selected_level == 5 ||
      preserve_unselected_early_stop_action) {
    decoded5 = execute_level56_slice(
        prep5, entries, 5, normalize_extrinsic, core_fn);
  }
  chase::DecoderCoreResult<SharedCoreLLR> decoded6{
      matrix::Matrix<float>(prep6.lin_matrix.rows(), prep6.lin_matrix.cols()),
      std::vector<bool>(prep6.lin_matrix.rows(), false)};
  if (selected_level == 0 || selected_level == 6 ||
      preserve_unselected_early_stop_action) {
    decoded6 = execute_level56_slice(
        prep6, entries, 6, normalize_extrinsic, core_fn);
  }
  for (auto& entry : entries) {
    if (!entry.needs_execution) {
      continue;
    }
    const auto& produced_rows =
        entry.source_level == 5 ? decoded5.produced_rows : decoded6.produced_rows;
    entry.produced = entry.source_local_row < produced_rows.size() &&
                     produced_rows[entry.source_local_row];
    entry.decode_status = entry.produced
                              ? Level56DecodeStatus::Produced
                              : Level56DecodeStatus::ActionFailed;
  }
  matrix::Matrix<LLR> tile_out5 = output_tile5 ? *output_tile5 : tile_in5;
  matrix::Matrix<LLR> tile_out6 = output_tile6 ? *output_tile6 : tile_in6;
  writeback_tile(prep5, decoded5, params5, tile_top5, false,
                 &tile_out5, nullptr);
  writeback_tile(prep6, decoded6, params6, tile_top6, true,
                 &tile_out6, last_tile_history_accum);

  Level56SharedResult<LLR> result;
  result.level5 = build_level56_tile_result(
      tile_in5, early5.raw, early5.effective, entries, 5, shared_invocation,
      params5, std::move(tile_out5));
  result.level6 = build_level56_tile_result(
      tile_in6, early6.raw, early6.effective, entries, 6, shared_invocation,
      params6, std::move(tile_out6));
  result.dispatch = std::move(entries);
  if (schedule_sample_ptr) {
    for (auto& code : schedule_sample.codes) {
      code.produced = result.dispatch[code.code_index].produced;
      code.decode_status = result.dispatch[code.code_index].decode_status;
    }
    schedule_sample.invocation = shared_invocation;
    result.has_schedule_sample = true;
    result.schedule_sample = std::move(schedule_sample);
  }
  return result;
}

}  // namespace detail
}  // namespace newcode
