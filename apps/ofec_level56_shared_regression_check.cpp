#include "ofec/detail/ofec_tile_impl.ipp"
#include "ofec/detail/ofec_level56_shared.ipp"

#include <algorithm>
#include <cmath>
#include <iostream>
#include <stdexcept>
#include <vector>

namespace {

using newcode::detail::Level56DispatchEntry;
using newcode::detail::Level56Eligibility;
using newcode::detail::Level56FinalAction;

std::vector<float> observed_core_betas;

void require(bool condition, const char* message) {
  if (!condition) {
    throw std::runtime_error(message);
  }
}

chase::DecoderCoreResult<float> controlled_soft_core(
    const matrix::Matrix<float>& lin_matrix,
    const matrix::Matrix<float>&,
    bool,
    const newcode::Params& p,
    const std::vector<bool>*,
    const std::vector<uint8_t>* mux_state) {
  observed_core_betas.push_back(p.beta);
  chase::DecoderCoreResult<float> result{
      matrix::Matrix<float>(lin_matrix.rows(), lin_matrix.cols()),
      std::vector<bool>(lin_matrix.rows(), false)};
  for (std::size_t row = 0; row < lin_matrix.rows(); ++row) {
    if (!mux_state || row >= mux_state->size() ||
        (*mux_state)[row] !=
            static_cast<uint8_t>(newcode::mux::StateTag::NeedSiso)) {
      continue;
    }
    result.produced_rows[row] = true;
    for (std::size_t col = 0; col < lin_matrix.cols(); ++col) {
      result.lout[row][col] = p.beta;
    }
  }
  return result;
}

std::vector<Level56DispatchEntry> make_entries(std::size_t rows_per_level) {
  std::vector<Level56DispatchEntry> entries;
  entries.reserve(rows_per_level * 2u);
  for (std::size_t level : {5u, 6u}) {
    for (std::size_t row = 0; row < rows_per_level; ++row) {
      entries.push_back(Level56DispatchEntry{
          .shared_row = entries.size(),
          .source_level = level,
          .source_local_row = row,
          .source_global_row = row,
          .early_stop_hit = false,
          .hybrid_class = newcode::detail::HybridRowClass::ParityOnly,
          .eligibility = Level56Eligibility::HisoOrSiso,
          .final_action = Level56FinalAction::Unscheduled,
          .assigned_core = -1,
      });
    }
  }
  return entries;
}

newcode::Params make_grouped_params() {
  newcode::Params p;
  p.LEVEL56_SHARED_ENABLE = true;
  p.LEVEL56_SHARED_HISO_ACTIVE = 8;
  p.LEVEL56_SHARED_SISO_ACTIVE = 8;
  p.LEVEL56_PRIORITY_MODE = newcode::Level56PriorityMode::Level5First;
  p.LEVEL56_SCHEDULE_MODE =
      newcode::Level56ScheduleMode::Group4LoadSortedMultiround;
  p.LEVEL56_SCHEDULE_OBSERVABILITY_ENABLE = true;
  p.LEVEL56_SINGLE_LEVEL_SELECT_ENABLE = false;
  p.HYBRID_CLASSIFIER_MODE =
      newcode::HybridClassifierMode::FriendS1S3WithS0Classifier;
  return p;
}

std::vector<Level56DispatchEntry> make_all_early_stop_grouped_entries() {
  auto entries = make_entries(32);
  for (auto& entry : entries) {
    entry.early_stop_hit = true;
    entry.hybrid_class = newcode::detail::HybridRowClass::None;
    entry.eligibility = Level56Eligibility::None;
    entry.final_action = Level56FinalAction::EarlyStopAction;
  }
  return entries;
}

void activate_grouped_row(
    std::vector<Level56DispatchEntry>* entries,
    std::size_t shared_row,
    newcode::detail::HybridRowClass row_class,
    Level56Eligibility eligibility = Level56Eligibility::HisoOrSiso) {
  auto& entry = entries->at(shared_row);
  entry.early_stop_hit = false;
  entry.hybrid_class = row_class;
  entry.eligibility = eligibility;
  entry.final_action = Level56FinalAction::Unscheduled;
}

void activate_group_prefix(
    std::vector<Level56DispatchEntry>* entries,
    std::size_t group,
    std::size_t count,
    newcode::detail::HybridRowClass row_class =
        newcode::detail::HybridRowClass::ParityOnly,
    Level56Eligibility eligibility = Level56Eligibility::HisoOrSiso) {
  for (std::size_t position = 0; position < count; ++position) {
    activate_grouped_row(entries, group * 4u + position, row_class,
                         eligibility);
  }
}

int grouped_entry_count(const std::vector<Level56DispatchEntry>& entries) {
  int max_slot = -1;
  for (const auto& entry : entries) {
    max_slot = std::max(max_slot, entry.assigned_entry_slot);
  }
  return max_slot + 1;
}

int group_for_slot(const std::vector<Level56DispatchEntry>& entries,
                   int slot) {
  int group = -1;
  for (const auto& entry : entries) {
    if (entry.assigned_entry_slot != slot) {
      continue;
    }
    const int entry_group = static_cast<int>(entry.shared_row / 4u);
    require(group == -1 || group == entry_group,
            "one entry slot must not select rows from multiple groups");
    group = entry_group;
  }
  return group;
}

void require_grouped_plan_invariants(
    const std::vector<Level56DispatchEntry>& entries) {
  require(grouped_entry_count(entries) <= 8,
          "grouped scheduler must use at most eight entry slots");
  for (const auto& entry : entries) {
    require(!(entry.planned_hiso && entry.planned_siso),
            "one row must not be reserved by both paths");
    if (entry.planned_hiso || entry.planned_siso) {
      require(entry.assigned_entry_slot >= 0 &&
                  entry.assigned_core == entry.assigned_entry_slot,
              "planned row core must equal its entry slot");
      require(!entry.remaining_for_schedule,
              "planned row must be removed from later rounds");
    }
  }
  for (int slot = 0; slot < grouped_entry_count(entries); ++slot) {
    require(group_for_slot(entries, slot) >= 0,
            "every consumed entry slot must select at least one row");
  }
}

std::size_t count_action(const std::vector<Level56DispatchEntry>& entries,
                         Level56FinalAction action) {
  return static_cast<std::size_t>(std::count_if(
      entries.begin(), entries.end(),
      [&](const auto& entry) { return entry.final_action == action; }));
}

void require_state_conservation(
    const std::vector<Level56DispatchEntry>& entries) {
  const std::size_t total =
      count_action(entries, Level56FinalAction::EarlyStopAction) +
      count_action(entries, Level56FinalAction::HisoDecode) +
      count_action(entries, Level56FinalAction::SisoDecode) +
      count_action(entries, Level56FinalAction::Unscheduled);
  require(total == entries.size(),
          "every shared row must have exactly one final action");
}

void check_full_siso_uses_no_hiso() {
  auto entries = make_entries(32);
  newcode::Params p;
  p.LEVEL56_SHARED_HISO_ACTIVE = 64;
  p.LEVEL56_SHARED_SISO_ACTIVE = 64;
  newcode::detail::schedule_level56_rows(&entries, p);
  newcode::detail::route_level56_g1(&entries, p);
  require(count_action(entries, Level56FinalAction::HisoDecode) == 0,
          "full SISO capacity must not reclaim HISO rows");
  require(count_action(entries, Level56FinalAction::SisoDecode) == 64,
          "full SISO capacity must schedule all rows on SISO");
  require_state_conservation(entries);
}

void check_hiso_only_handles_siso_overflow() {
  auto entries = make_entries(32);
  newcode::Params p;
  p.LEVEL56_SHARED_HISO_ACTIVE = 8;
  p.LEVEL56_SHARED_SISO_ACTIVE = 48;
  newcode::detail::schedule_level56_rows(&entries, p);
  newcode::detail::route_level56_g1(&entries, p);
  require(count_action(entries, Level56FinalAction::HisoDecode) == 8,
          "HISO reclaim must be capped by HISO capacity");
  require(count_action(entries, Level56FinalAction::SisoDecode) == 48,
          "SISO scheduling must use configured capacity");
  require(count_action(entries, Level56FinalAction::Unscheduled) == 8,
          "remaining overflow rows must be unscheduled");
  require_state_conservation(entries);
}

void check_siso_only_rows_never_use_hiso() {
  auto entries = make_entries(4);
  entries[0].hybrid_class = newcode::detail::HybridRowClass::HardFail;
  entries[0].eligibility = Level56Eligibility::SisoOnly;
  newcode::Params p;
  p.LEVEL56_SHARED_HISO_ACTIVE = 8;
  p.LEVEL56_SHARED_SISO_ACTIVE = 1;
  newcode::detail::schedule_level56_rows(&entries, p);
  newcode::detail::route_level56_g1(&entries, p);
  require(entries[0].final_action != Level56FinalAction::HisoDecode,
          "SisoOnly row must never be reclaimed to HISO");
  require_state_conservation(entries);
}

void check_early_stop_does_not_consume_shared_capacity() {
  auto entries = make_entries(4);
  entries[0].early_stop_hit = true;
  entries[0].eligibility = Level56Eligibility::None;
  entries[0].final_action = Level56FinalAction::EarlyStopAction;
  newcode::Params p;
  p.LEVEL56_SHARED_HISO_ACTIVE = 0;
  p.LEVEL56_SHARED_SISO_ACTIVE = 1;
  newcode::detail::schedule_level56_rows(&entries, p);
  newcode::detail::route_level56_g1(&entries, p);
  require(entries[0].final_action == Level56FinalAction::EarlyStopAction,
          "early-stop row must keep its action through shared scheduling");
  require(entries[0].assigned_core == -1,
          "early-stop row must not be routed to a shared core");
  require(count_action(entries, Level56FinalAction::SisoDecode) == 1,
          "early-stop row must not consume SISO capacity");
  require_state_conservation(entries);
}

void check_single_level_scheduler_bypasses_other_level() {
  auto entries = make_entries(4);
  newcode::Params p;
  p.LEVEL56_SHARED_HISO_ACTIVE = 1;
  p.LEVEL56_SHARED_SISO_ACTIVE = 1;
  newcode::detail::schedule_level56_rows(&entries, p, 6);
  newcode::detail::route_level56_g1(&entries, p, 6);

  for (const auto& entry : entries) {
    if (entry.source_level == 5) {
      require(entry.final_action == Level56FinalAction::Unscheduled,
              "single-level mode must bypass every unselected row");
      require(entry.assigned_core == -1,
              "bypassed Level 5 row must not be routed to a shared core");
    }
  }
  require(count_action(entries, Level56FinalAction::HisoDecode) == 1,
          "single-level mode must allocate HISO only to the selected level");
  require(count_action(entries, Level56FinalAction::SisoDecode) == 1,
          "single-level mode must allocate SISO only to the selected level");
  require_state_conservation(entries);
}

void check_single_level_preserves_unselected_early_stop_actions() {
  auto entries = make_entries(4);
  entries[0].early_stop_hit = true;
  entries[0].eligibility = Level56Eligibility::None;
  entries[0].final_action = Level56FinalAction::EarlyStopAction;

  newcode::Params p;
  p.LEVEL56_UNSELECTED_EARLY_STOP_ACTION_ENABLE = true;
  p.LEVEL56_SHARED_HISO_ACTIVE = 1;
  p.LEVEL56_SHARED_SISO_ACTIVE = 1;
  newcode::detail::schedule_level56_rows(&entries, p, 6);
  newcode::detail::route_level56_g1(&entries, p, 6);

  require(entries[0].final_action == Level56FinalAction::EarlyStopAction,
          "unselected early-stop hit must retain its action");
  require(entries[0].assigned_core == -1,
          "unselected early-stop action must not consume a shared core");
  for (std::size_t index = 1; index < 4; ++index) {
    require(entries[index].final_action == Level56FinalAction::Unscheduled,
            "unselected non-hit row must remain unscheduled");
  }
  require(count_action(entries, Level56FinalAction::HisoDecode) == 1,
          "preserved early-stop action must not consume HISO capacity");
  require(count_action(entries, Level56FinalAction::SisoDecode) == 1,
          "preserved early-stop action must not consume SISO capacity");
  require_state_conservation(entries);
}

void check_bypassed_slice_executes_only_preserved_early_stop_actions() {
  newcode::detail::TilePrepared<float> prep;
  prep.lin_matrix = matrix::Matrix<float>(2, newcode::Params::BCH_N);
  prep.lch_matrix = matrix::Matrix<float>(2, newcode::Params::BCH_N);
  prep.row_local_lookup = {0, 1};
  prep.row_global_lookup = {0, 1};
  prep.params_for_core.EARLY_STOP_ACTION_MODE = 1;
  prep.params_for_core.EARLY_STOP_ACTION_SIGN_BETA = 3.0f;

  newcode::TileEarlyStopResult early_stop;
  early_stop.row_passed_flags = {true, false};
  std::vector<Level56DispatchEntry> entries;
  newcode::detail::append_level56_bypassed_entries(
      prep, early_stop, true, 5, &entries);

  require(entries.size() == 2 && entries[0].early_stop_hit &&
              entries[0].final_action == Level56FinalAction::EarlyStopAction,
          "bypassed slice must preserve the detected early-stop hit");
  require(!entries[1].early_stop_hit &&
              entries[1].final_action == Level56FinalAction::Unscheduled,
          "bypassed slice must leave a non-hit row unscheduled");

  observed_core_betas.clear();
  const auto decoded = newcode::detail::execute_level56_slice(
      prep, entries, 5, false, &controlled_soft_core);
  require(decoded.produced_rows == std::vector<bool>({true, false}),
          "only the preserved early-stop row may produce output");
  require(observed_core_betas.empty(),
          "unselected early-stop execution must not call the SISO core");
}

void check_single_level_selection_uses_fewer_early_stops() {
  newcode::TileEarlyStopResult early5;
  newcode::TileEarlyStopResult early6;
  newcode::Params p;
  p.LEVEL56_SINGLE_LEVEL_SELECT_ENABLE = true;

  early5.rows_passed = 3;
  early6.rows_passed = 7;
  require(
      newcode::detail::select_level56_decode_level(early5, early6, p) == 5,
      "single-level selection must choose Level 5 when it has fewer early stops");

  early5.rows_passed = 9;
  early6.rows_passed = 2;
  require(
      newcode::detail::select_level56_decode_level(early5, early6, p) == 6,
      "single-level selection must choose Level 6 when it has fewer early stops");

  early5.rows_passed = 4;
  early6.rows_passed = 4;
  p.LEVEL56_PRIORITY_MODE = newcode::Level56PriorityMode::Level5First;
  require(newcode::detail::select_level56_decode_level(early5, early6, p) == 5,
          "Level5First tie-breaking must select Level 5");

  p.LEVEL56_PRIORITY_MODE = newcode::Level56PriorityMode::Level6First;
  require(newcode::detail::select_level56_decode_level(early5, early6, p) == 6,
          "Level6First tie-breaking must select Level 6");
}

void check_common_parameter_validation() {
  newcode::Params p;
  p.LEVEL56_SHARED_ENABLE = true;
  p.EARLY_STOP_ENABLE_LIST = {1, 1, 1, 1, 1, 1};
  p.HARD_TILE_LIST = {0, 0, 0, 0, 0, 0};
  p.HYBRID_ENABLE_LIST = {0, 0, 0, 0, 1, 1};
  p.HYBRID_CLASSIFIER_MODE =
      newcode::HybridClassifierMode::RepoFastClassifier;
  p.HYBRID_USE_FAST_CLASSIFIER = true;
  newcode::detail::validate_level56_shared_config(p);

  auto invalid_priority = p;
  invalid_priority.LEVEL56_PRIORITY_MODE =
      static_cast<newcode::Level56PriorityMode>(0);
  bool rejected = false;
  try {
    newcode::detail::validate_level56_shared_config(invalid_priority);
  } catch (const std::invalid_argument&) {
    rejected = true;
  }
  require(rejected, "removed priority value 0 must be rejected");

  auto grouped = p;
  grouped.LEVEL56_SHARED_HISO_ACTIVE = 8;
  grouped.LEVEL56_SHARED_SISO_ACTIVE = 8;
  grouped.LEVEL56_SCHEDULE_MODE =
      newcode::Level56ScheduleMode::Group4LoadSortedMultiround;
  grouped.HYBRID_CLASSIFIER_MODE =
      newcode::HybridClassifierMode::FriendS1S3WithS0Classifier;
  newcode::detail::validate_level56_shared_config(grouped);

  auto invalid_grouped_capacity = grouped;
  invalid_grouped_capacity.LEVEL56_SHARED_SISO_ACTIVE = 7;
  rejected = false;
  try {
    newcode::detail::validate_level56_shared_config(
        invalid_grouped_capacity);
  } catch (const std::invalid_argument&) {
    rejected = true;
  }
  require(rejected, "grouped mode must reject capacity other than 8/8");

  p.EARLY_STOP_ACTION_MODE_LIST = {1, 1, 1, 1, 1, 2};
  rejected = false;
  try {
    newcode::detail::validate_level56_shared_config(p);
  } catch (const std::invalid_argument&) {
    rejected = true;
  }
  require(rejected,
          "Level 5/6 common parameter mismatch must be rejected");
}

void check_per_level_siso_postprocessing() {
  auto make_prep = [](float beta, float alpha) {
    newcode::detail::TilePrepared<float> prep;
    prep.lin_matrix = matrix::Matrix<float>(2, newcode::Params::BCH_N);
    prep.lch_matrix = matrix::Matrix<float>(2, newcode::Params::BCH_N);
    prep.row_local_lookup = {0, 1};
    prep.row_global_lookup = {0, 1};
    prep.params_for_core.beta = beta;
    prep.params_for_core.ALPHA = alpha;
    return prep;
  };
  const auto prep5 = make_prep(2.0f, 0.5f);
  const auto prep6 = make_prep(3.0f, 0.25f);
  std::vector<Level56DispatchEntry> entries{
      {.source_level = 5,
       .source_local_row = 0,
       .eligibility = Level56Eligibility::SisoOnly,
       .final_action = Level56FinalAction::SisoDecode},
      {.source_level = 5,
       .source_local_row = 1,
       .eligibility = Level56Eligibility::SisoOnly,
       .final_action = Level56FinalAction::Unscheduled},
      {.source_level = 6,
       .source_local_row = 0,
       .eligibility = Level56Eligibility::SisoOnly,
       .final_action = Level56FinalAction::SisoDecode},
      {.source_level = 6,
       .source_local_row = 1,
       .eligibility = Level56Eligibility::SisoOnly,
       .final_action = Level56FinalAction::Unscheduled},
  };

  observed_core_betas.clear();
  const auto decoded5 = newcode::detail::execute_level56_slice(
      prep5, entries, 5, true, &controlled_soft_core);
  const auto decoded6 = newcode::detail::execute_level56_slice(
      prep6, entries, 6, true, &controlled_soft_core);
  require(observed_core_betas == std::vector<float>({2.0f, 3.0f}),
          "Level 5/6 SISO must use their source beta values");
  require(decoded5.produced_rows[0] && !decoded5.produced_rows[1] &&
              decoded6.produced_rows[0] && !decoded6.produced_rows[1],
          "only scheduled SISO rows may be normalized and produced");
  require(std::fabs(decoded5.lout[0][0] - 1.0f) < 1e-6f &&
              std::fabs(decoded6.lout[0][0] - 0.75f) < 1e-6f,
          "Level 5/6 SISO outputs must use their source alpha values");
}

void check_unscheduled_last_tile_history_passes_through_prior() {
  newcode::detail::TilePrepared<float> prep;
  prep.lin_matrix = matrix::Matrix<float>(1, newcode::Params::BCH_N);
  prep.lch_matrix = matrix::Matrix<float>(1, newcode::Params::BCH_N);
  prep.row_local_lookup = {0};
  prep.row_global_lookup = {352};
  chase::DecoderCoreResult<float> decoded{
      matrix::Matrix<float>(1, newcode::Params::BCH_N), {false}};
  matrix::Matrix<float> tile_out(352, 128);
  matrix::Matrix<float> history(704, 128);
  for (std::size_t row = 0; row < tile_out.rows(); ++row) {
    for (std::size_t col = 0; col < tile_out.cols(); ++col) {
      tile_out[row][col] = 3.0f;
    }
  }
  for (std::size_t row = 0; row < history.rows(); ++row) {
    std::fill(history[row].begin(), history[row].end(), 7.0f);
  }

  newcode::Params p;
  newcode::detail::writeback_tile(
      prep, decoded, p, 0, true, &tile_out, &history);
  std::size_t passthrough_count = 0;
  for (std::size_t row = 0; row < history.rows(); ++row) {
    for (std::size_t col = 0; col < history.cols(); ++col) {
      if (history[row][col] == 3.0f) {
        ++passthrough_count;
      } else {
        require(history[row][col] == 7.0f,
                "last-tile passthrough changed an unrelated history cell");
      }
    }
  }
  require(passthrough_count == 128,
          "an unscheduled last-tile row must pass all 128 prior values to history");
  for (std::size_t row = 0; row < tile_out.rows(); ++row) {
    for (std::size_t col = 0; col < tile_out.cols(); ++col) {
      require(tile_out[row][col] == 3.0f,
              "history passthrough must not modify an unscheduled tile row");
    }
  }

  matrix::Matrix<float> non_last_history(704, 128);
  for (std::size_t row = 0; row < non_last_history.rows(); ++row) {
    std::fill(non_last_history[row].begin(), non_last_history[row].end(),
              7.0f);
  }
  newcode::detail::writeback_tile(
      prep, decoded, p, 0, false, &tile_out, &non_last_history);
  for (std::size_t row = 0; row < non_last_history.rows(); ++row) {
    for (std::size_t col = 0; col < non_last_history.cols(); ++col) {
      require(non_last_history[row][col] == 7.0f,
              "a non-last tile must not update history");
    }
  }
}

void check_level_priority_modes() {
  for (const auto& [mode, expected_level] : {
           std::pair{newcode::Level56PriorityMode::Level5First, 5u},
           std::pair{newcode::Level56PriorityMode::Level6First, 6u}}) {
    auto entries = make_entries(4);
    newcode::Params p;
    p.LEVEL56_PRIORITY_MODE = mode;
    p.LEVEL56_SHARED_HISO_ACTIVE = 1;
    p.LEVEL56_SHARED_SISO_ACTIVE = 7;
    newcode::detail::schedule_level56_rows(&entries, p);
    newcode::detail::route_level56_g1(&entries, p);
    const auto selected = std::find_if(
        entries.begin(), entries.end(), [](const auto& entry) {
          return entry.final_action == Level56FinalAction::HisoDecode;
        });
    require(selected != entries.end() && selected->source_level == expected_level,
            "level-first mode selected the wrong source level");
    require(selected->source_local_row == 0,
            "level-first mode must use ascending source_local_row");
  }
}

void check_grouped_k0_uses_no_entries() {
  auto entries = make_all_early_stop_grouped_entries();
  const auto p = make_grouped_params();
  newcode::detail::schedule_level56_rows(&entries, p);
  newcode::detail::route_level56_g1(&entries, p);

  require(grouped_entry_count(entries) == 0,
          "K=0 must not consume an entry slot");
  require(count_action(entries, Level56FinalAction::EarlyStopAction) == 64,
          "K=0 must preserve all EarlyStopAction rows");
  require_grouped_plan_invariants(entries);
}

void check_grouped_k0_suppresses_early_stop_without_entries() {
  auto entries = make_all_early_stop_grouped_entries();
  auto p = make_grouped_params();
  p.LEVEL56_SCHEDULE_OBSERVABILITY_ENABLE = false;
  p.LEVEL56_EARLY_STOP_GROUP_UPDATE_MODE =
      newcode::Level56EarlyStopGroupUpdateMode::EnteredGroupsOnly;
  newcode::detail::schedule_level56_rows(&entries, p);

  require(grouped_entry_count(entries) == 0,
          "K=0 suppression mode must not consume an entry slot");
  require(count_action(entries, Level56FinalAction::EarlyStopAction) == 0 &&
              count_action(entries, Level56FinalAction::Unscheduled) == 64,
          "K=0 suppression mode must leave all early-stop rows unchanged");
  for (const auto& entry : entries) {
    require(entry.early_stop_hit && !entry.planned_hiso &&
                !entry.planned_siso && entry.assigned_entry_slot == -1 &&
                entry.assigned_core == -1,
            "suppressed K=0 early-stop row retained a planned route");
  }

  newcode::detail::TilePrepared<float> prep;
  prep.lin_matrix = matrix::Matrix<float>(32, newcode::Params::BCH_N);
  prep.lch_matrix = matrix::Matrix<float>(32, newcode::Params::BCH_N);
  observed_core_betas.clear();
  const auto decoded = newcode::detail::execute_level56_slice(
      prep, entries, 5, false, &controlled_soft_core);
  require(std::none_of(decoded.produced_rows.begin(),
                       decoded.produced_rows.end(),
                       [](bool produced) { return produced; }),
          "suppressed early-stop rows must not produce decoder output");
  require(observed_core_betas.empty(),
          "suppressed early-stop rows must not call the SISO core");

  require_grouped_plan_invariants(entries);
  require_state_conservation(entries);
}

void check_grouped_k_greater_than_8_selects_top_8() {
  auto entries = make_all_early_stop_grouped_entries();
  activate_group_prefix(&entries, 0, 4);
  for (std::size_t group = 1; group <= 8; ++group) {
    activate_group_prefix(&entries, group, 1);
  }
  const auto p = make_grouped_params();
  newcode::detail::schedule_level56_rows(&entries, p);
  newcode::detail::route_level56_g1(&entries, p);

  require(grouped_entry_count(entries) == 8,
          "K>8 must consume exactly eight entry slots");
  for (int slot = 0; slot < 8; ++slot) {
    require(group_for_slot(entries, slot) == slot,
            "K>8 load ties must prefer the smaller group index");
  }
  require(entries[32].final_action == Level56FinalAction::Unscheduled,
          "the ninth tied group must remain unscheduled");
  require(entries[32].source_level == 6,
          "the excluded tied group must be the later Level 6 group");
  require(entries[33].final_action == Level56FinalAction::EarlyStopAction,
          "default mode must preserve early-stop actions in an unentered group");
  require_grouped_plan_invariants(entries);
}

void check_grouped_k_greater_than_8_suppresses_unentered_group_early_stop() {
  auto entries = make_all_early_stop_grouped_entries();
  activate_group_prefix(&entries, 0, 4);
  for (std::size_t group = 1; group <= 8; ++group) {
    activate_group_prefix(&entries, group, 1);
  }
  auto p = make_grouped_params();
  p.LEVEL56_SCHEDULE_OBSERVABILITY_ENABLE = false;
  p.LEVEL56_EARLY_STOP_GROUP_UPDATE_MODE =
      newcode::Level56EarlyStopGroupUpdateMode::EnteredGroupsOnly;
  newcode::detail::schedule_level56_rows(&entries, p);
  newcode::detail::route_level56_g1(&entries, p);

  require(grouped_entry_count(entries) == 8,
          "suppression mode must not change the eight-entry schedule");
  for (int slot = 0; slot < 8; ++slot) {
    require(group_for_slot(entries, slot) == slot,
            "suppression mode changed grouped load ordering");
  }
  require(entries[5].early_stop_hit &&
              entries[5].final_action == Level56FinalAction::EarlyStopAction,
          "an entered group must preserve its early-stop actions");
  require(entries[32].final_action == Level56FinalAction::Unscheduled,
          "the unplanned non-early-stop row must remain unscheduled");
  for (std::size_t row = 33; row < 36; ++row) {
    require(entries[row].early_stop_hit &&
                entries[row].final_action == Level56FinalAction::Unscheduled &&
                !entries[row].planned_hiso && !entries[row].planned_siso &&
                entries[row].assigned_entry_slot == -1 &&
                entries[row].assigned_core == -1,
            "an unentered group's early-stop row was not suppressed");
  }
  require_grouped_plan_invariants(entries);
  require_state_conservation(entries);
}

void check_grouped_fill_idle_entries_k0_selects_first_8_groups() {
  auto entries = make_all_early_stop_grouped_entries();
  auto p = make_grouped_params();
  p.LEVEL56_EARLY_STOP_GROUP_UPDATE_MODE =
      newcode::Level56EarlyStopGroupUpdateMode::FillIdleEntries;
  newcode::Level56ScheduleSample sample;
  newcode::detail::schedule_level56_rows(&entries, p, 0, &sample);

  require(grouped_entry_count(entries) == 0,
          "K=0 fill mode must not create HISO/SISO plans");
  require(sample.total_group_entries == 8,
          "K=0 fill mode must consume eight supplemental group entries");
  require(sample.rounds.size() == 1 &&
              sample.rounds[0].selected_groups.size() == 8,
          "K=0 fill mode must record one eight-group supplemental round");
  for (std::size_t group = 0; group < 16; ++group) {
    const bool expected_entered = group < 8;
    require(sample.group_entry_counts[group] ==
                static_cast<std::size_t>(expected_entered),
            "K=0 fill mode selected the wrong supplemental group");
    for (std::size_t position = 0; position < 4; ++position) {
      const auto& entry = entries[group * 4 + position];
      require(entry.final_action ==
                  (expected_entered ? Level56FinalAction::EarlyStopAction
                                    : Level56FinalAction::Unscheduled),
              "K=0 fill mode applied the wrong early-stop action policy");
      require(!entry.planned_hiso && !entry.planned_siso &&
                  entry.assigned_entry_slot == -1 && entry.assigned_core == -1,
              "a supplemental early-stop entry must not reserve a core");
    }
  }

  auto entries_without_observability =
      make_all_early_stop_grouped_entries();
  p.LEVEL56_SCHEDULE_OBSERVABILITY_ENABLE = false;
  newcode::detail::schedule_level56_rows(&entries_without_observability, p);
  for (std::size_t index = 0; index < entries.size(); ++index) {
    require(entries_without_observability[index].final_action ==
                entries[index].final_action,
            "fill mode must not depend on schedule observability");
  }
  require_grouped_plan_invariants(entries);
  require_state_conservation(entries);
}

void check_grouped_fill_idle_entries_uses_ordinary_6_plus_2_supplemental() {
  auto entries = make_all_early_stop_grouped_entries();
  activate_group_prefix(&entries, 0, 3);
  activate_group_prefix(&entries, 1, 3);
  activate_group_prefix(&entries, 2, 2);
  activate_group_prefix(&entries, 3, 2);
  auto p = make_grouped_params();
  p.LEVEL56_EARLY_STOP_GROUP_UPDATE_MODE =
      newcode::Level56EarlyStopGroupUpdateMode::FillIdleEntries;
  newcode::Level56ScheduleSample sample;
  newcode::detail::schedule_level56_rows(&entries, p, 0, &sample);

  require(grouped_entry_count(entries) == 6,
          "fill mode changed the six-entry ordinary HISO/SISO schedule");
  require(sample.total_group_entries == 8,
          "fill mode must use the two idle entries after ordinary scheduling");
  require(sample.group_entry_counts[0] == 2 &&
              sample.group_entry_counts[1] == 2 &&
              sample.group_entry_counts[2] == 1 &&
              sample.group_entry_counts[3] == 1 &&
              sample.group_entry_counts[4] == 1 &&
              sample.group_entry_counts[5] == 1,
          "fill mode must select the lowest-index previously unentered groups");
  for (std::size_t row = 16; row < 24; ++row) {
    require(entries[row].early_stop_hit &&
                entries[row].final_action ==
                    Level56FinalAction::EarlyStopAction,
            "supplemental groups 4 and 5 must update early-stop rows");
  }
  for (std::size_t row = 24; row < 64; ++row) {
    require(entries[row].early_stop_hit &&
                entries[row].final_action == Level56FinalAction::Unscheduled,
            "groups after the two supplemental selections must not update");
  }
  require(sample.rounds.size() == 3 &&
              sample.rounds.back().selected_groups ==
                  std::vector<std::size_t>({4, 5}),
          "supplemental round must visit groups 4 and 5 in index order");
  require_grouped_plan_invariants(entries);
  require_state_conservation(entries);
}

void check_grouped_fill_idle_entries_does_not_change_full_schedule() {
  auto entered_only_entries = make_all_early_stop_grouped_entries();
  auto fill_entries = make_all_early_stop_grouped_entries();
  for (std::size_t group = 0; group < 8; ++group) {
    activate_group_prefix(&entered_only_entries, group, 1);
    activate_group_prefix(&fill_entries, group, 1);
  }

  auto entered_only_params = make_grouped_params();
  entered_only_params.LEVEL56_EARLY_STOP_GROUP_UPDATE_MODE =
      newcode::Level56EarlyStopGroupUpdateMode::EnteredGroupsOnly;
  auto fill_params = entered_only_params;
  fill_params.LEVEL56_EARLY_STOP_GROUP_UPDATE_MODE =
      newcode::Level56EarlyStopGroupUpdateMode::FillIdleEntries;
  newcode::Level56ScheduleSample fill_sample;
  newcode::detail::schedule_level56_rows(
      &entered_only_entries, entered_only_params);
  newcode::detail::schedule_level56_rows(
      &fill_entries, fill_params, 0, &fill_sample);

  require(fill_sample.total_group_entries == 8 &&
              fill_sample.rounds.size() == 1,
          "fill mode must not append supplemental entries after eight ordinary entries");
  for (std::size_t index = 0; index < fill_entries.size(); ++index) {
    require(fill_entries[index].final_action ==
                entered_only_entries[index].final_action &&
                fill_entries[index].assigned_entry_slot ==
                    entered_only_entries[index].assigned_entry_slot,
            "fill mode must match entered-only mode when all entries are used");
  }
  require_grouped_plan_invariants(fill_entries);
  require_state_conservation(fill_entries);
}

void check_grouped_k8_runs_once_per_group() {
  auto entries = make_all_early_stop_grouped_entries();
  for (std::size_t group = 0; group < 8; ++group) {
    activate_group_prefix(&entries, group, 4);
  }
  const auto p = make_grouped_params();
  newcode::detail::schedule_level56_rows(&entries, p);

  require(grouped_entry_count(entries) == 8,
          "K=8 must consume exactly eight entry slots");
  require(count_action(entries, Level56FinalAction::SisoDecode) == 8 &&
              count_action(entries, Level56FinalAction::HisoDecode) == 8,
          "K=8 must produce at most one SISO and HISO output per group");
  require(count_action(entries, Level56FinalAction::Unscheduled) == 16,
          "K=8 must not start a second scheduling round");
  require_grouped_plan_invariants(entries);
}

void check_grouped_k6_uses_6_plus_2_entries() {
  auto entries = make_all_early_stop_grouped_entries();
  for (std::size_t group = 0; group < 6; ++group) {
    activate_group_prefix(&entries, group, 4);
  }
  const auto p = make_grouped_params();
  newcode::detail::schedule_level56_rows(&entries, p);

  require(grouped_entry_count(entries) == 8,
          "K=6 must use six first-round and two second-round entries");
  require(group_for_slot(entries, 6) == 0 && group_for_slot(entries, 7) == 1,
          "second round must break equal remaining loads by group index");
  for (std::size_t row = 0; row < 8; ++row) {
    require(entries[row].final_action != Level56FinalAction::Unscheduled,
            "the two selected second-round groups must finish all four rows");
  }
  require_grouped_plan_invariants(entries);
}

void check_grouped_siso_only_multiround_reentry() {
  auto entries = make_all_early_stop_grouped_entries();
  activate_group_prefix(&entries, 0, 4,
                        newcode::detail::HybridRowClass::HardFail,
                        Level56Eligibility::SisoOnly);
  activate_group_prefix(&entries, 1, 4,
                        newcode::detail::HybridRowClass::HardFail,
                        Level56Eligibility::SisoOnly);
  const auto p = make_grouped_params();
  newcode::detail::schedule_level56_rows(&entries, p);

  require(grouped_entry_count(entries) == 8,
          "two all-SisoOnly groups must be allowed to re-enter four rounds");
  require(count_action(entries, Level56FinalAction::SisoDecode) == 8 &&
              count_action(entries, Level56FinalAction::HisoDecode) == 0,
          "SisoOnly rows must use only the SISO path");
  for (int slot = 0; slot < 8; ++slot) {
    require(group_for_slot(entries, slot) == slot % 2,
            "equal remaining SisoOnly loads must preserve group-index order");
  }
  require_grouped_plan_invariants(entries);
}

void check_grouped_intra_group_arbitration() {
  auto entries = make_all_early_stop_grouped_entries();
  activate_grouped_row(&entries, 0,
                       newcode::detail::HybridRowClass::ParityOnly);
  activate_grouped_row(&entries, 1,
                       newcode::detail::HybridRowClass::OneMain);
  activate_grouped_row(&entries, 2,
                       newcode::detail::HybridRowClass::HardFail,
                       Level56Eligibility::SisoOnly);
  activate_grouped_row(&entries, 3,
                       newcode::detail::HybridRowClass::TwoMain);
  const auto p = make_grouped_params();
  newcode::detail::schedule_level56_rows(&entries, p);

  require(entries[2].final_action == Level56FinalAction::SisoDecode,
          "SisoOnly HardFail must take SISO before flexible rows");
  require(entries[3].final_action == Level56FinalAction::HisoDecode,
          "TwoMain must be the highest-priority remaining flexible row");
  require(entries[1].final_action == Level56FinalAction::SisoDecode &&
              entries[0].final_action == Level56FinalAction::HisoDecode,
          "later round must schedule remaining flexible rows hardest first");
  require(entries[2].assigned_entry_slot == 0 &&
              entries[3].assigned_entry_slot == 0 &&
              entries[1].assigned_entry_slot == 1 &&
              entries[0].assigned_entry_slot == 1,
          "one group must re-enter with a new entry slot");
  require_grouped_plan_invariants(entries);
}

void check_grouped_single_flexible_row_uses_siso() {
  auto entries = make_all_early_stop_grouped_entries();
  activate_grouped_row(&entries, 36,
                       newcode::detail::HybridRowClass::TwoMain);
  const auto p = make_grouped_params();
  newcode::detail::schedule_level56_rows(&entries, p);

  require(entries[36].final_action == Level56FinalAction::SisoDecode &&
              entries[36].planned_siso && !entries[36].planned_hiso,
          "a single flexible row must select SISO before HISO");
  require(group_for_slot(entries, 0) == 9,
          "the selected row must remain inside its fixed Level 6 group");
  require_grouped_plan_invariants(entries);
}

void check_end_to_end_single_window() {
  newcode::Params p;
  p.CHASE_L = 2;
  p.CHASE_NTEST = 4;
  p.EARLY_STOP_ENABLE_LIST = {1, 1, 1, 1, 1, 1};
  p.EARLY_STOP_ACTION_SIGN_BETA_LIST = {1, 1, 1, 1, 1, 1};
  p.HARD_TILE_LIST = {0, 0, 0, 0, 0, 0};
  p.HYBRID_ENABLE_LIST = {0, 0, 0, 0, 1, 1};
  p.HYBRID_HARD_LLR_MAG_LIST = {4, 4, 4, 4, 4, 4};
  p.HYBRID_CLASSIFIER_MODE =
      newcode::HybridClassifierMode::RepoFastClassifier;
  p.HYBRID_USE_FAST_CLASSIFIER = true;
  p.ALPHA_LIST = {1, 1, 1, 1, 0.8f, 0.9f};
  p.beta_list = {1, 1, 1, 1, 1.2f, 1.4f};
  p.SISO_ACTIVE_LIST = {32, 32, 32, 32};
  p.HIHO_ACTIVE_LIST = {32, 32, 32, 32};
  p.LEVEL56_SHARED_ENABLE = true;
  p.LEVEL56_SHARED_HISO_ACTIVE = 64;
  p.LEVEL56_SHARED_SISO_ACTIVE = 64;

  matrix::Matrix<float> llr(p.win_height_rows(),
                            newcode::Params::NUM_SUBBLOCK_COLS *
                                newcode::Params::BITS_PER_SUBBLOCK_DIM);
  for (std::size_t row = 0; row < llr.rows(); ++row) {
    for (std::size_t col = 0; col < llr.cols(); ++col) {
      const int code = static_cast<int>((row * 17u + col * 13u) % 23u) - 11;
      llr[row][col] = static_cast<float>(code) * 0.35f +
                      ((row + col) & 1u ? 0.1f : -0.1f);
    }
  }

  std::vector<newcode::TileEarlyStopCounter> stats;
  const auto decoded = newcode::ofec_decode_llr_plain(
      llr, p, &stats, true, nullptr);
  require(decoded.rows() == llr.rows() && decoded.cols() == llr.cols(),
          "end-to-end shared decode changed the matrix shape");
  require(stats.size() == 6 && stats[4].total == 1 && stats[5].total == 1,
          "end-to-end shared decode did not process Level 5/6 once");
  for (std::size_t row = 0; row < decoded.rows(); ++row) {
    for (std::size_t col = 0; col < decoded.cols(); ++col) {
      require(std::isfinite(decoded[row][col]),
              "end-to-end shared decode produced a non-finite LLR");
    }
  }
}

void check_grouped_end_to_end_single_window() {
  newcode::Params p;
  p.CHASE_L = 2;
  p.CHASE_NTEST = 4;
  p.EARLY_STOP_ENABLE_LIST = {1, 1, 1, 1, 1, 1};
  p.EARLY_STOP_ACTION_SIGN_BETA_LIST = {1, 1, 1, 1, 1, 1};
  p.HARD_TILE_LIST = {0, 0, 0, 0, 0, 0};
  p.HYBRID_ENABLE_LIST = {0, 0, 0, 0, 1, 1};
  p.HYBRID_HARD_LLR_MAG_LIST = {4, 4, 4, 4, 4, 4};
  p.HYBRID_CLASSIFIER_MODE =
      newcode::HybridClassifierMode::FriendS1S3WithS0Classifier;
  p.HYBRID_USE_FAST_CLASSIFIER = true;
  p.ALPHA_LIST = {1, 1, 1, 1, 0.8f, 0.9f};
  p.beta_list = {1, 1, 1, 1, 1.2f, 1.4f};
  p.SISO_ACTIVE_LIST = {32, 32, 32, 32};
  p.HIHO_ACTIVE_LIST = {32, 32, 32, 32};
  p.LEVEL56_SHARED_ENABLE = true;
  p.LEVEL56_SHARED_HISO_ACTIVE = 8;
  p.LEVEL56_SHARED_SISO_ACTIVE = 8;
  p.LEVEL56_SCHEDULE_MODE =
      newcode::Level56ScheduleMode::Group4LoadSortedMultiround;
  p.LEVEL56_SCHEDULE_OBSERVABILITY_ENABLE = true;

  matrix::Matrix<float> llr(p.win_height_rows(),
                            newcode::Params::NUM_SUBBLOCK_COLS *
                                newcode::Params::BITS_PER_SUBBLOCK_DIM);
  for (std::size_t row = 0; row < llr.rows(); ++row) {
    for (std::size_t col = 0; col < llr.cols(); ++col) {
      const int code = static_cast<int>((row * 19u + col * 11u) % 29u) - 14;
      llr[row][col] = static_cast<float>(code) * 0.3f +
                      ((row + col) & 1u ? 0.05f : -0.05f);
    }
  }

  std::vector<newcode::TileEarlyStopCounter> stats;
  const auto decoded = newcode::ofec_decode_llr_plain(
      llr, p, &stats, true, nullptr);
  require(decoded.rows() == llr.rows() && decoded.cols() == llr.cols(),
          "grouped end-to-end decode changed the matrix shape");
  require(stats.size() == 6 && stats[4].total == 1 && stats[5].total == 1,
          "grouped end-to-end decode did not process Level 5/6 once");
  require(stats[4].level56_schedule_samples.size() == 1,
          "grouped end-to-end decode did not return a schedule sample");
  const auto& schedule = stats[4].level56_schedule_samples.front();
  require(schedule.codes.size() == 64,
          "grouped schedule sample must contain all 64 code rows");
  require(schedule.total_group_entries <= 8 &&
              schedule.planned_hiso_count <= 8 &&
              schedule.planned_siso_count <= 8,
          "grouped schedule sample exceeded the shared resource budget");
  if (!schedule.rounds.empty()) {
    require(schedule.rounds.back().used_entries_after ==
                schedule.total_group_entries,
            "grouped schedule rounds lost the final entry count");
  }
  for (std::size_t row = 0; row < decoded.rows(); ++row) {
    for (std::size_t col = 0; col < decoded.cols(); ++col) {
      require(std::isfinite(decoded[row][col]),
              "grouped end-to-end decode produced a non-finite LLR");
    }
  }
}

}  // namespace

int main() {
  try {
    check_full_siso_uses_no_hiso();
    check_hiso_only_handles_siso_overflow();
    check_siso_only_rows_never_use_hiso();
    check_early_stop_does_not_consume_shared_capacity();
    check_single_level_scheduler_bypasses_other_level();
    check_single_level_preserves_unselected_early_stop_actions();
    check_bypassed_slice_executes_only_preserved_early_stop_actions();
    check_single_level_selection_uses_fewer_early_stops();
    check_common_parameter_validation();
    check_per_level_siso_postprocessing();
    check_unscheduled_last_tile_history_passes_through_prior();
    check_level_priority_modes();
    check_grouped_k0_uses_no_entries();
    check_grouped_k0_suppresses_early_stop_without_entries();
    check_grouped_k_greater_than_8_selects_top_8();
    check_grouped_k_greater_than_8_suppresses_unentered_group_early_stop();
    check_grouped_fill_idle_entries_k0_selects_first_8_groups();
    check_grouped_fill_idle_entries_uses_ordinary_6_plus_2_supplemental();
    check_grouped_fill_idle_entries_does_not_change_full_schedule();
    check_grouped_k8_runs_once_per_group();
    check_grouped_k6_uses_6_plus_2_entries();
    check_grouped_siso_only_multiround_reentry();
    check_grouped_intra_group_arbitration();
    check_grouped_single_flexible_row_uses_siso();
    check_end_to_end_single_window();
    check_grouped_end_to_end_single_window();
    std::cout << "LEVEL56 shared scheduler regression checks passed\n";
    return 0;
  } catch (const std::exception& ex) {
    std::cerr << "LEVEL56 shared scheduler regression check failed: "
              << ex.what() << '\n';
    return 1;
  }
}
