#include "ofec/detail/ofec_tile_impl.ipp"
#include "ofec/detail/ofec_level56_shared.ipp"
#include "ofec/detail/ofec_level56_buffered_fifo.ipp"
#include "newcode/rx/ber/ber.hpp"

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

void check_twomain_parameterized_output() {
  std::vector<uint8_t> info(newcode::Params::BCH_K, 0u);
  for (std::size_t i = 0; i < info.size(); ++i) {
    info[i] = static_cast<uint8_t>(((i * 17u + 5u) % 11u) < 5u);
  }
  const auto transmitted = bch::bch_255_239_encode(info);
  auto received = transmitted;
  constexpr std::size_t kErrorA = 153;
  constexpr std::size_t kErrorB = 199;
  received[kErrorA] ^= 1u;
  received[kErrorB] ^= 1u;

  std::array<float, newcode::Params::BCH_N> lin{};
  for (std::size_t i = 0; i < lin.size(); ++i) {
    const float magnitude = 0.5f + static_cast<float>(i % 13u) * 0.125f;
    lin[i] = received[i] ? -magnitude : magnitude;
  }

  newcode::Params legacy;
  legacy.HYBRID_HARD_LLR_MAG = 99.0f;
  std::array<float, newcode::Params::BCH_N> legacy_lout{};
  std::array<uint8_t, newcode::Params::BCH_N> legacy_corrected{};
  newcode::detail::execute_hybrid_hard_class(
      newcode::detail::HybridRowClass::TwoMain, lin, legacy,
      &legacy_lout, &legacy_corrected);
  require(legacy_corrected == transmitted,
          "Legacy TwoMain HISO did not recover the transmitted codeword");
  for (std::size_t i = 0; i < legacy_lout.size(); ++i) {
    const float sign = transmitted[i] ? -1.0f : 1.0f;
    require(std::fabs(legacy_lout[i] - (sign * 99.0f - lin[i])) < 1e-6f,
            "Legacy TwoMain output changed after adding scheme 6");
  }

  auto uniform = legacy;
  uniform.TWOMAIN_HISO_OUTPUT_MODE =
      newcode::TwoMainHisoOutputMode::UnifiedParameterized;
  uniform.TWOMAIN_HISO_M2 = 32.0f;
  uniform.TWOMAIN_HISO_RHO_CORR = 1.0f;
  uniform.TWOMAIN_HISO_RHO_KEEP = 1.0f;
  std::array<float, newcode::Params::BCH_N> uniform_lout{};
  std::array<uint8_t, newcode::Params::BCH_N> uniform_corrected{};
  newcode::detail::execute_hybrid_hard_class(
      newcode::detail::HybridRowClass::TwoMain, lin, uniform,
      &uniform_lout, &uniform_corrected);
  require(uniform_corrected == legacy_corrected,
          "scheme 6-A changed BCH correction rather than only its output");
  for (std::size_t i = 0; i < uniform_lout.size(); ++i) {
    const float sign = transmitted[i] ? -1.0f : 1.0f;
    require(std::fabs(uniform_lout[i] - (sign * 32.0f - lin[i])) < 1e-6f,
            "scheme 6-A did not apply one uniform posterior magnitude");
  }

  auto differential = uniform;
  differential.TWOMAIN_HISO_RHO_CORR = 0.75f;
  differential.TWOMAIN_HISO_RHO_KEEP = 0.25f;
  std::array<float, newcode::Params::BCH_N> differential_lout{};
  newcode::detail::execute_hybrid_hard_class(
      newcode::detail::HybridRowClass::TwoMain, lin, differential,
      &differential_lout);
  for (std::size_t i = 0; i < differential_lout.size(); ++i) {
    const bool corrected = i == kErrorA || i == kErrorB;
    const float expected_mag = corrected ? 24.0f : 8.0f;
    const float sign = transmitted[i] ? -1.0f : 1.0f;
    require(std::fabs(differential_lout[i] -
                      (sign * expected_mag - lin[i])) < 1e-6f,
            "scheme 6-B/6-C did not distinguish corrected and kept positions");
  }

  // 方案六参数不得影响任何非 TwoMain 类别。
  auto parity_received = transmitted;
  parity_received[newcode::Params::BCH_OVERALL_IDX] ^= 1u;
  for (std::size_t i = 0; i < lin.size(); ++i) {
    lin[i] = parity_received[i] ? -1.25f : 1.25f;
  }
  std::array<float, newcode::Params::BCH_N> parity_legacy{};
  std::array<float, newcode::Params::BCH_N> parity_scheme6{};
  newcode::detail::execute_hybrid_hard_class(
      newcode::detail::HybridRowClass::ParityOnly, lin, legacy,
      &parity_legacy);
  newcode::detail::execute_hybrid_hard_class(
      newcode::detail::HybridRowClass::ParityOnly, lin, differential,
      &parity_scheme6);
  require(parity_legacy == parity_scheme6,
          "scheme 6 parameters affected a non-TwoMain HISO class");
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

void check_temporal_lookahead_three_windows() {
  newcode::Params p;
  p.CHASE_L = 2;
  p.CHASE_NTEST = 4;
  p.EARLY_STOP_ENABLE_LIST = {1, 1, 1, 1, 1, 1};
  p.EARLY_STOP_ACTION_SIGN_BETA_LIST = {1, 1, 1, 1, 1, 1};
  p.EARLY_STOP_BIND_GROUP_SIZE = 1;
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
  p.LEVEL56_TEMPORAL_LOOKAHEAD_ENABLE = true;
  p.LEVEL56_SHARED_HISO_ACTIVE = 8;
  p.LEVEL56_SHARED_SISO_ACTIVE = 8;
  p.LEVEL56_PRIORITY_MODE = newcode::Level56PriorityMode::Level5First;
  p.LEVEL56_SCHEDULE_MODE =
      newcode::Level56ScheduleMode::Group4LoadSortedMultiround;
  p.LEVEL56_EARLY_STOP_GROUP_UPDATE_MODE =
      newcode::Level56EarlyStopGroupUpdateMode::AllGroups;
  p.LEVEL56_SCHEDULE_OBSERVABILITY_ENABLE = true;

  const std::size_t rows = p.win_height_rows() + 2 * p.pop_push_rows();
  matrix::Matrix<float> llr(
      rows, newcode::Params::NUM_SUBBLOCK_COLS *
                newcode::Params::BITS_PER_SUBBLOCK_DIM);
  for (std::size_t row = 0; row < llr.rows(); ++row) {
    for (std::size_t col = 0; col < llr.cols(); ++col) {
      const int code = static_cast<int>((row * 31u + col * 17u) % 37u) - 18;
      llr[row][col] = static_cast<float>(code) * 0.22f +
                      ((row + col) & 1u ? 0.03f : -0.03f);
    }
  }

  std::vector<newcode::TileEarlyStopCounter> stats;
  const auto decoded = newcode::ofec_decode_llr_plain(
      llr, p, &stats, true, nullptr);
  require(decoded.rows() == llr.rows() && decoded.cols() == llr.cols(),
          "temporal decode changed the matrix shape");
  require(stats.size() == 6 && stats[4].total == 3 && stats[5].total == 3,
          "temporal decode did not process exactly three Level 5/6 batches");
  require(stats[4].level56_schedule_samples.size() == 3,
          "temporal decode did not export three schedule samples");

  const auto& samples = stats[4].level56_schedule_samples;
  require(samples.front().temporal_branch ==
              newcode::Level56TemporalBranch::NoHistory,
          "first temporal batch must use the no-history boundary branch");
  require(samples.back().temporal_branch ==
              newcode::Level56TemporalBranch::NoFuture,
          "last temporal batch must use the no-future boundary branch");
  require(samples[0].codes.size() == 128 &&
              samples[1].codes.size() == 192 &&
              samples[2].codes.size() == 128,
          "temporal boundaries must export 128/192/128 real code entries");
  for (std::size_t sample_index = 0; sample_index < samples.size();
       ++sample_index) {
    const auto& sample = samples[sample_index];
    require(sample.temporal_lookahead_enabled,
            "temporal sample did not record lookahead enablement");
    require(sample.total_group_entries <= 8 &&
                sample.temporal_t0_group_entries <= 8 &&
                sample.temporal_t1_group_entries <= 8,
            "temporal schedule exceeded the eight-entry budget");
    require(sample.temporal_t0_group_entries +
                sample.temporal_t1_group_entries <= 8,
            "temporal t=0/t=1 entries exceeded the shared budget");
    for (const auto& code : sample.codes) {
      require(code.code_index == code.time_index * 64 +
                                     (code.source_level == 5 ? 0 : 32) +
                                     code.source_local_row,
              "temporal code violated the fixed 192-code linear layout");
      require(code.time_offset == static_cast<int>(code.time_index) - 1,
              "temporal code has an invalid time offset");
      if (code.time_index == 0) {
        require(code.info_type ==
                    newcode::Level56TemporalInfoType::DecodeInfo,
                "t=0 must expose complete decode information");
      } else {
        require(code.info_type ==
                    newcode::Level56TemporalInfoType::EarlyStopInfo,
                "t=1/t=2 must expose EarlyStop information");
      }
      if (code.time_index == 1 && code.early_stop_hit) {
        require(code.final_action == 0,
                "AllGroups temporal mode suppressed an EarlyStopAction");
      }
      if (code.time_index == 2) {
        require(code.hybrid_class == 0 && code.resource_eligibility == 0 &&
                    code.final_action ==
                        static_cast<uint8_t>(Level56FinalAction::Unscheduled) &&
                    !code.produced && code.assigned_entry_slot == -1,
                "t=2 lookahead performed classification or decode work");
      }
    }
    const std::size_t first_time_index = sample_index == 0 ? 1 : 0;
    const std::size_t last_time_index = sample_index == 2 ? 1 : 2;
    require(sample.codes.front().time_index == first_time_index &&
                sample.codes.back().time_index == last_time_index,
            "temporal boundary exported a synthetic missing batch");
  }
  for (std::size_t row = 0; row < decoded.rows(); ++row) {
    for (std::size_t col = 0; col < decoded.cols(); ++col) {
      require(std::isfinite(decoded[row][col]),
              "temporal decode produced a non-finite LLR");
    }
  }
}

void check_temporal_rule_helpers() {
  require(newcode::detail::select_level56_temporal_branch(
              false, true, 0, 4, 4) ==
              newcode::Level56TemporalBranch::NoHistory,
          "missing history must select NoHistory");
  require(newcode::detail::select_level56_temporal_branch(
              true, false, 8, 0, 0) ==
              newcode::Level56TemporalBranch::NoFuture,
          "missing future must select NoFuture");
  require(newcode::detail::select_level56_temporal_branch(
              true, true, 3, 4, 4) ==
              newcode::Level56TemporalBranch::SupplementHistory,
          "strictly light temporal load must supplement history");
  require(newcode::detail::select_level56_temporal_branch(
              true, true, 3, 6, 7) ==
              newcode::Level56TemporalBranch::CurrentFirst,
          "equality at K1+K2=16-X must select current first");
  require(newcode::detail::select_level56_temporal_branch(
              true, true, 0, 0, 0) ==
              newcode::Level56TemporalBranch::CurrentFirst,
          "X=0 must process the current batch");
  require(newcode::detail::select_level56_temporal_branch(
              true, true, 3, 4, 4, 11) ==
              newcode::Level56TemporalBranch::CurrentFirst,
          "custom temporal threshold must switch to current first");
  require(newcode::detail::select_level56_temporal_branch(
              true, true, 3, 4, 4, 12) ==
              newcode::Level56TemporalBranch::SupplementHistory,
          "custom temporal threshold must allow history supplement");

  auto entries = make_entries(32);
  for (auto& entry : entries) {
    entry.early_stop_hit = true;
    entry.produced = false;
    entry.final_action = Level56FinalAction::EarlyStopAction;
  }
  entries[0].early_stop_hit = false;
  entries[0].final_action = Level56FinalAction::Unscheduled;
  entries[4].early_stop_hit = false;
  entries[4].final_action = Level56FinalAction::HisoDecode;
  entries[4].produced = true;
  require(newcode::detail::is_level56_pending_decode(entries[0]) &&
              !newcode::detail::is_level56_pending_decode(entries[1]) &&
              !newcode::detail::is_level56_pending_decode(entries[4]),
          "temporal pending predicate does not match the three-field rule");
  require(newcode::detail::level56_temporal_pending_group_count(entries) == 1,
          "X must count groups containing at least one pending code");
}

void check_temporal_two_stage_budget() {
  auto history = make_all_early_stop_grouped_entries();
  for (std::size_t group = 0; group < 3; ++group) {
    activate_grouped_row(&history, group * 4u,
                         newcode::detail::HybridRowClass::ParityOnly);
  }
  auto current = make_all_early_stop_grouped_entries();
  for (std::size_t group = 0; group < 16; ++group) {
    activate_grouped_row(&current, group * 4u,
                         newcode::detail::HybridRowClass::ParityOnly);
  }
  auto p = make_grouped_params();
  p.LEVEL56_EARLY_STOP_GROUP_UPDATE_MODE =
      newcode::Level56EarlyStopGroupUpdateMode::AllGroups;

  newcode::Level56ScheduleSample history_sample;
  newcode::detail::schedule_level56_rows(
      &history, p, 0, &history_sample, 3, 0, false);
  require(grouped_entry_count(history) == 3,
          "history stage must consume exactly used_t0 slots");
  require(history_sample.total_group_entries == 3,
          "history stage reported the wrong used_t0 budget");

  newcode::Level56ScheduleSample current_sample;
  newcode::detail::schedule_level56_rows(
      &current, p, 0, &current_sample, 5, 3, false);
  require(current_sample.total_group_entries == 8,
          "current stage must report the combined slot endpoint");
  for (const auto& entry : current) {
    if (entry.planned_hiso || entry.planned_siso) {
      require(entry.assigned_entry_slot >= 3 &&
                  entry.assigned_entry_slot < 8,
              "current stage used a history entry slot");
    }
  }
  require(grouped_entry_count(current) == 8,
          "two-stage schedule must use at most eight total slots");

  auto no_budget = make_all_early_stop_grouped_entries();
  activate_grouped_row(&no_budget, 0,
                       newcode::detail::HybridRowClass::ParityOnly);
  newcode::Level56ScheduleSample no_budget_sample;
  newcode::detail::schedule_level56_rows(
      &no_budget, p, 0, &no_budget_sample, 0, 8, false);
  require(no_budget_sample.total_group_entries == 8,
          "B=0 must preserve the consumed history slot endpoint");
  require(grouped_entry_count(no_budget) == 0,
          "B=0 must not plan ordinary HISO/SISO entries");
  require(no_budget[0].final_action == Level56FinalAction::Unscheduled,
          "B=0 must leave ordinary candidates unscheduled");
  for (std::size_t index = 1; index < no_budget.size(); ++index) {
    require(no_budget[index].final_action ==
                Level56FinalAction::EarlyStopAction,
            "AllGroups EarlyStopAction must not depend on B");
  }
  for (int slot = 3; slot < 8; ++slot) {
    require(group_for_slot(current, slot) >= 0,
            "every current-stage slot must select a group");
  }
}

void check_buffered_window_transition_formula() {
  using newcode::detail::level56_buffered_next_window_start;
  require(level56_buffered_next_window_start(32, 0, 32) == 30,
          "C_t=0 must move S_t back by two rows");
  require(level56_buffered_next_window_start(28, 1, 32) == 28,
          "C_t=1 must keep S_t unchanged");
  require(level56_buffered_next_window_start(28, 2, 32) == 30,
          "C_t=2 must recover S_t by two rows");
  require(level56_buffered_next_window_start(30, 2, 32) == 32,
          "S_t recovery must clip at R_buf");
  require(level56_buffered_next_window_start(0, 0, 32) == 0,
          "S_t must not move below zero");
  require(level56_buffered_next_window_start(0, 1, 32) == 0,
          "one completion at S_t=0 must only offset the fixed push");
  require(level56_buffered_next_window_start(0, 2, 32) == 2,
          "two completions at S_t=0 must begin recovery");
  require(level56_buffered_next_window_start(5, 0, 5) == 3 &&
              level56_buffered_next_window_start(3, 0, 5) == 1 &&
              level56_buffered_next_window_start(1, 0, 5) == 0,
          "odd R_buf must expose only complete two-row fallback steps");
}

void check_buffered_already_decoded_three_round_example() {
  auto entries = make_all_early_stop_grouped_entries();
  for (std::size_t group = 0; group < 12; ++group) {
    activate_group_prefix(&entries, group, 4);
  }
  auto p = make_grouped_params();
  p.LEVEL56_BUFFERED_FIFO_ENABLE = true;
  p.LEVEL56_EARLY_STOP_GROUP_UPDATE_MODE =
      newcode::Level56EarlyStopGroupUpdateMode::AllGroups;

  const auto complete_scheduled_round = [](auto* values) {
    for (auto& entry : *values) {
      if (entry.already_decoded) {
        continue;
      }
      if (entry.early_stop_hit) {
        entry.already_decoded = true;
        entry.pending = false;
        entry.produced = true;
        entry.writeback_complete = true;
        entry.decode_status = newcode::Level56DecodeStatus::Produced;
      } else if (entry.planned_hiso || entry.planned_siso) {
        entry.already_decoded = true;
        entry.pending = false;
        entry.produced = true;
        entry.writeback_complete = true;
        entry.decode_status = newcode::Level56DecodeStatus::Produced;
      } else {
        entry.pending = true;
        entry.produced = false;
      }
    }
  };
  const auto prepare_retry = [](auto* values) {
    for (auto& entry : *values) {
      entry.planned_hiso = false;
      entry.planned_siso = false;
      entry.assigned_entry_slot = -1;
      entry.assigned_core = -1;
      entry.needs_execution = false;
      entry.final_action = Level56FinalAction::Unscheduled;
      entry.remaining_for_schedule =
          !entry.already_decoded && entry.pending;
    }
  };
  const auto count_completed = [](const auto& values) {
    return static_cast<std::size_t>(std::count_if(
        values.begin(), values.end(),
        [](const auto& entry) { return entry.already_decoded; }));
  };

  newcode::detail::schedule_level56_rows(&entries, p);
  complete_scheduled_round(&entries);
  require(count_completed(entries) == 32 &&
              newcode::detail::level56_buffered_pending_count(entries) == 32,
          "t=0 must leave 32 completed and 32 pending codes");
  for (const auto& entry : entries) {
    if (entry.early_stop_hit) {
      require(entry.produced && entry.writeback_complete,
              "EarlyStop completion must produce its original action result");
    }
  }

  prepare_retry(&entries);
  newcode::detail::schedule_level56_rows(
      &entries, p, 0, nullptr, 8, 0, true);
  for (const auto& entry : entries) {
    if (entry.already_decoded) {
      require(!entry.planned_hiso && !entry.planned_siso &&
                  entry.assigned_entry_slot == -1,
              "AlreadyDecoded code re-entered the scheduler");
    }
  }
  complete_scheduled_round(&entries);
  require(count_completed(entries) == 48 &&
              newcode::detail::level56_buffered_pending_count(entries) == 16,
          "t=1 must schedule only pending codes and leave 16 pending");

  prepare_retry(&entries);
  newcode::detail::schedule_level56_rows(
      &entries, p, 0, nullptr, 8, 0, true);
  complete_scheduled_round(&entries);
  require(newcode::detail::level56_buffered_batch_complete(entries) &&
              newcode::detail::level56_buffered_pending_count(entries) == 0,
          "t=2 must complete the final 16 pending codes");
}

void check_buffered_retirement_markers() {
  auto full_early_stop = make_all_early_stop_grouped_entries();
  newcode::detail::mark_level56_full_early_stop_complete(
      &full_early_stop);
  require(newcode::detail::level56_buffered_batch_complete(full_early_stop),
          "full EarlyStop batch must retire as a completed batch");
  for (const auto& entry : full_early_stop) {
    require(entry.already_decoded && entry.produced && !entry.pending &&
                entry.writeback_complete &&
                entry.decode_status ==
                    newcode::Level56DecodeStatus::Produced,
            "full EarlyStop code has an invalid action/writeback completion state");
  }

  auto incomplete_writeback = full_early_stop;
  incomplete_writeback[7].writeback_complete = false;
  require(!newcode::detail::level56_buffered_batch_complete(
              incomplete_writeback),
          "batch must not retire before every required writeback completes");

  auto forced = make_entries(32);
  forced[0].already_decoded = true;
  forced[0].produced = true;
  forced[0].writeback_complete = true;
  for (std::size_t index = 1; index < forced.size(); ++index) {
    forced[index].pending = true;
  }
  newcode::detail::mark_level56_batch_forced_evicted(&forced);
  require(forced[0].already_decoded && !forced[0].forced_evicted,
          "forced eviction must preserve code that already completed");
  for (std::size_t index = 1; index < forced.size(); ++index) {
    require(!forced[index].already_decoded && forced[index].forced_evicted &&
                !forced[index].pending && !forced[index].produced &&
                forced[index].decode_status ==
                    newcode::Level56DecodeStatus::ForcedEvicted,
            "unfinished code did not enter the forced-evicted state");
  }
}

newcode::detail::Level56BufferedBatch<float> make_buffered_test_batch(
    std::size_t batch_id,
    bool full_early_stop) {
  newcode::detail::Level56BufferedBatch<float> batch;
  batch.batch_id = batch_id;
  batch.arrival_time = batch_id;
  batch.classified = true;
  batch.entries = full_early_stop
                      ? make_all_early_stop_grouped_entries()
                      : make_entries(32);
  return batch;
}

void mark_buffered_test_batch_complete(
    newcode::detail::Level56BufferedBatch<float>* batch) {
  for (auto& entry : batch->entries) {
    entry.already_decoded = true;
    entry.pending = false;
    entry.produced = true;
    entry.writeback_complete = true;
    entry.decode_status = newcode::Level56DecodeStatus::Produced;
  }
}

void check_buffered_service_fast_path_and_fifo_order() {
  newcode::detail::Level56BufferedFifoState<float> state;
  state.fifo.push_back(make_buffered_test_batch(0, true));
  state.fifo.push_back(make_buffered_test_batch(1, true));
  state.fifo.push_back(make_buffered_test_batch(2, false));
  state.fifo.push_back(make_buffered_test_batch(3, false));

  std::vector<std::size_t> classified;
  std::vector<std::size_t> ordinary;
  const auto outcome = newcode::detail::service_level56_buffered_fifo_time(
      &state, false,
      [&](auto* batch) { classified.push_back(batch->batch_id); },
      [&](auto* batch) { ordinary.push_back(batch->batch_id); },
      [](auto* batch) { mark_buffered_test_batch_complete(batch); });

  require(outcome.completed_batches == 2 &&
              outcome.full_early_stop_batches == 2 &&
              outcome.ordinary_service_used,
          "leading full-EarlyStop batches must retire before ordinary service");
  require(outcome.retirements.size() == 2 &&
              outcome.retirements[0].batch_id == 0 &&
              outcome.retirements[1].batch_id == 1,
          "full-EarlyStop retirement must preserve FIFO order");
  require(ordinary.size() == 1 && ordinary[0] == 2 &&
              state.fifo.size() == 2 && state.fifo.front().batch_id == 2,
          "one t must serve only the first ordinary FIFO head");
  require(classified == std::vector<std::size_t>({0, 1, 2}),
          "scheduler classified a future batch beyond the blocked head");

  newcode::detail::Level56BufferedFifoState<float> all_fast;
  for (std::size_t id = 0; id < 5; ++id) {
    all_fast.fifo.push_back(make_buffered_test_batch(id, true));
  }
  std::size_t fast_callbacks = 0;
  const auto all_fast_outcome =
      newcode::detail::service_level56_buffered_fifo_time(
          &all_fast, false, [](auto*) {}, [](auto*) {},
          [&](auto* batch) {
            mark_buffered_test_batch_complete(batch);
            ++fast_callbacks;
          });
  require(all_fast_outcome.completed_batches == 5 &&
              all_fast_outcome.full_early_stop_batches == 5 &&
              !all_fast_outcome.ordinary_service_used && all_fast.fifo.empty() &&
              fast_callbacks == 5,
          "one t must allow an unbounded run of full-EarlyStop heads");
}

void check_buffered_service_ordinary_budget_boundary() {
  newcode::detail::Level56BufferedFifoState<float> ordinary_then_ordinary;
  ordinary_then_ordinary.fifo.push_back(make_buffered_test_batch(0, false));
  ordinary_then_ordinary.fifo.push_back(make_buffered_test_batch(1, false));
  std::vector<std::size_t> ordinary_calls;
  std::vector<std::pair<std::size_t, int>> observed_shared_versions;
  int shared_version = 0;
  const auto ordinary_outcome =
      newcode::detail::service_level56_buffered_fifo_time(
          &ordinary_then_ordinary, false,
          [&](auto* batch) {
            observed_shared_versions.emplace_back(
                batch->batch_id, shared_version);
          },
          [&](auto* batch) {
            ordinary_calls.push_back(batch->batch_id);
            mark_buffered_test_batch_complete(batch);
            shared_version = 1;
          },
          [](auto* batch) { mark_buffered_test_batch_complete(batch); });
  require(ordinary_outcome.completed_batches == 1 &&
              ordinary_calls == std::vector<std::size_t>({0}) &&
              ordinary_then_ordinary.fifo.size() == 1 &&
              ordinary_then_ordinary.fifo.front().batch_id == 1,
          "ordinary completion must not transfer unused entries to B1");
  require(observed_shared_versions ==
              std::vector<std::pair<std::size_t, int>>({{0, 0}, {1, 1}}),
          "B1 must inspect shared SRAM only after B0 writeback completes");

  newcode::detail::Level56BufferedFifoState<float> ordinary_then_fast;
  ordinary_then_fast.fifo.push_back(make_buffered_test_batch(0, false));
  ordinary_then_fast.fifo.push_back(make_buffered_test_batch(1, true));
  ordinary_then_fast.fifo.push_back(make_buffered_test_batch(2, false));
  const auto mixed_outcome =
      newcode::detail::service_level56_buffered_fifo_time(
          &ordinary_then_fast, false, [](auto*) {},
          [&](auto* batch) { mark_buffered_test_batch_complete(batch); },
          [](auto* batch) { mark_buffered_test_batch_complete(batch); });
  require(mixed_outcome.completed_batches == 2 &&
              mixed_outcome.full_early_stop_batches == 1 &&
              mixed_outcome.retirements.size() == 2 &&
              mixed_outcome.retirements[0].reason ==
                  newcode::Level56BufferedBatchRetirement::Normal &&
              mixed_outcome.retirements[1].reason ==
                  newcode::Level56BufferedBatchRetirement::FullEarlyStop &&
              ordinary_then_fast.fifo.size() == 1 &&
              ordinary_then_fast.fifo.front().batch_id == 2,
          "ordinary completion must still allow the following full-EarlyStop path");
}

void check_buffered_service_forced_eviction_boundary() {
  newcode::detail::Level56BufferedFifoState<float> state;
  state.fifo.push_back(make_buffered_test_batch(0, false));
  state.fifo.push_back(make_buffered_test_batch(1, true));
  const auto outcome = newcode::detail::service_level56_buffered_fifo_time(
      &state, true, [](auto*) {}, [](auto*) {},
      [](auto* batch) { mark_buffered_test_batch_complete(batch); });
  require(outcome.completed_batches == 0 &&
              outcome.forced_evicted_batches == 1 &&
              outcome.forced_evicted_global_rows.size() == 64 &&
              outcome.retirements.size() == 1 &&
              outcome.retirements[0].batch_id == 0 &&
              outcome.retirements[0].reason ==
                  newcode::Level56BufferedBatchRetirement::ForcedEvicted &&
              state.fifo.size() == 1 && state.fifo.front().batch_id == 1,
          "S_t=0 must evict only the original unfinished head without C_t credit");

  newcode::detail::Level56BufferedFifoState<float> exposed_head;
  exposed_head.fifo.push_back(make_buffered_test_batch(0, true));
  exposed_head.fifo.push_back(make_buffered_test_batch(1, false));
  const auto exposed_outcome =
      newcode::detail::service_level56_buffered_fifo_time(
          &exposed_head, true, [](auto*) {}, [](auto*) {},
          [](auto* batch) { mark_buffered_test_batch_complete(batch); });
  require(exposed_outcome.completed_batches == 1 &&
              exposed_outcome.forced_evicted_batches == 0 &&
              exposed_head.fifo.size() == 1 &&
              exposed_head.fifo.front().batch_id == 1,
          "a newly exposed head must not be forced out in the same boundary t");
}

void check_buffered_configuration_validation() {
  auto p = make_grouped_params();
  p.HYBRID_ENABLE_LIST = {0, 0, 0, 0, 1, 1};
  p.LEVEL56_BUFFERED_FIFO_ENABLE = true;
  p.LEVEL56_BUFFER_ROWS = 4;
  p.LEVEL56_EARLY_STOP_GROUP_UPDATE_MODE =
      newcode::Level56EarlyStopGroupUpdateMode::AllGroups;
  newcode::detail::validate_level56_shared_config(p);

  newcode::detail::Level56BufferedFifoState<float> state;
  newcode::detail::initialize_level56_buffered_fifo_state(&state, p);
  require(state.initialized && state.window_start == 4,
          "buffered FIFO must initialize S_0 from the configured R_buf");

  auto temporal_conflict = p;
  temporal_conflict.LEVEL56_TEMPORAL_LOOKAHEAD_ENABLE = true;
  bool rejected = false;
  try {
    newcode::detail::validate_level56_shared_config(temporal_conflict);
  } catch (const std::invalid_argument&) {
    rejected = true;
  }
  require(rejected, "buffered FIFO must reject temporal lookahead");

  auto wrong_push = p;
  wrong_push.WINDOW_POP_PUSH = 1;
  rejected = false;
  try {
    newcode::detail::validate_level56_shared_config(wrong_push);
  } catch (const std::invalid_argument&) {
    rejected = true;
  }
  require(rejected, "buffered FIFO must require a two-row physical push");
}

newcode::Params make_buffered_pipeline_params() {
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
  p.LEVEL56_BUFFERED_FIFO_ENABLE = true;
  p.LEVEL56_BUFFER_ROWS = 32;
  p.LEVEL56_SHARED_HISO_ACTIVE = 8;
  p.LEVEL56_SHARED_SISO_ACTIVE = 8;
  p.LEVEL56_SCHEDULE_MODE =
      newcode::Level56ScheduleMode::Group4LoadSortedMultiround;
  p.LEVEL56_EARLY_STOP_GROUP_UPDATE_MODE =
      newcode::Level56EarlyStopGroupUpdateMode::AllGroups;
  p.LEVEL56_SCHEDULE_OBSERVABILITY_ENABLE = true;
  return p;
}

void check_buffered_end_to_end_full_early_stop() {
  auto p = make_buffered_pipeline_params();
  constexpr std::size_t kServiceTimes = 3;
  const std::size_t rows = p.win_height_rows() +
                           (kServiceTimes - 1u) * p.pop_push_rows();
  matrix::Matrix<float> llr(
      rows, newcode::Params::NUM_SUBBLOCK_COLS *
                newcode::Params::BITS_PER_SUBBLOCK_DIM);
  for (std::size_t row = 0; row < llr.rows(); ++row) {
    for (std::size_t col = 0; col < llr.cols(); ++col) {
      llr[row][col] = 8.0f;
    }
  }

  std::vector<newcode::TileEarlyStopCounter> stats;
  (void)newcode::ofec_decode_llr_plain(llr, p, &stats, true, nullptr);
  require(stats.size() == p.TILES_PER_WIN &&
              stats[4].level56_buffered_time_samples.size() == kServiceTimes,
          "full-EarlyStop integration did not emit one sample per t");
  require(stats[4].level56_schedule_samples.size() == kServiceTimes,
          "full-EarlyStop integration did not record its K=0 action pass");
  for (const auto& schedule : stats[4].level56_schedule_samples) {
    require(schedule.branch == newcode::Level56ScheduleBranch::K0 &&
                schedule.total_group_entries == 0 &&
                schedule.planned_hiso_count == 0 &&
                schedule.planned_siso_count == 0,
            "full-EarlyStop fast path consumed ordinary decode resources");
    require(std::all_of(
                schedule.codes.begin(), schedule.codes.end(),
                [](const auto& code) {
                  return code.early_stop_hit && code.produced &&
                         code.already_decoded && code.writeback_complete;
                }),
            "full-EarlyStop fast path skipped EarlyStopAction writeback");
  }
  const auto& samples = stats[4].level56_buffered_time_samples;
  for (std::size_t t = 0; t < samples.size(); ++t) {
    const auto& sample = samples[t];
    require(sample.service_time == t && sample.arrived_batch_id == t &&
                sample.fifo_depth_before == 1 &&
                sample.fifo_depth_after == 0 &&
                !sample.ordinary_service_used &&
                sample.completed_batches == 1 &&
                sample.full_early_stop_batches == 1 &&
                sample.forced_evicted_batches == 0 &&
                sample.window_start_before == p.LEVEL56_BUFFER_ROWS &&
                sample.window_start_after == p.LEVEL56_BUFFER_ROWS &&
                sample.retirements.size() == 1 &&
                sample.retirements.front().batch_id == t &&
                sample.retirements.front().reason ==
                    newcode::Level56BufferedBatchRetirement::FullEarlyStop,
            "real full-EarlyStop batch did not take the documented fast path");
  }
}

void check_buffered_end_to_end_multiple_windows() {
  auto p = make_buffered_pipeline_params();

  constexpr std::size_t kServiceTimes = 5;
  const std::size_t rows = p.win_height_rows() +
                           (kServiceTimes - 1u) * p.pop_push_rows();
  matrix::Matrix<float> llr(
      rows, newcode::Params::NUM_SUBBLOCK_COLS *
                newcode::Params::BITS_PER_SUBBLOCK_DIM);
  for (std::size_t row = 0; row < llr.rows(); ++row) {
    for (std::size_t col = 0; col < llr.cols(); ++col) {
      const int code =
          static_cast<int>((row * 19u + col * 11u) % 29u) - 14;
      llr[row][col] = static_cast<float>(code) * 0.3f +
                      ((row + col) & 1u ? 0.05f : -0.05f);
    }
  }

  std::vector<newcode::TileEarlyStopCounter> stats;
  const auto decoded = newcode::ofec_decode_llr_plain(
      llr, p, &stats, true, nullptr);
  require(decoded.rows() == llr.rows() && decoded.cols() == llr.cols(),
          "buffered end-to-end decode changed the matrix shape");
  require(stats.size() == p.TILES_PER_WIN &&
              stats[4].level56_buffered_time_samples.size() == kServiceTimes,
          "buffered end-to-end decode did not emit one sample per t");

  const auto& samples = stats[4].level56_buffered_time_samples;
  std::size_t previous_fifo_depth = 0;
  std::size_t next_retired_batch = 0;
  std::size_t ordinary_service_times = 0;
  std::size_t schedule_index = 0;
  bool observed_pending = false;
  bool observed_window_fallback = false;
  std::vector<std::vector<newcode::Level56ScheduleCodeSample>>
      prior_codes_by_batch;
  for (std::size_t t = 0; t < samples.size(); ++t) {
    const auto& sample = samples[t];
    require(sample.service_time == t && sample.arrived_batch_id == t,
            "buffered end-to-end arrival/service time sequence is invalid");
    require(sample.fifo_depth_before == previous_fifo_depth + 1u,
            "buffered end-to-end FIFO did not append exactly one batch");
    require(sample.had_head_before &&
                sample.head_batch_id_before == next_retired_batch,
            "buffered end-to-end service did not select the FIFO head");

    std::size_t normal_or_fast = 0;
    std::size_t forced = 0;
    for (const auto& retirement : sample.retirements) {
      require(retirement.batch_id == next_retired_batch,
              "buffered end-to-end retirement order is not FIFO");
      ++next_retired_batch;
      if (retirement.reason ==
          newcode::Level56BufferedBatchRetirement::ForcedEvicted) {
        ++forced;
      } else {
        ++normal_or_fast;
      }
    }
    require(normal_or_fast == sample.completed_batches &&
                forced == sample.forced_evicted_batches,
            "buffered end-to-end C_t included an invalid retirement type");
    require(sample.fifo_depth_after + sample.retirements.size() ==
                sample.fifo_depth_before,
            "buffered end-to-end FIFO depth violates arrival/retirement conservation");
    require(sample.window_start_after ==
                newcode::detail::level56_buffered_next_window_start(
                    sample.window_start_before, sample.completed_batches,
                    p.LEVEL56_BUFFER_ROWS),
            "buffered end-to-end S_t transition violates the documented formula");
    if (t > 0) {
      require(sample.window_start_before == samples[t - 1u].window_start_after,
              "buffered end-to-end S_t state was not carried to the next t");
    }
    if (sample.ordinary_service_used) {
      require(schedule_index < stats[4].level56_schedule_samples.size(),
              "buffered ordinary service has no linked schedule sample");
      const auto& schedule =
          stats[4].level56_schedule_samples[schedule_index++];
      require(schedule.buffered_fifo_enabled &&
                  schedule.buffered_service_time == sample.service_time &&
                  schedule.buffered_batch_id == sample.ordinary_batch_id &&
                  schedule.invocation ==
                      sample.ordinary_schedule_invocation,
              "buffered schedule sample is not linked to its t and batch");
      if (prior_codes_by_batch.size() <= schedule.buffered_batch_id) {
        prior_codes_by_batch.resize(schedule.buffered_batch_id + 1u);
      }
      const auto& prior_codes =
          prior_codes_by_batch[schedule.buffered_batch_id];
      for (std::size_t code_index = 0; code_index < schedule.codes.size();
           ++code_index) {
        const auto& code = schedule.codes[code_index];
        require(code.code_index == code_index &&
                    code.group_index == code_index / 4u &&
                    code.position_in_group == code_index % 4u,
                "buffered retry changed the fixed 64-code/group layout");
        require(code.source_global_row >= code.source_local_row,
                "buffered schedule exported an invalid global-row mapping");
        if (code.already_decoded) {
          require(code.writeback_complete && !code.pending &&
                      !code.forced_evicted,
                  "AlreadyDecoded code lacks terminal writeback state");
        }
        if (code.pending) {
          require(!code.already_decoded && !code.writeback_complete &&
                      !code.forced_evicted,
                  "pending code overlaps a terminal buffered state");
        }
        if (!prior_codes.empty()) {
          const auto& prior = prior_codes[code_index];
          require(code.code_index == prior.code_index &&
                      code.source_level == prior.source_level &&
                      code.source_local_row == prior.source_local_row &&
                      code.source_global_row == prior.source_global_row &&
                      code.group_index == prior.group_index &&
                      code.position_in_group == prior.position_in_group &&
                      code.early_stop_hit == prior.early_stop_hit &&
                      code.hybrid_class == prior.hybrid_class &&
                      code.resource_eligibility ==
                          prior.resource_eligibility,
                  "pending retry refreshed frozen classification or row mapping");
          if (prior.already_decoded) {
            require(code.already_decoded && code.writeback_complete &&
                        !code.planned_hiso && !code.planned_siso &&
                        code.assigned_entry_slot == -1 &&
                        code.assigned_core == -1 && !code.pending,
                    "AlreadyDecoded code re-entered a later ordinary schedule");
          }
        }
      }
      prior_codes_by_batch[schedule.buffered_batch_id] = schedule.codes;
      ++ordinary_service_times;
    }
    observed_pending = observed_pending || sample.pending_after > 0;
    observed_window_fallback =
        observed_window_fallback ||
        sample.window_start_after < p.LEVEL56_BUFFER_ROWS;
    previous_fifo_depth = sample.fifo_depth_after;
  }

  require(stats[4].level56_schedule_samples.size() ==
              ordinary_service_times &&
              schedule_index == ordinary_service_times,
          "buffered end-to-end ordinary service/sample count mismatch");
  for (const auto& schedule : stats[4].level56_schedule_samples) {
    require(schedule.codes.size() == 64 &&
                schedule.total_group_entries <= 8 &&
                schedule.planned_hiso_count <= 8 &&
                schedule.planned_siso_count <= 8,
            "buffered end-to-end ordinary service exceeded eight entries");
  }
  require(samples.front().pending_after > 0 && observed_pending &&
              observed_window_fallback,
          "buffered end-to-end input did not exercise pending window fallback");
  for (std::size_t row = 0; row < decoded.rows(); ++row) {
    for (std::size_t col = 0; col < decoded.cols(); ++col) {
      require(std::isfinite(decoded[row][col]),
              "buffered end-to-end decode produced a non-finite LLR");
    }
  }

  auto boundary_params = p;
  boundary_params.LEVEL56_BUFFER_ROWS = 0;
  matrix::Matrix<float> boundary_llr(
      boundary_params.win_height_rows(), llr.cols());
  for (std::size_t row = 0; row < boundary_llr.rows(); ++row) {
    for (std::size_t col = 0; col < boundary_llr.cols(); ++col) {
      boundary_llr[row][col] = llr[row][col];
    }
  }
  std::vector<newcode::TileEarlyStopCounter> boundary_stats;
  (void)newcode::ofec_decode_llr_plain(
      boundary_llr, boundary_params, &boundary_stats, true, nullptr);
  require(boundary_stats.size() == boundary_params.TILES_PER_WIN &&
              boundary_stats[4].level56_buffered_time_samples.size() == 1,
          "buffered boundary decode did not emit its t sample");
  const auto& boundary =
      boundary_stats[4].level56_buffered_time_samples.front();
  require(boundary.window_start_before == 0 &&
              boundary.window_start_after == 0 &&
              boundary.completed_batches == 0 &&
              boundary.forced_evicted_batches == 1 &&
              boundary.fifo_depth_after == 0 &&
              boundary.retirements.size() == 1 &&
              boundary.retirements.front().reason ==
                  newcode::Level56BufferedBatchRetirement::ForcedEvicted &&
              boundary.forced_evicted_global_rows.size() ==
                  samples.front().pending_after,
          "buffered S_t=0 boundary did not force only unfinished codes");
}

void check_full_frame_ber_keeps_all_evaluation_bits() {
  newcode::Params p;
  p.TILES_PER_WIN = 2;
  p.TILE_HEIGHT_BR = 1;
  p.TILE_OVERLAP_BR = 0;

  constexpr std::size_t kInfoBitsPerRow =
      newcode::Params::BCH_K -
      newcode::Params::NUM_SUBBLOCK_COLS *
          newcode::Params::BITS_PER_SUBBLOCK_DIM;
  const std::size_t window_bits = p.win_height_rows() * kInfoBitsPerRow;
  const std::size_t tile_bits = p.tile_height_rows() * kInfoBitsPerRow;
  // The fixed BER interval excludes four leading and two trailing windows.
  // Keep two complete windows in the middle and place both errors there.
  const std::size_t bit_count = 8u * window_bits;
  const std::size_t first_error = 4u * window_bits + 10u;
  const std::size_t second_error =
      4u * window_bits + tile_bits + 10u;

  std::vector<uint8_t> reference(bit_count, 0u);
  std::vector<uint8_t> received(bit_count, 0u);
  received[first_error] = 1u;
  received[second_error] = 1u;

  std::vector<std::size_t> error_positions;
  const auto aggregate =
      newcode::compute_ber(reference, received, p, &error_positions);
  require(aggregate.errors == 2 && aggregate.total == 2u * window_bits &&
              error_positions ==
                  std::vector<std::size_t>({first_error, second_error}),
          "full-frame post-FEC BER must keep every bit in the standard evaluation interval");

  const auto windows =
      newcode::compute_ber_per_window(reference, received, p);
  require(windows.size() == 8 && windows[4].errors == 2 &&
              windows[4].total == window_bits,
          "window post-FEC BER must not exclude forced-evicted positions");

  const auto tiles =
      newcode::compute_ber_per_tile_window(reference, received, p);
  require(tiles.size() == 16 && tiles[8].errors == 1 &&
              tiles[8].total == tile_bits &&
              tiles[9].errors == 1 && tiles[9].total == tile_bits,
          "tile post-FEC BER must not exclude forced-evicted positions");
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
    check_twomain_parameterized_output();
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
    check_temporal_rule_helpers();
    check_temporal_two_stage_budget();
    check_temporal_lookahead_three_windows();
    check_buffered_window_transition_formula();
    check_buffered_already_decoded_three_round_example();
    check_buffered_retirement_markers();
    check_buffered_service_fast_path_and_fifo_order();
    check_buffered_service_ordinary_budget_boundary();
    check_buffered_service_forced_eviction_boundary();
    check_buffered_configuration_validation();
    check_buffered_end_to_end_full_early_stop();
    check_buffered_end_to_end_multiple_windows();
    check_full_frame_ber_keeps_all_evaluation_bits();
    std::cout << "LEVEL56 shared scheduler regression checks passed\n";
    return 0;
  } catch (const std::exception& ex) {
    std::cerr << "LEVEL56 shared scheduler regression check failed: "
              << ex.what() << '\n';
    return 1;
  }
}
