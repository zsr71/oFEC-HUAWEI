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

void check_unscheduled_history_is_unchanged() {
  newcode::detail::TilePrepared<float> prep;
  prep.lin_matrix = matrix::Matrix<float>(1, newcode::Params::BCH_N);
  prep.lch_matrix = matrix::Matrix<float>(1, newcode::Params::BCH_N);
  prep.row_local_lookup = {0};
  prep.row_global_lookup = {352};
  chase::DecoderCoreResult<float> decoded{
      matrix::Matrix<float>(1, newcode::Params::BCH_N), {false}};
  matrix::Matrix<float> tile_out(352, 128);
  matrix::Matrix<float> history(704, 128);
  for (std::size_t row = 0; row < history.rows(); ++row) {
    for (std::size_t col = 0; col < history.cols(); ++col) {
      history[row][col] = 7.0f;
    }
  }

  newcode::Params p;
  newcode::detail::writeback_tile(
      prep, decoded, p, 0, true, &tile_out, &history);
  for (std::size_t row = 0; row < history.rows(); ++row) {
    for (std::size_t col = 0; col < history.cols(); ++col) {
      require(history[row][col] == 7.0f,
              "produced=false row must not update last-tile history");
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

}  // namespace

int main() {
  try {
    check_full_siso_uses_no_hiso();
    check_hiso_only_handles_siso_overflow();
    check_siso_only_rows_never_use_hiso();
    check_early_stop_does_not_consume_shared_capacity();
    check_single_level_scheduler_bypasses_other_level();
    check_single_level_selection_uses_fewer_early_stops();
    check_common_parameter_validation();
    check_per_level_siso_postprocessing();
    check_unscheduled_history_is_unchanged();
    check_level_priority_modes();
    check_end_to_end_single_window();
    std::cout << "LEVEL56 shared scheduler regression checks passed\n";
    return 0;
  } catch (const std::exception& ex) {
    std::cerr << "LEVEL56 shared scheduler regression check failed: "
              << ex.what() << '\n';
    return 1;
  }
}
