#include "newcode/ofec/mux/mux_state_builder.hpp"

namespace newcode::mux {

std::vector<uint8_t> build_state_from_early_stop(
    const TileEarlyStopResult& stats) {
  std::vector<uint8_t> state(stats.row_passed_flags.size(),
                             static_cast<uint8_t>(StateTag::NeedSiso));
  for (std::size_t i = 0; i < stats.row_passed_flags.size(); ++i) {
    if (stats.row_passed_flags[i]) {
      state[i] = static_cast<uint8_t>(StateTag::EarlyStopped);
    }
  }
  return state;
}

StateCount count_state012(const std::vector<uint8_t>& state) {
  StateCount cnt{};
  for (const uint8_t tag : state) {
    if (tag == static_cast<uint8_t>(StateTag::NeedSiso)) {
      ++cnt.n0_need_siso;
    } else if (tag == static_cast<uint8_t>(StateTag::EarlyStopped)) {
      ++cnt.n1_early_stopped;
    } else if (tag == static_cast<uint8_t>(StateTag::Unscheduled)) {
      ++cnt.n2_unscheduled;
    }
  }
  return cnt;
}

}  // namespace newcode::mux

