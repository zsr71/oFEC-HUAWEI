#include "ofec_sweep_detail.hpp"

namespace ofec_sweep {
namespace detail {

std::vector<float> build_ebn0_values(const SweepParameterConfig& config) {
  if (!config.ebn0_candidates.empty()) {
    return config.ebn0_candidates;
  }

  if (config.ebn0_points <= 0) {
    return {newcode::DEFAULT_EBN0_DB};
  }

  if (config.ebn0_points == 1) {
    return {config.ebn0_start};
  }

  std::vector<float> values;
  values.reserve(static_cast<std::size_t>(config.ebn0_points));
  const float step = (config.ebn0_end - config.ebn0_start) /
                     static_cast<float>(config.ebn0_points - 1);
  for (int i = 0; i < config.ebn0_points; ++i) {
    values.push_back(config.ebn0_start + step * static_cast<float>(i));
  }
  return values;
}

}  // namespace detail
}  // namespace ofec_sweep
