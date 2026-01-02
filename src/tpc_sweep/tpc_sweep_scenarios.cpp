#include "tpc_sweep_detail.hpp"

#include <iomanip>
#include <limits>
#include <random>
#include <sstream>
#include <unordered_set>

namespace tpc_sweep {
namespace detail {

std::vector<int> generate_random_seeds(int count) {
  std::vector<int> seeds;
  if (count <= 0) {
    return seeds;
  }

  std::random_device rd;
  std::mt19937 rng(rd());
  std::uniform_int_distribution<int> dist(1, std::numeric_limits<int>::max());
  std::unordered_set<int> seen;
  seen.reserve(static_cast<std::size_t>(count));
  seeds.reserve(static_cast<std::size_t>(count));

  while (seeds.size() < static_cast<std::size_t>(count)) {
    int candidate = dist(rng);
    if (seen.insert(candidate).second) {
      seeds.push_back(candidate);
    }
  }
  return seeds;
}

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

std::vector<SweepScenario> build_scenarios(const SweepParameterConfig& /*config*/,
                                           const std::vector<float>& ebn0_candidates,
                                           const std::vector<int>& bitgen_seeds,
                                           const std::vector<int>& channel_seeds) {
  std::vector<SweepScenario> scenarios;
  for (float ebn0_db : ebn0_candidates) {
    for (int bitgen_seed : bitgen_seeds) {
      for (int channel_seed : channel_seeds) {
        SweepScenario scenario;
        scenario.bitgen_seed = bitgen_seed;
        scenario.channel_seed = channel_seed;
        scenario.ebn0_db = ebn0_db;
        std::ostringstream oss;
        oss << std::fixed << std::setprecision(3)
            << "EbN0_" << ebn0_db
            << "_bitSeed" << bitgen_seed
            << "_chanSeed" << channel_seed;
        scenario.name = oss.str();
        scenarios.push_back(std::move(scenario));
      }
    }
  }
  return scenarios;
}

}  // namespace detail
}  // namespace tpc_sweep
