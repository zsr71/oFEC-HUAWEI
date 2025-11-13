#include "ofec_sweep_detail.hpp"

#include <algorithm>
#include <cmath>
#include <iomanip>
#include <random>
#include <sstream>
#include <unordered_set>

namespace ofec_sweep {
namespace detail {

std::vector<float> generate_sequence(float start, float step, std::size_t length) {
  std::vector<float> seq(length, start);
  for (std::size_t i = 0; i < length; ++i) {
    seq[i] = start + step * static_cast<float>(i);
  }
  return seq;
}

float infer_step(const std::vector<float>& values) {
  return values.size() >= 2 ? values[1] - values[0] : 0.0f;
}

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

newcode::PipelineConfig make_pipeline_config(const SweepParameterConfig& config) {
  newcode::PipelineConfig cfg;
  cfg.decoder_name = config.decoder_name;
  cfg.interleaver_name = config.interleaver_name;
  cfg.normalize_extrinsic = config.normalize_extrinsic;
  cfg.bits_per_symbol = config.bits_per_symbol;
  cfg.quiet = config.quiet_pipeline;
  return cfg;
}

std::vector<SweepScenario> build_scenarios(const SweepParameterConfig& config,
                                           const std::vector<float>& ebn0_candidates,
                                           const std::vector<int>& bitgen_seeds,
                                           const std::vector<int>& channel_seeds) {
  const auto& base_params = config.base_params;
  std::vector<SweepScenario> scenarios;

  SweepScenario baseline;
  baseline.name = "baseline";
  baseline.alpha_list = base_params.ALPHA_LIST;
  baseline.beta_list = base_params.beta_list;
  if (!baseline.alpha_list.empty()) {
    baseline.alpha_start = baseline.alpha_list.front();
    baseline.alpha_step = infer_step(baseline.alpha_list);
  }
  if (!baseline.beta_list.empty()) {
    baseline.beta_start = baseline.beta_list.front();
    baseline.beta_step = infer_step(baseline.beta_list);
  }
  baseline.chase_L = base_params.CHASE_L;
  baseline.chase_n_test = 1 << base_params.CHASE_L;
  baseline.bitgen_seed = base_params.BITGEN_SEED;
  baseline.channel_seed = base_params.CHANNEL_SEED;
  baseline.ebn0_db = !ebn0_candidates.empty()
                       ? ebn0_candidates.front()
                       : newcode::DEFAULT_EBN0_DB;
  scenarios.push_back(baseline);

  auto matches_baseline = [&](const SweepScenario& candidate) -> bool {
    const SweepScenario& base = scenarios.front();
    return candidate.alpha_list == base.alpha_list &&
           candidate.beta_list == base.beta_list &&
           candidate.chase_L == base.chase_L &&
           candidate.bitgen_seed == base.bitgen_seed &&
           candidate.channel_seed == base.channel_seed &&
           std::fabs(candidate.ebn0_db - base.ebn0_db) < 1e-6f;
  };

  const std::size_t len = base_params.TILES_PER_WIN;
  if (config.explicit_patterns.empty()) {
    for (float alpha_start : config.alpha_start_candidates) {
      for (float alpha_step : config.alpha_step_candidates) {
        for (float beta_start : config.beta_start_candidates) {
          for (float beta_step : config.beta_step_candidates) {
            for (int chase_L : config.chase_l_candidates) {
              for (float ebn0_db : ebn0_candidates) {
                for (int bitgen_seed : bitgen_seeds) {
                  for (int channel_seed : channel_seeds) {
                    SweepScenario scenario;
                    scenario.alpha_start = alpha_start;
                    scenario.alpha_step = alpha_step;
                    scenario.beta_start = beta_start;
                    scenario.beta_step = beta_step;
                    scenario.chase_L = chase_L;
                    scenario.chase_n_test = 1 << chase_L;
                    scenario.bitgen_seed = bitgen_seed;
                    scenario.channel_seed = channel_seed;
                    scenario.ebn0_db = ebn0_db;
                    scenario.alpha_list = generate_sequence(alpha_start, alpha_step, len);
                    scenario.beta_list = generate_sequence(beta_start, beta_step, len);

                    if (matches_baseline(scenario)) {
                      continue;
                    }

                    std::ostringstream oss;
                    oss << std::fixed << std::setprecision(3)
                        << "alphaS" << alpha_start << "_d" << alpha_step
                        << "_betaS" << beta_start << "_d" << beta_step
                        << "_chL" << chase_L
                        << "_EbN0_" << ebn0_db
                        << "_bitSeed" << bitgen_seed
                        << "_chanSeed" << channel_seed;
                    scenario.name = oss.str();
                    scenarios.push_back(std::move(scenario));
                  }
                }
              }
            }
          }
        }
      }
    }
  }

  for (std::size_t idx = 0; idx < config.explicit_patterns.size(); ++idx) {
    const auto& pattern = config.explicit_patterns[idx];
    if (pattern.alpha_list.size() != len || pattern.beta_list.size() != len) {
      continue;
    }

    const std::string base_label = pattern.label.empty()
                                     ? ("explicit" + std::to_string(idx))
                                     : pattern.label;

    for (int chase_L : config.chase_l_candidates) {
      for (float ebn0_db : ebn0_candidates) {
        for (int bitgen_seed : bitgen_seeds) {
          for (int channel_seed : channel_seeds) {
            SweepScenario scenario;
            scenario.alpha_list = pattern.alpha_list;
            scenario.beta_list = pattern.beta_list;
            if (!scenario.alpha_list.empty()) {
              scenario.alpha_start = scenario.alpha_list.front();
              scenario.alpha_step = infer_step(scenario.alpha_list);
            }
            if (!scenario.beta_list.empty()) {
              scenario.beta_start = scenario.beta_list.front();
              scenario.beta_step = infer_step(scenario.beta_list);
            }
            scenario.chase_L = chase_L;
            scenario.chase_n_test = 1 << chase_L;
            scenario.bitgen_seed = bitgen_seed;
            scenario.channel_seed = channel_seed;
            scenario.ebn0_db = ebn0_db;

            if (matches_baseline(scenario)) {
              continue;
            }

            std::ostringstream oss;
            oss << base_label
                << "_chL" << chase_L
                << "_EbN0_" << std::fixed << std::setprecision(3) << ebn0_db
                << "_bitSeed" << bitgen_seed
                << "_chanSeed" << channel_seed;
            scenario.name = oss.str();
            scenarios.push_back(std::move(scenario));
          }
        }
      }
    }
  }

  return scenarios;
}

}  // namespace detail
}  // namespace ofec_sweep
