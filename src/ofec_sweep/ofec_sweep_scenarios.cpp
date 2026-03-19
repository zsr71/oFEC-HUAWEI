#include "ofec_sweep_detail.hpp"

#include <algorithm>
#include <cmath>
#include <iomanip>
#include <sstream>

namespace ofec_sweep {
namespace detail {
namespace {

template <typename T>
std::vector<T> choose_candidates(const std::vector<T>& provided, const T& fallback) {
  if (!provided.empty()) {
    return provided;
  }
  return {fallback};
}

std::vector<float> default_action_beta_list(const SweepParameterConfig& config,
                                            const SweepScenario& scenario,
                                            std::size_t tiles) {
  if (std::isfinite(config.early_stop_action_sign_beta_fill)) {
    return std::vector<float>(tiles, config.early_stop_action_sign_beta_fill);
  }
  return scenario.beta_list;
}

void finalize_beta_lists(SweepScenario& scenario) {
  if (!scenario.beta_list.empty()) {
    scenario.beta_start = scenario.beta_list.front();
    scenario.beta_step = infer_step(scenario.beta_list);
  }
  if (!scenario.early_stop_action_sign_beta_list.empty()) {
    scenario.early_stop_beta_start = scenario.early_stop_action_sign_beta_list.front();
    scenario.early_stop_beta_step = infer_step(scenario.early_stop_action_sign_beta_list);
  } else {
    scenario.early_stop_beta_start = scenario.beta_start;
    scenario.early_stop_beta_step = scenario.beta_step;
  }
}

void append_common_name(std::ostringstream& oss, const SweepScenario& scenario) {
  oss << "_chN" << scenario.chase_n_test;
  if (scenario.decoder_name == "chase_topk_pruned") {
    oss << "_topk" << scenario.chase_topk_keep;
  }
  if (scenario.decoder_name == "chase_group_minima") {
    oss << "_grpBits" << scenario.chase_group_minima_bits;
  }
  oss << "_cond" << scenario.early_stop_condition_mode
      << "_act" << scenario.early_stop_action_mode;
  if (scenario.early_stop_condition_mode == 1) {
    oss << "_v1bch" << (scenario.early_stop_cond_v1_require_bch ? 1 : 0)
        << "_v1ov" << (scenario.early_stop_cond_v1_require_overall ? 1 : 0);
  } else if (scenario.early_stop_condition_mode == 2) {
    oss << "_v2thr" << std::fixed << std::setprecision(3)
        << scenario.early_stop_v2_llr_abs_threshold
        << "_v2unr" << scenario.early_stop_v2_max_unreliable_bits
        << "_v2ov" << (scenario.early_stop_cond_v2_include_overall ? 1 : 0)
        << std::defaultfloat;
  }
  if (scenario.early_stop_action_mode == 1) {
    oss << "_esBetaS" << std::fixed << std::setprecision(3)
        << scenario.early_stop_beta_start
        << "_d" << scenario.early_stop_beta_step
        << std::defaultfloat;
  } else if (scenario.early_stop_action_mode == 3) {
    oss << "_hardMag" << std::fixed << std::setprecision(3)
        << scenario.early_stop_action_hard_llr_mag
        << std::defaultfloat;
  }
}

template <typename ScenarioBuilder>
void enumerate_common_axes(const SweepParameterConfig& config,
                           const std::vector<std::string>& decoder_names,
                           const std::vector<int>& chase_l_values,
                           const std::vector<int>& chase_n_test_values,
                           const std::vector<int>& condition_modes,
                           const std::vector<int>& action_modes,
                           const std::vector<float>& ebn0_candidates,
                           const std::vector<int>& bitgen_seeds,
                           const std::vector<int>& channel_seeds,
                           ScenarioBuilder&& build_one) {
  for (const auto& decoder_name : decoder_names) {
    const std::vector<int> topk_keep_values =
        (decoder_name == "chase_topk_pruned")
            ? choose_candidates(config.chase_topk_keep_candidates,
                                config.chase_topk_keep)
            : std::vector<int>{config.chase_topk_keep};
    const std::vector<int> group_minima_bits_values =
        (decoder_name == "chase_group_minima")
            ? choose_candidates(config.chase_group_minima_bits_candidates,
                                config.chase_group_minima_bits)
            : std::vector<int>{config.chase_group_minima_bits};

    for (int chase_L : chase_l_values) {
      for (int chase_n_test : chase_n_test_values) {
        for (int group_minima_bits : group_minima_bits_values) {
          for (int condition_mode : condition_modes) {
            const auto v1_require_bch_values =
                (condition_mode == 1)
                    ? choose_candidates(
                          config.early_stop_cond_v1_require_bch_candidates,
                          config.early_stop_cond_v1_require_bch)
                    : std::vector<bool>{config.early_stop_cond_v1_require_bch};
            const auto v1_require_overall_values =
                (condition_mode == 1)
                    ? choose_candidates(
                          config.early_stop_cond_v1_require_overall_candidates,
                          config.early_stop_cond_v1_require_overall)
                    : std::vector<bool>{config.early_stop_cond_v1_require_overall};
            const auto v2_threshold_values =
                (condition_mode == 2)
                    ? choose_candidates(
                          config.early_stop_v2_llr_abs_threshold_candidates,
                          config.early_stop_v2_llr_abs_threshold)
                    : std::vector<float>{config.early_stop_v2_llr_abs_threshold};
            const auto v2_unreliable_values =
                (condition_mode == 2)
                    ? choose_candidates(
                          config.early_stop_v2_max_unreliable_bits_candidates,
                          config.early_stop_v2_max_unreliable_bits)
                    : std::vector<int>{config.early_stop_v2_max_unreliable_bits};

            for (int action_mode : action_modes) {
              const auto hard_llr_mag_values =
                  (action_mode == 3)
                      ? choose_candidates(
                            config.early_stop_action_hard_llr_mag_candidates,
                            config.early_stop_action_hard_llr_mag)
                      : std::vector<float>{config.early_stop_action_hard_llr_mag};

              for (int topk_keep : topk_keep_values) {
                for (bool v1_require_bch : v1_require_bch_values) {
                  for (bool v1_require_overall : v1_require_overall_values) {
                    for (float v2_threshold : v2_threshold_values) {
                      for (int v2_unreliable : v2_unreliable_values) {
                        for (float hard_llr_mag : hard_llr_mag_values) {
                          for (float ebn0_db : ebn0_candidates) {
                            for (int bitgen_seed : bitgen_seeds) {
                              for (int channel_seed : channel_seeds) {
                                build_one(decoder_name,
                                          chase_L,
                                          chase_n_test,
                                          group_minima_bits,
                                          condition_mode,
                                          action_mode,
                                          v1_require_bch,
                                          v1_require_overall,
                                          v2_threshold,
                                          v2_unreliable,
                                          hard_llr_mag,
                                          topk_keep,
                                          ebn0_db,
                                          bitgen_seed,
                                          channel_seed);
                              }
                            }
                          }
                        }
                      }
                    }
                  }
                }
              }
            }
          }
        }
      }
    }
  }
}

}  // namespace

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

  const std::size_t len = base_params.TILES_PER_WIN;
  const std::vector<std::string> decoder_names =
      choose_candidates(config.decoder_name_candidates, config.decoder_name);
  const std::vector<int> chase_l_values =
      choose_candidates(config.chase_l_candidates, base_params.CHASE_L);
  const std::vector<int> chase_n_test_values =
      choose_candidates(config.chase_n_test_candidates, config.chase_n_test);
  const std::vector<int> condition_modes =
      choose_candidates(config.early_stop_condition_candidates,
                        config.early_stop_condition_mode);
  const std::vector<int> action_modes =
      choose_candidates(config.early_stop_action_candidates,
                        config.early_stop_action_mode);

  if (config.explicit_patterns.empty()) {
    for (float alpha_start : config.alpha_start_candidates) {
      for (float alpha_step : config.alpha_step_candidates) {
        for (float beta_start : config.beta_start_candidates) {
          for (float beta_step : config.beta_step_candidates) {
            enumerate_common_axes(
                config,
                decoder_names,
                chase_l_values,
                chase_n_test_values,
                condition_modes,
                action_modes,
                ebn0_candidates,
                bitgen_seeds,
                channel_seeds,
                [&](const std::string& decoder_name,
                    int chase_L,
                    int chase_n_test,
                    int group_minima_bits,
                    int condition_mode,
                    int action_mode,
                    bool v1_require_bch,
                    bool v1_require_overall,
                    float v2_threshold,
                    int v2_unreliable,
                    float hard_llr_mag,
                    int topk_keep,
                    float ebn0_db,
                    int bitgen_seed,
                    int channel_seed) {
                  SweepScenario scenario;
                  scenario.decoder_name = decoder_name;
                  scenario.alpha_start = alpha_start;
                  scenario.alpha_step = alpha_step;
                  scenario.beta_start = beta_start;
                  scenario.beta_step = beta_step;
                  scenario.chase_L = chase_L;
                  scenario.chase_n_test = chase_n_test;
                  scenario.chase_topk_keep = topk_keep;
                  scenario.chase_group_minima_bits = group_minima_bits;
                  scenario.early_stop_condition_mode = condition_mode;
                  scenario.early_stop_action_mode = action_mode;
                  scenario.early_stop_cond_v1_require_bch = v1_require_bch;
                  scenario.early_stop_cond_v1_require_overall = v1_require_overall;
                  scenario.early_stop_v2_llr_abs_threshold = v2_threshold;
                  scenario.early_stop_v2_max_unreliable_bits = v2_unreliable;
                  scenario.early_stop_cond_v2_include_overall =
                      config.early_stop_cond_v2_include_overall;
                  scenario.early_stop_action_hard_llr_mag = hard_llr_mag;
                  scenario.bitgen_seed = bitgen_seed;
                  scenario.channel_seed = channel_seed;
                  scenario.ebn0_db = ebn0_db;
                  scenario.alpha_list =
                      generate_sequence(alpha_start, alpha_step, len);
                  scenario.beta_list =
                      generate_sequence(beta_start, beta_step, len);

                  if (action_mode == 1) {
                    const auto early_stop_beta_start_values =
                        choose_candidates(
                            config.early_stop_action_beta_start_candidates,
                            beta_start);
                    const auto early_stop_beta_step_values =
                        choose_candidates(
                            config.early_stop_action_beta_step_candidates,
                            beta_step);
                    for (float es_beta_start : early_stop_beta_start_values) {
                      for (float es_beta_step : early_stop_beta_step_values) {
                        SweepScenario expanded = scenario;
                        expanded.early_stop_action_sign_beta_list =
                            generate_sequence(es_beta_start, es_beta_step, len);
                        finalize_beta_lists(expanded);

                        std::ostringstream oss;
                        oss << decoder_name << "_"
                            << std::fixed << std::setprecision(3)
                            << "alphaS" << alpha_start << "_d" << alpha_step
                            << "_betaS" << beta_start << "_d" << beta_step
                            << "_chL" << chase_L;
                        append_common_name(oss, expanded);
                        oss << "_EbN0_" << ebn0_db
                            << "_bitSeed" << bitgen_seed
                            << "_chanSeed" << channel_seed;
                        expanded.name = oss.str();
                        scenarios.push_back(std::move(expanded));
                      }
                    }
                  } else {
                    scenario.early_stop_action_sign_beta_list =
                        default_action_beta_list(config, scenario, len);
                    finalize_beta_lists(scenario);

                    std::ostringstream oss;
                    oss << decoder_name << "_"
                        << std::fixed << std::setprecision(3)
                        << "alphaS" << alpha_start << "_d" << alpha_step
                        << "_betaS" << beta_start << "_d" << beta_step
                        << "_chL" << chase_L;
                    append_common_name(oss, scenario);
                    oss << "_EbN0_" << ebn0_db
                        << "_bitSeed" << bitgen_seed
                        << "_chanSeed" << channel_seed;
                    scenario.name = oss.str();
                    scenarios.push_back(std::move(scenario));
                  }
                });
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

    const std::string base_label =
        pattern.label.empty() ? ("explicit" + std::to_string(idx)) : pattern.label;

    enumerate_common_axes(
        config,
        decoder_names,
        chase_l_values,
        chase_n_test_values,
        condition_modes,
        action_modes,
        ebn0_candidates,
        bitgen_seeds,
        channel_seeds,
        [&](const std::string& decoder_name,
            int chase_L,
            int chase_n_test,
            int group_minima_bits,
            int condition_mode,
            int action_mode,
            bool v1_require_bch,
            bool v1_require_overall,
            float v2_threshold,
            int v2_unreliable,
            float hard_llr_mag,
            int topk_keep,
            float ebn0_db,
            int bitgen_seed,
            int channel_seed) {
          SweepScenario scenario;
          scenario.decoder_name = decoder_name;
          scenario.alpha_list = pattern.alpha_list;
          scenario.beta_list = pattern.beta_list;
          scenario.chase_L = chase_L;
          scenario.chase_n_test = chase_n_test;
          scenario.chase_topk_keep = topk_keep;
          scenario.chase_group_minima_bits = group_minima_bits;
          scenario.early_stop_condition_mode = condition_mode;
          scenario.early_stop_action_mode = action_mode;
          scenario.early_stop_cond_v1_require_bch = v1_require_bch;
          scenario.early_stop_cond_v1_require_overall = v1_require_overall;
          scenario.early_stop_v2_llr_abs_threshold = v2_threshold;
          scenario.early_stop_v2_max_unreliable_bits = v2_unreliable;
          scenario.early_stop_cond_v2_include_overall =
              config.early_stop_cond_v2_include_overall;
          scenario.early_stop_action_hard_llr_mag = hard_llr_mag;
          scenario.bitgen_seed = bitgen_seed;
          scenario.channel_seed = channel_seed;
          scenario.ebn0_db = ebn0_db;

          if (action_mode == 1 &&
              pattern.early_stop_action_sign_beta_list.size() == len) {
            scenario.early_stop_action_sign_beta_list =
                pattern.early_stop_action_sign_beta_list;
          } else {
            scenario.early_stop_action_sign_beta_list =
                default_action_beta_list(config, scenario, len);
          }

          if (!scenario.alpha_list.empty()) {
            scenario.alpha_start = scenario.alpha_list.front();
            scenario.alpha_step = infer_step(scenario.alpha_list);
          }
          finalize_beta_lists(scenario);

          std::ostringstream oss;
          oss << decoder_name << "_" << base_label
              << "_chL" << chase_L;
          append_common_name(oss, scenario);
          oss << "_EbN0_" << std::fixed << std::setprecision(3) << ebn0_db
              << "_bitSeed" << bitgen_seed
              << "_chanSeed" << channel_seed;
          scenario.name = oss.str();
          scenarios.push_back(std::move(scenario));
        });
  }

  return scenarios;
}

}  // namespace detail
}  // namespace ofec_sweep
