#include <iostream>
#include <limits>
#include <sstream>
#include <string>
#include <vector>

#include "newcode/params.hpp"
#include "newcode/pipeline_runner.hpp"

using namespace newcode;

struct SweepScenario {
  std::string name;
  std::vector<float> alpha_list;
  std::vector<float> beta_list;
};

static std::vector<SweepScenario> build_scenarios(const Params& base_params,
                                                  const std::vector<float>& alpha_candidates,
                                                  const std::vector<float>& beta_candidates)
{
  std::vector<SweepScenario> scenarios;
  scenarios.push_back({"baseline", base_params.ALPHA_LIST, base_params.beta_list});

  for (float alpha_val : alpha_candidates) {
    for (float beta_val : beta_candidates) {
      SweepScenario sc;
      std::ostringstream oss;
      oss << "alpha" << alpha_val << "_beta" << beta_val;
      sc.name = oss.str();
      sc.alpha_list.assign(base_params.TILES_PER_WIN, alpha_val);
      sc.beta_list.assign(base_params.TILES_PER_WIN, beta_val);
      scenarios.push_back(std::move(sc));
    }
  }
  return scenarios;
}

int main()
{
  Params base_params;

  // 在此调整搜索范围：每个候选值会在所有 tile 上使用同一数值
  //const std::vector<float> alpha_candidates = {0.001f,0.01f,0.05f,0.07f,0.1f, 0.2f, 0.3f,0.4f, 0.5f,1.0f,10.0f};
  //const std::vector<float> beta_candidates   = {0.1f,1.0f,5.0f,10.0f,25.0f,50.0f, 100.0f, 200.0f, 300.0f,500.0f};
  const std::vector<float> alpha_candidates = {0.2f};
  const std::vector<float> beta_candidates   = {25.0f};
  auto scenarios = build_scenarios(base_params, alpha_candidates, beta_candidates);

  double best_post_ber = std::numeric_limits<double>::infinity();
  std::size_t best_index = static_cast<std::size_t>(-1);
  PipelineResult best_result{};
  std::vector<std::string> scenario_summaries;
  scenario_summaries.reserve(scenarios.size());

  for (std::size_t idx = 0; idx < scenarios.size(); ++idx) {
    const auto& scenario = scenarios[idx];
    if (scenario.alpha_list.size() != base_params.TILES_PER_WIN ||
        scenario.beta_list.size() != base_params.TILES_PER_WIN) {
      std::cerr << "[WARN] Scenario '" << scenario.name
                << "' skipped due to list size mismatch (expected "
                << base_params.TILES_PER_WIN << ")\n";
      continue;
    }

    Params params = base_params;
    params.ALPHA_LIST = scenario.alpha_list;
    params.beta_list  = scenario.beta_list;
    if (!params.ALPHA_LIST.empty()) params.ALPHA = params.ALPHA_LIST.front();
    if (!params.beta_list.empty())  params.beta  = params.beta_list.front();

    auto result = run_pipeline(params, scenario.name);

    std::ostringstream summary;
    summary << "[SUMMARY] " << scenario.name
            << " Pre-FEC BER=" << result.pre_fec.ber
            << " (errs=" << result.pre_fec.errors << "/" << result.pre_fec.total << ")"
            << " | Post-FEC BER=" << result.post_fec.ber
            << " (errs=" << result.post_fec.errors << "/" << result.post_fec.total << ")";
    scenario_summaries.push_back(summary.str());

    if (result.post_fec.total > 0 && result.post_fec.ber < best_post_ber) {
      best_post_ber = result.post_fec.ber;
      best_index = idx;
      best_result = result;
    }
  }

  if (best_index == static_cast<std::size_t>(-1)) {
    std::cerr << "[ERROR] No valid scenarios evaluated.\n";
    return 1;
  }

  std::cout << "\n[SUMMARY] All scenarios:\n";
  for (const auto& line : scenario_summaries) {
    std::cout << "  " << line << '\n';
  }

  const auto& best_scenario = scenarios[best_index];
  std::cout << "\n[RESULT] Best scenario: " << best_scenario.name
            << " with Post-FEC BER=" << best_result.post_fec.ber
            << " (errs=" << best_result.post_fec.errors << "/" << best_result.post_fec.total << ")\n";

  std::cout << "[RESULT] Best ALPHA_LIST: ";
  for (size_t i = 0; i < best_scenario.alpha_list.size(); ++i) {
    std::cout << best_scenario.alpha_list[i]
              << (i + 1 < best_scenario.alpha_list.size() ? ", " : "\n");
  }

  std::cout << "[RESULT] Best beta_list: ";
  for (size_t i = 0; i < best_scenario.beta_list.size(); ++i) {
    std::cout << best_scenario.beta_list[i]
              << (i + 1 < best_scenario.beta_list.size() ? ", " : "\n");
  }

  return 0;
}







