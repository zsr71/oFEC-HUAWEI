#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>

#include "newcode/ofec_sweep_runner.hpp"
#include "ofec_sweep_detail.hpp"

namespace {

constexpr size_t kTilesPerWindow = 4;
constexpr float kEvalEbN0 = 3.07f;
constexpr size_t kStage1Bits = 3 * 110 * 16 * 111;   // fast coarse sweep
constexpr size_t kStage2Bits = 6 * 110 * 16 * 111;   // slow fine sweep
constexpr float kKeepRatio = 0.20f;      // keep top 20%

static constexpr const char* kInterleaverName = "identity";
static constexpr const char* kDecoderName = "plain";
static constexpr unsigned    kBitsPerSymbol = 1;
static constexpr bool        kNormalizeExtrinsic = false;
static constexpr bool        kGenerateRandomBits = true;
static constexpr bool        kNormalizeKnownPrefixTail = false;


// Default scan grids (can be tweaked before calling build_round0_shapes).
const std::vector<float> kAlphaLowGrid      = {0.30f, 0.40f, 0.50f};
const std::vector<float> kAlphaHighGrid     = {0.50f, 0.60f, 0.70f};
const std::vector<float> kBetaLowGrid       = {0.60f, 0.70f, 0.80f};
const std::vector<float> kBetaHighGrid      = {0.80f, 0.90f, 1.00f};
const std::vector<float> kGammaAlphaGrid    = {0.7f, 1.0f, 1.5f};
const std::vector<float> kGammaBetaGrid     = {0.7f, 1.0f, 1.5f};

struct Shape {
  float alpha_low;
  float alpha_high;
  float gamma_alpha;
  float beta_low;
  float beta_high;
  float gamma_beta;
  std::string tag;
};

std::vector<float> build_sequence(float low, float high, float gamma, size_t count) {
  std::vector<float> seq(count, low);
  if (count <= 1) return seq;
  const float span = std::max(0.0f, high - low);
  for (size_t i = 0; i < count; ++i) {
    float t = (count == 1) ? 0.0f : static_cast<float>(i) / static_cast<float>(count - 1);
    float shaped = (gamma == 1.0f) ? t : std::pow(std::clamp(t, 0.0f, 1.0f), gamma);
    seq[i] = low + span * shaped;
  }
  return seq;
}

// neighbor_values was used by a removed refinement stage; no longer needed.

void ensure_sweep2_csv_header(const std::string& csv_path) {
  std::ifstream fin(csv_path);
  if (fin.good() && fin.peek() != std::ifstream::traits_type::eof()) {
    return;
  }
  std::ofstream fout(csv_path, std::ios::out | std::ios::app);
  fout << "timestamp,run_id,stage,num_bits,label,alpha_start,beta_start,alpha_step,beta_step,"
          "chase_L,chase_n_test,bitgen_seed,channel_seed,ebn0_db,alpha_list,beta_list,"
          "pre_ber,pre_errs,pre_total,post_ber,post_errs,post_total\n";
}

void write_csv_row(std::ofstream& csv,
                   const std::string& run_id,
                   const std::string& stage_tag,
                   size_t num_bits,
                   const ofec_sweep::detail::SweepScenario& scenario,
                   const newcode::PipelineResult& result) {
  const float alpha_step = ofec_sweep::detail::infer_step(scenario.alpha_list);
  const float beta_step = ofec_sweep::detail::infer_step(scenario.beta_list);
  csv << ofec_sweep::detail::now_stamp() << ","
      << run_id << ","
      << stage_tag << ","
      << num_bits << ","
      << scenario.name << ","
      << scenario.alpha_start << ","
      << scenario.beta_start << ","
      << alpha_step << ","
      << beta_step << ","
      << scenario.chase_L << ","
      << scenario.chase_n_test << ","
      << scenario.bitgen_seed << ","
      << scenario.channel_seed << ","
      << scenario.ebn0_db << ","
      << '"' << ofec_sweep::detail::join_vec(scenario.alpha_list, '|', 6) << "\","
      << '"' << ofec_sweep::detail::join_vec(scenario.beta_list, '|', 6) << "\","
      << result.pre_fec.ber << ","
      << result.pre_fec.errors << ","
      << result.pre_fec.total << ","
      << result.post_fec.ber << ","
      << result.post_fec.errors << ","
      << result.post_fec.total << "\n";
  csv.flush();
}

ofec_sweep::ExplicitAlphaBetaPattern shape_to_pattern(const Shape& shape,
                                                      const std::string& phase,
                                                      size_t ordinal) {
  auto alpha_seq = build_sequence(shape.alpha_low, shape.alpha_high,
                                  shape.gamma_alpha, kTilesPerWindow);
  auto beta_seq = build_sequence(shape.beta_low, shape.beta_high,
                                 shape.gamma_beta, kTilesPerWindow);

  std::ostringstream oss;
  oss << "sweep2_" << phase << "_" << ordinal
      << "_a(" << shape.alpha_low << "," << shape.alpha_high << "," << shape.gamma_alpha
      << ")_b(" << shape.beta_low << "," << shape.beta_high << "," << shape.gamma_beta << ")"
      << shape.tag;

  ofec_sweep::ExplicitAlphaBetaPattern pattern;
  pattern.label = oss.str();
  pattern.alpha_list = alpha_seq;
  pattern.beta_list = beta_seq;
  return pattern;
}

std::vector<Shape> build_round0_shapes(const std::vector<float>& alpha_low_grid,
                                       const std::vector<float>& alpha_high_grid,
                                       const std::vector<float>& gamma_alpha_grid,
                                       const std::vector<float>& beta_low_grid,
                                       const std::vector<float>& beta_high_grid,
                                       const std::vector<float>& gamma_beta_grid) {
  std::vector<Shape> shapes;
  for (float aL : alpha_low_grid) {
    for (float aH : alpha_high_grid) {
      if (aH - aL < 0.05f) continue;  // ensure meaningful span
      for (float bL : beta_low_grid) {
        for (float bH : beta_high_grid) {
          if (bH - bL < 0.05f) continue;
          for (float gA : gamma_alpha_grid) {
            for (float gB : gamma_beta_grid) {
              Shape shape;
              shape.alpha_low = aL;
              shape.alpha_high = aH;
              shape.gamma_alpha = gA;
              shape.beta_low = bL;
              shape.beta_high = bH;
              shape.gamma_beta = gB;
              shape.tag = "_r0";
              shapes.push_back(shape);
            }
          }
        }
      }
    }
  }
  return shapes;
}

// build_round1_shapes was used to refine endpoints and gammas; removed per request.

// build_micro_patterns removed per request.

std::vector<ofec_sweep::ExplicitAlphaBetaPattern> build_schedule() {
  std::vector<ofec_sweep::ExplicitAlphaBetaPattern> patterns;
  auto round0 = build_round0_shapes(kAlphaLowGrid,
                                    kAlphaHighGrid,
                                    kGammaAlphaGrid,
                                    kBetaLowGrid,
                                    kBetaHighGrid,
                                    kGammaBetaGrid);
  for (size_t i = 0; i < round0.size(); ++i) {
    patterns.push_back(shape_to_pattern(round0[i], "r0", i));
  }
  return patterns;
}

ofec_sweep::SweepParameterConfig build_base_config() {
  ofec_sweep::SweepParameterConfig config;
  config.base_params.TILES_PER_WIN = kTilesPerWindow;
  config.base_params.ALPHA_LIST.assign(kTilesPerWindow, 0.3f);
  config.base_params.beta_list.assign(kTilesPerWindow, 0.6f);
  config.base_params.HARD_TILE_LIST.assign(kTilesPerWindow, 0);
  config.base_params.BITGEN_RANDOM_BITS = kGenerateRandomBits;
  config.base_params.NORMALIZE_KNOWN_PREFIX_TAIL = kNormalizeKnownPrefixTail;

  config.interleaver_name = kInterleaverName;
  config.decoder_name = kDecoderName;
  config.bits_per_symbol = kBitsPerSymbol;
  config.normalize_extrinsic = kNormalizeExtrinsic;
  config.generate_random_bits = kGenerateRandomBits;
  config.normalize_known_prefix_tail = kNormalizeKnownPrefixTail;

  config.chase_l_candidates = {config.base_params.CHASE_L};
  config.alpha_start_candidates.clear();
  config.alpha_step_candidates.clear();
  config.beta_start_candidates.clear();
  config.beta_step_candidates.clear();
  config.bitgen_seed_count = 1;
  config.channel_seed_count = 1;
  config.ebn0_start = kEvalEbN0;
  config.ebn0_end = kEvalEbN0;
  config.ebn0_points = 1;
  config.explicit_patterns.clear();
  config.base_params.debug_trace = {};
  return config;
}

struct Evaluation {
  ofec_sweep::ExplicitAlphaBetaPattern pattern;
  newcode::PipelineResult result;
};

ofec_sweep::detail::SweepScenario make_scenario(
    const ofec_sweep::ExplicitAlphaBetaPattern& pattern,
    const ofec_sweep::SweepParameterConfig& config) {
  ofec_sweep::detail::SweepScenario sc;
  sc.name = pattern.label;
  sc.alpha_list = pattern.alpha_list;
  sc.beta_list = pattern.beta_list;
  if (!sc.alpha_list.empty()) {
    sc.alpha_start = sc.alpha_list.front();
  }
  if (!sc.beta_list.empty()) {
    sc.beta_start = sc.beta_list.front();
  }
  sc.chase_L = config.base_params.CHASE_L;
  sc.chase_n_test = 1 << sc.chase_L;
  sc.bitgen_seed = config.base_params.BITGEN_SEED;
  sc.channel_seed = config.base_params.CHANNEL_SEED;
  sc.ebn0_db = kEvalEbN0;
  return sc;
}

std::vector<Evaluation> run_stage(const std::vector<ofec_sweep::ExplicitAlphaBetaPattern>& patterns,
                                  const ofec_sweep::SweepParameterConfig& config_template,
                                  size_t num_bits,
                                  const std::string& stage_tag,
                                  const std::string& run_id,
                                  ofec_sweep::detail::DualOut& out,
                                  std::ofstream& csv) {
  if (patterns.empty()) return {};
  ofec_sweep::SweepParameterConfig config = config_template;
  newcode::Params base_params = config.base_params;
  base_params.NUM_INFO_BITS = num_bits;
  newcode::PipelineConfig pipeline_cfg = ofec_sweep::detail::make_pipeline_config(config);

  out << "[INFO] Stage " << stage_tag << " evaluating " << patterns.size()
      << " patterns with NUM_INFO_BITS=" << num_bits << "\n";

  std::vector<Evaluation> evaluations;
  evaluations.reserve(patterns.size());

  for (const auto& pattern : patterns) {
    auto scenario = make_scenario(pattern, config);
    newcode::Params params = base_params;
    params.ALPHA_LIST = scenario.alpha_list;
    params.beta_list = scenario.beta_list;
    if (!params.ALPHA_LIST.empty()) params.ALPHA = params.ALPHA_LIST.front();
    if (!params.beta_list.empty()) params.beta = params.beta_list.front();
    params.CHASE_L = scenario.chase_L;
    params.CHASE_NTEST = scenario.chase_n_test;
    params.BITGEN_SEED = scenario.bitgen_seed;
    params.CHANNEL_SEED = scenario.channel_seed;

    auto result = newcode::run_pipeline(params, pipeline_cfg,
                                        scenario.name + "_" + stage_tag,
                                        scenario.ebn0_db);

    out << "[RESULT-" << stage_tag << "] " << scenario.name
        << " Post-BER=" << result.post_fec.ber
        << " (errs=" << result.post_fec.errors << "/" << result.post_fec.total << ")\n";

    write_csv_row(csv, run_id, stage_tag, num_bits, scenario, result);

    evaluations.push_back(Evaluation{pattern, result});
  }
  return evaluations;
}

std::vector<ofec_sweep::ExplicitAlphaBetaPattern> select_top_patterns(
    std::vector<Evaluation>& evals) {
  if (evals.empty()) return {};
  std::sort(evals.begin(), evals.end(),
            [](const Evaluation& a, const Evaluation& b) {
              return a.result.post_fec.ber < b.result.post_fec.ber;
            });
  const size_t keep = std::max<size_t>(1, static_cast<size_t>(std::ceil(evals.size() * kKeepRatio)));
  std::vector<ofec_sweep::ExplicitAlphaBetaPattern> top;
  top.reserve(keep);
  for (size_t i = 0; i < keep && i < evals.size(); ++i) {
    top.push_back(evals[i].pattern);
  }
  return top;
}

void print_final_summary(const Evaluation& best,
                         ofec_sweep::detail::DualOut& out) {
  out << "\n[SUMMARY] Best pattern: " << best.pattern.label
      << " | Post-BER=" << best.result.post_fec.ber
      << " (errs=" << best.result.post_fec.errors << "/"
      << best.result.post_fec.total << ")\n";
  out << "  Alphas: ";
  for (size_t i = 0; i < best.pattern.alpha_list.size(); ++i) {
    out << best.pattern.alpha_list[i]
        << (i + 1 < best.pattern.alpha_list.size() ? ", " : "\n");
  }
  out << "  Betas : ";
  for (size_t i = 0; i < best.pattern.beta_list.size(); ++i) {
    out << best.pattern.beta_list[i]
        << (i + 1 < best.pattern.beta_list.size() ? ", " : "\n");
  }
}

}  // namespace

int main() {
  const std::filesystem::path data_dir = "data";
  ofec_sweep::detail::ensure_dir(data_dir);
  const std::string run_id = ofec_sweep::detail::now_stamp();
  const std::string log_path = (data_dir / ("run_" + run_id + "_sweep2.log")).string();
  ofec_sweep::detail::DualOut log(std::cout, log_path);
  const std::string csv_path =
      (data_dir / ("ofec_sweep2_results_" + run_id + ".csv")).string();
  ensure_sweep2_csv_header(csv_path);
  std::ofstream csv(csv_path, std::ios::out | std::ios::app);
  csv.setf(std::ios::fixed);
  csv << std::setprecision(8);

  auto full_patterns = build_schedule();
  log << "[INFO] Stage 1 candidate count: " << full_patterns.size() << "\n";
  auto config = build_base_config();

  auto stage1_results = run_stage(full_patterns, config, kStage1Bits, "stage1",
                                  run_id, log, csv);
  if (stage1_results.empty()) {
    log << "[ERROR] Stage 1 produced no results\n";
    return 1;
  }

  auto top_patterns = select_top_patterns(stage1_results);
  log << "[INFO] Stage 2 candidate count: " << top_patterns.size() << "\n";
  auto stage2_results = run_stage(top_patterns, config, kStage2Bits, "stage2",
                                  run_id, log, csv);
  if (stage2_results.empty()) {
    log << "[ERROR] Stage 2 produced no results\n";
    return 1;
  }

  auto best_it = std::min_element(stage2_results.begin(), stage2_results.end(),
                                  [](const Evaluation& a, const Evaluation& b) {
                                    return a.result.post_fec.ber < b.result.post_fec.ber;
                                  });
  print_final_summary(*best_it, log);
  return 0;
}
