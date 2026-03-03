#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <numeric>
#include <sstream>
#include <string>
#include <vector>

#include "newcode/io/ensure_dir.hpp"
#include "newcode/utils/linspace.hpp"
#include "newcode/utils/now_stamp.hpp"
#include "newcode/ofec_sweep_runner.hpp"
#include "ofec_sweep_detail.hpp"

namespace {

constexpr size_t kTilesPerWindow = 2;
constexpr float  kEvalEbN0       = 3.07f;
constexpr size_t kStage1Bits     = 4 * 110 * 16 * 111;
constexpr size_t kStage2Bits     = 16 * 110 * 16 * 111;
constexpr float  kKeepRatio      = 0.20f;

static constexpr const char* kInterleaverName          = "identity";
static constexpr const char* kDecoderName             = "plain";
static constexpr unsigned    kBitsPerSymbol           = 1;
static constexpr bool        kNormalizeExtrinsic      = true;
static constexpr bool        kGenerateRandomBits      = true;
static constexpr bool        kNormalizeKnownPrefixTail = true;
static constexpr std::size_t kLlrBits                 = 16;
static constexpr float       kQuantClipRatio          = 0.0f;
const std::vector<int> kSisoActiveList                = {32, 30};

const std::vector<float> kAlphaLowGrid   = utils::linspace(0.00f, 1.50f, 3);
const std::vector<float> kAlphaHighGrid  = utils::linspace(0.00f, 1.50f, 3);
const std::vector<float> kBetaLowGrid    = utils::linspace(0.00f, 1.50f, 3);
const std::vector<float> kBetaHighGrid   = utils::linspace(0.00f, 1.50f, 3);
const std::vector<float> kGammaAlphaGrid = utils::linspace(0.70f, 1.30f, 1);
const std::vector<float> kGammaBetaGrid  = utils::linspace(0.70f, 1.30f, 1);

struct Shape {
  float alpha_low;
  float alpha_high;
  float gamma_alpha;
  float beta_low;
  float beta_high;
  float gamma_beta;
  std::string tag = "_r0";
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
      << ")_b(" << shape.beta_low << "," << shape.beta_high << "," << shape.gamma_beta << ")";

  ofec_sweep::ExplicitAlphaBetaPattern pattern;
  pattern.label = oss.str();
  pattern.alpha_list = alpha_seq;
  pattern.beta_list = beta_seq;
  return pattern;
}

std::vector<Shape> build_shapes() {
  std::vector<Shape> shapes;
  for (float aL : kAlphaLowGrid) {
    for (float aH : kAlphaHighGrid) {
      if (aH - aL < 0.05f) continue;
      for (float bL : kBetaLowGrid) {
        for (float bH : kBetaHighGrid) {
          if (bH - bL < 0.05f) continue;
          for (float gA : kGammaAlphaGrid) {
            for (float gB : kGammaBetaGrid) {
              Shape s;
              s.alpha_low = aL;
              s.alpha_high = aH;
              s.gamma_alpha = gA;
              s.beta_low = bL;
              s.beta_high = bH;
              s.gamma_beta = gB;
              shapes.push_back(s);
            }
          }
        }
      }
    }
  }
  return shapes;
}

std::vector<ofec_sweep::ExplicitAlphaBetaPattern> build_schedule(const std::vector<Shape>& shapes) {
  std::vector<ofec_sweep::ExplicitAlphaBetaPattern> patterns;
  patterns.reserve(shapes.size());
  for (size_t i = 0; i < shapes.size(); ++i) {
    patterns.push_back(shape_to_pattern(shapes[i], "r0", i));
  }
  return patterns;
}

ofec_sweep::SweepParameterConfig build_base_config() {
  ofec_sweep::SweepParameterConfig config;
  config.base_params.TILES_PER_WIN = kTilesPerWindow;
  config.base_params.ALPHA_LIST.assign(kTilesPerWindow, 0.3f);
  config.base_params.beta_list.assign(kTilesPerWindow, 0.6f);
  config.base_params.HARD_TILE_LIST.assign(kTilesPerWindow, 0);
  config.base_params.SISO_ACTIVE_LIST = kSisoActiveList;
  config.base_params.BITGEN_RANDOM_BITS = kGenerateRandomBits;
  config.base_params.NORMALIZE_KNOWN_PREFIX_TAIL = kNormalizeKnownPrefixTail;
  config.base_params.LLR_CLIP_RATIO = kQuantClipRatio;
  config.base_params.LLR_BITS = kLlrBits;

  config.interleaver_name = kInterleaverName;
  config.decoder_name = kDecoderName;
  config.bits_per_symbol = kBitsPerSymbol;
  config.normalize_extrinsic = kNormalizeExtrinsic;
  config.generate_random_bits = kGenerateRandomBits;
  config.normalize_known_prefix_tail = kNormalizeKnownPrefixTail;
  config.quant_clip_ratio = kQuantClipRatio;
  config.siso_active_list = kSisoActiveList;

  config.chase_l_candidates = {config.base_params.CHASE_L};
  config.bitgen_seed_count = 1;
  config.channel_seed_count = 1;
  config.ebn0_start = kEvalEbN0;
  config.ebn0_end = kEvalEbN0;
  config.ebn0_points = 1;

  config.base_params.debug_trace = {};
  return config;
}

struct Evaluation {
  ofec_sweep::ExplicitAlphaBetaPattern pattern;
  newcode::PipelineResult result;
  Shape shape;
};

ofec_sweep::detail::SweepScenario make_scenario(
    const ofec_sweep::ExplicitAlphaBetaPattern& pattern,
    const Shape& shape,
    const ofec_sweep::SweepParameterConfig& config) {
  ofec_sweep::detail::SweepScenario sc;
  sc.name = pattern.label;
  sc.alpha_list = pattern.alpha_list;
  sc.beta_list = pattern.beta_list;
  if (!sc.alpha_list.empty()) sc.alpha_start = sc.alpha_list.front();
  if (!sc.beta_list.empty()) sc.beta_start = sc.beta_list.front();
  sc.alpha_low = shape.alpha_low;
  sc.alpha_high = shape.alpha_high;
  sc.gamma_alpha = shape.gamma_alpha;
  sc.beta_low = shape.beta_low;
  sc.beta_high = shape.beta_high;
  sc.gamma_beta = shape.gamma_beta;
  sc.chase_L = config.base_params.CHASE_L;
  sc.chase_n_test = 1 << sc.chase_L;
  sc.bitgen_seed = config.base_params.BITGEN_SEED;
  sc.channel_seed = config.base_params.CHANNEL_SEED;
  sc.ebn0_db = kEvalEbN0;
  return sc;
}

std::vector<Evaluation> run_stage(const std::vector<ofec_sweep::ExplicitAlphaBetaPattern>& patterns,
                                  const std::vector<Shape>& shapes,
                                  const ofec_sweep::SweepParameterConfig& config_template,
                                  size_t num_bits,
                                  const std::string& stage_tag,
                                  const std::string& run_id,
                                  ofec_sweep::detail::DualOut& log,
                                  std::ofstream& csv) {
  if (patterns.empty()) return {};
  ofec_sweep::SweepParameterConfig config = config_template;
  newcode::Params base_params = config.base_params;
  base_params.NUM_INFO_BITS = num_bits;
  config.base_params = base_params;

  log << "[INFO] Stage " << stage_tag << " evaluating " << patterns.size()
      << " patterns with NUM_INFO_BITS=" << num_bits << "\n";

  std::vector<ofec_sweep::detail::SweepScenario> scenarios;
  scenarios.reserve(patterns.size());
  for (size_t idx = 0; idx < patterns.size(); ++idx) {
    scenarios.push_back(make_scenario(patterns[idx], shapes[idx], config));
  }

  auto scenario_outputs = ofec_sweep::detail::run_scenarios_parallel(
      scenarios, config, /*max_workers_hint=*/0, stage_tag, &log);
  if (scenario_outputs.empty()) {
    return {};
  }

  std::vector<Evaluation> evaluations;
  evaluations.reserve(scenario_outputs.size());

  for (const auto& output : scenario_outputs) {
    const size_t idx = output.idx;
    const auto& scenario = scenarios[idx];
    const auto& result = output.result;

    log << "[RESULT-" << stage_tag << "] " << scenario.name
        << " Post-BER=" << result.post_fec.ber
        << " (errs=" << result.post_fec.errors << "/" << result.post_fec.total << ")\n";

    ofec_sweep::detail::write_csv_row(csv,
                                      utils::now_stamp(),
                                      run_id,
                                      stage_tag,
                                      num_bits,
                                      scenario,
                                      result,
                                      ofec_sweep::detail::CsvFormat::Extended);

    evaluations.push_back(Evaluation{patterns[idx], result, shapes[idx]});
  }
  return evaluations;
}

std::vector<ofec_sweep::ExplicitAlphaBetaPattern> select_top_patterns(
    std::vector<Evaluation>& evals,
    std::vector<Shape>& shapes) {
  if (evals.empty()) return {};
  std::vector<size_t> indices(evals.size());
  std::iota(indices.begin(), indices.end(), 0);
  std::sort(indices.begin(), indices.end(),
            [&](size_t a, size_t b) {
              return evals[a].result.post_fec.ber < evals[b].result.post_fec.ber;
            });
  const size_t keep = std::max<size_t>(1, static_cast<size_t>(std::ceil(evals.size() * kKeepRatio)));
  std::vector<ofec_sweep::ExplicitAlphaBetaPattern> top_patterns;
  std::vector<Shape> top_shapes;
  top_patterns.reserve(keep);
  top_shapes.reserve(keep);
  for (size_t i = 0; i < keep && i < indices.size(); ++i) {
    top_patterns.push_back(evals[indices[i]].pattern);
    top_shapes.push_back(evals[indices[i]].shape);
  }
  shapes.swap(top_shapes);
  return top_patterns;
}

void print_final_summary(const Evaluation& best,
                         ofec_sweep::detail::DualOut& log) {
  log << "\n[SUMMARY] Best pattern: " << best.pattern.label
      << " | Post-BER=" << best.result.post_fec.ber
      << " (errs=" << best.result.post_fec.errors << "/"
      << best.result.post_fec.total << ")\n";
  log << "  Alphas: ";
  for (size_t i = 0; i < best.pattern.alpha_list.size(); ++i) {
    log << best.pattern.alpha_list[i]
        << (i + 1 < best.pattern.alpha_list.size() ? ", " : "\n");
  }
  log << "  Betas : ";
  for (size_t i = 0; i < best.pattern.beta_list.size(); ++i) {
    log << best.pattern.beta_list[i]
        << (i + 1 < best.pattern.beta_list.size() ? ", " : "\n");
  }
}

}  // namespace

int main() {
  const std::filesystem::path data_dir = "data";
  io::ensure_dir(data_dir);
  const std::string run_id = utils::now_stamp();
  const std::string log_path = (data_dir / ("run_" + run_id + "_sweep2.log")).string();
  ofec_sweep::detail::DualOut log(std::cout, log_path);
  const std::string csv_path =
      (data_dir / ("ofec_sweep2_results_" + run_id + ".csv")).string();
  ofec_sweep::detail::ensure_csv_header_v2(csv_path);
  std::ofstream csv(csv_path, std::ios::out | std::ios::app);
  csv.setf(std::ios::fixed);
  csv << std::setprecision(8);

  auto shapes = build_shapes();
  auto patterns = build_schedule(shapes);
  log << "[INFO] Stage 1 candidate count: " << patterns.size() << "\n";
  auto config = build_base_config();
  log << "[INFO] Stage 1 workers: "
      << ofec_sweep::detail::resolve_worker_count(config) << "\n";

  auto stage1_results = run_stage(patterns, shapes, config, kStage1Bits,
                                  "stage1", run_id, log, csv);
  if (stage1_results.empty()) {
    log << "[ERROR] Stage 1 produced no results\n";
    return 1;
  }

  auto top_patterns = select_top_patterns(stage1_results, shapes);
  log << "[INFO] Stage 2 candidate count: " << top_patterns.size() << "\n";
  auto stage2_results = run_stage(top_patterns, shapes, config, kStage2Bits,
                                  "stage2", run_id, log, csv);
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
