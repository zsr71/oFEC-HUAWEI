#include <algorithm>
#include <cmath>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

#include "new_float_only/io/ensure_dir.hpp"
#include "new_float_only/sweep_task_runner.hpp"
#include "new_float_only/utils/now_stamp.hpp"

namespace {

constexpr const char* kLabel = "alpha_beta_sweep_float";
constexpr float kEbN0Db = 3.07f;
constexpr std::size_t kTrialCount = 1;
constexpr std::size_t kStage1Bits = 16 * 132 * 16 * 111;
constexpr std::size_t kStage2Bits = 32 * 132 * 16 * 111;
constexpr std::size_t kStage2KeepCount = 8;
constexpr int kBitgenSeedBase = 1521867291;
constexpr int kChannelSeedBase = 998258255;
constexpr unsigned kBitsPerSymbol = 1;
constexpr bool kGenerateRandomBits = true;
constexpr bool kNormalizeExtrinsic = true;
constexpr unsigned kMaxWorkersOverride = 0;
constexpr bool kQuietPipeline = true;
constexpr bool kQuietLogs = false;

constexpr bool kNormalizeKnownPrefixTail = false;
constexpr int kChaseL = 6;

// alpha 序列起点候选：决定某个 pattern 里第一个 tile 使用的 alpha 下界。
const std::vector<float> kAlphaLowGrid = {0.1f, 0.2f, 0.3f, 0.4f};
// alpha 序列终点候选：决定某个 pattern 里最后一个 tile 使用的 alpha 上界。
const std::vector<float> kAlphaHighGrid = {0.6f, 0.7f, 0.8f, 0.9f, 1.0f, 1.1f};
// beta 序列起点候选：决定某个 pattern 里第一个 tile 使用的 beta 下界。
const std::vector<float> kBetaLowGrid = {0.1f, 0.2f, 0.3f, 0.4f};
// beta 序列终点候选：决定某个 pattern 里最后一个 tile 使用的 beta 上界。
const std::vector<float> kBetaHighGrid = {0.6f, 0.7f, 0.8f, 0.9f, 1.0f, 1.1f};
// gamma_alpha 控制 alpha 序列从 low 到 high 的弯曲程度；1.0 表示线性。
const std::vector<float> kGammaAlphaGrid = {1.0f};
// gamma_beta 控制 beta 序列从 low 到 high 的弯曲程度；1.0 表示线性。
const std::vector<float> kGammaBetaGrid = {1.0f};

struct Shape {
  float alpha_low = 0.0f;
  float alpha_high = 0.0f;
  float gamma_alpha = 1.0f;
  float beta_low = 0.0f;
  float beta_high = 0.0f;
  float gamma_beta = 1.0f;
};

struct PatternCandidate {
  Shape shape;
  std::string pattern_label;
  std::vector<float> alpha_list;
  std::vector<float> beta_list;
};

struct StageEvaluation {
  PatternCandidate candidate;
  new_float_only::SweepPatternSummary summary;
};

std::vector<float> build_sequence(float low, float high, float gamma, std::size_t count) {
  std::vector<float> seq(count, low);
  if (count <= 1) {
    return seq;
  }

  const float span = std::max(0.0f, high - low);
  for (std::size_t i = 0; i < count; ++i) {
    const float t = static_cast<float>(i) / static_cast<float>(count - 1);
    const float shaped =
        (gamma == 1.0f) ? t : std::pow(std::clamp(t, 0.0f, 1.0f), gamma);
    seq[i] = low + span * shaped;
  }
  return seq;
}

std::string format_float(float value) {
  std::ostringstream oss;
  oss << std::fixed << std::setprecision(2) << value;
  return oss.str();
}

std::string format_sequence(const std::vector<float>& seq) {
  std::ostringstream oss;
  for (std::size_t i = 0; i < seq.size(); ++i) {
    if (i) {
      oss << ';';
    }
    oss << std::fixed << std::setprecision(3) << seq[i];
  }
  return oss.str();
}

std::vector<Shape> build_shapes() {
  std::vector<Shape> shapes;
  for (float alpha_low : kAlphaLowGrid) {
    for (float alpha_high : kAlphaHighGrid) {
      if (alpha_high - alpha_low < 0.05f) {
        continue;
      }
      for (float beta_low : kBetaLowGrid) {
        for (float beta_high : kBetaHighGrid) {
          if (beta_high - beta_low < 0.05f) {
            continue;
          }
          for (float gamma_alpha : kGammaAlphaGrid) {
            for (float gamma_beta : kGammaBetaGrid) {
              shapes.push_back(Shape{
                  .alpha_low = alpha_low,
                  .alpha_high = alpha_high,
                  .gamma_alpha = gamma_alpha,
                  .beta_low = beta_low,
                  .beta_high = beta_high,
                  .gamma_beta = gamma_beta,
              });
            }
          }
        }
      }
    }
  }
  return shapes;
}

std::vector<PatternCandidate> build_candidates(const new_float_only::Params& params) {
  const std::vector<Shape> shapes = build_shapes();
  std::vector<PatternCandidate> candidates;
  candidates.reserve(shapes.size());

  for (const Shape& shape : shapes) {
    PatternCandidate candidate;
    candidate.shape = shape;
    candidate.pattern_label =
        "a(" + format_float(shape.alpha_low) + "," + format_float(shape.alpha_high) + "," +
        format_float(shape.gamma_alpha) + ")" +
        "_b(" + format_float(shape.beta_low) + "," + format_float(shape.beta_high) + "," +
        format_float(shape.gamma_beta) + ")";
    candidate.alpha_list =
        build_sequence(shape.alpha_low, shape.alpha_high, shape.gamma_alpha, params.TILES_PER_WIN);
    candidate.beta_list =
        build_sequence(shape.beta_low, shape.beta_high, shape.gamma_beta, params.TILES_PER_WIN);
    candidates.push_back(candidate);
  }

  return candidates;
}

std::vector<new_float_only::SweepPattern> materialize_patterns(
    const std::vector<PatternCandidate>& candidates) {
  std::vector<new_float_only::SweepPattern> patterns;
  patterns.reserve(candidates.size());

  for (std::size_t i = 0; i < candidates.size(); ++i) {
    new_float_only::SweepPattern pattern;
    pattern.pattern_index = i;
    pattern.pattern_label = candidates[i].pattern_label;
    pattern.alpha_list = candidates[i].alpha_list;
    pattern.beta_list = candidates[i].beta_list;
    patterns.push_back(pattern);
  }
  return patterns;
}

void write_summary_csv_header(std::ofstream& csv) {
  csv << "run_id,label,stage,rank,pattern_index,pattern_label,ebn0_db,trial_count,max_workers,"
         "num_info_bits,"
         "alpha_low,alpha_high,gamma_alpha,beta_low,beta_high,gamma_beta,"
         "alpha_list,beta_list,"
         "pre_ber,pre_errors,pre_total,pre_frame_error_trials,"
         "post_ber,post_errors,post_total,post_frame_error_trials\n";
}

void write_summary_csv_row(std::ofstream& csv,
                           const std::string& run_id,
                           const std::string& stage_name,
                           std::size_t rank,
                           std::size_t num_info_bits,
                           unsigned max_workers,
                           const new_float_only::SweepTaskRunnerConfig& config,
                           const StageEvaluation& evaluation) {
  const Shape& shape = evaluation.candidate.shape;
  const auto& summary = evaluation.summary;

  csv << run_id << ","
      << '"' << config.label << '"' << ","
      << stage_name << ","
      << rank << ","
      << summary.pattern_index << ","
      << '"' << summary.pattern_label << '"' << ","
      << config.ebn0_db << ","
      << summary.trials_completed << ","
      << max_workers << ","
      << num_info_bits << ","
      << shape.alpha_low << ","
      << shape.alpha_high << ","
      << shape.gamma_alpha << ","
      << shape.beta_low << ","
      << shape.beta_high << ","
      << shape.gamma_beta << ","
      << '"' << format_sequence(summary.alpha_list) << '"' << ","
      << '"' << format_sequence(summary.beta_list) << '"' << ","
      << summary.pre_fec.ber << ","
      << summary.pre_fec.errors << ","
      << summary.pre_fec.total << ","
      << summary.pre_frame_error_trials << ","
      << summary.post_fec.ber << ","
      << summary.post_fec.errors << ","
      << summary.post_fec.total << ","
      << summary.post_frame_error_trials << "\n";
}

new_float_only::SweepTaskRunnerConfig build_stage_config(
    const new_float_only::Params& params,
    const std::string& stage_name) {
  new_float_only::SweepTaskRunnerConfig config;
  config.label = std::string(kLabel) + "_" + stage_name;
  config.ebn0_db = kEbN0Db;
  config.bits_per_symbol = kBitsPerSymbol;
  config.generate_random_bits = kGenerateRandomBits;
  config.normalize_extrinsic = kNormalizeExtrinsic;
  config.base_decoder = params;
  config.max_workers_override = kMaxWorkersOverride;
  config.quiet_pipeline = kQuietPipeline;
  config.quiet_logs = kQuietLogs;
  return config;
}

std::vector<StageEvaluation> run_stage(const std::string& stage_name,
                                       std::size_t num_info_bits,
                                       const new_float_only::DecoderConfig& decoder_template,
                                       const std::vector<PatternCandidate>& candidates,
                                       const std::vector<new_float_only::SweepSeedPair>& seed_schedule,
                                       std::ofstream& csv,
                                       const std::string& run_id) {
  if (candidates.empty()) {
    return {};
  }

  new_float_only::DecoderConfig decoder = decoder_template;
  decoder.NUM_INFO_BITS = num_info_bits;
  const new_float_only::Params params =
      new_float_only::normalize_sweep_decoder_config(decoder);
  const std::vector<new_float_only::SweepPattern> patterns =
      materialize_patterns(candidates);
  const std::vector<new_float_only::SweepTask> tasks =
      new_float_only::build_sweep_tasks(patterns, seed_schedule);
  const new_float_only::SweepTaskRunnerConfig config =
      build_stage_config(params, stage_name);
  const unsigned max_workers = new_float_only::resolve_sweep_worker_count(config);

  std::cout << "[APP] " << stage_name
            << ": patterns=" << patterns.size()
            << ", shared_seeds=" << seed_schedule.size()
            << ", tasks=" << tasks.size()
            << ", num_info_bits=" << num_info_bits
            << ", workers=" << max_workers << "\n";

  const std::vector<new_float_only::SweepTaskResult> task_results =
      new_float_only::run_sweep_tasks(config, tasks);
  std::vector<new_float_only::SweepPatternSummary> summaries =
      new_float_only::aggregate_sweep_task_results(patterns, task_results);

  std::sort(summaries.begin(),
            summaries.end(),
            [](const new_float_only::SweepPatternSummary& lhs,
               const new_float_only::SweepPatternSummary& rhs) {
              if (lhs.post_fec.ber != rhs.post_fec.ber) {
                return lhs.post_fec.ber < rhs.post_fec.ber;
              }
              return lhs.pattern_index < rhs.pattern_index;
            });

  std::vector<StageEvaluation> evaluations;
  evaluations.reserve(summaries.size());
  for (const auto& summary : summaries) {
    StageEvaluation evaluation;
    evaluation.candidate = candidates.at(summary.pattern_index);
    evaluation.summary = summary;
    evaluations.push_back(evaluation);
  }

  for (std::size_t rank = 0; rank < evaluations.size(); ++rank) {
    write_summary_csv_row(csv,
                          run_id,
                          stage_name,
                          rank,
                          num_info_bits,
                          max_workers,
                          config,
                          evaluations[rank]);
  }

  return evaluations;
}

std::vector<PatternCandidate> select_top_candidates(
    const std::vector<StageEvaluation>& evaluations,
    std::size_t keep_count) {
  if (evaluations.empty()) {
    return {};
  }

  const std::size_t keep = std::min(keep_count, evaluations.size());
  std::vector<PatternCandidate> top_candidates;
  top_candidates.reserve(keep);
  for (std::size_t i = 0; i < keep; ++i) {
    top_candidates.push_back(evaluations[i].candidate);
  }
  return top_candidates;
}

void print_stage_results(const std::string& stage_name,
                         const std::vector<StageEvaluation>& evaluations) {
  for (const auto& evaluation : evaluations) {
    const auto& summary = evaluation.summary;
    std::cout << "[APP] " << stage_name
              << " pattern#" << summary.pattern_index
              << " " << summary.pattern_label
              << " post-BER=" << summary.post_fec.ber
              << " (" << summary.post_fec.errors << "/" << summary.post_fec.total << ")"
              << ", pre-BER=" << summary.pre_fec.ber
              << " (" << summary.pre_fec.errors << "/" << summary.pre_fec.total << ")\n";
  }
}

}  // namespace

int main() {
  new_float_only::DecoderConfig decoder;
  decoder.NUM_INFO_BITS = kStage2Bits;
  decoder.NORMALIZE_KNOWN_PREFIX_TAIL = kNormalizeKnownPrefixTail;
  decoder.CHASE_L = kChaseL;
  decoder.CHASE_NTEST = 1 << decoder.CHASE_L;
  decoder.HARD_DECODE_DEFAULT = false;
  decoder.HARD_TILE_LIST = {0, 0, 0, 0};
  decoder.DUMP_WORK_LLR = false;
  decoder.debug_trace.enable = false;

  try {
    const new_float_only::Params params =
        new_float_only::normalize_sweep_decoder_config(decoder);
    const std::vector<PatternCandidate> candidates = build_candidates(params);
    const std::vector<new_float_only::SweepSeedPair> seed_schedule =
        new_float_only::build_sweep_seed_schedule(kTrialCount,
                                                  kBitgenSeedBase,
                                                  kChannelSeedBase);

    const std::filesystem::path data_dir = "data";
    io::ensure_dir(data_dir);
    const std::string run_id = utils::now_stamp();
    const std::string csv_path =
        (data_dir / ("ofec_alpha_beta_sweep_float_" + run_id + ".csv")).string();
    std::ofstream csv(csv_path, std::ios::out | std::ios::trunc);
    if (!csv.is_open()) {
      throw std::runtime_error("Failed to open summary CSV: " + csv_path);
    }
    write_summary_csv_header(csv);

    const std::vector<StageEvaluation> stage1_evaluations =
        run_stage("stage1", kStage1Bits, decoder, candidates, seed_schedule, csv, run_id);
    if (stage1_evaluations.empty()) {
      throw std::runtime_error("Stage1 produced no evaluations");
    }

    const std::vector<PatternCandidate> stage2_candidates =
        select_top_candidates(stage1_evaluations, kStage2KeepCount);
    std::cout << "[APP] stage1 keep top " << stage2_candidates.size()
              << " / " << stage1_evaluations.size() << " patterns for stage2\n";

    const std::vector<StageEvaluation> stage2_evaluations =
        run_stage("stage2", kStage2Bits, decoder, stage2_candidates, seed_schedule, csv, run_id);
    if (stage2_evaluations.empty()) {
      throw std::runtime_error("Stage2 produced no evaluations");
    }

    print_stage_results("stage2", stage2_evaluations);
    std::cout << "[APP] summary CSV: " << csv_path << "\n";
    return 0;
  } catch (const std::exception& ex) {
    std::cerr << "[APP] failed: " << ex.what() << "\n";
    return 2;
  }
}
