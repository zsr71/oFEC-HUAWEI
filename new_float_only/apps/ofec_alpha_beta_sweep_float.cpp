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

// 下面这组常量控制扫描网格生成方式：只需要指定起始值、终点值和总点数。
// 程序会自动用线性插值生成 alpha/beta/gamma 的候选取值。
constexpr float kAlphaLowStart = 0.1f;
constexpr float kAlphaLowEnd = 0.4f;
constexpr std::size_t kAlphaLowCount = 4;

constexpr float kAlphaHighStart = 0.6f;
constexpr float kAlphaHighEnd = 1.1f;
constexpr std::size_t kAlphaHighCount = 6;

constexpr float kBetaLowStart = 0.1f;
constexpr float kBetaLowEnd = 0.4f;
constexpr std::size_t kBetaLowCount = 4;

constexpr float kBetaHighStart = 0.6f;
constexpr float kBetaHighEnd = 1.1f;
constexpr std::size_t kBetaHighCount = 6;

constexpr float kGammaAlphaStart = 1.0f;
constexpr float kGammaAlphaEnd = 1.0f;
constexpr std::size_t kGammaAlphaCount = 1;

constexpr float kGammaBetaStart = 1.0f;
constexpr float kGammaBetaEnd = 1.0f;
constexpr std::size_t kGammaBetaCount = 1;

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

std::vector<float> build_grid(float start, float end, std::size_t count) {
  // 功能：根据起始值、终点值和总点数生成一个线性扫描网格。
  // 输入：start/end 定义区间两端，count 定义需要多少个候选点。
  // 输出：返回长度为 count 的网格；当 count<=1 时仅返回 {start}。
  std::vector<float> grid(count == 0 ? 1 : count, start);
  if (count <= 1) {
    return grid;
  }

  for (std::size_t i = 0; i < count; ++i) {
    const float t = static_cast<float>(i) / static_cast<float>(count - 1);
    grid[i] = start + (end - start) * t;
  }
  return grid;
}

std::vector<float> build_sequence(float low, float high, float gamma, std::size_t count) {
  // 功能：根据 low/high/gamma 生成一个按 tile 排列的参数序列。
  // 输入：low/high 为序列两端，gamma 为形状参数，count 为需要生成的元素个数。
  // 输出：返回长度为 count 的 alpha 或 beta 序列；当 gamma=1 时退化为线性插值。
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
  // 功能：把单个浮点参数格式化成固定两位小数的字符串。
  // 输入：value 为要显示的 alpha/beta/gamma 数值。
  // 输出：返回用于 pattern label 的短字符串。
  std::ostringstream oss;
  oss << std::fixed << std::setprecision(2) << value;
  return oss.str();
}

std::string format_sequence(const std::vector<float>& seq) {
  // 功能：把一个 tile 参数序列转成可写入 CSV 的文本。
  // 输入：seq 为某个 pattern 的 alpha_list 或 beta_list。
  // 输出：返回以 ';' 拼接的字符串，便于后续复核每个 tile 的取值。
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
  // 功能：把 low/high/gamma 网格展开成所有候选 shape。
  // 输入：读取文件顶部定义好的 alpha/beta/gamma 网格常量，不额外接收参数。
  // 输出：返回所有合法 shape 组合；过近的 low/high 组合会被过滤掉。
  const std::vector<float> alpha_low_grid =
      build_grid(kAlphaLowStart, kAlphaLowEnd, kAlphaLowCount);
  const std::vector<float> alpha_high_grid =
      build_grid(kAlphaHighStart, kAlphaHighEnd, kAlphaHighCount);
  const std::vector<float> beta_low_grid =
      build_grid(kBetaLowStart, kBetaLowEnd, kBetaLowCount);
  const std::vector<float> beta_high_grid =
      build_grid(kBetaHighStart, kBetaHighEnd, kBetaHighCount);
  const std::vector<float> gamma_alpha_grid =
      build_grid(kGammaAlphaStart, kGammaAlphaEnd, kGammaAlphaCount);
  const std::vector<float> gamma_beta_grid =
      build_grid(kGammaBetaStart, kGammaBetaEnd, kGammaBetaCount);

  std::vector<Shape> shapes;
  for (float alpha_low : alpha_low_grid) {
    for (float alpha_high : alpha_high_grid) {
      if (alpha_high - alpha_low < 0.05f) {
        continue;
      }
      for (float beta_low : beta_low_grid) {
        for (float beta_high : beta_high_grid) {
          if (beta_high - beta_low < 0.05f) {
            continue;
          }
          for (float gamma_alpha : gamma_alpha_grid) {
            for (float gamma_beta : gamma_beta_grid) {
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
  // 功能：把 shape 进一步具体化成可运行的 pattern 候选。
  // 输入：params 主要提供 TILES_PER_WIN，用于决定 alpha/beta 序列长度。
  // 输出：返回带 label、alpha_list、beta_list 的候选集合，供 stage1/stage2 直接使用。
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
  // 功能：把应用层的 PatternCandidate 转成公共 scheduler 使用的 SweepPattern。
  // 输入：candidates 为当前阶段准备评估的一组参数模式。
  // 输出：返回可直接送入 build_sweep_tasks() 的 pattern 列表。
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
  // 功能：写 summary CSV 的表头。
  // 输入：csv 为已经打开的输出文件流。
  // 输出：无返回值；副作用是向文件写入列名，覆盖 stage/rank/shape/BER 等字段。
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
  // 功能：把单个 pattern 在某个 stage 的结果写成一行 CSV。
  // 输入：run_id/stage_name/rank 标识本轮扫描位置；config 给出公共运行参数；
  // evaluation 提供 shape、alpha/beta 序列和聚合后的 pre/post BER。
  // 输出：无返回值；副作用是向 csv 追加一行可复盘记录。
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
  // 功能：组装某个阶段共用的 task runner 配置。
  // 输入：params 是已经规范化的 decoder 参数；stage_name 用于区分 stage1/stage2 标签。
  // 输出：返回传给 run_sweep_tasks() 的配置对象。
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
  // 功能：执行一次完整阶段评估，包括构造 task、并行运行、聚合排序和写 CSV。
  // 输入：stage_name 标识当前阶段；num_info_bits 控制该阶段用多少比特；
  // decoder_template 提供基础解码参数；candidates 是待评估 pattern；
  // seed_schedule 是所有 pattern 共享的 seed；csv/run_id 用于结果落盘。
  // 输出：返回按 post-BER 从优到劣排序后的 StageEvaluation 列表。
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
  // 功能：从 stage1 的排序结果里截取前 keep_count 个模式，送入 stage2 精扫。
  // 输入：evaluations 必须已经按 post-BER 排好序；keep_count 为保留数量上限。
  // 输出：返回进入下一阶段的 PatternCandidate 子集。
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
  // 功能：把某个阶段的 pattern 排名打印到控制台。
  // 输入：stage_name 用于标记 stage1/stage2；evaluations 为该阶段的排序结果。
  // 输出：无返回值；副作用是打印每个 pattern 的 pre/post BER。
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
  // 功能：构造默认扫描配置，执行 stage1 快扫和 stage2 精扫，并生成 summary CSV。
  // 输入：不接收命令行参数，所有实验参数都来自本文件顶部常量。
  // 输出：成功时返回 0，并在 data/ 下写出 CSV；失败时返回 2 并打印异常信息。
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
