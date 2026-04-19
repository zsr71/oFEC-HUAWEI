#include <algorithm>
#include <cmath>
#include <filesystem>
#include <fstream>
#include <future>
#include <iomanip>
#include <iostream>
#include <limits>
#include <optional>
#include <sstream>
#include <string>
#include <thread>
#include <utility>
#include <vector>

#include "mux_bypass_edges.hpp"
#include "newcode/io/dualwriter.hpp"
#include "newcode/ofec_single_runner.hpp"

namespace {

// ======== 用户可调参数区域 ========

// 输出参数
static constexpr const char* kLabelPrefix = "ber_window_probe"; // 运行标签前缀，会进入日志名和内部任务名
static constexpr const char* kOutputDir = "data/ber_window_probe"; // 本 app 的输出目录
static constexpr bool kQuietPipeline = true; // true=压低单次 pipeline 的控制台输出，避免 probe 时日志过多
static constexpr unsigned kMaxParallelEbN0 = 0;  // Eb/N0 并行度上限；0=自动使用硬件并发
static constexpr unsigned kMaxParallelSeeds = 0; // 每个 Eb/N0 内部的 seed 并行度上限；0=自动使用硬件并发

// 测试点与 seed
static const std::vector<float> kEbN0List = {3.052051f}; // 需要测试的 Eb/N0 列表
static constexpr int kBitgenSeedBase = 20260319;         // 基础比特种子；实际每个 seed 在此基础上递增
static constexpr int kChannelSeedBase = 3182026;         // 基础信道种子；实际每个 seed 在此基础上递增
static constexpr int kSeedCount = 10;                    // 每个 Eb/N0 点重复运行的 seed 数量

// 发射端参数
static constexpr int         kChaseL_override    = 6;    // Chase L，-1 表示沿用 Params 默认值
static constexpr bool        kGenerateRandomBits = true; // true=发送随机信息比特，false=发送全 0 比特

// 信道参数
static constexpr unsigned    kBitsPerSymbol      = 1;    // 每个调制符号携带的比特数：1=BPSK，偶数=QAM

// 早停参数
static constexpr bool        kEnableEarlyStop              = true;   // 早停总开关
static const std::vector<int> kEarlyStopEnableList         = {0,0,0,0,1,1}; // 按 tile 覆盖早停总开关：0=关，非 0=开
static constexpr int         kEarlyStopConditionMode       = 1;      // 早停条件模式：1=v1，2=v2
static const std::vector<int> kEarlyStopConditionModeList  = {};     // 按 tile 覆盖条件模式；空表示沿用全局值
static constexpr int         kEarlyStopActionMode          = 4;      // 早停动作模式：1~6
static const std::vector<int> kEarlyStopActionModeList     = {};     // 按 tile 覆盖动作模式；空表示沿用全局值
static constexpr int         kEarlyStopBindGroupSize       = 1;      // 条件1的绑定组大小：1=逐 row，4=四个绑定
static const std::vector<int> kEarlyStopBindGroupSizeList  = {};     // 按 tile 覆盖绑定组大小；空表示沿用全局值
static constexpr bool        kEarlyStopCondV1RequireBch    = true;   // 条件1是否要求 BCH syndrome 为 0
static constexpr bool        kEarlyStopCondV1RequireOverall = true;  // 条件1是否要求 overall parity 通过
static constexpr float       kEarlyStopV2LlrAbsThreshold   = 26.0f;  // 条件2里把 bit 视为不可靠的 |LLR| 阈值
static constexpr int         kEarlyStopV2MaxUnreliableBits = 25;     // 条件2允许的不可靠 bit 数上限
static constexpr bool        kEarlyStopCondV2IncludeOverall = true;  // 条件2统计不可靠 bit 时是否把 overall parity bit 算进去

// LLR 量化相关参数
static constexpr std::size_t kLlrBits        = 6;      // LLR 位宽：16=浮点，2~15=qfloat
static constexpr float       kQuantClipRatio = 0.5f;   // 动态 clip 比例，0 表示禁用自适应 clip

// 解码主参数
static constexpr const char* kInterleaverName          = "identity"; // 交织器名称，identity 表示不交织
static constexpr const char* kDecoderName              = "chase_baseline"; // decoder 名称
static constexpr int         kChaseNTestOverride       = -1;         // Chase NTEST，-1 表示默认按 2^L
static constexpr int         kChaseTopkKeep            = 24;         // top-k/pruned decoder 保留候选数
static constexpr int         kChaseGroupMinimaBits     = 3;          // group-minima decoder 的分组 bit 数
static constexpr bool        kNormalizeExtrinsic       = false;      // 是否对 decoder 输出 extrinsic 做归一化
static constexpr bool        kNormalizeKnownPrefixTail = false;      // 是否对 known-prefix 后的尾部 LLR 做归一化
static constexpr float kEarlyStopActionResidualDivisor = 1.0f;       // 动作2里 residual 的除数
static constexpr float kEarlyStopActionHardLlrMag      = 1.0f;       // 动作3里输出的固定 |LLR| 幅度

static const std::vector<float> kAlphaExplicit = {
    0.428571f, 0.447738f, 0.482782f, 0.528162f, 0.581902f, 0.642857f // 每个 tile 的 alpha 显式列表
};
static const std::vector<float> kBetaExplicit = {
    2.857143f, 6.179301f, 12.253626f, 20.119585f, 29.434408f, 40.000000f // 每个 tile 的普通 beta 显式列表
};
static const std::vector<float> kEarlyStopActionBetaExplicit = {
    2.857143f, 6.179301f, 12.253626f, 20.119585f, 29.434408f, 40.000000f // 每个 tile 的 early-stop 专用 beta 显式列表
};
static const std::vector<int> kSisoActiveList = {32, 32, 32, 32, 32, 32}; // 每个 tile 允许参与 SISO 的行数预算
static constexpr int  kMuxGroupG          = 1;                             // MUX 分组粒度，1 表示全局池化
static constexpr int  kMuxSchedulingMode  = 0;                             // MUX 调度模式：0=legacy，1=按 early-stop 细节排序
static constexpr int  kMuxPriorityRule    = 0;                             // 新 MUX 的优先级规则：0=更差优先，1=更接近通过优先
static constexpr bool kMuxEnableReconfig  = false;                         // true 表示启用重配置版 MUX 调度
static constexpr int  kMuxBypassScheme    = 1;                             // 旁路边集合方案编号：1=scheme1，2=scheme2

// ==================================

struct AggregatedWindowStats {
  std::size_t window_idx = 0;
  std::size_t agg_pre_errs = 0;
  std::size_t agg_pre_bits = 0;
  std::size_t agg_post_errs = 0;
  std::size_t agg_post_bits = 0;
};

struct CandidateStartStats {
  std::size_t start_window_k = 0;
  std::size_t post_errs_from_k = 0;
  std::size_t post_bits_from_k = 0;
  double post_ber_from_k = 0.0;
};

struct SeedTask {
  int seed_index = 0;
  int bitgen_seed = 0;
  int channel_seed = 0;
  float ebn0_db = 0.0f;
  std::string label;
  newcode::Params params;
};

struct SeedRunResult {
  int seed_index = 0;
  int bitgen_seed = 0;
  int channel_seed = 0;
  newcode::PipelineResult pipeline_result;
};

struct RawWindowRow {
  float ebn0_db = 0.0f;
  int seed_index = 0;
  int bitgen_seed = 0;
  int channel_seed = 0;
  std::size_t window_idx = 0;
  std::size_t pre_errs = 0;
  std::size_t pre_bits = 0;
  double pre_ber = 0.0;
  std::size_t post_errs = 0;
  std::size_t post_bits = 0;
  double post_ber = 0.0;
};

struct EbN0RunResult {
  float ebn0_db = 0.0f;
  std::vector<AggregatedWindowStats> aggregated;
  std::vector<CandidateStartStats> candidates;
  std::vector<RawWindowRow> raw_rows;
};

struct EbN0TaskGroup {
  float ebn0_db = 0.0f;
  std::vector<SeedTask> seed_tasks;
};

std::string format_float(float value) {
  std::ostringstream oss;
  oss << std::fixed << std::setprecision(6) << value;
  return oss.str();
}

void ensure_aggregated_size(std::vector<AggregatedWindowStats>& stats,
                            std::size_t size) {
  while (stats.size() < size) {
    AggregatedWindowStats entry;
    entry.window_idx = stats.size();
    stats.push_back(entry);
  }
}

unsigned resolved_parallel_ebn0() {
  if (kMaxParallelEbN0 > 0) {
    return kMaxParallelEbN0;
  }
  const unsigned hc = std::thread::hardware_concurrency();
  if (hc == 0) {
    return 1u;
  }
  return std::min<unsigned>(hc, static_cast<unsigned>(std::max<std::size_t>(1, kEbN0List.size())));
}

unsigned resolved_parallel_seeds_for_ebn0(unsigned ebn0_parallelism) {
  if (kMaxParallelSeeds > 0) {
    return kMaxParallelSeeds;
  }
  const unsigned hc = std::thread::hardware_concurrency();
  if (hc == 0) {
    return 1u;
  }
  return std::max(1u, hc / std::max(1u, ebn0_parallelism));
}

void build_candidates(const std::vector<AggregatedWindowStats>& aggregated,
                      std::vector<CandidateStartStats>& candidates_out) {
  candidates_out.clear();
  if (aggregated.empty()) {
    return;
  }

  const std::size_t num_windows = aggregated.size();
  candidates_out.reserve(num_windows);
  for (std::size_t k = 0; k < num_windows; ++k) {
    CandidateStartStats entry;
    entry.start_window_k = k;
    for (std::size_t w = k; w < num_windows; ++w) {
      entry.post_errs_from_k += aggregated[w].agg_post_errs;
      entry.post_bits_from_k += aggregated[w].agg_post_bits;
    }
    entry.post_ber_from_k =
        (entry.post_bits_from_k == 0)
            ? 0.0
            : static_cast<double>(entry.post_errs_from_k) /
                  static_cast<double>(entry.post_bits_from_k);
    candidates_out.push_back(entry);
  }
}

void write_raw_header(std::ofstream& out) {
  out << "ebn0_db,seed_index,bitgen_seed,channel_seed,window_idx,"
         "pre_errs,pre_bits,pre_ber,post_errs,post_bits,post_ber\n";
  out.flush();
}

void write_aggregated_header(std::ofstream& out) {
  out << "ebn0_db,window_idx,"
         "agg_pre_errs,agg_pre_bits,agg_pre_ber,"
         "agg_post_errs,agg_post_bits,agg_post_ber\n";
  out.flush();
}

void write_candidates_header(std::ofstream& out) {
  out << "ebn0_db,start_window_k,post_errs_from_k,post_bits_from_k,post_ber_from_k\n";
  out.flush();
}

ofec_single::Config make_base_config() {
  const auto& selected_mux_bypass_edges =
      app_mux::bypass_edges_for_scheme(kMuxBypassScheme);

  ofec_single::Config config{
      .label = kLabelPrefix,
      .ebn0_db = kEbN0List.front(),
      .chaseL_override = kChaseL_override,
      .chase_n_test_override = kChaseNTestOverride,
      .chase_topk_keep = kChaseTopkKeep,
      .chase_group_minima_bits = kChaseGroupMinimaBits,
      .normalize_extrinsic = kNormalizeExtrinsic,
      .bits_per_symbol = kBitsPerSymbol,
      .bitgen_seed = kBitgenSeedBase,
      .channel_seed = kChannelSeedBase,
      .enable_early_stop = kEnableEarlyStop,
      .early_stop_enable_list = kEarlyStopEnableList,
      .early_stop_condition_mode = kEarlyStopConditionMode,
      .early_stop_condition_mode_list = kEarlyStopConditionModeList,
      .early_stop_action_mode = kEarlyStopActionMode,
      .early_stop_action_mode_list = kEarlyStopActionModeList,
      .early_stop_bind_group_size = kEarlyStopBindGroupSize,
      .early_stop_bind_group_size_list = kEarlyStopBindGroupSizeList,
      .early_stop_cond_v1_require_bch = kEarlyStopCondV1RequireBch,
      .early_stop_cond_v1_require_overall = kEarlyStopCondV1RequireOverall,
      .early_stop_v2_llr_abs_threshold = kEarlyStopV2LlrAbsThreshold,
      .early_stop_v2_max_unreliable_bits = kEarlyStopV2MaxUnreliableBits,
      .early_stop_cond_v2_include_overall = kEarlyStopCondV2IncludeOverall,
      .alpha_explicit = kAlphaExplicit,
      .beta_explicit = kBetaExplicit,
      .early_stop_action_sign_beta_explicit = kEarlyStopActionBetaExplicit,
      .early_stop_action_residual_divisor = kEarlyStopActionResidualDivisor,
      .early_stop_action_hard_llr_mag = kEarlyStopActionHardLlrMag,
      .siso_active_list = kSisoActiveList,
      .mux_group_g = kMuxGroupG,
      .mux_scheduling_mode = kMuxSchedulingMode,
      .mux_early_stop_priority_rule = kMuxPriorityRule,
      .mux_enable_reconfig = kMuxEnableReconfig,
      .mux_extra_bypass_edges = selected_mux_bypass_edges,
      .interleaver_name = kInterleaverName,
      .decoder_name = kDecoderName,
      .generate_random_bits = kGenerateRandomBits,
      .normalize_known_prefix_tail = kNormalizeKnownPrefixTail,
      .quant_clip_ratio = kQuantClipRatio,
      .llr_bits = kLlrBits,
      .dump_quantized_llr = false,
      .quantized_llr_output_path = {},
      .dump_work_llr = false,
      .work_llr_output_path = {},
      .dump_tile_early_stop_samples = false,
      .tile_early_stop_samples_output_path = {},
      .debug_trace = {},
  };
  config.debug_trace.enable = false;
  return config;
}

std::vector<SeedRunResult> run_seed_tasks(std::vector<SeedTask> tasks,
                                          const newcode::PipelineConfig& pipeline_cfg,
                                          unsigned max_parallel_seeds) {
  std::vector<SeedRunResult> seed_results;
  seed_results.reserve(tasks.size());
  std::vector<std::future<SeedRunResult>> pending;

  auto launch_task = [&](SeedTask task) {
    return std::async(std::launch::async, [task = std::move(task), pipeline_cfg]() mutable {
      SeedRunResult out;
      out.seed_index = task.seed_index;
      out.bitgen_seed = task.bitgen_seed;
      out.channel_seed = task.channel_seed;
      out.pipeline_result =
          newcode::run_pipeline(task.params, pipeline_cfg, task.label, task.ebn0_db);
      return out;
    });
  };

  for (auto& task : tasks) {
    pending.push_back(launch_task(std::move(task)));
    if (pending.size() >= max_parallel_seeds) {
      seed_results.push_back(pending.front().get());
      pending.erase(pending.begin());
    }
  }
  for (auto& future : pending) {
    seed_results.push_back(future.get());
  }
  std::sort(seed_results.begin(), seed_results.end(),
            [](const SeedRunResult& lhs, const SeedRunResult& rhs) {
              return lhs.seed_index < rhs.seed_index;
            });
  return seed_results;
}

EbN0RunResult run_ebn0_probe(EbN0TaskGroup task_group,
                             const newcode::PipelineConfig& pipeline_cfg,
                             unsigned seed_parallelism) {
  std::vector<SeedRunResult> seed_results =
      run_seed_tasks(std::move(task_group.seed_tasks), pipeline_cfg,
                     std::max(1u, seed_parallelism));

  EbN0RunResult out;
  out.ebn0_db = task_group.ebn0_db;

  for (const auto& seed_result : seed_results) {
    const newcode::PipelineResult& result = seed_result.pipeline_result;
    const std::size_t max_windows =
        std::max(result.pre_fec_windows.size(), result.post_fec_windows.size());
    ensure_aggregated_size(out.aggregated, max_windows);

    for (std::size_t w = 0; w < max_windows; ++w) {
      const auto pre = (w < result.pre_fec_windows.size())
                           ? result.pre_fec_windows[w]
                           : newcode::WindowBerStats{w, 0, 0, 0.0};
      const auto post = (w < result.post_fec_windows.size())
                            ? result.post_fec_windows[w]
                            : newcode::WindowBerStats{w, 0, 0, 0.0};

      out.aggregated[w].agg_pre_errs += pre.errors;
      out.aggregated[w].agg_pre_bits += pre.total;
      out.aggregated[w].agg_post_errs += post.errors;
      out.aggregated[w].agg_post_bits += post.total;

      out.raw_rows.push_back(RawWindowRow{
          .ebn0_db = task_group.ebn0_db,
          .seed_index = seed_result.seed_index,
          .bitgen_seed = seed_result.bitgen_seed,
          .channel_seed = seed_result.channel_seed,
          .window_idx = w,
          .pre_errs = pre.errors,
          .pre_bits = pre.total,
          .pre_ber = pre.ber,
          .post_errs = post.errors,
          .post_bits = post.total,
          .post_ber = post.ber,
      });
    }
  }

  build_candidates(out.aggregated, out.candidates);
  return out;
}

} // namespace

int main() {
  const std::filesystem::path output_dir = kOutputDir;
  std::filesystem::create_directories(output_dir);

  std::ofstream raw_csv(output_dir / "per_seed_per_window.csv");
  std::ofstream aggregated_csv(output_dir / "aggregated_window_ber.csv");
  std::ofstream candidates_csv(output_dir / "start_window_candidates.csv");
  std::ofstream log_file(output_dir / "probe.log");

  write_raw_header(raw_csv);
  write_aggregated_header(aggregated_csv);
  write_candidates_header(candidates_csv);

  io::DualWriter log(log_file);
  log << "[INFO] BER window probe started\n";
  log << "[INFO] seed_count=" << kSeedCount
      << ", max_parallel_ebn0=" << resolved_parallel_ebn0()
      << ", max_parallel_seeds=" << resolved_parallel_seeds_for_ebn0(resolved_parallel_ebn0())
      << '\n';

  ofec_single::Config base_cfg = make_base_config();
  newcode::PipelineConfig pipeline_cfg = ofec_single::detail::build_pipeline_config(base_cfg);
  pipeline_cfg.quiet = kQuietPipeline;

  std::vector<EbN0TaskGroup> ebn0_task_groups;
  ebn0_task_groups.reserve(kEbN0List.size());
  for (float ebn0_db : kEbN0List) {
    log << "[INFO] prepare Eb/N0=" << ebn0_db << " dB\n";
    EbN0TaskGroup group;
    group.ebn0_db = ebn0_db;
    group.seed_tasks.reserve(kSeedCount);
    for (int seed_index = 0; seed_index < kSeedCount; ++seed_index) {
      ofec_single::Config cfg = base_cfg;
      cfg.ebn0_db = ebn0_db;
      cfg.bitgen_seed = kBitgenSeedBase + seed_index;
      cfg.channel_seed = kChannelSeedBase + seed_index;
      cfg.label = std::string(kLabelPrefix) + "_e" + format_float(ebn0_db) +
                  "_s" + std::to_string(seed_index);
      auto params_opt = ofec_single::detail::build_params(cfg, log);
      if (!params_opt.has_value()) {
        log << "[ERROR] failed to build params for Eb/N0=" << ebn0_db
            << ", seed_index=" << seed_index << '\n';
        return 2;
      }
      group.seed_tasks.push_back(SeedTask{
          .seed_index = seed_index,
          .bitgen_seed = cfg.bitgen_seed,
          .channel_seed = cfg.channel_seed,
          .ebn0_db = ebn0_db,
          .label = cfg.label,
          .params = std::move(*params_opt),
      });
    }
    ebn0_task_groups.push_back(std::move(group));
  }

  const unsigned ebn0_parallel = std::max(1u, resolved_parallel_ebn0());
  const unsigned seed_parallel = std::max(1u, resolved_parallel_seeds_for_ebn0(ebn0_parallel));
  std::vector<EbN0RunResult> ebn0_results;
  ebn0_results.reserve(kEbN0List.size());
  std::vector<std::future<EbN0RunResult>> pending_ebn0;

  for (auto& task_group : ebn0_task_groups) {
    log << "[INFO] schedule Eb/N0=" << task_group.ebn0_db
        << " dB, seed_parallelism=" << seed_parallel << '\n';
    pending_ebn0.push_back(std::async(std::launch::async, [task_group = std::move(task_group), pipeline_cfg, seed_parallel]() mutable {
      return run_ebn0_probe(std::move(task_group), pipeline_cfg, seed_parallel);
    }));
    if (pending_ebn0.size() >= ebn0_parallel) {
      ebn0_results.push_back(pending_ebn0.front().get());
      pending_ebn0.erase(pending_ebn0.begin());
    }
  }
  for (auto& future : pending_ebn0) {
    ebn0_results.push_back(future.get());
  }
  std::sort(ebn0_results.begin(), ebn0_results.end(),
            [](const EbN0RunResult& lhs, const EbN0RunResult& rhs) {
              return lhs.ebn0_db < rhs.ebn0_db;
            });

  for (const auto& ebn0_result : ebn0_results) {
    for (const auto& row : ebn0_result.raw_rows) {
      raw_csv << std::fixed << std::setprecision(6)
              << row.ebn0_db << ','
              << row.seed_index << ','
              << row.bitgen_seed << ','
              << row.channel_seed << ','
              << row.window_idx << ','
              << row.pre_errs << ','
              << row.pre_bits << ','
              << std::setprecision(12) << row.pre_ber << ','
              << row.post_errs << ','
              << row.post_bits << ','
              << row.post_ber << '\n';
    }

    for (const auto& entry : ebn0_result.aggregated) {
      const double agg_pre_ber =
          (entry.agg_pre_bits == 0)
              ? 0.0
              : static_cast<double>(entry.agg_pre_errs) /
                    static_cast<double>(entry.agg_pre_bits);
      const double agg_post_ber =
          (entry.agg_post_bits == 0)
              ? 0.0
              : static_cast<double>(entry.agg_post_errs) /
                    static_cast<double>(entry.agg_post_bits);

      aggregated_csv << std::fixed << std::setprecision(6)
                     << ebn0_result.ebn0_db << ','
                     << entry.window_idx << ','
                     << entry.agg_pre_errs << ','
                     << entry.agg_pre_bits << ','
                     << std::setprecision(12) << agg_pre_ber << ','
                     << entry.agg_post_errs << ','
                     << entry.agg_post_bits << ','
                     << agg_post_ber << '\n';
    }

    for (const auto& candidate : ebn0_result.candidates) {
      candidates_csv << std::fixed << std::setprecision(6)
                     << ebn0_result.ebn0_db << ','
                     << candidate.start_window_k << ','
                     << candidate.post_errs_from_k << ','
                     << candidate.post_bits_from_k << ','
                     << std::setprecision(12) << candidate.post_ber_from_k << '\n';
    }

    log << "[INFO] Eb/N0=" << ebn0_result.ebn0_db
        << " candidates_generated=" << ebn0_result.candidates.size() << '\n';
  }

  log << "[INFO] outputs saved under " << output_dir.string() << "\n";
  return 0;
}
