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
static const std::vector<float> kEbN0List = {3.05f}; // 需要测试的 Eb/N0 列表
static constexpr int kBitgenSeedBase = 20260319;         // 基础比特种子；实际每个 seed 在此基础上递增
static constexpr int kChannelSeedBase = 3182026;         // 基础信道种子；实际每个 seed 在此基础上递增
static constexpr int kSeedCount = 1;                    // 每个 Eb/N0 点重复运行的 seed 数量

// 发射端参数
static constexpr bool        kGenerateRandomBits = true; // true=发送随机信息比特，false=发送全 0 比特

// 信道参数
static constexpr unsigned    kBitsPerSymbol      = 1;    // 每个调制符号携带的比特数：1=BPSK，偶数=QAM

// 早停参数
static constexpr bool        kEnableEarlyStop              = true;   // 早停总开关
static const std::vector<int> kEarlyStopEnableList         = {0,0,0,0,1,1}; // 按 tile 覆盖早停总开关：0=关，非 0=开
static constexpr int         kEarlyStopConditionMode       = 1;      // 早停条件模式：1=v1，2=v2
static const std::vector<int> kEarlyStopConditionModeList  = {};     // 按 tile 覆盖条件模式；空表示沿用全局值
static constexpr int         kEarlyStopActionMode          = 1;      // 早停动作模式：1~6
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
static constexpr int         kChaseL_override          = 6;           // Chase L，-1 表示沿用 Params 默认值
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
    99.857143f, 99.179301f, 99.253626f, 99.119585f, 99.434408f, 99.000000f // 每个 tile 的 early-stop 专用 beta 显式列表
};
static const std::vector<int> kSisoActiveList = {32, 32, 32, 32, 16, 4}; // 每个 tile 允许参与 SISO 的行数预算
static constexpr int  kMuxGroupG          = 1;                             // MUX 分组粒度，1 表示全局池化
static constexpr int  kMuxSchedulingMode  = 0;                             // MUX 调度模式：0=legacy，1=按 early-stop 细节排序
static constexpr int  kMuxPriorityRule    = 0;                             // 新 MUX 的优先级规则：0=更差优先，1=更接近通过优先
static constexpr bool kMuxEnableReconfig  = false;                         // true 表示启用重配置版 MUX 调度
static constexpr int  kMuxBypassScheme    = 1;                             // 旁路边集合方案编号：1=scheme1，2=scheme2
static constexpr bool kHybridEnable       = true;                         // true=方案三软硬混合前置分流开关
static const std::vector<int> kHybridEnableList = {0,0,0,0,1,1};                      // 按 tile 覆盖 hybrid 开关：空=沿用 kHybridEnable
static constexpr float kHybridHardLlrMag = 99.0f;                        // hybrid hard-finish 默认输出 |LLR| 幅度
static const std::vector<float> kHybridHardLlrMagList = {};              // 按 tile 覆盖 hybrid hard-finish |LLR| 幅度；空=沿用 kHybridHardLlrMag
static constexpr newcode::HybridClassifierMode kHybridClassifierMode =
    newcode::HybridClassifierMode::FriendS1S3Classifier;                       // LegacyHardDecode / RepoFastClassifier / FriendS1S3Classifier / FriendS1S3WithS0Classifier
static constexpr newcode::HybridSisoBackfillMode kHybridSisoBackfillMode =
    newcode::HybridSisoBackfillMode::TwoErrorOnly;                                 // Disabled / TwoErrorOnly / OneAndTwoErrorPriority
static constexpr bool kHybridNormalizeSoftOnly = false;                    // true=只归一化 soft rows，false=保持当前兼容行为

// ==================================

// 单个 window 在多 seed 聚合后的统计量。
// 这里保存的是“同一个 Eb/N0 点下，所有 seed 累加后”的 per-window 误差和比特数。
struct AggregatedWindowStats {
  std::size_t window_idx = 0;
  std::size_t agg_pre_errs = 0;
  std::size_t agg_pre_bits = 0;
  std::size_t agg_post_errs = 0;
  std::size_t agg_post_bits = 0;
};

// 单个 seed 的完整任务输入。
// 一个 seed 对应一组独立的 bitgen/channel 随机种子，以及其专属的 Params。
struct SeedTask {
  int seed_index = 0;
  int bitgen_seed = 0;
  int channel_seed = 0;
  float ebn0_db = 0.0f;
  std::string label;
  newcode::Params params;
};

// 单个 seed 跑完后的结果。
// 这里只保存 pipeline 的原始输出，后续再由外层统一做聚合和导出。
struct SeedRunResult {
  int seed_index = 0;
  int bitgen_seed = 0;
  int channel_seed = 0;
  newcode::PipelineResult pipeline_result;
};

// 逐 window / tile 的原始导出行。
// 这份数据是给 Matlab / Python 画时域曲线或做窗口 / tile 分布分析用的。
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

// 单个 seed 的逐 tile early-stop 样本行。
// 这份数据用于和 per_seed_per_tile.csv 对齐，观察 early-stop 命中数在时域上的变化。
struct TileEarlyStopSampleRow {
  float ebn0_db = 0.0f;
  int seed_index = 0;
  int bitgen_seed = 0;
  int channel_seed = 0;
  std::size_t invocation = 0;
  std::size_t tile_index = 0;
  std::size_t rows_total = 0;
  std::size_t rows_passed = 0;
  std::size_t rows_hard_finish = 0;
  std::size_t rows_need_siso_before_mux = 0;
  std::size_t rows_unscheduled = 0;
};

// 某个 Eb/N0 点的完整结果。
// 包含：
// - 聚合后的 per-window 统计
// - 聚合后的 per-tile 统计
// - 原始逐 seed、逐 window / tile 行
struct EbN0RunResult {
  float ebn0_db = 0.0f;
  newcode::BerStats pre_fec{};
  newcode::BerStats pre_fec_quantized_hard{};
  newcode::BerStats post_fec{};
  bool has_pre_fec_quantized_hard = false;
  std::vector<AggregatedWindowStats> aggregated;
  std::vector<AggregatedWindowStats> tile_aggregated;
  std::vector<RawWindowRow> raw_rows;
  std::vector<RawWindowRow> raw_tile_rows;
  std::vector<TileEarlyStopSampleRow> tile_early_stop_rows;
};

// 一个 Eb/N0 下的任务组。
// 组内包含多个 seed 任务，供 seed 并行调度使用。
struct EbN0TaskGroup {
  float ebn0_db = 0.0f;
  std::vector<SeedTask> seed_tasks;
};

// 把浮点数格式化成稳定的十进制字符串。
// 用于日志和任务名，避免默认格式在不同平台/编译器上产生差异。
std::string format_float(float value) {
  std::ostringstream oss;
  oss << std::fixed << std::setprecision(6) << value;
  return oss.str();
}

// 确保聚合窗口数组长度足够。
// 如果当前结果里还没有某个 window 的统计，就先补一个空槽位。
void ensure_aggregated_size(std::vector<AggregatedWindowStats>& stats,
                            std::size_t size) {
  while (stats.size() < size) {
    AggregatedWindowStats entry;
    entry.window_idx = stats.size();
    stats.push_back(entry);
  }
}

// 计算外层 Eb/N0 的并行度。
// 若用户显式给了上限，就优先用上限；否则根据硬件并发和 Eb/N0 点数自动裁剪。
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

// 计算单个 Eb/N0 任务内可用的 seed 并行度。
// 如果外层已经并行跑了多个 Eb/N0，就会把硬件并发再均摊给 seed。
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

// 输出一个 Eb/N0 点的整体 BER 汇总，风格接近 ofec_single。
void log_ebn0_result_summary(io::DualWriter& log,
                             const EbN0RunResult& result) {
  log << "[RESULT] Eb/N0=" << result.ebn0_db << " dB"
      << " | Pre-FEC BER=" << result.pre_fec.ber
      << " (errs=" << result.pre_fec.errors
      << "/" << result.pre_fec.total << ")";
  if (result.has_pre_fec_quantized_hard) {
    log << " | Pre-FEC BER (quantized hard)="
        << result.pre_fec_quantized_hard.ber
        << " (errs=" << result.pre_fec_quantized_hard.errors
        << "/" << result.pre_fec_quantized_hard.total << ")";
  }
  log << " | Post-FEC BER=" << result.post_fec.ber
      << " (errs=" << result.post_fec.errors
      << "/" << result.post_fec.total << ")\n";

  log << "[DETAIL] Eb/N0=" << result.ebn0_db
      << " seeds=" << kSeedCount
      << " windows=" << result.aggregated.size()
      << " tiles=" << result.tile_aggregated.size()
      << " raw_window_rows=" << result.raw_rows.size()
      << " raw_tile_rows=" << result.raw_tile_rows.size()
      << " early_stop_sample_rows=" << result.tile_early_stop_rows.size()
      << "\n";
}

// 输出所有 Eb/N0 点累加后的总汇总。
void log_probe_total_summary(io::DualWriter& log,
                             const std::vector<EbN0RunResult>& results) {
  newcode::BerStats pre{};
  newcode::BerStats pre_quantized{};
  newcode::BerStats post{};
  bool has_pre_quantized = false;

  for (const auto& result : results) {
    pre.errors += result.pre_fec.errors;
    pre.total += result.pre_fec.total;
    if (result.has_pre_fec_quantized_hard) {
      has_pre_quantized = true;
      pre_quantized.errors += result.pre_fec_quantized_hard.errors;
      pre_quantized.total += result.pre_fec_quantized_hard.total;
    }
    post.errors += result.post_fec.errors;
    post.total += result.post_fec.total;
  }

  pre.ber = (pre.total == 0)
                ? 0.0
                : static_cast<double>(pre.errors) /
                      static_cast<double>(pre.total);
  pre_quantized.ber =
      (pre_quantized.total == 0)
          ? 0.0
          : static_cast<double>(pre_quantized.errors) /
                static_cast<double>(pre_quantized.total);
  post.ber = (post.total == 0)
                 ? 0.0
                 : static_cast<double>(post.errors) /
                       static_cast<double>(post.total);

  log << "[SUMMARY] Total Eb/N0 points=" << results.size()
      << " | Pre-FEC BER=" << pre.ber
      << " (errs=" << pre.errors << "/" << pre.total << ")";
  if (has_pre_quantized) {
    log << " | Pre-FEC BER (quantized hard)=" << pre_quantized.ber
        << " (errs=" << pre_quantized.errors
        << "/" << pre_quantized.total << ")";
  }
  log << " | Post-FEC BER=" << post.ber
      << " (errs=" << post.errors << "/" << post.total << ")\n";
}

// 写逐 seed、逐 window 的原始 CSV 表头。
// 这份表主要用于时域曲线、窗口分布和 Matlab 事后分析。
void write_raw_header(std::ofstream& out) {
  out << "ebn0_db,seed_index,bitgen_seed,channel_seed,window_idx,"
         "pre_errs,pre_bits,pre_ber,post_errs,post_bits,post_ber\n";
  out.flush();
}

// 写逐 seed、逐 tile 的原始 CSV 表头。
void write_raw_tile_header(std::ofstream& out) {
  out << "ebn0_db,seed_index,bitgen_seed,channel_seed,tile_idx,"
         "pre_errs,pre_bits,pre_ber,post_errs,post_bits,post_ber\n";
  out.flush();
}

// 写 tile early-stop 样本 CSV 表头。
void write_tile_early_stop_samples_header(std::ofstream& out) {
  out << "ebn0_db,seed_index,bitgen_seed,channel_seed,invocation,tile_index,"
         "rows_total,rows_passed,rows_hard_finish,rows_need_siso_before_mux,"
         "rows_unscheduled\n";
  out.flush();
}

// 写按 window 聚合后的 CSV 表头。
// 这份表用于观察同一 Eb/N0 下，各 window 的总体 BER 走势。
void write_aggregated_header(std::ofstream& out) {
  out << "ebn0_db,window_idx,"
         "agg_pre_errs,agg_pre_bits,agg_pre_ber,"
         "agg_post_errs,agg_post_bits,agg_post_ber\n";
  out.flush();
}

// 写按 tile 聚合后的 CSV 表头。
void write_aggregated_tile_header(std::ofstream& out) {
  out << "ebn0_db,tile_idx,"
         "agg_pre_errs,agg_pre_bits,agg_pre_ber,"
         "agg_post_errs,agg_post_bits,agg_post_ber\n";
  out.flush();
}

// 聚合一个 seed 的整体 BER 到当前 Eb/N0 结果。
void accumulate_overall_ber(EbN0RunResult& out,
                            const newcode::PipelineResult& result) {
  out.pre_fec.errors += result.pre_fec.errors;
  out.pre_fec.total += result.pre_fec.total;
  out.pre_fec.ber = (out.pre_fec.total == 0)
                        ? 0.0
                        : static_cast<double>(out.pre_fec.errors) /
                              static_cast<double>(out.pre_fec.total);

  if (result.has_pre_fec_quantized_hard) {
    out.has_pre_fec_quantized_hard = true;
    out.pre_fec_quantized_hard.errors += result.pre_fec_quantized_hard.errors;
    out.pre_fec_quantized_hard.total += result.pre_fec_quantized_hard.total;
    out.pre_fec_quantized_hard.ber =
        (out.pre_fec_quantized_hard.total == 0)
            ? 0.0
            : static_cast<double>(out.pre_fec_quantized_hard.errors) /
                  static_cast<double>(out.pre_fec_quantized_hard.total);
  }

  out.post_fec.errors += result.post_fec.errors;
  out.post_fec.total += result.post_fec.total;
  out.post_fec.ber = (out.post_fec.total == 0)
                         ? 0.0
                         : static_cast<double>(out.post_fec.errors) /
                               static_cast<double>(out.post_fec.total);
}

// 构造 probe 用的基础 ofec_single 配置。
// 这里复用单跑入口的参数接线，但把 seed / EbN0 / 输出等改成 probe 需要的形式。
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
      .hybrid_enable = kHybridEnable,
      .hybrid_enable_list = kHybridEnableList,
      .hybrid_hard_llr_mag = kHybridHardLlrMag,
      .hybrid_hard_llr_mag_list = kHybridHardLlrMagList,
      .hybrid_classifier_mode = kHybridClassifierMode,
      .hybrid_siso_backfill_mode = kHybridSisoBackfillMode,
      .hybrid_normalize_soft_only = kHybridNormalizeSoftOnly,
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

// 在一个 Eb/N0 点内部并行跑多个 seed。
// 每个 seed 完整执行一次 pipeline，然后统一回收结果并按 seed_index 排序。
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
  //限制的是 seed 任务并行数
  for (auto& task : tasks) {
    pending.push_back(launch_task(std::move(task)));
    if (pending.size() >= max_parallel_seeds) {
      seed_results.push_back(pending.front().get());
      pending.erase(pending.begin());
    }
  }
  //收集剩余 seed 异步任务结果，并按 seed 编号排序
  for (auto& future : pending) {
    seed_results.push_back(future.get());
  }
  std::sort(seed_results.begin(), seed_results.end(),
            [](const SeedRunResult& lhs, const SeedRunResult& rhs) {
              return lhs.seed_index < rhs.seed_index;
            });
  return seed_results;
}

// 跑完一个 Eb/N0 点的完整 probe。
// 内部先并行跑 seed，再把所有 seed 的 per-window 结果做聚合和候选构造。
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
    accumulate_overall_ber(out, result);

    const std::size_t max_windows =
        std::max(result.pre_fec_windows.size(), result.post_fec_windows.size());
    ensure_aggregated_size(out.aggregated, max_windows);
    const std::size_t max_tiles =
        std::max(result.pre_fec_tile_windows.size(), result.post_fec_tile_windows.size());
    ensure_aggregated_size(out.tile_aggregated, max_tiles);

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

    for (std::size_t t = 0; t < max_tiles; ++t) {
      const auto pre = (t < result.pre_fec_tile_windows.size())
                           ? result.pre_fec_tile_windows[t]
                           : newcode::WindowBerStats{t, 0, 0, 0.0};
      const auto post = (t < result.post_fec_tile_windows.size())
                            ? result.post_fec_tile_windows[t]
                            : newcode::WindowBerStats{t, 0, 0, 0.0};

      out.tile_aggregated[t].agg_pre_errs += pre.errors;
      out.tile_aggregated[t].agg_pre_bits += pre.total;
      out.tile_aggregated[t].agg_post_errs += post.errors;
      out.tile_aggregated[t].agg_post_bits += post.total;

      out.raw_tile_rows.push_back(RawWindowRow{
          .ebn0_db = task_group.ebn0_db,
          .seed_index = seed_result.seed_index,
          .bitgen_seed = seed_result.bitgen_seed,
          .channel_seed = seed_result.channel_seed,
          .window_idx = t,
          .pre_errs = pre.errors,
          .pre_bits = pre.total,
          .pre_ber = pre.ber,
          .post_errs = post.errors,
          .post_bits = post.total,
          .post_ber = post.ber,
      });
    }

    for (const auto& sample : result.tile_early_stop_samples) {
      out.tile_early_stop_rows.push_back(TileEarlyStopSampleRow{
          .ebn0_db = task_group.ebn0_db,
          .seed_index = seed_result.seed_index,
          .bitgen_seed = seed_result.bitgen_seed,
          .channel_seed = seed_result.channel_seed,
          .invocation = sample.invocation,
          .tile_index = sample.tile_index,
          .rows_total = sample.rows_total,
          .rows_passed = sample.rows_passed,
          .rows_hard_finish = sample.rows_hard_finish,
          .rows_need_siso_before_mux = sample.rows_need_siso_before_mux,
          .rows_unscheduled = sample.rows_unscheduled,
      });
    }
  }
  return out;
}

} // namespace

// 程序入口：
// 1) 准备输出目录和 CSV/日志文件
// 2) 构造每个 Eb/N0、每个 seed 的任务
// 3) 先按 Eb/N0 并行，再在每个 Eb/N0 内按 seed 并行
// 4) 输出逐 window 原始结果、聚合结果和候选起点统计
int main() {
  // 1) 准备输出目录，确保后续 CSV 和日志能够落盘。
  const std::filesystem::path output_dir = kOutputDir;
  std::filesystem::create_directories(output_dir);

  // 2) 打开原始逐 seed CSV、tile CSV、聚合 CSV 和候选起点 CSV、日志文件。
  std::ofstream raw_csv(output_dir / "per_seed_per_window.csv");
  std::ofstream raw_tile_csv(output_dir / "per_seed_per_tile.csv");
  std::ofstream tile_early_stop_samples_csv(output_dir / "per_seed_tile_early_stop_samples.csv");
  std::ofstream aggregated_csv(output_dir / "aggregated_window_ber.csv");
  std::ofstream aggregated_tile_csv(output_dir / "aggregated_tile_ber.csv");
  std::ofstream log_file(output_dir / "probe.log");

  // 3) 写各个 CSV 的表头。
  write_raw_header(raw_csv);
  write_raw_tile_header(raw_tile_csv);
  write_tile_early_stop_samples_header(tile_early_stop_samples_csv);
  write_aggregated_header(aggregated_csv);
  write_aggregated_tile_header(aggregated_tile_csv);

  // 4) 打印 probe 启动时的基本信息，包括 seed 数和并行配置。
  io::DualWriter log(log_file);
  log << "[INFO] BER window probe started\n";
  log << "[INFO] seed_count=" << kSeedCount
      << ", max_parallel_ebn0=" << resolved_parallel_ebn0()
      << ", max_parallel_seeds=" << resolved_parallel_seeds_for_ebn0(resolved_parallel_ebn0())
      << '\n';

  // 5) 先构造一份基础 ofec_single 配置，再转换成 pipeline 配置。
  ofec_single::Config base_cfg = make_base_config();
  newcode::PipelineConfig pipeline_cfg = ofec_single::detail::build_pipeline_config(base_cfg);
  pipeline_cfg.quiet = kQuietPipeline;

  // 6) 按 Eb/N0 组织任务组；每个 Eb/N0 下再展开多个 seed 任务。
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

  // 7) 先按 Eb/N0 并行，再在每个 Eb/N0 内按 seed 并行跑完整 probe。
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
    }));//移交给一个异步线程执行
    //并行限流代码
    if (pending_ebn0.size() >= ebn0_parallel) {
      ebn0_results.push_back(pending_ebn0.front().get());
      pending_ebn0.erase(pending_ebn0.begin());
    }
  }
  for (auto& future : pending_ebn0) {
    ebn0_results.push_back(future.get());
  }//收集所有剩余异步任务结果
  //然后按 Eb/N0 从低到高排列结果
  std::sort(ebn0_results.begin(), ebn0_results.end(),
            [](const EbN0RunResult& lhs, const EbN0RunResult& rhs) {
              return lhs.ebn0_db < rhs.ebn0_db;
            });

  // 8) 把逐 seed、逐 window / tile 的原始结果和聚合结果写回 CSV。
  for (const auto& ebn0_result : ebn0_results) {
    log_ebn0_result_summary(log, ebn0_result);

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

    for (const auto& row : ebn0_result.raw_tile_rows) {
      raw_tile_csv << std::fixed << std::setprecision(6)
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

    for (const auto& sample : ebn0_result.tile_early_stop_rows) {
      tile_early_stop_samples_csv << std::fixed << std::setprecision(6)
                                  << sample.ebn0_db << ','
                                  << sample.seed_index << ','
                                  << sample.bitgen_seed << ','
                                  << sample.channel_seed << ','
                                  << sample.invocation << ','
                                  << sample.tile_index << ','
                                  << sample.rows_total << ','
                                  << sample.rows_passed << ','
                                  << sample.rows_hard_finish << ','
                                  << sample.rows_need_siso_before_mux << ','
                                  << sample.rows_unscheduled << '\n';
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

    for (const auto& entry : ebn0_result.tile_aggregated) {
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

      aggregated_tile_csv << std::fixed << std::setprecision(6)
                          << ebn0_result.ebn0_db << ','
                          << entry.window_idx << ','
                          << entry.agg_pre_errs << ','
                          << entry.agg_pre_bits << ','
                          << std::setprecision(12) << agg_pre_ber << ','
                          << entry.agg_post_errs << ','
                          << entry.agg_post_bits << ','
                          << agg_post_ber << '\n';
    }
  }
  log_probe_total_summary(log, ebn0_results);

  // 9) 收尾日志，提示输出已经写到哪个目录。
  log << "[INFO] outputs saved under " << output_dir.string() << "\n";
  return 0;
}
