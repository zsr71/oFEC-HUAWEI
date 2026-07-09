#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstdint>
#include <filesystem>
#include <fstream>
#include <future>
#include <iomanip>
#include <iostream>
#include <limits>
#include <sstream>
#include <stdexcept>
#include <string>
#include <thread>
#include <utility>
#include <vector>

#include "mux_bypass_edges.hpp"
#include "newcode/io/ensure_dir.hpp"
#include "newcode/ofec_single_runner.hpp"
#include "newcode/two_stream_shared_runner.hpp"
#include "newcode/utils/now_stamp.hpp"
#include "ofec_sweep_detail.hpp"

namespace {

// ======== 用户可调参数区域 ========
// 这个 app 的定位是：
// - 复用现有 two_stream_shared 主流程
// - 复用 sweep3 的 low-BER chunk 聚合框架
// - 第一版只关心 A/B 两路 BER 曲线，不聚合 shared 可观测性

// 运行标签
static constexpr const char* kLabel = "two_native_stream_shared_sweep"; // 运行标签：日志文件、CSV 文件和 chunk 内部 label 的统一前缀

// 发端参数：交织 / 比特源 / 调制入口
static constexpr const char* kInterleaverName = "identity"; // 交织器名称：identity=不交织；后续若要对齐别的入口，可改成仓库中已有的其他交织器名
static constexpr unsigned kBitsPerSymbol = 1;               // 每个调制符号携带的比特数：1=BPSK；偶数时通常表示对应阶数的 QAM
static constexpr bool kGenerateRandomBits = true;           // true=每个 chunk 发送随机信息比特；false=发送全 0，比对极简回归时可用
static constexpr bool kNormalizeKnownPrefixTail = false;    // 是否对 known-prefix 后的尾部 LLR 做归一化；false=保持当前 two-stream 单点入口口径

// 信道参数：噪声强度 / 扫描范围
static constexpr float kEbN0Start = 3.00f;                 // 扫描起始 Eb/N0（dB）；通常从单点 app 已经能跑通的附近起扫
static constexpr float kEbN0End = 3.20f;                   // 扫描结束 Eb/N0（dB）；第一版先做窄区间，方便和单点结果逐点对照
static constexpr int kEbN0Points = 5;                      // Eb/N0 采样点数；含首尾端点，5 表示把 [start,end] 均匀切成 5 个点

// 量化参数：只影响接收端 LLR 量化口径
static constexpr std::size_t kLlrBits = 6;                 // LLR 位宽：16=浮点直通；2~15=qfloat 量化；当前与 two-stream 单点入口保持 6 bit
static constexpr float kQuantClipRatio = 0.5f;             // 动态 clip 比例：0=禁用自适应 clip；非 0 表示按全部 LLR 样本分位口径计算共享 clip

// 双流各自的基础随机种子
static constexpr int kBitgenSeedA = 20260320;              // Stream A 的基础信息比特种子；每个 chunk 会基于它再派生 chunk 级 seed
static constexpr int kChannelSeedA = 3182027;              // Stream A 的基础信道种子；每个 chunk 会基于它再派生 chunk 级 seed
static constexpr int kBitgenSeedB = 20260319;              // Stream B 的基础信息比特种子；与 A 分开是为了保持两路原生独立
static constexpr int kChannelSeedB = 3182026;              // Stream B 的基础信道种子；与 A 分开是为了避免两路信道相关

// 低 BER 聚合参数：每个 Eb/N0 点会拆成多个 chunk，直到满足停止条件
static constexpr std::size_t kChunkNumInfoBits =
    8u * 132u * 16u * 111u;                                // 每个 chunk、每一路生成的信息比特数；值越大，每个 chunk 更重，调度粒度更粗
static constexpr std::size_t kTargetPostErrors = 50;       // 每一路累计到这么多 post-FEC 错误后即可停止该点；值越大，BER 估计更稳但耗时更长
static constexpr std::size_t kMaxPostFecTotalBits = 2e8;   // 每一路允许累计比较的最大 post-FEC 比特数；用于防止极低 BER 点无限跑下去
static constexpr unsigned kMaxTotalWorkers = 0;            // 全局同时运行的 chunk 数上限；0=自动按机器并发数/NTHREADS 决定
static constexpr unsigned kMaxInflightChunksPerPoint = 16; // 单个 Eb/N0 点默认最多挂起多少个 chunk；0 不建议在这里使用，当前保持显式上限
static constexpr bool kEnableDynamicInflightPerPoint = true; // true=剩余点变少时自动提高单点并发，尽量吃满 worker
static constexpr unsigned kMaxDynamicInflightChunksPerPoint = 0; // 动态单点并发封顶；0=不额外封顶，最多提升到全局 worker 数
static constexpr bool kEnableZeroErrorUpperBound = true;   // true=当某一路 post-FEC 仍为 0 错时，用上置信界做提前停止
static constexpr double kTargetBerUpperBound = 1e-8;       // 零错上界目标；若上界已经低于该目标，可接受“零错但非 BER=0”的提前停止
static constexpr double kConfidenceLevel = 0.95;           // 零错上界使用的置信水平；0.95 表示采用 95% 置信界

// 解码参数：decoder 选择、Chase 参数、交织口径
static constexpr int kChaseL = 6;                          // Chase L：选择多少个最不可靠位置；常见值如 4/5/6
static constexpr int kChaseNTestOverride = -1;            // Chase 测试序列数量：-1=按 2^L 自动生成；非负值可强制覆盖测试 pattern 数
static constexpr int kChaseTopkKeep = 8;                  // chase_topk_pruned 中保留参与外信息计算的 Top-K 候选数；当前 baseline 虽不用 topk，也保持参数齐全
static constexpr int kChaseGroupMinimaBits = 3;           // chase_group_minima 的分组 bit 数；当前 baseline 虽不用 group_minima，也保持参数齐全
static constexpr const char* kDecoderName =
    "two_stream_shared_chase_baseline";                   // 当前 two-stream shared 主流程入口名；底层会映射到 shared runner 支持的 core
static constexpr bool kNormalizeExtrinsic = false;         // 是否对 Chase 输出的 extrinsic 做归一化；false=保持当前单点入口口径

// Early-stop 参数：总开关 -> 条件 -> 动作 -> 细节参数
static constexpr bool kEnableEarlyStop = true;             // early-stop 总开关：true=启用；false=完全关闭，所有行都继续进入后续 hybrid/MUX/soft path
static const std::vector<int> kEarlyStopEnableList = {1, 1, 1, 1, 1, 1}; // 按 tile 覆盖 early-stop 开关：1=开，0=关；长度应与 TILES_PER_WIN 一致
static constexpr int kEarlyStopConditionMode = 1;          // early-stop 条件模式：1=v1（BCH/overall 条件）；2=v2（LLR 阈值 + 不可靠位数）
static const std::vector<int> kEarlyStopConditionModeList = {}; // 按 tile 覆盖条件模式；空=所有 tile 都沿用 kEarlyStopConditionMode
static constexpr int kEarlyStopActionMode = 1;             // 命中 early-stop 后的动作模式：1~8；完整语义见单点 app 注释
static const std::vector<int> kEarlyStopActionModeList = {}; // 按 tile 覆盖动作模式；空=所有 tile 都沿用 kEarlyStopActionMode
static constexpr int kEarlyStopBindGroupSize = 1;          // 条件1专用组绑定大小：1=逐 row 判；4=四个 row 绑定后整体通过才 early-stop
static const std::vector<int> kEarlyStopBindGroupSizeList = {}; // 按 tile 覆盖组绑定大小；空=沿用全局默认值
static constexpr bool kEarlyStopCondV1RequireBch = true;   // 条件1里是否要求 BCH syndrome 为 0；true=要求 BCH 通过
static constexpr bool kEarlyStopCondV1RequireOverall = true; // 条件1里是否要求 overall parity 一致；true=要求 overall 也通过
static constexpr float kEarlyStopV2LlrAbsThreshold = 0.5f; // 条件2里“不可靠位”的 |LLR| 阈值；值越大，越容易把 bit 判成不可靠
static constexpr int kEarlyStopV2MaxUnreliableBits = 8;    // 条件2允许的不可靠 bit 数上限；超过该数则认为这行不能 early-stop
static constexpr bool kEarlyStopCondV2IncludeOverall = true; // 条件2统计不可靠位时是否把 overall parity bit 也纳入
static const std::vector<float> kEarlyStopActionBetaExplicit = {}; // 每个 tile 的 early-stop 动作 beta 显式列表；空=沿用普通 beta 配置
static constexpr float kEarlyStopActionResidualDivisor = 1.0f; // 动作2里 residual 的除数；当前 sweep 只是显式保留，与单点入口对齐
static constexpr float kEarlyStopActionHardLlrMag = 1.0f; // 动作3里硬解成功后输出的固定 |LLR| 幅度；当前 sweep 只是显式保留，与单点入口对齐

// alpha / beta / shared SISO 预算
static const std::vector<float> kAlphaExplicit = {
    0.428571f, 0.447738f, 0.482782f, 0.528162f, 0.581902f, 0.642857f}; // 每个 tile 的 alpha 显式列表；长度需与 TILES_PER_WIN 一致
static const std::vector<float> kBetaExplicit = {
    2.857143f, 6.179301f, 12.253626f, 20.119585f, 29.434408f, 40.0f}; // 每个 tile 的 Chase/fallback beta 显式列表；长度需与 TILES_PER_WIN 一致
static const std::vector<int> kSisoActiveList = {64, 64, 64, 32, 16, 8}; // 每个 tile 在 shared 64-row 域里可参与 SISO 的预算；前 3 个 tile 全开，后 3 个 tile 逐步收紧
static const std::vector<int> kHiHoActiveList = {64, 64, 64, 64, 64, 64}; // 每个 tile 在 shared 64-row 域里可参与 HIHO 硬解码的预算

// MUX 调度参数
static constexpr int kMuxGroupG = 1;                    // MUX 分组数：1=全局池化；>1 时会把共享 code 域均分成多个组分别切预算
static constexpr int kMuxSchedulingMode = 0;            // MUX 调度模式：0=legacy 顺序裁剪；1=按 early-stop 细节排序后再裁剪
static constexpr int kMuxPriorityRule = 0;              // MUX 优先级规则：0=harder_first；1=near_threshold_first；仅在 scheduling_mode=1 时真正起作用
static constexpr bool kMuxEnableReconfig = false;       // true=启用 staged reconfig 调度；false=直接用 grouped budget 裁剪
static const std::vector<newcode::mux::MuxEdge> kMuxExtraBypassEdges =
    app_mux::bypass_edges_for_scheme(1);                // reconfig 模式使用的额外旁路边集合；当前 mode=false 时只是保持参数区完整

// Hybrid 参数
static constexpr bool kHybridEnable = true;             // hybrid 总开关：true=允许做前置硬分类/硬完成；false=全部保留到 soft path
static const std::vector<int> kHybridEnableList = {0, 0, 0, 1, 1, 1}; // 按 tile 覆盖 hybrid 开关：当前前 3 个 tile 关，后 3 个 tile 开
static constexpr float kHybridHardLlrMag = 99.0f;       // hybrid hard-finish 默认输出 |LLR| 幅度；通常取较大值以体现“强硬判”含义
static const std::vector<float> kHybridHardLlrMagList = {}; // 按 tile 覆盖 hybrid hard-finish |LLR| 幅度；空=所有 tile 都沿用全局默认值
static constexpr newcode::HybridClassifierMode kHybridClassifierMode =
    newcode::HybridClassifierMode::FriendS1S3WithS0Classifier; // hybrid 分类器模式：当前使用 FriendS1S3WithS0Classifier，与单点入口保持一致
static constexpr newcode::HybridSisoBackfillMode kHybridSisoBackfillMode =
    newcode::HybridSisoBackfillMode::OneAndTwoErrorPriority; // hybrid SISO 回填模式：Disabled / TwoErrorOnly / OneAndTwoErrorPriority / ParityOneAndTwoErrorPriority
static constexpr bool kHybridNormalizeSoftOnly = false;  // true=只归一化 soft rows；false=保持当前单点入口兼容口径

// Debug 参数
static constexpr bool kQuietConsole = false;            // true=减少控制台输出；false=保留每个 chunk 的 BER 输出，便于第一版排查
static constexpr bool kDecoderTraceEnable = false;      // decoder trace 总开关；true 时会把 trace 配置灌进 two-stream shared 主流程
static constexpr long kDecoderTraceRow = -1;            // 追踪目标的全局 row；-1=不指定具体 row
static constexpr long kDecoderTraceCol = -1;            // 追踪目标的全局 col；-1=不指定具体 col
static constexpr bool kDecoderTraceLogRead = false;     // 是否打印 tile 读取映射；追查 merged/slice back 时有用
static constexpr bool kDecoderTraceLogWrite = false;    // 是否打印 tile 写回映射；追查 writeback 时有用
static constexpr bool kDecoderTraceLogMismatch = false; // 是否打印同坐标写回不一致告警；仅在细查内部行为时打开

// ==================================

// ======== 结果结构 / 运行状态 ========
enum class StopReason {
  TargetPostErrors,
  MaxPostFecBits,
  ZeroErrorUpperBound,
  NoData
};

struct StreamAggregate {
  newcode::BerStats pre_fec{};
  newcode::BerStats pre_fec_quantized_hard{};
  newcode::BerStats post_fec{};
  bool has_pre_fec_quantized_hard = false;
  bool post_ber_is_upper_bound = false;
  double post_ber_upper_bound = std::numeric_limits<double>::quiet_NaN();
  StopReason stop_reason = StopReason::NoData;
};

struct ChunkResult {
  std::size_t chunk_index = 0;
  int bitgen_seed_a = 0;
  int channel_seed_a = 0;
  int bitgen_seed_b = 0;
  int channel_seed_b = 0;
  newcode::two_stream_shared::Result result;
};

struct AggregatedPointResult {
  float ebn0_db = std::numeric_limits<float>::quiet_NaN();
  StreamAggregate stream_a;
  StreamAggregate stream_b;
  std::size_t chunks_completed = 0;
  double elapsed_seconds = 0.0;
};

struct PointState {
  AggregatedPointResult point;
  std::size_t next_chunk_index = 0;
  unsigned inflight_chunks = 0;
  bool started = false;
  bool launch_stopped = false;
  bool completed = false;
  std::chrono::steady_clock::time_point start_time{};
};

struct ActiveChunk {
  std::size_t point_index = 0;
  std::size_t chunk_index = 0;
  std::future<ChunkResult> future;
};

// ======== 基础格式 / 字符串辅助 ========
std::string stop_reason_to_string(StopReason reason) {
  switch (reason) {
    case StopReason::TargetPostErrors:
      return "target_post_errors";
    case StopReason::MaxPostFecBits:
      return "max_post_fec_total_bits";
    case StopReason::ZeroErrorUpperBound:
      return "zero_error_upper_bound";
    case StopReason::NoData:
    default:
      return "no_data";
  }
}

std::string join_vec(const std::vector<float>& values, char sep, int precision) {
  return ofec_sweep::detail::join_vec(values, sep, precision);
}

std::string join_int_vec(const std::vector<int>& values, char sep) {
  std::ostringstream oss;
  for (std::size_t i = 0; i < values.size(); ++i) {
    if (i > 0) {
      oss << sep;
    }
    oss << values[i];
  }
  return oss.str();
}

const char* hybrid_classifier_mode_name(newcode::HybridClassifierMode mode) {
  switch (mode) {
    case newcode::HybridClassifierMode::LegacyHardDecode:
      return "legacy_hard_decode";
    case newcode::HybridClassifierMode::RepoFastClassifier:
      return "repo_fast_classifier";
    case newcode::HybridClassifierMode::FriendS1S3Classifier:
      return "friend_s1s3_classifier";
    case newcode::HybridClassifierMode::FriendS1S3WithS0Classifier:
      return "friend_s1s3_with_s0_classifier";
  }
  return "unknown";
}

const char* hybrid_siso_backfill_mode_name(newcode::HybridSisoBackfillMode mode) {
  switch (mode) {
    case newcode::HybridSisoBackfillMode::Disabled:
      return "disabled";
    case newcode::HybridSisoBackfillMode::TwoErrorOnly:
      return "two_error_only";
    case newcode::HybridSisoBackfillMode::OneAndTwoErrorPriority:
      return "one_and_two_error_priority";
    case newcode::HybridSisoBackfillMode::ParityOneAndTwoErrorPriority:
      return "parity_one_and_two_error_priority";
  }
  return "unknown";
}

std::string format_duration(std::chrono::duration<double> duration) {
  if (duration.count() < 0.0) {
    duration = std::chrono::duration<double>(0.0);
  }
  using Seconds = std::chrono::seconds;
  const auto total_seconds = std::chrono::duration_cast<Seconds>(duration).count();
  const auto hours = total_seconds / 3600;
  const auto minutes = (total_seconds % 3600) / 60;
  const auto seconds = total_seconds % 60;

  std::ostringstream oss;
  oss << std::setfill('0');
  if (hours > 0) {
    oss << hours << ':' << std::setw(2) << minutes << ':' << std::setw(2) << seconds;
  } else {
    oss << minutes << ':' << std::setw(2) << seconds;
  }
  return oss.str();
}

// ======== chunk runner：seed 派生 / 单 chunk 执行 ========
std::uint32_t mix_u32(std::uint32_t x) {
  x ^= x >> 16;
  x *= 0x7feb352dU;
  x ^= x >> 15;
  x *= 0x846ca68bU;
  x ^= x >> 16;
  return x;
}

int derive_seed(int base, std::size_t chunk_index, std::uint32_t salt) {
  const std::uint32_t raw =
      mix_u32(static_cast<std::uint32_t>(base) ^
              static_cast<std::uint32_t>(chunk_index * 0x9e3779b9ULL) ^
              salt);
  const int positive = static_cast<int>(raw & 0x7fffffffU);
  return positive == 0 ? 1 : positive;
}

newcode::two_stream_shared::Config build_chunk_config(float ebn0_db,
                                                      std::size_t chunk_index) {
  ofec_single::Config base_cfg{};
  base_cfg.label = kLabel;
  base_cfg.ebn0_db = ebn0_db;
  base_cfg.chaseL_override = kChaseL;
  base_cfg.chase_n_test_override = kChaseNTestOverride;
  base_cfg.chase_topk_keep = kChaseTopkKeep;
  base_cfg.chase_group_minima_bits = kChaseGroupMinimaBits;
  base_cfg.normalize_extrinsic = kNormalizeExtrinsic;
  base_cfg.bits_per_symbol = kBitsPerSymbol;
  base_cfg.bitgen_seed = kBitgenSeedA;
  base_cfg.channel_seed = kChannelSeedA;
  base_cfg.enable_early_stop = kEnableEarlyStop;
  base_cfg.early_stop_enable_list = kEarlyStopEnableList;
  base_cfg.early_stop_condition_mode = kEarlyStopConditionMode;
  base_cfg.early_stop_condition_mode_list = kEarlyStopConditionModeList;
  base_cfg.early_stop_action_mode = kEarlyStopActionMode;
  base_cfg.early_stop_action_mode_list = kEarlyStopActionModeList;
  base_cfg.early_stop_bind_group_size = kEarlyStopBindGroupSize;
  base_cfg.early_stop_bind_group_size_list = kEarlyStopBindGroupSizeList;
  base_cfg.early_stop_cond_v1_require_bch = kEarlyStopCondV1RequireBch;
  base_cfg.early_stop_cond_v1_require_overall = kEarlyStopCondV1RequireOverall;
  base_cfg.early_stop_v2_llr_abs_threshold = kEarlyStopV2LlrAbsThreshold;
  base_cfg.early_stop_v2_max_unreliable_bits = kEarlyStopV2MaxUnreliableBits;
  base_cfg.early_stop_cond_v2_include_overall = kEarlyStopCondV2IncludeOverall;
  base_cfg.alpha_explicit = kAlphaExplicit;
  base_cfg.beta_explicit = kBetaExplicit;
  base_cfg.early_stop_action_sign_beta_explicit = kEarlyStopActionBetaExplicit;
  base_cfg.early_stop_action_residual_divisor =
      kEarlyStopActionResidualDivisor;
  base_cfg.early_stop_action_hard_llr_mag = kEarlyStopActionHardLlrMag;
  base_cfg.siso_active_list = kSisoActiveList;
  base_cfg.hiho_active_list = kHiHoActiveList;
  base_cfg.mux_group_g = kMuxGroupG;
  base_cfg.mux_scheduling_mode = kMuxSchedulingMode;
  base_cfg.mux_early_stop_priority_rule = kMuxPriorityRule;
  base_cfg.mux_enable_reconfig = kMuxEnableReconfig;
  base_cfg.mux_extra_bypass_edges = kMuxExtraBypassEdges;
  base_cfg.hybrid_enable = kHybridEnable;
  base_cfg.hybrid_enable_list = kHybridEnableList;
  base_cfg.hybrid_hard_llr_mag = kHybridHardLlrMag;
  base_cfg.hybrid_hard_llr_mag_list = kHybridHardLlrMagList;
  base_cfg.hybrid_classifier_mode = kHybridClassifierMode;
  base_cfg.hybrid_siso_backfill_mode = kHybridSisoBackfillMode;
  base_cfg.hybrid_normalize_soft_only = kHybridNormalizeSoftOnly;
  base_cfg.interleaver_name = kInterleaverName;
  base_cfg.decoder_name = kDecoderName;
  base_cfg.generate_random_bits = kGenerateRandomBits;
  base_cfg.normalize_known_prefix_tail = kNormalizeKnownPrefixTail;
  base_cfg.quant_clip_ratio = kQuantClipRatio;
  base_cfg.llr_bits = kLlrBits;
  base_cfg.debug_trace = newcode::Params::DebugTraceConfig{
      .enable = kDecoderTraceEnable,
      .log_read_mapping = kDecoderTraceLogRead,
      .log_write_mapping = kDecoderTraceLogWrite,
      .log_mismatch = kDecoderTraceLogMismatch,
      .log_chase_detail = false,
      .dump_chase_csv = false,
      .row = kDecoderTraceRow,
      .col = kDecoderTraceCol,
      .chase_decoder_row = -1,
      .chase_decoder_col = -1,
      .chase_tile_index = -1,
      .chase_invocation = -1,
      .chase_csv_dir = {},
      .chase_expected_bits = {},
      .chase_candidate_s1 = {},
      .chase_candidate_s3 = {},
      .chase_candidate_good = {},
      .chase_candidate_corrected_errors = {},
      .targets = {},
      .active_chase_entries = {},
  };

  std::ofstream null_file;
  io::DualWriter log(null_file);
  auto params = ofec_single::detail::build_params(base_cfg, log);
  if (!params.has_value()) {
    throw std::runtime_error(
        "failed to build params for two-stream shared sweep chunk");
  }
  params->NUM_INFO_BITS = kChunkNumInfoBits;

  newcode::PipelineConfig pipeline =
      ofec_single::detail::build_pipeline_config(base_cfg);
  pipeline.quiet = kQuietConsole;

  const int bitgen_seed_a = derive_seed(kBitgenSeedA, chunk_index, 0x13579bdfU);
  const int channel_seed_a =
      derive_seed(kChannelSeedA, chunk_index, 0x2468ace0U);
  const int bitgen_seed_b = derive_seed(kBitgenSeedB, chunk_index, 0x10293847U);
  const int channel_seed_b =
      derive_seed(kChannelSeedB, chunk_index, 0x56473829U);

  return newcode::two_stream_shared::Config{
      .params = *params,
      .pipeline = pipeline,
      .stream_a =
          {
              .label = std::string(kLabel) + "_A_chunk" +
                       std::to_string(chunk_index),
              .bitgen_seed = bitgen_seed_a,
              .channel_seed = channel_seed_a,
              .ebn0_db = ebn0_db,
          },
      .stream_b =
          {
              .label = std::string(kLabel) + "_B_chunk" +
                       std::to_string(chunk_index),
              .bitgen_seed = bitgen_seed_b,
              .channel_seed = channel_seed_b,
              .ebn0_db = ebn0_db,
          },
  };
}

ChunkResult run_chunk(float ebn0_db, std::size_t chunk_index) {
  auto config = build_chunk_config(ebn0_db, chunk_index);
  ChunkResult chunk;
  chunk.chunk_index = chunk_index;
  chunk.bitgen_seed_a = config.stream_a.bitgen_seed;
  chunk.channel_seed_a = config.stream_a.channel_seed;
  chunk.bitgen_seed_b = config.stream_b.bitgen_seed;
  chunk.channel_seed_b = config.stream_b.channel_seed;
  chunk.result = newcode::two_stream_shared::run_two_stream_shared(config);
  return chunk;
}

// ======== 点内聚合 / 停止条件 ========
double zero_error_upper_bound(std::size_t total_bits, double confidence_level) {
  if (total_bits == 0) {
    return std::numeric_limits<double>::infinity();
  }
  const double clamped = std::clamp(confidence_level, 1e-12, 1.0 - 1e-12);
  return -std::log(1.0 - clamped) / static_cast<double>(total_bits);
}

StopReason evaluate_stop_reason(const newcode::BerStats& post_fec) {
  if (post_fec.errors >= kTargetPostErrors) {
    return StopReason::TargetPostErrors;
  }
  if (post_fec.total >= kMaxPostFecTotalBits) {
    return StopReason::MaxPostFecBits;
  }
  if (kEnableZeroErrorUpperBound && post_fec.errors == 0 && post_fec.total > 0) {
    const double upper = zero_error_upper_bound(post_fec.total, kConfidenceLevel);
    if (upper <= kTargetBerUpperBound) {
      return StopReason::ZeroErrorUpperBound;
    }
  }
  return StopReason::NoData;
}

void update_upper_bound_fields(StreamAggregate* stream) {
  if (stream->post_fec.errors == 0 && stream->post_fec.total > 0) {
    stream->post_ber_upper_bound =
        zero_error_upper_bound(stream->post_fec.total, kConfidenceLevel);
    stream->post_ber_is_upper_bound =
        (stream->stop_reason == StopReason::ZeroErrorUpperBound);
  } else {
    stream->post_ber_upper_bound = std::numeric_limits<double>::quiet_NaN();
    stream->post_ber_is_upper_bound = false;
  }
}

void accumulate_stream(StreamAggregate* agg, const newcode::PipelineResult& result) {
  agg->pre_fec.errors += result.pre_fec.errors;
  agg->pre_fec.total += result.pre_fec.total;
  agg->pre_fec.ber =
      agg->pre_fec.total == 0
          ? 0.0
          : static_cast<double>(agg->pre_fec.errors) /
                static_cast<double>(agg->pre_fec.total);

  if (result.has_pre_fec_quantized_hard) {
    agg->has_pre_fec_quantized_hard = true;
    agg->pre_fec_quantized_hard.errors += result.pre_fec_quantized_hard.errors;
    agg->pre_fec_quantized_hard.total += result.pre_fec_quantized_hard.total;
    agg->pre_fec_quantized_hard.ber =
        agg->pre_fec_quantized_hard.total == 0
            ? 0.0
            : static_cast<double>(agg->pre_fec_quantized_hard.errors) /
                  static_cast<double>(agg->pre_fec_quantized_hard.total);
  }

  agg->post_fec.errors += result.post_fec.errors;
  agg->post_fec.total += result.post_fec.total;
  agg->post_fec.ber =
      agg->post_fec.total == 0
          ? 0.0
          : static_cast<double>(agg->post_fec.errors) /
                static_cast<double>(agg->post_fec.total);
}

void accumulate_chunk(AggregatedPointResult* point, const ChunkResult& chunk) {
  accumulate_stream(&point->stream_a, chunk.result.stream_a);
  accumulate_stream(&point->stream_b, chunk.result.stream_b);
  ++point->chunks_completed;
}

bool point_is_stopped(const AggregatedPointResult& point) {
  return point.stream_a.stop_reason != StopReason::NoData &&
         point.stream_b.stop_reason != StopReason::NoData;
}

void refresh_stop_state(AggregatedPointResult* point) {
  point->stream_a.stop_reason = evaluate_stop_reason(point->stream_a.post_fec);
  point->stream_b.stop_reason = evaluate_stop_reason(point->stream_b.post_fec);
  update_upper_bound_fields(&point->stream_a);
  update_upper_bound_fields(&point->stream_b);
}

void log_point_progress(const PointState& state,
                        bool force,
                        ofec_sweep::detail::DualOut& log) {
  const auto& point = state.point;
  if (!force && point.chunks_completed > 0 && point.chunks_completed % 5 != 0 &&
      point.chunks_completed > 3) {
    return;
  }

  std::ostringstream oss;
  oss << "[POINT] Eb/N0=" << point.ebn0_db
      << " chunks=" << point.chunks_completed
      << " inflight=" << state.inflight_chunks
      << " | A post=" << point.stream_a.post_fec.errors << "/"
      << point.stream_a.post_fec.total
      << " ber=" << std::scientific << point.stream_a.post_fec.ber
      << std::defaultfloat
      << " | B post=" << point.stream_b.post_fec.errors << "/"
      << point.stream_b.post_fec.total
      << " ber=" << std::scientific << point.stream_b.post_fec.ber
      << std::defaultfloat;
  if (state.started) {
    const auto elapsed =
        std::chrono::duration<double>(std::chrono::steady_clock::now() -
                                      state.start_time);
    oss << " elapsed=" << format_duration(elapsed);
  }
  log << oss.str() << "\n";
}

void maybe_mark_stop(PointState* state, ofec_sweep::detail::DualOut& log) {
  if (state->launch_stopped) {
    return;
  }
  refresh_stop_state(&state->point);
  if (!point_is_stopped(state->point)) {
    return;
  }
  state->launch_stopped = true;
  log_point_progress(*state, true, log);
  log << "[INFO] stop reason for Eb/N0=" << state->point.ebn0_db
      << ": A=" << stop_reason_to_string(state->point.stream_a.stop_reason)
      << ", B=" << stop_reason_to_string(state->point.stream_b.stop_reason)
      << "\n";
}

void finalize_point_if_ready(PointState* state,
                             ofec_sweep::detail::DualOut& log,
                             std::size_t* completed_count) {
  if (state->completed || !state->launch_stopped || state->inflight_chunks != 0) {
    return;
  }

  refresh_stop_state(&state->point);
  if (state->started) {
    state->point.elapsed_seconds =
        std::chrono::duration<double>(std::chrono::steady_clock::now() -
                                      state->start_time)
            .count();
  }

  log << "[SUMMARY] Eb/N0=" << state->point.ebn0_db
      << " | A Post-FEC BER=" << state->point.stream_a.post_fec.ber
      << " (errs=" << state->point.stream_a.post_fec.errors << "/"
      << state->point.stream_a.post_fec.total << ")"
      << " | B Post-FEC BER=" << state->point.stream_b.post_fec.ber
      << " (errs=" << state->point.stream_b.post_fec.errors << "/"
      << state->point.stream_b.post_fec.total << ")"
      << " | chunks=" << state->point.chunks_completed
      << " | stopA=" << stop_reason_to_string(state->point.stream_a.stop_reason)
      << " | stopB=" << stop_reason_to_string(state->point.stream_b.stop_reason)
      << " | elapsed=" << format_duration(
             std::chrono::duration<double>(state->point.elapsed_seconds))
      << "\n";

  state->completed = true;
  ++(*completed_count);
}

// ======== 全局 scheduler ========
std::vector<AggregatedPointResult> run_low_ber_points_global(
    const std::vector<float>& ebn0_values,
    unsigned total_workers,
    unsigned base_inflight_per_point,
    unsigned dynamic_inflight_cap,
    ofec_sweep::detail::DualOut& log) {
  std::vector<PointState> states;
  states.reserve(ebn0_values.size());
  for (float ebn0_db : ebn0_values) {
    PointState state;
    state.point.ebn0_db = ebn0_db;
    states.push_back(std::move(state));
  }

  std::vector<ActiveChunk> active;
  active.reserve(total_workers);
  std::size_t completed_count = 0;
  std::size_t next_launch_point = 0;
  unsigned last_logged_effective_limit = 0;

  auto count_launchable_points = [&]() {
    std::size_t count = 0;
    for (const auto& state : states) {
      if (!state.completed && !state.launch_stopped) {
        ++count;
      }
    }
    return count;
  };

  auto effective_inflight_limit = [&]() {
    unsigned limit = base_inflight_per_point;
    if (kEnableDynamicInflightPerPoint) {
      const std::size_t active_points = count_launchable_points();
      if (active_points > 0) {
        const unsigned dynamic_limit =
            static_cast<unsigned>((total_workers + active_points - 1) /
                                  active_points);
        limit = std::max(limit, dynamic_limit);
      }
    }
    if (dynamic_inflight_cap > 0) {
      limit = std::min(limit, dynamic_inflight_cap);
    }
    return std::max(1u, std::min(total_workers, limit));
  };

  auto can_launch_for = [&](std::size_t point_index) {
    const PointState& state = states[point_index];
    return !state.completed && !state.launch_stopped &&
           state.inflight_chunks < effective_inflight_limit();
  };

  auto launch_one = [&](std::size_t point_index) {
    PointState& state = states[point_index];
    if (!state.started) {
      state.started = true;
      state.start_time = std::chrono::steady_clock::now();
      log << "\n[RUN] Eb/N0=" << state.point.ebn0_db << "\n";
    }

    const std::size_t chunk_index = state.next_chunk_index++;
    const float ebn0_db = state.point.ebn0_db;
    active.push_back(ActiveChunk{
        point_index,
        chunk_index,
        std::async(std::launch::async, [ebn0_db, chunk_index]() {
          return run_chunk(ebn0_db, chunk_index);
        })});
    ++state.inflight_chunks;
  };

  auto fill_workers = [&]() {
    bool launched_any = false;
    while (active.size() < total_workers && completed_count < states.size()) {
      const unsigned current_limit = effective_inflight_limit();
      if (current_limit != last_logged_effective_limit) {
        last_logged_effective_limit = current_limit;
        log << "[INFO] effective max inflight chunks per point = "
            << current_limit
            << " (launchable_points=" << count_launchable_points()
            << ")\n";
      }

      bool launched_this_round = false;
      for (std::size_t offset = 0; offset < states.size(); ++offset) {
        const std::size_t point_index =
            (next_launch_point + offset) % states.size();
        if (!can_launch_for(point_index)) {
          continue;
        }

        launch_one(point_index);
        next_launch_point = (point_index + 1) % states.size();
        launched_this_round = true;
        launched_any = true;
        break;
      }

      if (!launched_this_round) {
        break;
      }
    }
    return launched_any;
  };

  fill_workers();

  while (completed_count < states.size()) {
    bool consumed = false;
    for (auto it = active.begin(); it != active.end(); ++it) {
      if (it->future.wait_for(std::chrono::milliseconds(0)) !=
          std::future_status::ready) {
        continue;
      }

      const std::size_t point_index = it->point_index;
      ChunkResult chunk = it->future.get();
      active.erase(it);

      PointState& state = states[point_index];
      if (state.inflight_chunks > 0) {
        --state.inflight_chunks;
      }
      accumulate_chunk(&state.point, chunk);
      refresh_stop_state(&state.point);
      log_point_progress(state, false, log);
      maybe_mark_stop(&state, log);
      finalize_point_if_ready(&state, log, &completed_count);
      fill_workers();

      consumed = true;
      break;
    }

    if (consumed) {
      continue;
    }

    if (active.empty()) {
      log << "[ERROR] global scheduler has no active chunks before all points "
             "completed\n";
      break;
    }
    std::this_thread::sleep_for(std::chrono::milliseconds(20));
  }

  std::vector<AggregatedPointResult> results;
  results.reserve(states.size());
  for (const auto& state : states) {
    results.push_back(state.point);
  }
  return results;
}

// ======== CSV writer ========
void ensure_csv_header_two_stream_shared_sweep(const std::string& csv_path) {
  std::ifstream fin(csv_path);
  if (fin.good() && fin.peek() != std::ifstream::traits_type::eof()) {
    return;
  }

  std::ofstream fout(csv_path, std::ios::out | std::ios::app);
  fout << "timestamp,run_id,ebn0_db,"
          "pre_ber_a,pre_errs_a,pre_total_a,pre_quant_ber_a,pre_quant_errs_a,pre_quant_total_a,post_ber_a,post_errs_a,post_total_a,"
          "pre_ber_b,pre_errs_b,pre_total_b,pre_quant_ber_b,pre_quant_errs_b,pre_quant_total_b,post_ber_b,post_errs_b,post_total_b,"
          "chunk_num_info_bits,chunks_completed,target_post_errors,max_post_fec_total_bits,confidence_level,"
          "stop_reason_a,post_ber_is_upper_bound_a,post_ber_upper_bound_a,"
          "stop_reason_b,post_ber_is_upper_bound_b,post_ber_upper_bound_b,"
          "elapsed_seconds,"
          "bitgen_seed_a_base,channel_seed_a_base,bitgen_seed_b_base,channel_seed_b_base,"
          "alpha_list,beta_list,siso_active_list,"
          "chase_L,chase_n_test,chase_topk_keep,chase_group_minima_bits,"
          "mux_group_g,mux_scheduling_mode,mux_early_stop_priority_rule,mux_enable_reconfig,"
          "hybrid_enable,hybrid_enable_list,hybrid_classifier_mode,hybrid_siso_backfill_mode,hybrid_normalize_soft_only,"
          "early_stop_enable,early_stop_enable_list,early_stop_condition_mode,early_stop_action_mode,early_stop_bind_group_size,"
          "early_stop_cond_v1_require_bch,early_stop_cond_v1_require_overall,"
          "early_stop_v2_llr_abs_threshold,early_stop_v2_max_unreliable_bits,early_stop_cond_v2_include_overall\n";
}

void write_csv_row_two_stream_shared_sweep(std::ostream& csv,
                                           const std::string& timestamp,
                                           const std::string& run_id,
                                           const AggregatedPointResult& point) {
  csv << timestamp << ','
      << run_id << ','
      << point.ebn0_db << ','
      << point.stream_a.pre_fec.ber << ','
      << point.stream_a.pre_fec.errors << ','
      << point.stream_a.pre_fec.total << ',';
  if (point.stream_a.has_pre_fec_quantized_hard) {
    csv << point.stream_a.pre_fec_quantized_hard.ber << ','
        << point.stream_a.pre_fec_quantized_hard.errors << ','
        << point.stream_a.pre_fec_quantized_hard.total << ',';
  } else {
    csv << ",,,";
  }
  csv << point.stream_a.post_fec.ber << ','
      << point.stream_a.post_fec.errors << ','
      << point.stream_a.post_fec.total << ','
      << point.stream_b.pre_fec.ber << ','
      << point.stream_b.pre_fec.errors << ','
      << point.stream_b.pre_fec.total << ',';
  if (point.stream_b.has_pre_fec_quantized_hard) {
    csv << point.stream_b.pre_fec_quantized_hard.ber << ','
        << point.stream_b.pre_fec_quantized_hard.errors << ','
        << point.stream_b.pre_fec_quantized_hard.total << ',';
  } else {
    csv << ",,,";
  }
  csv << point.stream_b.post_fec.ber << ','
      << point.stream_b.post_fec.errors << ','
      << point.stream_b.post_fec.total << ','
      << kChunkNumInfoBits << ','
      << point.chunks_completed << ','
      << kTargetPostErrors << ','
      << kMaxPostFecTotalBits << ','
      << kConfidenceLevel << ','
      << stop_reason_to_string(point.stream_a.stop_reason) << ','
      << (point.stream_a.post_ber_is_upper_bound ? 1 : 0) << ',';
  if (std::isnan(point.stream_a.post_ber_upper_bound)) {
    csv << ',';
  } else {
    csv << point.stream_a.post_ber_upper_bound << ',';
  }
  csv << stop_reason_to_string(point.stream_b.stop_reason) << ','
      << (point.stream_b.post_ber_is_upper_bound ? 1 : 0) << ',';
  if (std::isnan(point.stream_b.post_ber_upper_bound)) {
    csv << ',';
  } else {
    csv << point.stream_b.post_ber_upper_bound << ',';
  }
  csv << point.elapsed_seconds << ','
      << kBitgenSeedA << ','
      << kChannelSeedA << ','
      << kBitgenSeedB << ','
      << kChannelSeedB << ','
      << '"' << join_vec(kAlphaExplicit, '|', 6) << "\","
      << '"' << join_vec(kBetaExplicit, '|', 6) << "\","
      << '"' << join_int_vec(kSisoActiveList, '|') << "\","
      << kChaseL << ','
      << (kChaseNTestOverride < 0 ? (1 << kChaseL) : kChaseNTestOverride)
      << ','
      << kChaseTopkKeep << ','
      << kChaseGroupMinimaBits << ','
      << kMuxGroupG << ','
      << kMuxSchedulingMode << ','
      << kMuxPriorityRule << ','
      << (kMuxEnableReconfig ? 1 : 0) << ','
      << (kHybridEnable ? 1 : 0) << ','
      << '"' << join_int_vec(kHybridEnableList, '|') << "\","
      << hybrid_classifier_mode_name(kHybridClassifierMode) << ','
      << hybrid_siso_backfill_mode_name(kHybridSisoBackfillMode) << ','
      << (kHybridNormalizeSoftOnly ? 1 : 0) << ','
      << (kEnableEarlyStop ? 1 : 0) << ','
      << '"' << join_int_vec(kEarlyStopEnableList, '|') << "\","
      << kEarlyStopConditionMode << ','
      << kEarlyStopActionMode << ','
      << kEarlyStopBindGroupSize << ','
      << (kEarlyStopCondV1RequireBch ? 1 : 0) << ','
      << (kEarlyStopCondV1RequireOverall ? 1 : 0) << ','
      << kEarlyStopV2LlrAbsThreshold << ','
      << kEarlyStopV2MaxUnreliableBits << ','
      << (kEarlyStopCondV2IncludeOverall ? 1 : 0) << '\n';
}

void write_results_csv_two_stream_shared_sweep(
    std::ostream& csv,
    const std::string& run_id,
    const std::vector<AggregatedPointResult>& results) {
  for (const auto& point : results) {
    write_csv_row_two_stream_shared_sweep(csv, utils::now_stamp(), run_id,
                                          point);
  }
}

// ======== main 辅助函数 ========
struct SweepOutputFiles {
  std::filesystem::path data_dir = "data";
  std::string run_id;
  std::string log_path;
  std::string csv_path;
};

struct WorkerLaunchConfig {
  unsigned available_workers = 0;
  unsigned total_workers = 0;
  unsigned base_inflight_per_point = 0;
  unsigned dynamic_inflight_cap = 0;
};

std::vector<float> build_ebn0_values() {
  ofec_sweep::SweepParameterConfig config;
  config.ebn0_start = kEbN0Start;
  config.ebn0_end = kEbN0End;
  config.ebn0_points = kEbN0Points;
  return ofec_sweep::detail::build_ebn0_values(config);
}

SweepOutputFiles build_output_files() {
  SweepOutputFiles files;
  io::ensure_dir(files.data_dir);
  files.run_id = utils::now_stamp();
  files.log_path =
      (files.data_dir / ("run_" + files.run_id + "_two_stream_shared_sweep.log"))
          .string();
  files.csv_path =
      (files.data_dir /
       ("ofec_two_stream_shared_sweep_results_" + files.run_id + ".csv"))
          .string();
  return files;
}

std::ofstream open_csv_writer_two_stream_shared_sweep(
    const std::string& csv_path) {
  ensure_csv_header_two_stream_shared_sweep(csv_path);
  std::ofstream csv(csv_path, std::ios::out | std::ios::app);
  csv.setf(std::ios::fixed);
  csv << std::setprecision(10);
  return csv;
}

WorkerLaunchConfig resolve_worker_launch_config() {
  ofec_sweep::SweepParameterConfig worker_cfg;
  WorkerLaunchConfig config;
  config.available_workers =
      ofec_sweep::detail::resolve_worker_count(worker_cfg);
  config.total_workers =
      (kMaxTotalWorkers == 0)
          ? config.available_workers
          : std::max(1u,
                     std::min(config.available_workers, kMaxTotalWorkers));
  config.base_inflight_per_point =
      (kMaxInflightChunksPerPoint == 0)
          ? config.total_workers
          : std::max(1u,
                     std::min(config.total_workers, kMaxInflightChunksPerPoint));
  config.dynamic_inflight_cap =
      (kMaxDynamicInflightChunksPerPoint == 0)
          ? config.total_workers
          : std::max(1u, std::min(config.total_workers,
                                  kMaxDynamicInflightChunksPerPoint));
  return config;
}

void log_sweep_plan(const std::vector<float>& ebn0_values,
                    const WorkerLaunchConfig& workers,
                    ofec_sweep::detail::DualOut& out) {
  out << "[INFO] total Eb/N0 points = " << ebn0_values.size() << "\n";
  out << "[INFO] total workers = " << workers.total_workers << " (available="
      << workers.available_workers << ")\n";
  out << "[INFO] base max inflight chunks per point = "
      << workers.base_inflight_per_point << "\n";
  out << "[INFO] dynamic inflight per point = "
      << (kEnableDynamicInflightPerPoint ? "enabled" : "disabled")
      << " (cap=" << workers.dynamic_inflight_cap << ")\n";
  out << "[INFO] chunk_num_info_bits = " << kChunkNumInfoBits << "\n";
  out << "[INFO] target_post_errors = " << kTargetPostErrors << "\n";
  out << "[INFO] max_post_fec_total_bits = " << kMaxPostFecTotalBits << "\n";
  if (kEnableZeroErrorUpperBound) {
    out << "[INFO] zero-error upper-bound stop enabled, target="
        << std::scientific << kTargetBerUpperBound << std::defaultfloat
        << ", confidence=" << kConfidenceLevel << "\n";
  }
}

void sort_results_by_ebn0(std::vector<AggregatedPointResult>* results) {
  std::sort(results->begin(), results->end(),
            [](const AggregatedPointResult& lhs,
               const AggregatedPointResult& rhs) {
              return lhs.ebn0_db < rhs.ebn0_db;
            });
}

}  // namespace

int main() {
  const SweepOutputFiles files = build_output_files();
  ofec_sweep::detail::DualOut out(std::cout, files.log_path, !kQuietConsole);
  std::ofstream csv = open_csv_writer_two_stream_shared_sweep(files.csv_path);

  const std::vector<float> ebn0_values = build_ebn0_values();
  if (ebn0_values.empty()) {
    std::cerr << "[ERROR] no Eb/N0 points generated for two-stream shared sweep\n";
    return 1;
  }

  const WorkerLaunchConfig workers = resolve_worker_launch_config();
  log_sweep_plan(ebn0_values, workers, out);

  auto results = run_low_ber_points_global(ebn0_values,
                                           workers.total_workers,
                                           workers.base_inflight_per_point,
                                           workers.dynamic_inflight_cap,
                                           out);
  sort_results_by_ebn0(&results);
  write_results_csv_two_stream_shared_sweep(csv, files.run_id, results);
  csv.flush();
  return 0;
}
