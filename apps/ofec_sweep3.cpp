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
#include "newcode/ofec_sweep_runner.hpp"
#include "newcode/utils/now_stamp.hpp"
#include "ofec_sweep_detail.hpp"

namespace {

// ======== 用户可调参数区域 ========

// 发端参数：交织 / 比特源 / 调制入口
static constexpr const char* kInterleaverName = "identity"; // 交织器名称，identity 表示不交织
static constexpr unsigned kBitsPerSymbol = 1;               // 每个调制符号携带的比特数：1=BPSK，偶数=QAM
static constexpr bool kGenerateRandomBits = true;           // true=发送随机信息比特，false=发送全 0 比特
static constexpr int kBitgenSeed = 20260319;                // 基础比特种子；每个 chunk 会在此基础上派生

// 信道参数：噪声强度 / 信道随机性
static constexpr float kEbN0Start = 3.00f;   // 扫描起始 Eb/N0（dB）
static constexpr float kEbN0End = 3.20f;     // 扫描结束 Eb/N0（dB）
static constexpr int kEbN0Points = 20;        // Eb/N0 采样点数；含首尾端点
static constexpr int kChannelSeed = 3192026; // 基础信道种子；每个 chunk 会在此基础上派生

// 量化参数：只影响 LLR 量化口径
static constexpr std::size_t kLlrBits = 6;      // LLR 位宽：16=浮点，2~15=qfloat
static constexpr float kQuantClipRatio = 0.5f;  // 动态 clip 比例，0 表示禁用自适应 clip

// 低 BER 聚合参数：每点多 chunk 聚合，直到达到停止条件
static constexpr std::size_t kTilesPerWindow = 6;               // 本 app 使用的 TILES_PER_WIN；需与显式 alpha/beta 列表长度一致
static constexpr std::size_t kChunkNumInfoBits = 8 * 132 * 16 * 111; // 每个 Monte Carlo chunk 的输入信息比特数
static constexpr std::size_t kTargetPostErrors = 50;            // 单个 Eb/N0 点累计到这么多 post-FEC 错误后即可停止
static constexpr std::size_t kMaxPostFecTotalBits = 2e8;        // 单个 Eb/N0 点允许累计比较的最大 post-FEC 比特数
static constexpr unsigned kMaxTotalWorkers = 0;                 // 全局同时运行 chunk 数上限；0=自动使用 NTHREADS/机器可用 worker
static constexpr unsigned kMaxInflightChunksPerPoint = 16;      // 单个 Eb/N0 点的基础挂起 chunk 上限；0=不额外限制
static constexpr bool kEnableDynamicInflightPerPoint = true;    // true=剩余 Eb/N0 点变少时动态提高单点挂起上限
static constexpr unsigned kMaxDynamicInflightChunksPerPoint = 0; // 动态单点挂起上限封顶；0=不封顶，最多到全局 worker
static constexpr bool kEnableZeroErrorUpperBound = true;        // true=零错时使用上置信界提前停止
static constexpr double kTargetBerUpperBound = 1e-8;            // 零错上界目标：若上界已低于此值则提前停止
static constexpr double kConfidenceLevel = 0.95;                // 零错上界使用的置信水平，例如 0.95 表示 95%

// 解码参数：decoder 选择、Chase 参数、alpha/beta、MUX 调度
static constexpr const char* kDecoderName = "chase_baseline";        // 默认 decoder 名称
static const std::vector<const char*> kDecoderNameCandidates = {};   // decoder 扫描候选；空表示不扫 decoder 维度
static constexpr bool kNormalizeExtrinsic = false;                   // 是否对 decoder 输出 extrinsic 做归一化
static constexpr bool kNormalizeKnownPrefixTail = false;             // 是否对 known-prefix 后的尾部 LLR 做归一化
static const std::vector<int> kChaseLCandidates = {6};              // Chase L 扫描候选
static constexpr int kChaseNTest = 64;                              // 默认 Chase NTEST；若不单独扫则等于实际使用值
static const std::vector<int> kChaseNTestCandidates = {};           // Chase NTEST 扫描候选；空表示不单独扫描
static constexpr int kChaseTopkKeep = 8;                            // top-k/pruned decoder 保留候选数
static const std::vector<int> kChaseTopkKeepCandidates = {};        // top-k 保留数扫描候选
static constexpr int kChaseGroupMinimaBits = 4;                     // group-minima decoder 的分组 bit 数
static const std::vector<int> kChaseGroupMinimaBitsCandidates = {4}; // group-minima 分组 bit 数扫描候选
static const std::vector<int> kSisoActiveList = {32, 32, 32, 32, 32, 32};   // 每个 tile 的 SISO 预算
static constexpr int kMuxGroupG = 1;                                // MUX 分组粒度；1=全局池化
static constexpr int kMuxSchedulingMode = 0;                        // MUX 调度模式：0=legacy，1=按 early-stop 细节排序
static const std::vector<int> kMuxSchedulingModeCandidates = {};    // MUX 调度模式扫描候选
static constexpr int kMuxPriorityRule = 0;                          // 新 MUX 的优先级规则：0=更差优先，1=更接近通过优先
static const std::vector<int> kMuxPriorityRuleCandidates = {};      // MUX 优先级规则扫描候选
static constexpr bool kMuxEnableReconfig = false;                   // 是否启用 reconfig MUX 调度
static constexpr int kMuxBypassScheme = 3;                          // 旁路边方案编号：仅在 reconfig 打开时真正参与调度
static constexpr bool kHybridEnable = false;                        // true=启用方案三软硬混合前置分流
static const std::vector<int> kHybridEnableList = {};               // 按 tile 覆盖 hybrid 开关：空=沿用 kHybridEnable
static constexpr newcode::HybridClassifierMode kHybridClassifierMode =
    newcode::HybridClassifierMode::LegacyHardDecode;                // LegacyHardDecode / RepoFastClassifier / FriendS1S3Classifier / FriendS1S3WithS0Classifier
static constexpr bool kHybridNormalizeSoftOnly = false;             // true=只归一化 soft rows，false=保持兼容行为
static const std::vector<ofec_sweep::ExplicitAlphaBetaPattern> kExplicitAlphaBetaSets = {
    {"custom_label",                                             // 该组显式 alpha/beta 的标签，会进入场景名
     {0.428571,0.447738,0.482782,0.528162,0.581902,0.642857},               // 每个 tile 的 alpha 显式列表
     {2.857143,6.179301,12.253626,20.119585,29.434408,40.000000},            // 每个 tile 的普通 beta 显式列表
     {99,99,99,99,99,99}},                                                        // 每个 tile 的 early-stop 专用 beta 显式列表；空表示后续按默认规则回退
};

// 早停参数：总开关 -> 条件 -> 条件细参 -> 动作 -> 动作细参
static constexpr bool kEnableEarlyStop = false;                      // 早停总开关
static const std::vector<int> kEarlyStopEnableList = {};            // 按 tile 覆盖早停开关；空表示全部沿用总开关
static constexpr int kEarlyStopConditionMode = 1;                   // 早停条件模式：1=v1，2=v2
static const std::vector<int> kEarlyStopConditionCandidates = {};   // 早停条件模式扫描候选
static constexpr int kEarlyStopActionMode = 1;                      // 早停动作模式：1~6
static constexpr int kEarlyStopBindGroupSize = 1;                   // 条件1的组绑定大小；1=逐 row，4=四个绑定
static const std::vector<int> kEarlyStopActionCandidates = {};      // 早停动作模式扫描候选
static constexpr bool kEarlyStopCondV1RequireBch = true;            // 条件1是否要求 BCH syndrome 全 0
static const std::vector<bool> kEarlyStopCondV1RequireBchCandidates = {}; // 条件1 BCH 开关扫描候选
static constexpr bool kEarlyStopCondV1RequireOverall = true;        // 条件1是否要求 overall parity 通过
static const std::vector<bool> kEarlyStopCondV1RequireOverallCandidates = {}; // 条件1 overall 开关扫描候选
static constexpr float kEarlyStopV2LlrAbsThreshold = 0.5f;          // 条件2判定“不可靠位”的 |LLR| 阈值
static const std::vector<float> kEarlyStopV2LlrAbsThresholdCandidates = {}; // 条件2阈值扫描候选
static constexpr int kEarlyStopV2MaxUnreliableBits = 8;             // 条件2允许的不可靠 bit 数上限
static const std::vector<int> kEarlyStopV2MaxUnreliableBitsCandidates = {}; // 条件2 bit 数上限扫描候选
static constexpr bool kEarlyStopCondV2IncludeOverall = true;        // 条件2是否把 overall parity bit 一起纳入统计
static constexpr float kEarlyStopActionResidualDivisor = 0.4f;      // 动作2里 residual 的除数
static constexpr float kEarlyStopActionHardLlrMag = 1.0f;           // 动作3里输出的固定 |LLR| 幅度
static const std::vector<float> kEarlyStopActionBetaStartCandidates = {}; // early-stop 专用 beta 起点扫描候选
static const std::vector<float> kEarlyStopActionBetaStepCandidates = {};  // early-stop 专用 beta 步长扫描候选
static const std::vector<float> kEarlyStopActionHardLlrMagCandidates = {}; // 动作3 hard LLR 幅度扫描候选

// Debug 参数：控制日志与 decoder trace
static constexpr bool kQuietConsole = false;           // true=减少控制台输出
static constexpr bool kDecoderTraceEnable = false;     // decoder trace 总开关
static constexpr long kDecoderTraceRow = -1;           // 追踪目标的全局 row；-1 表示不指定
static constexpr long kDecoderTraceCol = -1;           // 追踪目标的全局 col；-1 表示不指定
static constexpr bool kDecoderTraceLogRead = false;    // 是否打印 tile 读取映射
static constexpr bool kDecoderTraceLogWrite = false;   // 是否打印 tile 写回映射
static constexpr bool kDecoderTraceLogMismatch = false; // 是否打印同坐标写回不一致告警

// ==================================

enum class StopReason {
  TargetPostErrors,
  MaxPostFecBits,
  ZeroErrorUpperBound,
  NoData
};

struct ChunkResult {
  std::size_t chunk_index = 0;
  int bitgen_seed = 0;
  int channel_seed = 0;
  newcode::PipelineResult result;
};

struct AggregatedPointResult {
  ofec_sweep::detail::SweepScenario scenario;
  newcode::BerStats pre_fec{};
  newcode::BerStats pre_fec_quantized_hard{};
  newcode::BerStats post_fec{};
  bool has_pre_fec_quantized_hard = false;
  std::size_t chunks_completed = 0;
  bool post_ber_is_upper_bound = false;
  double post_ber_upper_bound = std::numeric_limits<double>::quiet_NaN();
  StopReason stop_reason = StopReason::NoData;
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

void ensure_csv_header_sweep3(const std::string& csv_path) {
  std::ifstream fin(csv_path);
  if (fin.good() && fin.peek() != std::ifstream::traits_type::eof()) {
    return;
  }

  std::ofstream fout(csv_path, std::ios::out | std::ios::app);
  fout << "timestamp,run_id,scenario,decoder_name,ebn0_db,"
          "pre_ber,pre_errs,pre_total,post_ber,post_errs,post_total,"
          "chunk_num_info_bits,chunks_completed,target_post_errors,max_post_fec_total_bits,"
          "confidence_level,post_ber_is_upper_bound,post_ber_upper_bound,stop_reason,elapsed_seconds,"
          "bitgen_seed_base,channel_seed_base,"
          "alpha_list,beta_list,early_stop_beta_list,"
          "chase_L,chase_n_test,chase_topk_keep,chase_group_minima_bits,"
          "mux_group_g,mux_scheduling_mode,mux_early_stop_priority_rule,mux_bypass_scheme,"
          "hybrid_enable,hybrid_enable_list,hybrid_classifier_mode,hybrid_normalize_soft_only,"
          "early_stop_condition_mode,early_stop_action_mode,early_stop_bind_group_size,"
          "early_stop_cond_v1_require_bch,early_stop_cond_v1_require_overall,"
          "early_stop_v2_llr_abs_threshold,early_stop_v2_max_unreliable_bits,early_stop_cond_v2_include_overall\n";
}

void write_csv_row_sweep3(std::ostream& csv,
                          const std::string& timestamp,
                          const std::string& run_id,
                          const AggregatedPointResult& point) {
  csv << timestamp << ','
      << run_id << ','
      << point.scenario.name << ','
      << point.scenario.decoder_name << ','
      << point.scenario.ebn0_db << ','
      << point.pre_fec.ber << ','
      << point.pre_fec.errors << ','
      << point.pre_fec.total << ','
      << point.post_fec.ber << ','
      << point.post_fec.errors << ','
      << point.post_fec.total << ','
      << kChunkNumInfoBits << ','
      << point.chunks_completed << ','
      << kTargetPostErrors << ','
      << kMaxPostFecTotalBits << ','
      << kConfidenceLevel << ','
      << (point.post_ber_is_upper_bound ? 1 : 0) << ',';
  if (std::isnan(point.post_ber_upper_bound)) {
    csv << ',';
  } else {
    csv << point.post_ber_upper_bound << ',';
  }
  csv << stop_reason_to_string(point.stop_reason) << ','
      << point.elapsed_seconds << ','
      << point.scenario.bitgen_seed << ','
      << point.scenario.channel_seed << ','
      << '"' << join_vec(point.scenario.alpha_list, '|', 6) << "\","
      << '"' << join_vec(point.scenario.beta_list, '|', 6) << "\","
      << '"' << join_vec(point.scenario.early_stop_action_sign_beta_list, '|', 6) << "\","
      << point.scenario.chase_L << ','
      << point.scenario.chase_n_test << ','
      << point.scenario.chase_topk_keep << ','
      << point.scenario.chase_group_minima_bits << ','
      << point.scenario.mux_group_g << ','
      << point.scenario.mux_scheduling_mode << ','
      << point.scenario.mux_early_stop_priority_rule << ','
      << point.scenario.mux_bypass_scheme << ','
      << (kHybridEnable ? 1 : 0) << ','
      << '"' << join_int_vec(kHybridEnableList, '|') << "\","
      << hybrid_classifier_mode_name(kHybridClassifierMode) << ','
      << (kHybridNormalizeSoftOnly ? 1 : 0) << ','
      << point.scenario.early_stop_condition_mode << ','
      << point.scenario.early_stop_action_mode << ','
      << point.scenario.early_stop_bind_group_size << ','
      << (point.scenario.early_stop_cond_v1_require_bch ? 1 : 0) << ','
      << (point.scenario.early_stop_cond_v1_require_overall ? 1 : 0) << ','
      << point.scenario.early_stop_v2_llr_abs_threshold << ','
      << point.scenario.early_stop_v2_max_unreliable_bits << ','
      << (point.scenario.early_stop_cond_v2_include_overall ? 1 : 0) << '\n';
}

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

double zero_error_upper_bound(std::size_t total_bits, double confidence_level) {
  if (total_bits == 0) {
    return std::numeric_limits<double>::infinity();
  }
  const double clamped = std::clamp(confidence_level, 1e-12, 1.0 - 1e-12);
  return -std::log(1.0 - clamped) / static_cast<double>(total_bits);
}

StopReason evaluate_stop_reason(const AggregatedPointResult& point) {
  if (point.post_fec.errors >= kTargetPostErrors) {
    return StopReason::TargetPostErrors;
  }
  if (point.post_fec.total >= kMaxPostFecTotalBits) {
    return StopReason::MaxPostFecBits;
  }
  if (kEnableZeroErrorUpperBound && point.post_fec.errors == 0 &&
      point.post_fec.total > 0) {
    const double upper =
        zero_error_upper_bound(point.post_fec.total, kConfidenceLevel);
    if (upper <= kTargetBerUpperBound) {
      return StopReason::ZeroErrorUpperBound;
    }
  }
  return StopReason::NoData;
}

void update_upper_bound_fields(AggregatedPointResult& point) {
  if (point.post_fec.errors == 0 && point.post_fec.total > 0) {
    point.post_ber_upper_bound =
        zero_error_upper_bound(point.post_fec.total, kConfidenceLevel);
    point.post_ber_is_upper_bound =
        (point.stop_reason == StopReason::ZeroErrorUpperBound);
  } else {
    point.post_ber_upper_bound = std::numeric_limits<double>::quiet_NaN();
    point.post_ber_is_upper_bound = false;
  }
}

void accumulate_chunk(AggregatedPointResult& point, const ChunkResult& chunk) {
  point.pre_fec.errors += chunk.result.pre_fec.errors;
  point.pre_fec.total += chunk.result.pre_fec.total;
  point.pre_fec.ber =
      point.pre_fec.total == 0
          ? 0.0
          : static_cast<double>(point.pre_fec.errors) /
                static_cast<double>(point.pre_fec.total);

  if (chunk.result.has_pre_fec_quantized_hard) {
    point.has_pre_fec_quantized_hard = true;
    point.pre_fec_quantized_hard.errors +=
        chunk.result.pre_fec_quantized_hard.errors;
    point.pre_fec_quantized_hard.total +=
        chunk.result.pre_fec_quantized_hard.total;
    point.pre_fec_quantized_hard.ber =
        point.pre_fec_quantized_hard.total == 0
            ? 0.0
            : static_cast<double>(point.pre_fec_quantized_hard.errors) /
                  static_cast<double>(point.pre_fec_quantized_hard.total);
  }

  point.post_fec.errors += chunk.result.post_fec.errors;
  point.post_fec.total += chunk.result.post_fec.total;
  point.post_fec.ber =
      point.post_fec.total == 0
          ? 0.0
          : static_cast<double>(point.post_fec.errors) /
                static_cast<double>(point.post_fec.total);

  ++point.chunks_completed;
}

newcode::Params make_params_for_scenario(const ofec_sweep::detail::SweepScenario& scenario,
                                         const ofec_sweep::SweepParameterConfig& config) {
  newcode::Params params = config.base_params;
  params.BITGEN_RANDOM_BITS = config.generate_random_bits;
  params.NORMALIZE_KNOWN_PREFIX_TAIL = config.normalize_known_prefix_tail;
  params.ALPHA_LIST = scenario.alpha_list;
  params.beta_list = scenario.beta_list;
  params.EARLY_STOP_ACTION_SIGN_BETA_LIST =
      scenario.early_stop_action_sign_beta_list;
  if (!params.ALPHA_LIST.empty()) {
    params.ALPHA = params.ALPHA_LIST.front();
  }
  if (!params.beta_list.empty()) {
    params.beta = params.beta_list.front();
  }
  if (!params.EARLY_STOP_ACTION_SIGN_BETA_LIST.empty()) {
    params.EARLY_STOP_ACTION_SIGN_BETA =
        params.EARLY_STOP_ACTION_SIGN_BETA_LIST.front();
  }
  params.CHASE_L = scenario.chase_L;
  params.CHASE_NTEST = scenario.chase_n_test;
  params.CHASE_TOPK_KEEP = scenario.chase_topk_keep;
  params.CHASE_GROUP_MINIMA_BITS = scenario.chase_group_minima_bits;
  params.MUX_GROUP_G = scenario.mux_group_g;
  params.MUX_SCHEDULING_MODE = scenario.mux_scheduling_mode;
  params.MUX_EARLY_STOP_PRIORITY_RULE =
      scenario.mux_early_stop_priority_rule;
  params.EARLY_STOP_CONDITION_MODE = scenario.early_stop_condition_mode;
  params.EARLY_STOP_ACTION_MODE = scenario.early_stop_action_mode;
  params.EARLY_STOP_BIND_GROUP_SIZE = scenario.early_stop_bind_group_size;
  params.EARLY_STOP_COND_V1_REQUIRE_BCH =
      scenario.early_stop_cond_v1_require_bch;
  params.EARLY_STOP_COND_V1_REQUIRE_OVERALL =
      scenario.early_stop_cond_v1_require_overall;
  params.EARLY_STOP_V2_LLR_ABS_THRESHOLD =
      scenario.early_stop_v2_llr_abs_threshold;
  params.EARLY_STOP_V2_MAX_UNRELIABLE_BITS =
      scenario.early_stop_v2_max_unreliable_bits;
  params.EARLY_STOP_COND_V2_INCLUDE_OVERALL =
      scenario.early_stop_cond_v2_include_overall;
  params.EARLY_STOP_ACTION_HARD_LLR_MAG =
      scenario.early_stop_action_hard_llr_mag;
  return params;
}

ChunkResult run_chunk(const ofec_sweep::detail::SweepScenario& scenario,
                      const ofec_sweep::SweepParameterConfig& config,
                      std::size_t chunk_index) {
  ChunkResult chunk;
  chunk.chunk_index = chunk_index;
  chunk.bitgen_seed = derive_seed(scenario.bitgen_seed, chunk_index, 0x13579bdfU);
  chunk.channel_seed = derive_seed(scenario.channel_seed, chunk_index, 0x2468ace0U);

  newcode::Params params = make_params_for_scenario(scenario, config);
  params.NUM_INFO_BITS = kChunkNumInfoBits;
  params.BITGEN_SEED = chunk.bitgen_seed;
  params.CHANNEL_SEED = chunk.channel_seed;

  newcode::PipelineConfig pipeline_cfg = ofec_sweep::detail::make_pipeline_config(config);
  pipeline_cfg.decoder_name = scenario.decoder_name;

  const std::string label =
      scenario.name + "_chunk" + std::to_string(chunk_index);
  chunk.result =
      newcode::run_pipeline(params, pipeline_cfg, label, scenario.ebn0_db);
  return chunk;
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
  oss << "[POINT] " << point.scenario.name
      << " chunks=" << point.chunks_completed
      << " inflight=" << state.inflight_chunks
      << " post=" << point.post_fec.errors << "/" << point.post_fec.total;
  if (point.post_fec.total > 0) {
    oss << " ber=" << std::scientific << point.post_fec.ber << std::defaultfloat;
  }
  if (point.post_fec.errors == 0 && point.post_fec.total > 0) {
    oss << " ub≈" << std::scientific << point.post_ber_upper_bound
        << std::defaultfloat;
  }
  if (state.started) {
    const auto elapsed =
        std::chrono::duration<double>(std::chrono::steady_clock::now() -
                                      state.start_time);
    oss << " elapsed=" << format_duration(elapsed);
  }
  log << oss.str() << "\n";
}

bool mark_stop_if_needed(PointState& state, ofec_sweep::detail::DualOut& log) {
  if (state.launch_stopped) {
    return false;
  }

  state.point.stop_reason = evaluate_stop_reason(state.point);
  update_upper_bound_fields(state.point);
  if (state.point.stop_reason == StopReason::NoData) {
    return false;
  }

  state.launch_stopped = true;
  log_point_progress(state, true, log);
  log << "[INFO] stop reason for " << state.point.scenario.name << ": "
      << stop_reason_to_string(state.point.stop_reason) << "\n";
  return true;
}

void finalize_point_if_ready(PointState& state,
                             ofec_sweep::detail::DualOut& log,
                             std::size_t& completed_count) {
  if (state.completed || !state.launch_stopped || state.inflight_chunks != 0) {
    return;
  }

  if (state.point.stop_reason == StopReason::NoData) {
    state.point.stop_reason = evaluate_stop_reason(state.point);
    update_upper_bound_fields(state.point);
  } else {
    update_upper_bound_fields(state.point);
  }
  if (state.started) {
    state.point.elapsed_seconds =
        std::chrono::duration<double>(std::chrono::steady_clock::now() -
                                      state.start_time)
            .count();
  }

  log << "[SUMMARY] " << state.point.scenario.name
      << " Post-FEC BER=" << state.point.post_fec.ber
      << " (errs=" << state.point.post_fec.errors << "/"
      << state.point.post_fec.total << ")";
  if (state.point.post_fec.errors == 0 &&
      !std::isnan(state.point.post_ber_upper_bound)) {
    log << " upper_bound<=" << state.point.post_ber_upper_bound;
  }
  log << " | chunks=" << state.point.chunks_completed
      << " | stop=" << stop_reason_to_string(state.point.stop_reason)
      << " | elapsed=" << format_duration(
             std::chrono::duration<double>(state.point.elapsed_seconds))
      << "\n";

  state.completed = true;
  ++completed_count;
}

std::vector<AggregatedPointResult> run_low_ber_scenarios_global(
    const std::vector<ofec_sweep::detail::SweepScenario>& scenarios,
    const ofec_sweep::SweepParameterConfig& config,
    unsigned total_workers,
    unsigned base_inflight_per_point,
    unsigned dynamic_inflight_cap,
    ofec_sweep::detail::DualOut& log) {
  std::vector<PointState> states;
  states.reserve(scenarios.size());
  for (const auto& scenario : scenarios) {
    PointState state;
    state.point.scenario = scenario;
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
      log << "\n[RUN] " << state.point.scenario.name << "\n";
    }

    const std::size_t chunk_index = state.next_chunk_index++;
    const auto scenario = state.point.scenario;
    active.push_back(ActiveChunk{
        point_index,
        chunk_index,
        std::async(std::launch::async, [scenario, config, chunk_index]() {
          return run_chunk(scenario, config, chunk_index);
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
      accumulate_chunk(state.point, chunk);
      update_upper_bound_fields(state.point);
      log_point_progress(state, false, log);
      mark_stop_if_needed(state, log);
      finalize_point_if_ready(state, log, completed_count);
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

ofec_sweep::SweepParameterConfig build_config() {
  const auto& selected_mux_bypass_edges =
      app_mux::bypass_edges_for_scheme(kMuxBypassScheme);

  ofec_sweep::SweepParameterConfig config;
  config.base_params.TILES_PER_WIN = kTilesPerWindow;
  config.base_params.HARD_TILE_LIST.assign(kTilesPerWindow, 0);
  config.base_params.BITGEN_RANDOM_BITS = kGenerateRandomBits;
  config.base_params.NORMALIZE_KNOWN_PREFIX_TAIL = kNormalizeKnownPrefixTail;
  config.base_params.ENABLE_EARLY_STOP = kEnableEarlyStop;
  config.base_params.LLR_CLIP_RATIO = kQuantClipRatio;
  config.base_params.LLR_BITS = kLlrBits;
  config.base_params.SISO_ACTIVE_LIST = kSisoActiveList;
  config.base_params.HYBRID_ENABLE = kHybridEnable;
  config.base_params.HYBRID_ENABLE_LIST = kHybridEnableList;
  config.base_params.HYBRID_CLASSIFIER_MODE = kHybridClassifierMode;
  config.base_params.HYBRID_USE_FAST_CLASSIFIER =
      kHybridClassifierMode != newcode::HybridClassifierMode::LegacyHardDecode;
  config.base_params.HYBRID_NORMALIZE_SOFT_ONLY = kHybridNormalizeSoftOnly;
  if (!kHybridEnableList.empty() &&
      kHybridEnableList.size() != kTilesPerWindow) {
    throw std::invalid_argument(
        "kHybridEnableList length must equal kTilesPerWindow when non-empty");
  }

  config.interleaver_name = kInterleaverName;
  config.decoder_name = kDecoderName;
  for (const char* decoder_name : kDecoderNameCandidates) {
    if (decoder_name && decoder_name[0] != '\0') {
      config.decoder_name_candidates.emplace_back(decoder_name);
    }
  }
  config.bits_per_symbol = kBitsPerSymbol;
  config.normalize_extrinsic = kNormalizeExtrinsic;
  config.generate_random_bits = kGenerateRandomBits;
  config.normalize_known_prefix_tail = kNormalizeKnownPrefixTail;
  config.quant_clip_ratio = kQuantClipRatio;
  config.quiet_pipeline = kQuietConsole;
  config.quiet_logs = kQuietConsole;

  config.enable_early_stop = kEnableEarlyStop;
  config.early_stop_enable_list = kEarlyStopEnableList;
  config.early_stop_condition_mode = kEarlyStopConditionMode;
  config.early_stop_action_mode = kEarlyStopActionMode;
  config.early_stop_bind_group_size = kEarlyStopBindGroupSize;
  config.early_stop_cond_v1_require_bch = kEarlyStopCondV1RequireBch;
  config.early_stop_cond_v1_require_overall = kEarlyStopCondV1RequireOverall;
  config.early_stop_v2_llr_abs_threshold = kEarlyStopV2LlrAbsThreshold;
  config.early_stop_v2_max_unreliable_bits = kEarlyStopV2MaxUnreliableBits;
  config.early_stop_cond_v2_include_overall = kEarlyStopCondV2IncludeOverall;
  config.early_stop_action_residual_divisor = kEarlyStopActionResidualDivisor;
  config.early_stop_action_hard_llr_mag = kEarlyStopActionHardLlrMag;
  config.chase_n_test = kChaseNTest;
  config.chase_topk_keep = kChaseTopkKeep;
  config.chase_group_minima_bits = kChaseGroupMinimaBits;

  config.chase_n_test_candidates = kChaseNTestCandidates;
  config.chase_topk_keep_candidates = kChaseTopkKeepCandidates;
  config.chase_group_minima_bits_candidates = kChaseGroupMinimaBitsCandidates;
  config.early_stop_action_beta_start_candidates =
      kEarlyStopActionBetaStartCandidates;
  config.early_stop_action_beta_step_candidates =
      kEarlyStopActionBetaStepCandidates;
  config.early_stop_action_hard_llr_mag_candidates =
      kEarlyStopActionHardLlrMagCandidates;
  config.early_stop_condition_candidates = kEarlyStopConditionCandidates;
  config.early_stop_action_candidates = kEarlyStopActionCandidates;
  config.early_stop_cond_v1_require_bch_candidates =
      kEarlyStopCondV1RequireBchCandidates;
  config.early_stop_cond_v1_require_overall_candidates =
      kEarlyStopCondV1RequireOverallCandidates;
  config.early_stop_v2_llr_abs_threshold_candidates =
      kEarlyStopV2LlrAbsThresholdCandidates;
  config.early_stop_v2_max_unreliable_bits_candidates =
      kEarlyStopV2MaxUnreliableBitsCandidates;

  config.siso_active_list = kSisoActiveList;
  config.mux_group_g = kMuxGroupG;
  config.mux_scheduling_mode = kMuxSchedulingMode;
  config.mux_scheduling_mode_candidates = kMuxSchedulingModeCandidates;
  config.mux_early_stop_priority_rule = kMuxPriorityRule;
  config.mux_early_stop_priority_rule_candidates = kMuxPriorityRuleCandidates;
  config.mux_enable_reconfig = kMuxEnableReconfig;
  config.mux_bypass_scheme = kMuxBypassScheme;
  config.mux_extra_bypass_edges = selected_mux_bypass_edges;
  config.chase_l_candidates = kChaseLCandidates;

  config.base_params.BITGEN_SEED = kBitgenSeed;
  config.base_params.CHANNEL_SEED = kChannelSeed;
  config.bitgen_seed_candidates = {kBitgenSeed};
  config.channel_seed_candidates = {kChannelSeed};
  config.bitgen_seed_count = 1;
  config.channel_seed_count = 1;

  config.ebn0_start = kEbN0Start;
  config.ebn0_end = kEbN0End;
  config.ebn0_points = kEbN0Points;

  config.explicit_patterns = kExplicitAlphaBetaSets;

  config.base_params.debug_trace = newcode::Params::DebugTraceConfig{
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

  config.base_params.CHASE_NTEST = kChaseNTest;
  config.base_params.CHASE_TOPK_KEEP = kChaseTopkKeep;
  config.base_params.CHASE_GROUP_MINIMA_BITS = kChaseGroupMinimaBits;
  config.base_params.MUX_GROUP_G = kMuxGroupG;
  config.base_params.MUX_SCHEDULING_MODE = kMuxSchedulingMode;
  config.base_params.MUX_EARLY_STOP_PRIORITY_RULE = kMuxPriorityRule;
  config.base_params.MUX_ENABLE_RECONFIG = kMuxEnableReconfig;
  config.base_params.MUX_EXTRA_BYPASS_EDGES = selected_mux_bypass_edges;
  config.base_params.HYBRID_ENABLE = kHybridEnable;
  config.base_params.HYBRID_ENABLE_LIST = kHybridEnableList;
  config.base_params.HYBRID_CLASSIFIER_MODE = kHybridClassifierMode;
  config.base_params.HYBRID_USE_FAST_CLASSIFIER =
      kHybridClassifierMode != newcode::HybridClassifierMode::LegacyHardDecode;
  config.base_params.HYBRID_NORMALIZE_SOFT_ONLY = kHybridNormalizeSoftOnly;
  config.base_params.EARLY_STOP_ENABLE_LIST = kEarlyStopEnableList;
  config.base_params.EARLY_STOP_CONDITION_MODE = kEarlyStopConditionMode;
  config.base_params.EARLY_STOP_ACTION_MODE = kEarlyStopActionMode;
  config.base_params.EARLY_STOP_BIND_GROUP_SIZE = kEarlyStopBindGroupSize;
  config.base_params.EARLY_STOP_COND_V1_REQUIRE_BCH =
      kEarlyStopCondV1RequireBch;
  config.base_params.EARLY_STOP_COND_V1_REQUIRE_OVERALL =
      kEarlyStopCondV1RequireOverall;
  config.base_params.EARLY_STOP_V2_LLR_ABS_THRESHOLD =
      kEarlyStopV2LlrAbsThreshold;
  config.base_params.EARLY_STOP_V2_MAX_UNRELIABLE_BITS =
      kEarlyStopV2MaxUnreliableBits;
  config.base_params.EARLY_STOP_COND_V2_INCLUDE_OVERALL =
      kEarlyStopCondV2IncludeOverall;
  config.base_params.EARLY_STOP_ACTION_RESIDUAL_DIVISOR =
      kEarlyStopActionResidualDivisor;
  config.base_params.EARLY_STOP_ACTION_HARD_LLR_MAG =
      kEarlyStopActionHardLlrMag;

  return config;
}

}  // namespace

int main() {
  const std::filesystem::path data_dir = "data";
  io::ensure_dir(data_dir);
  const std::string run_id = utils::now_stamp();
  const std::string log_path =
      (data_dir / ("run_" + run_id + "_sweep3.log")).string();
  ofec_sweep::detail::DualOut out(std::cout, log_path, !kQuietConsole);

  const std::string csv_path =
      (data_dir / ("ofec_sweep3_results_" + run_id + ".csv")).string();
  ensure_csv_header_sweep3(csv_path);
  std::ofstream csv(csv_path, std::ios::out | std::ios::app);
  csv.setf(std::ios::fixed);
  csv << std::setprecision(10);

  ofec_sweep::SweepParameterConfig config = build_config();
  const std::vector<float> ebn0_values =
      ofec_sweep::detail::build_ebn0_values(config);
  const std::vector<int> bitgen_seeds =
      config.bitgen_seed_candidates.empty()
          ? std::vector<int>{config.base_params.BITGEN_SEED}
          : config.bitgen_seed_candidates;
  const std::vector<int> channel_seeds =
      config.channel_seed_candidates.empty()
          ? std::vector<int>{config.base_params.CHANNEL_SEED}
          : config.channel_seed_candidates;
  auto scenarios = ofec_sweep::detail::build_scenarios(
      config, ebn0_values, bitgen_seeds, channel_seeds);

  if (scenarios.empty()) {
    std::cerr << "[ERROR] no scenarios generated for ofec_sweep3\n";
    return 1;
  }

  const unsigned available_workers =
      ofec_sweep::detail::resolve_worker_count(config);
  const unsigned total_workers =
      (kMaxTotalWorkers == 0)
          ? available_workers
          : std::max(1u, std::min(available_workers, kMaxTotalWorkers));
  const unsigned base_inflight_per_point =
      (kMaxInflightChunksPerPoint == 0)
          ? total_workers
          : std::max(1u, std::min(total_workers, kMaxInflightChunksPerPoint));
  const unsigned dynamic_inflight_cap =
      (kMaxDynamicInflightChunksPerPoint == 0)
          ? total_workers
          : std::max(1u,
                     std::min(total_workers,
                              kMaxDynamicInflightChunksPerPoint));

  out << "[INFO] total scenarios = " << scenarios.size() << "\n";
  out << "[INFO] total workers = " << total_workers << " (available="
      << available_workers << ")\n";
  out << "[INFO] base max inflight chunks per point = "
      << base_inflight_per_point << "\n";
  out << "[INFO] dynamic inflight per point = "
      << (kEnableDynamicInflightPerPoint ? "enabled" : "disabled")
      << " (cap=" << dynamic_inflight_cap << ")\n";
  out << "[INFO] chunk_num_info_bits = " << kChunkNumInfoBits << "\n";
  out << "[INFO] target_post_errors = " << kTargetPostErrors << "\n";
  out << "[INFO] max_post_fec_total_bits = " << kMaxPostFecTotalBits << "\n";
  if (kEnableZeroErrorUpperBound) {
    out << "[INFO] zero-error upper-bound stop enabled, target="
        << std::scientific << kTargetBerUpperBound << std::defaultfloat
        << ", confidence=" << kConfidenceLevel << "\n";
  }

  auto results = run_low_ber_scenarios_global(scenarios,
                                              config,
                                              total_workers,
                                              base_inflight_per_point,
                                              dynamic_inflight_cap,
                                              out);
  std::sort(results.begin(),
            results.end(),
            [](const AggregatedPointResult& lhs,
               const AggregatedPointResult& rhs) {
              if (lhs.scenario.ebn0_db != rhs.scenario.ebn0_db) {
                return lhs.scenario.ebn0_db < rhs.scenario.ebn0_db;
              }
              return lhs.scenario.name < rhs.scenario.name;
            });
  for (const auto& point : results) {
    write_csv_row_sweep3(csv, utils::now_stamp(), run_id, point);
  }
  csv.flush();

  return 0;
}
