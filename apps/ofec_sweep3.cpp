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

// 基础发射与信道参数
static constexpr const char* kInterleaverName = "identity";  // 交织器名称，identity 表示不交织
static constexpr unsigned kBitsPerSymbol = 1;                  // 每个调制符号携带的比特数：1=BPSK，偶数=QAM
static constexpr bool kGenerateRandomBits = true;              // true=发送随机信息比特，false=发送全 0 比特
static constexpr int kBitgenSeed = 618911;                     // 基础比特种子；每个 chunk 会在此基础上派生
static constexpr float kEbN0Start = 3.05f;                     // 扫描起始 Eb/N0（dB）
static constexpr float kEbN0End = 3.13f;                       // 扫描结束 Eb/N0（dB）
static constexpr int kEbN0Points = 9;                          // Eb/N0 采样点数；含首尾端点
static constexpr int kChannelSeed = 1701;                      // 基础信道种子；每个 chunk 会在此基础上派生
static constexpr std::size_t kLlrBits = 6;                     // LLR 位宽：16=浮点，2~15=qfloat
static constexpr float kQuantClipRatio = 0.5f;                 // 动态 clip 比例，0 表示禁用自适应 clip

// Monte Carlo 聚合与并发控制
static constexpr std::size_t kTilesPerWindow = 6;               // 需与显式 alpha/beta 列表长度一致；新 Level 5/6 共享模式要求为 6
static constexpr std::size_t kChunkNumInfoBits = 16 * 132 * 16 * 111; // 每个 Monte Carlo chunk 的输入信息比特数
static constexpr std::size_t kTargetPostErrors = 5000;          // 单个 Eb/N0 点累计到这么多 post-FEC 错误后即可停止
static constexpr std::size_t kMaxPostFecTotalBits = 50e8;       // 单个 Eb/N0 点允许累计比较的最大 post-FEC 比特数
static constexpr unsigned kMaxTotalWorkers = 0;                 // 全局同时运行 chunk 数上限；0=自动使用 NTHREADS/机器可用 worker
static constexpr unsigned kMaxInflightChunksPerPoint = 16;      // 单个 Eb/N0 点的基础挂起 chunk 上限；0=不额外限制
static constexpr bool kEnableDynamicInflightPerPoint = true;    // true=剩余 Eb/N0 点变少时动态提高单点挂起上限
static constexpr unsigned kMaxDynamicInflightChunksPerPoint = 0; // 动态单点挂起上限封顶；0=不封顶，最多到全局 worker
static constexpr double kProgressReportIntervalSeconds = 15.0;  // 进度快照输出周期；<=0 表示关闭周期性汇总
static constexpr bool kEnableZeroErrorUpperBound = false;       // true=零错时使用上置信界提前停止
static constexpr double kTargetBerUpperBound = 1e-8;             // 零错上界目标：若上界已低于此值则提前停止
static constexpr double kConfidenceLevel = 0.95;                 // 零错上界使用的置信水平，例如 0.95 表示 95%

// Decoder 与 Chase 扫描参数
static constexpr const char* kDecoderName = "chase_baseline";       // 默认 decoder 名称
static const std::vector<const char*> kDecoderNameCandidates = {};  // decoder 扫描候选；空表示不扫 decoder 维度
static constexpr bool kNormalizeExtrinsic = false;                  // 是否对 decoder 输出 extrinsic 做归一化
static constexpr bool kNormalizeKnownPrefixTail = false;            // 是否对 known-prefix 后的尾部 LLR 做归一化
static const std::vector<int> kChaseLCandidates = {6};              // Chase L 扫描候选
static constexpr int kChaseNTest = 64;                              // 默认 Chase NTEST；若不单独扫则等于实际使用值
static const std::vector<int> kChaseNTestCandidates = {};           // Chase NTEST 扫描候选；空表示不单独扫描
static constexpr int kChaseTopkKeep = 24;                           // 与阶段三冻结基线一致；chase_baseline 下仅作为公共参数留档
static const std::vector<int> kChaseTopkKeepCandidates = {};        // top-k 保留数扫描候选
static constexpr int kChaseGroupMinimaBits = 4;                     // group-minima decoder 的分组 bit 数
static const std::vector<int> kChaseGroupMinimaBitsCandidates = {4}; // group-minima 分组 bit 数扫描候选

// Level 1-4 常规 MUX 参数；Level 5/6 共享模式开启后，末两项不参与独立 MUX 分配
static const std::vector<int> kSisoActiveList = {32, 32, 32, 32, 8, 4}; // 每个 tile 的 SISO 预算
static const std::vector<int> kHiHoActiveList = {32, 32, 32, 32, 8, 4}; // 每个 tile 的 HISO 预算
static constexpr int kMuxGroupG = 1;                                // MUX 分组粒度；1=全局池化
static constexpr int kMuxSchedulingMode = 0;                        // 0=legacy，1=按 early-stop 细节排序
static const std::vector<int> kMuxSchedulingModeCandidates = {};    // MUX 调度模式扫描候选
static constexpr int kMuxPriorityRule = 0;                          // 0=更差优先，1=更接近通过优先
static const std::vector<int> kMuxPriorityRuleCandidates = {};      // MUX 优先级规则扫描候选
static constexpr bool kMuxEnableReconfig = false;                   // 是否启用 reconfig MUX 调度
static constexpr int kMuxBypassScheme = 1;                          // 旁路边方案编号；仅在 reconfig 打开时参与调度

// Hybrid 前置分类与 hard-finish 参数
static constexpr bool kHybridEnable = true;                         // true=启用方案三软硬混合前置分流
static const std::vector<int> kHybridEnableList = {0, 0, 0, 0, 1, 1}; // 按 tile 覆盖 hybrid 开关；空=沿用 kHybridEnable
static constexpr float kHybridHardLlrMag = 99.0f;                   // hybrid hard-finish 默认输出 |LLR| 幅度
static const std::vector<float> kHybridHardLlrMagList = {};         // 按 tile 覆盖 hybrid hard-finish |LLR| 幅度
static constexpr newcode::HybridClassifierMode kHybridClassifierMode =
    newcode::HybridClassifierMode::FriendS1S3WithS0Classifier;      // LegacyHardDecode / RepoFastClassifier / FriendS1S3Classifier / FriendS1S3WithS0Classifier
static constexpr newcode::HybridSisoBackfillMode kHybridSisoBackfillMode =
    newcode::HybridSisoBackfillMode::ParityOneAndTwoErrorPriority;  // Disabled / TwoErrorOnly / OneAndTwoErrorPriority / ParityOneAndTwoErrorPriority
static constexpr bool kHybridNormalizeSoftOnly = false;             // true=只归一化 soft rows，false=保持兼容行为

// 阶段三 TwoMain 策略预设。一次进程只选择其中一种策略，并只展开
// Eb/N0 维度；三种策略分别启动三次，避免在同一任务中混合算法配置。
struct TwoMainSweepPolicy {
  const char* label;
  newcode::TwoMainHisoOutputMode output_mode;
  float m2;
  float rho_corr;
  float rho_keep;
  unsigned hiso_class_mask;
};

enum class TwoMainSweepSelection {
  Legacy,
  Scheme1NoTwoMainHiso,
  Scheme6AM2_48,
};

static constexpr TwoMainSweepSelection kDefaultTwoMainSweepSelection =
    TwoMainSweepSelection::Legacy;

constexpr TwoMainSweepPolicy selected_twomain_policy(
    TwoMainSweepSelection selection) {
  switch (selection) {
    case TwoMainSweepSelection::Legacy:
      return {"legacy", newcode::TwoMainHisoOutputMode::Legacy,
              99.0f, 1.0f, 1.0f, 0x0fu};
    case TwoMainSweepSelection::Scheme1NoTwoMainHiso:
      return {"scheme1_no_twomain_hiso",
              newcode::TwoMainHisoOutputMode::Legacy,
              99.0f, 1.0f, 1.0f, 0x07u};
    case TwoMainSweepSelection::Scheme6AM2_48:
      return {"scheme6a_m2_48",
              newcode::TwoMainHisoOutputMode::UnifiedParameterized,
              48.0f, 1.0f, 1.0f, 0x0fu};
  }
  return {"legacy", newcode::TwoMainHisoOutputMode::Legacy,
          99.0f, 1.0f, 1.0f, 0x0fu};
}

TwoMainSweepSelection parse_twomain_policy_selection(
    const std::string& value) {
  if (value == "legacy") {
    return TwoMainSweepSelection::Legacy;
  }
  if (value == "scheme1_no_twomain_hiso") {
    return TwoMainSweepSelection::Scheme1NoTwoMainHiso;
  }
  if (value == "scheme6a_m2_48") {
    return TwoMainSweepSelection::Scheme6AM2_48;
  }
  throw std::invalid_argument(
      "unknown TwoMain policy '" + value +
      "'; expected legacy, scheme1_no_twomain_hiso, or scheme6a_m2_48");
}

// Level 5/6 共享：四行分组、负载排序与多轮 MUX 调度
static constexpr bool kLevel56SharedEnable = true;                  // true=第五/六级共享 HISO/SISO
static constexpr bool kLevel56TemporalLookaheadEnable = false;      // 旧三时刻 lookahead；与 buffered FIFO 互斥
static constexpr int kLevel56TemporalGroupLoadThreshold = 16;       // temporal 分支阈值：X+K1+K2 小于该值时补解 t=0
static constexpr bool kLevel56BufferedFifoEnable = true;            // 新方案：完整 64-code batch FIFO
static constexpr std::size_t kLevel56BufferRows = 16000;            // 阶段三参考条件：足够长的 R_buf，用于隔离 TwoMain 动作
static constexpr bool kLevel56DrainAtFrameEnd = true;               // 每个 Monte Carlo chunk 末尾排空 FIFO
static constexpr int kLevel56SharedHisoActive = 8;                  // 新分组方案固定使用 8 个共享 HISO entry slot
static constexpr int kLevel56SharedSisoActive = 8;                  // 新分组方案固定使用 8 个共享 SISO entry slot
static constexpr std::size_t kLevel56Group4MaxEntries = 8;          // 与阶段三冻结基线一致
static constexpr newcode::Level56PriorityMode kLevel56PriorityMode =
    newcode::Level56PriorityMode::Level5First; // Level5First=组负载相同时 Level 5 优先；Level6First=Level 6 优先
static constexpr newcode::Level56ScheduleMode kLevel56ScheduleMode =
    newcode::Level56ScheduleMode::Group4LoadSortedMultiround; // GlobalPriority=旧全局优先级；Group4LoadSortedMultiround=16 个四行组按负载排序并多轮调度
static constexpr newcode::Level56EarlyStopGroupUpdateMode
    kLevel56EarlyStopGroupUpdateMode =
        newcode::Level56EarlyStopGroupUpdateMode::AllGroups; // temporal 方案固定使用 AllGroups
static constexpr bool kLevel56SingleLevelSelectEnable = false; // 新分组多轮模式必须关闭；true=按 early-stop 命中数动态只解一级
static constexpr bool kLevel56UnselectedEarlyStopActionEnable = true; // 与阶段三冻结基线一致；single-level selection 关闭时不改变动作
// 只读观测：用于在服务器扫描结果中同时核查 FIFO 闭环和 TwoMain 动作。
// 不参与分类、调度、译码或写回决策。
static constexpr bool kLevel56ScheduleObservabilityEnable = true;

// 早停参数：总开关 -> 条件 -> 条件细参 -> 动作 -> 动作细参
static constexpr bool kEnableEarlyStop = true;                      // 早停总开关
static const std::vector<int> kEarlyStopEnableList = {1, 1, 1, 1, 1, 1}; // 按 tile 覆盖早停开关；空表示全部沿用总开关
static constexpr int kEarlyStopConditionMode = 1;                   // 早停条件模式：1=v1，2=v2
static const std::vector<int> kEarlyStopConditionCandidates = {};   // 早停条件模式扫描候选
static constexpr int kEarlyStopActionMode = 7;                      // 早停动作模式：1~8
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

// 显式 alpha/beta 场景；每个列表长度必须等于 kTilesPerWindow
static const std::vector<ofec_sweep::ExplicitAlphaBetaPattern> kExplicitAlphaBetaSets = {
    {"custom_label",  // 该组显式 alpha/beta 的标签，会进入场景名
     {0.428571, 0.447738, 0.482782, 0.528162, 0.581902, 0.642857},
     {2.857143, 6.179301, 12.253626, 20.119585, 29.434408, 40.000000},
     {99.857143, 99.179301, 99.253626, 99.119585, 99, 99}},
};

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
  std::size_t level56_completed_batches = 0;
  std::size_t level56_forced_evicted_batches = 0;
  std::size_t level56_max_fifo_depth = 0;
  std::size_t twomain_hiso_actions = 0;
  std::size_t twomain_siso_actions = 0;
  std::size_t twomain_pending_actions = 0;
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

const char* level56_early_stop_group_update_mode_name(
    newcode::Level56EarlyStopGroupUpdateMode mode) {
  switch (mode) {
    case newcode::Level56EarlyStopGroupUpdateMode::AllGroups:
      return "all_groups";
    case newcode::Level56EarlyStopGroupUpdateMode::EnteredGroupsOnly:
      return "entered_groups_only";
    case newcode::Level56EarlyStopGroupUpdateMode::FillIdleEntries:
      return "fill_idle_entries";
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

const char* twomain_hiso_output_mode_name(
    newcode::TwoMainHisoOutputMode mode) {
  switch (mode) {
    case newcode::TwoMainHisoOutputMode::Legacy:
      return "legacy";
    case newcode::TwoMainHisoOutputMode::UnifiedParameterized:
      return "unified_parameterized";
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

std::string format_eta_seconds(double seconds) {
  if (!std::isfinite(seconds) || seconds < 0.0) {
    return "--:--";
  }
  return format_duration(std::chrono::duration<double>(seconds));
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
          "level56_completed_batches,level56_forced_evicted_batches,level56_max_fifo_depth,twomain_hiso_actions,twomain_siso_actions,twomain_pending_actions,"
          "bitgen_seed_base,channel_seed_base,"
          "alpha_list,beta_list,early_stop_beta_list,"
          "chase_L,chase_n_test,chase_topk_keep,chase_group_minima_bits,"
          "mux_group_g,mux_scheduling_mode,mux_early_stop_priority_rule,mux_bypass_scheme,"
          "hybrid_enable,hybrid_enable_list,hybrid_classifier_mode,hybrid_siso_backfill_mode,hybrid_normalize_soft_only,"
          "twomain_policy,twomain_hiso_output_mode,twomain_hiso_m2,twomain_hiso_rho_corr,twomain_hiso_rho_keep,level56_hiso_allowed_class_mask,"
          "level56_shared_enable,level56_temporal_lookahead_enable,level56_buffered_fifo_enable,level56_buffer_rows,level56_drain_at_frame_end,level56_group4_max_entries,level56_shared_hiso_active,level56_shared_siso_active,level56_priority_mode,level56_schedule_mode,level56_early_stop_group_update_mode,"
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
      << point.level56_completed_batches << ','
      << point.level56_forced_evicted_batches << ','
      << point.level56_max_fifo_depth << ','
      << point.twomain_hiso_actions << ','
      << point.twomain_siso_actions << ','
      << point.twomain_pending_actions << ','
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
      << hybrid_siso_backfill_mode_name(kHybridSisoBackfillMode) << ','
      << (kHybridNormalizeSoftOnly ? 1 : 0) << ','
      << point.scenario.twomain_policy_label << ','
      << twomain_hiso_output_mode_name(
             point.scenario.twomain_hiso_output_mode)
      << ','
      << point.scenario.twomain_hiso_m2 << ','
      << point.scenario.twomain_hiso_rho_corr << ','
      << point.scenario.twomain_hiso_rho_keep << ','
      << point.scenario.level56_hiso_allowed_class_mask << ','
      << (kLevel56SharedEnable ? 1 : 0) << ','
      << (kLevel56TemporalLookaheadEnable ? 1 : 0) << ','
      << (kLevel56BufferedFifoEnable ? 1 : 0) << ','
      << kLevel56BufferRows << ','
      << (kLevel56DrainAtFrameEnd ? 1 : 0) << ','
      << kLevel56Group4MaxEntries << ','
      << point.scenario.level56_shared_hiso_active << ','
      << point.scenario.level56_shared_siso_active << ','
      << static_cast<int>(kLevel56PriorityMode) << ','
      << static_cast<int>(kLevel56ScheduleMode) << ','
      << level56_early_stop_group_update_mode_name(
             kLevel56EarlyStopGroupUpdateMode)
      << ','
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

  for (const auto& sample : chunk.result.level56_buffered_time_samples) {
    point.level56_completed_batches += sample.completed_batches;
    point.level56_forced_evicted_batches += sample.forced_evicted_batches;
    point.level56_max_fifo_depth = std::max(
        point.level56_max_fifo_depth,
        std::max(sample.fifo_depth_before, sample.fifo_depth_after));
  }

  constexpr uint8_t kTwoMainClass =
      static_cast<uint8_t>(newcode::detail::HybridRowClass::TwoMain);
  constexpr uint8_t kHisoAction = 1;
  constexpr uint8_t kSisoAction = 2;
  constexpr uint8_t kPendingAction = 3;
  for (const auto& sample : chunk.result.level56_schedule_samples) {
    for (const auto& code : sample.codes) {
      if (code.hybrid_class != kTwoMainClass) {
        continue;
      }
      if (code.final_action == kHisoAction) {
        ++point.twomain_hiso_actions;
      } else if (code.final_action == kSisoAction) {
        ++point.twomain_siso_actions;
      } else if (code.final_action == kPendingAction) {
        ++point.twomain_pending_actions;
      }
    }
  }
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
  params.TWOMAIN_HISO_OUTPUT_MODE = scenario.twomain_hiso_output_mode;
  params.TWOMAIN_HISO_M2 = scenario.twomain_hiso_m2;
  params.TWOMAIN_HISO_RHO_CORR = scenario.twomain_hiso_rho_corr;
  params.TWOMAIN_HISO_RHO_KEEP = scenario.twomain_hiso_rho_keep;
  params.LEVEL56_HISO_ALLOWED_CLASS_MASK =
      scenario.level56_hiso_allowed_class_mask;
  params.LEVEL56_SHARED_HISO_ACTIVE = scenario.level56_shared_hiso_active;
  params.LEVEL56_SHARED_SISO_ACTIVE = scenario.level56_shared_siso_active;
  return params;
}

std::vector<ofec_sweep::detail::SweepScenario> apply_twomain_policy(
    const std::vector<ofec_sweep::detail::SweepScenario>& base_scenarios,
    const TwoMainSweepPolicy& policy) {
  std::vector<ofec_sweep::detail::SweepScenario> scenarios;
  scenarios.reserve(base_scenarios.size());

  for (const auto& base : base_scenarios) {
    ofec_sweep::detail::SweepScenario scenario = base;
    scenario.twomain_policy_label = policy.label;
    scenario.twomain_hiso_output_mode = policy.output_mode;
    scenario.twomain_hiso_m2 = policy.m2;
    scenario.twomain_hiso_rho_corr = policy.rho_corr;
    scenario.twomain_hiso_rho_keep = policy.rho_keep;
    scenario.level56_hiso_allowed_class_mask = policy.hiso_class_mask;
    scenario.level56_shared_hiso_active = kLevel56SharedHisoActive;
    scenario.level56_shared_siso_active = kLevel56SharedSisoActive;
    scenario.name = std::string(policy.label) + "_" + base.name;
    scenarios.push_back(std::move(scenario));
  }
  return scenarios;
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

bool mark_stop_if_needed(PointState& state, ofec_sweep::detail::DualOut&) {
  if (state.launch_stopped) {
    return false;
  }

  state.point.stop_reason = evaluate_stop_reason(state.point);
  update_upper_bound_fields(state.point);
  if (state.point.stop_reason == StopReason::NoData) {
    return false;
  }

  state.launch_stopped = true;
  return true;
}

void finalize_point_if_ready(PointState& state,
                             ofec_sweep::detail::DualOut&,
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

  state.completed = true;
  ++completed_count;
}

double point_progress_fraction(const PointState& state) {
  if (state.completed) {
    return 1.0;
  }

  double progress = 0.0;
  if (kTargetPostErrors > 0) {
    progress = std::max(
        progress,
        static_cast<double>(state.point.post_fec.errors) /
            static_cast<double>(kTargetPostErrors));
  }
  if (kMaxPostFecTotalBits > 0) {
    progress = std::max(
        progress,
        static_cast<double>(state.point.post_fec.total) /
            static_cast<double>(kMaxPostFecTotalBits));
  }
  if (kEnableZeroErrorUpperBound && state.point.post_fec.errors == 0) {
    const double bits_needed_for_zero_error_target =
        -std::log(1.0 - std::clamp(kConfidenceLevel, 1e-12, 1.0 - 1e-12)) /
        kTargetBerUpperBound;
    if (bits_needed_for_zero_error_target > 0.0 &&
        std::isfinite(bits_needed_for_zero_error_target)) {
      progress = std::max(
          progress,
          static_cast<double>(state.point.post_fec.total) /
              bits_needed_for_zero_error_target);
    }
  }

  return std::clamp(progress, 0.0, 0.999);
}

double point_elapsed_seconds(const PointState& state,
                             std::chrono::steady_clock::time_point now) {
  if (!state.started) {
    return 0.0;
  }
  return std::chrono::duration<double>(now - state.start_time).count();
}

double point_eta_seconds(const PointState& state,
                         std::chrono::steady_clock::time_point now) {
  if (state.completed) {
    return 0.0;
  }
  const double progress = point_progress_fraction(state);
  if (progress <= 0.0) {
    return std::numeric_limits<double>::infinity();
  }
  const double elapsed = point_elapsed_seconds(state, now);
  return elapsed * (1.0 / progress - 1.0);
}

void log_progress_snapshot(const std::vector<PointState>& states,
                           std::size_t completed_count,
                           std::size_t active_chunk_count,
                           unsigned total_workers,
                           std::chrono::steady_clock::time_point start_time,
                           ofec_sweep::detail::DualOut& log) {
  const auto now = std::chrono::steady_clock::now();
  const double total_elapsed_seconds =
      std::chrono::duration<double>(now - start_time).count();

  double overall_progress_sum = 0.0;
  for (const auto& state : states) {
    overall_progress_sum += point_progress_fraction(state);
  }
  const double overall_progress =
      states.empty() ? 1.0 : overall_progress_sum / static_cast<double>(states.size());
  const double overall_eta =
      overall_progress > 0.0
          ? total_elapsed_seconds * (1.0 / overall_progress - 1.0)
          : std::numeric_limits<double>::infinity();

  log << "\n[PROGRESS] completed=" << completed_count << "/" << states.size()
      << " active_chunks=" << active_chunk_count << "/" << total_workers
      << " overall≈" << std::fixed << std::setprecision(1)
      << (overall_progress * 100.0) << "%"
      << " elapsed=" << format_duration(std::chrono::duration<double>(total_elapsed_seconds))
      << " eta≈" << format_eta_seconds(overall_eta) << std::defaultfloat << "\n";

  for (const auto& state : states) {
    const char* status = "PENDING";
    if (state.completed) {
      status = "DONE";
    } else if (state.launch_stopped) {
      status = "DRAIN";
    } else if (state.started) {
      status = "RUN";
    }

    log << "  [" << status << "]"
        << " Eb/N0=" << std::fixed << std::setprecision(3)
        << state.point.scenario.ebn0_db;

    if (state.completed) {
      log << " eta=00:00";
    } else if (state.started) {
      log << " eta≈" << format_eta_seconds(point_eta_seconds(state, now));
    } else {
      log << " eta=--:--";
    }
    log << "\n";
  }
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
  const auto scheduler_start_time = std::chrono::steady_clock::now();
  auto next_progress_report_time =
      scheduler_start_time + std::chrono::duration<double>(kProgressReportIntervalSeconds);

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
    if (kProgressReportIntervalSeconds > 0.0 &&
        std::chrono::steady_clock::now() >= next_progress_report_time) {
      log_progress_snapshot(states,
                            completed_count,
                            active.size(),
                            total_workers,
                            scheduler_start_time,
                            log);
      next_progress_report_time =
          std::chrono::steady_clock::now() +
          std::chrono::duration<double>(kProgressReportIntervalSeconds);
    }

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
  config.base_params.HIHO_ACTIVE_LIST = kHiHoActiveList;
  config.base_params.HYBRID_ENABLE = kHybridEnable;
  config.base_params.HYBRID_ENABLE_LIST = kHybridEnableList;
  config.base_params.HYBRID_HARD_LLR_MAG = kHybridHardLlrMag;
  config.base_params.HYBRID_HARD_LLR_MAG_LIST = kHybridHardLlrMagList;
  config.base_params.HYBRID_CLASSIFIER_MODE = kHybridClassifierMode;
  config.base_params.HYBRID_SISO_BACKFILL_MODE = kHybridSisoBackfillMode;
  config.base_params.HYBRID_USE_FAST_CLASSIFIER =
      kHybridClassifierMode != newcode::HybridClassifierMode::LegacyHardDecode;
  config.base_params.HYBRID_NORMALIZE_SOFT_ONLY = kHybridNormalizeSoftOnly;
  config.base_params.LEVEL56_SHARED_ENABLE = kLevel56SharedEnable;
  config.base_params.LEVEL56_TEMPORAL_LOOKAHEAD_ENABLE =
      kLevel56TemporalLookaheadEnable;
  config.base_params.LEVEL56_TEMPORAL_GROUP_LOAD_THRESHOLD =
      kLevel56TemporalGroupLoadThreshold;
  config.base_params.LEVEL56_BUFFERED_FIFO_ENABLE =
      kLevel56BufferedFifoEnable;
  config.base_params.LEVEL56_BUFFER_ROWS = kLevel56BufferRows;
  config.base_params.LEVEL56_BUFFERED_FIFO_DRAIN_AT_FRAME_END =
      kLevel56DrainAtFrameEnd;
  config.base_params.LEVEL56_GROUP4_MAX_ENTRIES =
      kLevel56Group4MaxEntries;
  config.base_params.LEVEL56_SHARED_HISO_ACTIVE = kLevel56SharedHisoActive;
  config.base_params.LEVEL56_SHARED_SISO_ACTIVE = kLevel56SharedSisoActive;
  config.base_params.LEVEL56_PRIORITY_MODE = kLevel56PriorityMode;
  config.base_params.LEVEL56_SCHEDULE_MODE = kLevel56ScheduleMode;
  config.base_params.LEVEL56_EARLY_STOP_GROUP_UPDATE_MODE =
      kLevel56EarlyStopGroupUpdateMode;
  config.base_params.LEVEL56_SINGLE_LEVEL_SELECT_ENABLE =
      kLevel56SingleLevelSelectEnable;
  config.base_params.LEVEL56_UNSELECTED_EARLY_STOP_ACTION_ENABLE =
      kLevel56UnselectedEarlyStopActionEnable;
  config.base_params.LEVEL56_SCHEDULE_OBSERVABILITY_ENABLE =
      kLevel56ScheduleObservabilityEnable;
  if (!kHybridEnableList.empty() &&
      kHybridEnableList.size() != kTilesPerWindow) {
    throw std::invalid_argument(
        "kHybridEnableList length must equal kTilesPerWindow when non-empty");
  }
  if (!kHybridHardLlrMagList.empty() &&
      kHybridHardLlrMagList.size() != kTilesPerWindow) {
    throw std::invalid_argument(
        "kHybridHardLlrMagList length must equal kTilesPerWindow when non-empty");
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
  config.quiet_pipeline = true;
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
  config.hiho_active_list = kHiHoActiveList;
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
  config.base_params.HYBRID_HARD_LLR_MAG = kHybridHardLlrMag;
  config.base_params.HYBRID_HARD_LLR_MAG_LIST = kHybridHardLlrMagList;
  config.base_params.HYBRID_CLASSIFIER_MODE = kHybridClassifierMode;
  config.base_params.HYBRID_SISO_BACKFILL_MODE = kHybridSisoBackfillMode;
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

int main(int argc, char** argv) {
  TwoMainSweepSelection selection = kDefaultTwoMainSweepSelection;
  bool dry_run = false;
  try {
    for (int index = 1; index < argc; ++index) {
      const std::string arg = argv[index];
      if (arg == "--dry-run") {
        dry_run = true;
      } else if (arg == "--policy") {
        if (index + 1 >= argc) {
          throw std::invalid_argument("--policy requires a value");
        }
        selection = parse_twomain_policy_selection(argv[++index]);
      } else if (arg.rfind("--policy=", 0) == 0) {
        selection = parse_twomain_policy_selection(arg.substr(9));
      } else {
        throw std::invalid_argument("unknown argument '" + arg + "'");
      }
    }
  } catch (const std::invalid_argument& error) {
    std::cerr << "[ERROR] " << error.what() << "\n"
              << "Usage: " << argv[0]
              << " [--policy legacy|scheme1_no_twomain_hiso|scheme6a_m2_48]"
                 " [--dry-run]\n";
    return 2;
  }
  const TwoMainSweepPolicy policy = selected_twomain_policy(selection);

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
  const auto base_scenarios = ofec_sweep::detail::build_scenarios(
      config, ebn0_values, bitgen_seeds, channel_seeds);
  auto scenarios = apply_twomain_policy(base_scenarios, policy);

  if (scenarios.empty()) {
    std::cerr << "[ERROR] no scenarios generated for ofec_sweep3\n";
    return 1;
  }

  if (dry_run) {
    std::cout << "[DRY-RUN] ofec_sweep3 will not start Monte Carlo work\n"
              << "[DRY-RUN] base scenarios=" << base_scenarios.size()
              << ", total points=" << scenarios.size() << "\n"
              << "[DRY-RUN] Eb/N0=" << kEbN0Start << ".." << kEbN0End
              << " (" << kEbN0Points << " points), bitgen/channel seed="
              << kBitgenSeed << '/' << kChannelSeed << "\n"
              << "[DRY-RUN] chunk bits=" << kChunkNumInfoBits
              << ", target post errors=" << kTargetPostErrors
              << ", max post bits=" << kMaxPostFecTotalBits << "\n"
              << "[DRY-RUN] Level56 FIFO R_buf=" << kLevel56BufferRows
              << ", drain=" << (kLevel56DrainAtFrameEnd ? "on" : "off")
              << ", HISO/SISO=" << kLevel56SharedHisoActive << '/'
              << kLevel56SharedSisoActive << "\n"
              << "[DRY-RUN] policy=" << policy.label
              << ", output_mode="
              << twomain_hiso_output_mode_name(policy.output_mode)
              << ", M2=" << policy.m2
              << ", rho_corr/keep=" << policy.rho_corr << '/'
              << policy.rho_keep
              << ", HISO class mask=" << policy.hiso_class_mask << "\n";
    return 0;
  }

  const std::filesystem::path data_dir = "data";
  io::ensure_dir(data_dir);
  const std::string run_id = utils::now_stamp();
  const std::string log_path =
      (data_dir /
       ("run_" + run_id + "_sweep3_" + policy.label + ".log"))
          .string();
  ofec_sweep::detail::DualOut out(std::cout, log_path, !kQuietConsole);

  const std::string csv_path =
      (data_dir /
       ("ofec_sweep3_results_" + run_id + "_" + policy.label + ".csv"))
          .string();
  ensure_csv_header_sweep3(csv_path);
  std::ofstream csv(csv_path, std::ios::out | std::ios::app);
  csv.setf(std::ios::fixed);
  csv << std::setprecision(10);

  out << "[CONFIG] TwoMain policy=" << policy.label
      << " output_mode=" << twomain_hiso_output_mode_name(policy.output_mode)
      << " M2=" << policy.m2
      << " rho_corr/keep=" << policy.rho_corr << '/' << policy.rho_keep
      << " HISO_class_mask=" << policy.hiso_class_mask
      << " points=" << scenarios.size() << "\n";

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
