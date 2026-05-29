#include <fstream>
#include <iomanip>
#include <iostream>
#include <vector>

#include "mux_bypass_edges.hpp"
#include "newcode/ofec_single_runner.hpp"
#include "newcode/two_stream_shared_runner.hpp"

namespace {

// ======== 用户可改区域 ========
// two-stream shared decoder 的主验证入口。
// 这里尽量保持和 ofec_ber_window_probe.cpp 一样的“常量区 + main 接线”结构，
// 方便后续逐项对齐单流程配置。

// 运行标签
static constexpr const char* kLabel = "two_native_stream_shared_bchhard"; // 运行标签：日志和结果输出的名字前缀

// 信道 / 前端参数
static constexpr float kEbN0Db = 3.05f;                // 双流共用的信道 Eb/N0，单位 dB
static constexpr unsigned kBitsPerSymbol = 1;          // 每个调制符号携带的比特数：1=BPSK，偶数=QAM
static constexpr std::size_t kNumInfoBits =
    64u * 132u * 16u * 111u;                          // 每一路独立数据流生成的信息比特总数
static constexpr std::size_t kLlrBits = 6;             // LLR 位宽：16=浮点，2~15=qfloat
static constexpr float kQuantClipRatio = 0.5f;         // 动态 clip 比例，0 表示禁用自适应 clip

// 双流各自的随机种子
static constexpr int kBitgenSeedA = 20260320;          // Stream A 的信息比特随机种子
static constexpr int kChannelSeedA = 3182027;          // Stream A 的信道噪声随机种子
static constexpr int kBitgenSeedB = 20260319;          // Stream B 的信息比特随机种子
static constexpr int kChannelSeedB = 3182026;          // Stream B 的信道噪声随机种子

// Chase / decoder core 参数
static constexpr int kChaseL = 6;                      // Chase L：选择多少个最不可靠位置
static constexpr int kChaseNTestOverride = -1;        // Chase 测试序列数量，-1 表示默认按 2^L 生成
static constexpr int kChaseTopkKeep = 8;              // chase_topk_pruned 中保留参与外信息计算的 Top-K 候选数
static constexpr int kChaseGroupMinimaBits = 3;       // chase_group_minima 按前多少个 test-pattern 位分组
static constexpr const char* kInterleaverName = "identity"; // 交织器名称，identity 表示不交织
static constexpr const char* kDecoderName = "two_stream_shared_chase_baseline"; // 解码器名称：当前 two-stream shared 主流程入口；底层可映射到 chase_baseline / chase_overall_parity_search / chase_topk_pruned / chase_global_pair / chase_group_minima
static constexpr bool kNormalizeExtrinsic = false;    // 是否对 Chase 输出的 extrinsic 做归一化
static constexpr bool kGenerateRandomBits = true;     // true=发送随机信息比特，false=发送全 0 比特
static constexpr bool kNormalizeKnownPrefixTail = false; // 是否对 known-prefix 之后的尾部 LLR 做归一化

// Early-stop 参数
static constexpr bool kEnableEarlyStop = true;       // true=启用 early-stop，false=完全关闭
static const std::vector<int> kEarlyStopEnableList = {1, 1, 1, 1, 1, 1}; // 按 tile 覆盖 early-stop 总开关：0=关，非 0=开；空表示全部沿用 kEnableEarlyStop
static constexpr int kEarlyStopConditionMode = 1;     // early-stop 条件模式：1=v1，2=v2
static const std::vector<int> kEarlyStopConditionModeList = {}; // 按 tile 覆盖条件模式；空表示全部沿用 kEarlyStopConditionMode
static constexpr int kEarlyStopActionMode = 6;        // 命中 early-stop 后的动作模式：1=sign beta，2=residual only，3=hard-decode sign LLR，4=sign beta pre-div alpha，5=residual pre-div alpha plus sign beta，6=sign beta without BCH hard-decode
static const std::vector<int> kEarlyStopActionModeList = {}; // 按 tile 覆盖动作模式；空表示全部沿用 kEarlyStopActionMode
static constexpr int kEarlyStopBindGroupSize = 1;     // 条件1专用：多少个 row 绑定为一组；1=逐 row，4=整组都通过才 early-stop
static const std::vector<int> kEarlyStopBindGroupSizeList = {}; // 按 tile 覆盖绑定组大小；空表示全部沿用 kEarlyStopBindGroupSize
static constexpr bool kEarlyStopCondV1RequireBch = true; // 条件1里是否要求 BCH syndrome 为 0
static constexpr bool kEarlyStopCondV1RequireOverall = true; // 条件1里是否要求 overall parity 一致
static constexpr float kEarlyStopV2LlrAbsThreshold = 0.5f; // 条件2里把 bit 视为“不可靠”的 |LLR| 阈值
static constexpr int kEarlyStopV2MaxUnreliableBits = 8; // 条件2里允许的不可靠 bit 数上限
static constexpr bool kEarlyStopCondV2IncludeOverall = true; // 条件2统计不可靠 bit 时是否把 overall bit 算进去
static const std::vector<float> kEarlyStopActionBetaExplicit = { 99.857143,99.179301,99.253626,99.119585,99.434408,99.000000}; // 每个 tile 的 early-stop 动作 beta 显式列表；空表示沿用 beta 配置
static constexpr float kEarlyStopActionResidualDivisor = 1.0f; // 动作2里 residual 的除数
static constexpr float kEarlyStopActionHardLlrMag = 1.0f; // 动作3里硬解成功后输出的固定 |LLR| 幅度

// alpha / beta / shared SISO 预算
static const std::vector<float> kAlphaExplicit = {
    0.428571f, 0.447738f, 0.482782f, 0.528162f, 0.581902f, 0.642857f}; // 每个 tile 的 alpha 显式列表
static const std::vector<float> kBetaExplicit = {
    2.857143f, 6.179301f, 12.253626f, 20.119585f, 29.434408f, 40.0f}; // 每个 tile 的 Chase/fallback beta 显式列表
static const std::vector<int> kSisoActiveList = {64, 64, 64, 32, 16, 8}; // 每个 tile 在 shared 64-row 域里可参与 SISO 的预算

// MUX 调度参数
static constexpr int kMuxGroupG = 1;                  // MUX 分组数：1=全局池化，>1=按组平均切预算
static constexpr int kMuxSchedulingMode = 0;          // MUX 调度模式：0=legacy 顺序裁剪，1=按 early-stop 细节排序
static constexpr int kMuxPriorityRule = 0;            // MUX 优先级规则：0=harder_first，1=near_threshold_first
static constexpr bool kMuxEnableReconfig = false;     // true=启用 staged reconfig 调度，false=直接预算裁剪
static const std::vector<newcode::mux::MuxEdge> kMuxExtraBypassEdges =
    app_mux::bypass_edges_for_scheme(1);              // reconfig 模式使用的额外旁路边集合；当前 scheme 1

// Hybrid 参数
static constexpr bool kHybridEnable = true;          // true=启用 hybrid prepass，false=全部保留到 soft path
static const std::vector<int> kHybridEnableList = {0,0,0,1,1,1}; // 按 tile 覆盖 hybrid 开关：0=关，非 0=开；空表示全部沿用 kHybridEnable
static constexpr float kHybridHardLlrMag = 99.0f;     // hybrid hard-finish 默认输出 |LLR| 幅度
static const std::vector<float> kHybridHardLlrMagList = {}; // 按 tile 覆盖 hybrid hard-finish |LLR| 幅度；空表示全部沿用 kHybridHardLlrMag
static constexpr newcode::HybridClassifierMode kHybridClassifierMode =
    newcode::HybridClassifierMode::FriendS1S3WithS0Classifier;  // hybrid 分类器模式：LegacyHardDecode / RepoFastClassifier / FriendS1S3Classifier / FriendS1S3WithS0Classifier
static constexpr newcode::HybridSisoBackfillMode kHybridSisoBackfillMode =
    newcode::HybridSisoBackfillMode::OneAndTwoErrorPriority;        // hybrid SISO 回填模式：Disabled / TwoErrorOnly / OneAndTwoErrorPriority
static constexpr bool kHybridNormalizeSoftOnly = false; // true=只归一化 soft rows，false=对所有 produced rows 保持兼容行为

}  // namespace

void print_quant_stats(
    const char* label,
    const newcode::two_stream_shared::StreamQuantizationStats& stats) {
  if (stats.total_values == 0) {
    std::cout << "[OBS] " << label << " quantization: bypassed\n";
    return;
  }

  std::cout << "[OBS] " << label << " quantization: clip=" << stats.quant_clip
            << ", saturated=" << stats.saturated_values << "/"
            << stats.total_values << " (ratio=" << stats.saturation_ratio
            << ")\n";
}

void print_shared_core_stats(
    const newcode::two_stream_shared::SharedCoreAggregateStats& stats) {
  std::cout << "[OBS] shared core: batches=" << stats.total_batches
            << ", rows=" << stats.total_rows
            << ", produced=" << stats.produced_rows
            << ", failed=" << stats.failed_rows << "\n";

  if (stats.produced_rows_per_stream.size() > 0 ||
      stats.failed_rows_per_stream.size() > 0) {
    const std::size_t produced_a =
        stats.produced_rows_per_stream.size() > 0
            ? stats.produced_rows_per_stream[0]
            : 0;
    const std::size_t failed_a =
        stats.failed_rows_per_stream.size() > 0
            ? stats.failed_rows_per_stream[0]
            : 0;
    std::cout << "[OBS] shared core stream A: produced=" << produced_a
              << ", failed=" << failed_a << "\n";
  }

  if (stats.produced_rows_per_stream.size() > 1 ||
      stats.failed_rows_per_stream.size() > 1) {
    const std::size_t produced_b =
        stats.produced_rows_per_stream.size() > 1
            ? stats.produced_rows_per_stream[1]
            : 0;
    const std::size_t failed_b =
        stats.failed_rows_per_stream.size() > 1
            ? stats.failed_rows_per_stream[1]
            : 0;
    std::cout << "[OBS] shared core stream B: produced=" << produced_b
              << ", failed=" << failed_b << "\n";
  }
}

int main() {
  ofec_single::Config base_cfg{};
  base_cfg.label = kLabel;
  base_cfg.ebn0_db = kEbN0Db;
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

  std::ofstream null_file;
  io::DualWriter log(null_file);
  auto params = ofec_single::detail::build_params(base_cfg, log);
  if (!params.has_value()) {
    std::cerr << "[ERROR] failed to build params for two-stream shared run\n";
    return 2;
  }
  params->NUM_INFO_BITS = kNumInfoBits;

  newcode::PipelineConfig pipeline = ofec_single::detail::build_pipeline_config(base_cfg);
  pipeline.quiet = false;

  newcode::two_stream_shared::Config config{
      .params = *params,
      .pipeline = pipeline,
      .stream_a =
          {
              .label = std::string(kLabel) + "_A",
              .bitgen_seed = kBitgenSeedA,
              .channel_seed = kChannelSeedA,
              .ebn0_db = kEbN0Db,
          },
      .stream_b =
          {
              .label = std::string(kLabel) + "_B",
              .bitgen_seed = kBitgenSeedB,
              .channel_seed = kChannelSeedB,
              .ebn0_db = kEbN0Db,
          },
  };

  auto result = newcode::two_stream_shared::run_two_stream_shared(config);

  std::cout << std::setprecision(8);
  std::cout << "[OBS] shared quant clip="
            << result.observability.shared_quant_clip << "\n";
  print_quant_stats("stream A", result.observability.stream_a_quantization);
  print_quant_stats("stream B", result.observability.stream_b_quantization);
  print_shared_core_stats(result.observability.shared_core);

  std::cout << "[SUMMARY] Stream A post-FEC BER=" << result.stream_a.post_fec.ber
            << " (" << result.stream_a.post_fec.errors << "/"
            << result.stream_a.post_fec.total << ")\n";
  std::cout << "[SUMMARY] Stream B post-FEC BER=" << result.stream_b.post_fec.ber
            << " (" << result.stream_b.post_fec.errors << "/"
            << result.stream_b.post_fec.total << ")\n";
  return 0;
}
