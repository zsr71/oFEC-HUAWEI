#include <vector>

#include "mux_bypass_edges.hpp"
#include "newcode/ofec_sweep_runner.hpp"

// ======== 用户可调参数区域 ========

// 发端参数：交织 / 比特源 / 调制入口
static constexpr const char* kInterleaverName               = "identity"; // 交织器名称，identity 表示不改变顺序
static constexpr unsigned    kBitsPerSymbol                 = 1;          // 每个调制符号携带的比特数：1=BPSK，偶数=QAM
static constexpr bool        kGenerateRandomBits            = true;       // true=随机信息比特，false=全 0 比特
static constexpr int         kBitgenSeed                    = 20260319; // 顶层固定 bitgen seed；本 app 中优先级最高，会覆盖随机 SeedCount 路径
static constexpr int         kBitgenSeedCount               = 1;          // 自动生成的 bitgen seed 数量

// 信道参数：噪声强度 / 信道随机性
static constexpr float       kEbN0Start                     = 3.07f;      // 扫描起始 Eb/N0
static constexpr float       kEbN0End                       = 3.37f;      // 扫描结束 Eb/N0
static constexpr int         kEbN0Points                    = 30;          // Eb/N0 采样点数
static constexpr int         kChannelSeed                   = 3192026;  // 顶层固定 channel seed；本 app 中优先级最高，会覆盖随机 SeedCount 路径
static constexpr int         kChannelSeedCount              = 1;          // 自动生成的 channel seed 数量

// 量化参数：只影响 LLR 量化口径
static constexpr std::size_t kLlrBits                       = 6;          // LLR 位宽：16=浮点，2~15=qfloat
static constexpr float       kQuantClipRatio                = 0.5f;       // 动态 clip 比例，0 表示禁用

// 解码参数：decoder 选择、Chase 参数、alpha/beta、MUX 调度
// 解码器名称：
// chase_baseline=逐 bit 搜索 Cplus/Cminus 的基线 Chase
// chase_topk_pruned=Top-K 裁剪版（pruned=裁剪，只保留前 K 个 good 候选）
// chase_global_pair=全局固定一对 best/second 来算所有 bit 的可靠度
// chase_group_minima=分组组内最优版（minima=每组里度量最小/score 最大的代表）
// chase_overall_parity_search=把 overall parity 也纳入搜索的变体
static constexpr const char* kDecoderName                   = "chase_baseline";
static const std::vector<const char*> kDecoderNameCandidates = {};         // decoder_name 扫描候选，空表示只跑固定解码器；非空时会优先按这里展开多个解码方法
static constexpr bool        kNormalizeExtrinsic            = false;      // 是否对 extrinsic 做归一化
static constexpr bool        kNormalizeKnownPrefixTail      = false;      // 是否对 known-prefix 之后的尾部做归一化
static const std::vector<int> kChaseLCandidates             = {6};        // Chase L 候选列表
static constexpr int         kChaseNTest                    = 64;         // Chase 测试序列数量固定值
static const std::vector<int> kChaseNTestCandidates         = {};         // Chase 测试序列数量扫描候选，空表示沿用固定值

static constexpr int         kChaseTopkKeep                 = 8;          // chase_topk_pruned 的默认 Top-K 保留数；只对“裁剪版”解码器生效
static const std::vector<int> kChaseTopkKeepCandidates      = {};         // chase_topk_pruned 的 Top-K 扫描候选，空表示沿用固定值

static constexpr int         kChaseGroupMinimaBits          = 4;          // chase_group_minima 按前多少个 test-pattern 位分组；取 3 时对应 2^3=8 个组
static const std::vector<int> kChaseGroupMinimaBitsCandidates = {4};       // chase_group_minima 分组位数扫描候选，空表示沿用固定值

static const std::vector<int> kSisoActiveList               = {32, 32, 32, 32}; // 每个 tile 的 SISO 行数预算
static constexpr int         kMuxGroupG                     = 1;          // MUX 分组粒度，1 表示全局池化
static constexpr int         kMuxSchedulingMode             = 0;          // MUX 调度模式：0=legacy，1=按 early-stop 细节排序
static const std::vector<int> kMuxSchedulingModeCandidates  = {};         // MUX 调度模式扫描候选，空表示沿用固定值
static constexpr int         kMuxPriorityRule               = 0;          // 新 MUX 的优先级规则：0=更差优先，1=更接近通过优先
static const std::vector<int> kMuxPriorityRuleCandidates    = {};         // 新 MUX 优先级规则扫描候选，空表示沿用固定值
static constexpr bool        kMuxEnableReconfig             = false;      // 是否启用重配置版 MUX 调度
static constexpr int         kMuxBypassScheme               = 3;          // 旁路边集合方案编号：1=scheme1，2=scheme2
static const std::vector<ofec_sweep::ExplicitAlphaBetaPattern> kExplicitAlphaBetaSets = { // 显式给出 alpha/beta/early-stop beta 列表的方案集合
 {"custom_label",
  {0.342857,0.387439,0.435806,0.485714},
  {8.571428,10.037715,16.865997,31.428572},
  {}},
};

// 早停参数：总开关 -> 条件 -> 条件细参 -> 动作 -> 动作细参
static constexpr bool        kEnableEarlyStop               = false;      // 是否启用 early-stop 总开关
static const std::vector<int> kEarlyStopEnableList          = {};         // 按 tile 覆盖 early-stop 总开关：0=关，非 0=开；空表示所有 tile 沿用 kEnableEarlyStop
static constexpr int         kEarlyStopConditionMode        = 1;          // 早停条件编号：1=v1，2=v2
static const std::vector<int> kEarlyStopConditionModeList   = {};         // 按 tile 覆盖早停条件模式；空表示所有 tile 沿用 kEarlyStopConditionMode
static const std::vector<int> kEarlyStopConditionCandidates = {};    // 早停条件候选列表
static constexpr int         kEarlyStopActionMode           = 1;          // 早停命中后的动作编号：1=sign beta，2=residual only，3=硬解成功后直接输出 ±hard_mag，4=sign beta 后预除 alpha，5=residual 预除 alpha 后再加 sign beta，6=直接输出 ±early-stop beta
static const std::vector<int> kEarlyStopActionModeList      = {};         // 按 tile 覆盖早停动作模式；空表示所有 tile 沿用 kEarlyStopActionMode
static constexpr int         kEarlyStopBindGroupSize        = 1;          // 条件1专用的组绑定大小：1=逐 row；4=每 4 个 row 都通过才整体 early-stop
static const std::vector<int> kEarlyStopBindGroupSizeList   = {};         // 按 tile 覆盖条件1绑定组大小；空表示所有 tile 沿用 kEarlyStopBindGroupSize
static const std::vector<int> kEarlyStopActionCandidates    = {}; // 早停动作候选列表

static constexpr bool        kEarlyStopCondV1RequireBch     = true;       // 条件1里是否要求 BCH syndrome 为 0
static const std::vector<bool> kEarlyStopCondV1RequireBchCandidates = {}; // 条件1里是否要求 BCH 的候选，空表示沿用固定值
static constexpr bool        kEarlyStopCondV1RequireOverall = true;       // 条件1里是否要求 overall parity 一致
static const std::vector<bool> kEarlyStopCondV1RequireOverallCandidates = {}; // 条件1里是否要求 overall 的候选，空表示沿用固定值

static constexpr float       kEarlyStopV2LlrAbsThreshold    = 0.5f;       // 条件2中判不可靠 bit 的 |LLR| 阈值
static const std::vector<float> kEarlyStopV2LlrAbsThresholdCandidates = {}; // 条件2阈值候选，空表示沿用固定值
static constexpr int         kEarlyStopV2MaxUnreliableBits  = 8;          // 条件2允许的不可靠 bit 数上限
static const std::vector<int> kEarlyStopV2MaxUnreliableBitsCandidates = {}; // 条件2不可靠 bit 上限候选，空表示沿用固定值
static constexpr bool        kEarlyStopCondV2IncludeOverall = true;       // 条件2统计时是否把 overall bit 纳入

static constexpr float       kEarlyStopActionResidualDivisor = 0.4f;      // 动作2里 residual 的除数
static constexpr float       kEarlyStopActionHardLlrMag     = 1.0f;       // 动作3里硬解成功后输出的固定 |LLR| 幅度

static const std::vector<float> kEarlyStopActionBetaStartCandidates = {}; // early-stop 动作 beta 起点候选，空表示不单独扫描
static const std::vector<float> kEarlyStopActionBetaStepCandidates = {};  // early-stop 动作 beta 步进候选，空表示不单独扫描

static const std::vector<float> kEarlyStopActionHardLlrMagCandidates = {}; // 动作3的 |hard_llr_mag| 候选，空表示沿用固定值

// Debug 参数：控制日志与 decoder trace
static constexpr bool        kQuietConsole                  = false;      // true=减少控制台打印，false=保留详细日志
static constexpr bool        kDecoderTraceEnable            = false;      // 是否启用 decoder trace
static constexpr long        kDecoderTraceRow               = -1;         // 需要跟踪的全局行号，-1 表示不指定
static constexpr long        kDecoderTraceCol               = -1;         // 需要跟踪的全局列号，-1 表示不指定
static constexpr bool        kDecoderTraceLogRead           = false;      // 是否打印 tile 读取映射
static constexpr bool        kDecoderTraceLogWrite          = false;      // 是否打印 tile 写回映射
static constexpr bool        kDecoderTraceLogMismatch       = false;      // 是否打印写回冲突告警

// ==================================

int main() {
  const auto& selected_mux_bypass_edges =
      app_mux::bypass_edges_for_scheme(kMuxBypassScheme);
  ofec_sweep::SweepParameterConfig config;

  // 基础链路配置
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

  // 早停固定配置
  config.enable_early_stop = kEnableEarlyStop;
  config.early_stop_enable_list = kEarlyStopEnableList;
  config.early_stop_condition_mode = kEarlyStopConditionMode;
  config.early_stop_condition_mode_list = kEarlyStopConditionModeList;
  config.early_stop_action_mode = kEarlyStopActionMode;
  config.early_stop_action_mode_list = kEarlyStopActionModeList;
  config.early_stop_bind_group_size = kEarlyStopBindGroupSize;
  config.early_stop_bind_group_size_list = kEarlyStopBindGroupSizeList;
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

  // 扫描候选：early-stop 动作 beta
  config.early_stop_action_beta_start_candidates = kEarlyStopActionBetaStartCandidates;
  config.early_stop_action_beta_step_candidates = kEarlyStopActionBetaStepCandidates;
  config.early_stop_action_hard_llr_mag_candidates =
      kEarlyStopActionHardLlrMagCandidates;

  // 扫描候选：early-stop 条件 / 动作模式
  config.early_stop_condition_candidates = kEarlyStopConditionCandidates;
  config.early_stop_action_candidates = kEarlyStopActionCandidates;

  // 扫描候选：条件 1 / 条件 2 细参数
  config.early_stop_cond_v1_require_bch_candidates =
      kEarlyStopCondV1RequireBchCandidates;
  config.early_stop_cond_v1_require_overall_candidates =
      kEarlyStopCondV1RequireOverallCandidates;
  config.early_stop_v2_llr_abs_threshold_candidates =
      kEarlyStopV2LlrAbsThresholdCandidates;
  config.early_stop_v2_max_unreliable_bits_candidates =
      kEarlyStopV2MaxUnreliableBitsCandidates;

  // MUX / 调度配置
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

  // Monte Carlo seed 配置
  config.base_params.BITGEN_SEED = kBitgenSeed;
  config.base_params.CHANNEL_SEED = kChannelSeed;
  config.bitgen_seed_candidates = {kBitgenSeed};
  config.channel_seed_candidates = {kChannelSeed};
  config.bitgen_seed_count = kBitgenSeedCount;
  config.channel_seed_count = kChannelSeedCount;

  // Eb/N0 扫描区间
  config.ebn0_start = kEbN0Start;
  config.ebn0_end = kEbN0End;
  config.ebn0_points = kEbN0Points;

  // 显式列表方案
  config.explicit_patterns = kExplicitAlphaBetaSets;

  // Debug trace 基础配置
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

  // 直接落到 Params 的基础字段
  config.base_params.BITGEN_RANDOM_BITS = kGenerateRandomBits;
  config.base_params.NORMALIZE_KNOWN_PREFIX_TAIL = kNormalizeKnownPrefixTail;
  config.base_params.ENABLE_EARLY_STOP = kEnableEarlyStop;
  config.base_params.LLR_CLIP_RATIO = kQuantClipRatio;
  config.base_params.LLR_BITS = kLlrBits;
  config.base_params.SISO_ACTIVE_LIST = kSisoActiveList;

  return ofec_sweep::run_sweep(config);
}
