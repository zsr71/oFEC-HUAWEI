#include <vector>

#include "mux_bypass_edges.hpp"
#include "newcode/utils/linspace.hpp"
#include "newcode/ofec_sweep_runner.hpp"

// ======== 用户可调参数区域 ========

// 基础运行入口
static constexpr const char* kInterleaverName               = "identity"; // 交织器名称，identity 表示不改变顺序
static constexpr const char* kDecoderName                   = "plain";    // 解码器名称，当前 sweep 默认扫 plain

// 信道与量化口径
static constexpr unsigned    kBitsPerSymbol                 = 1;          // 每个调制符号携带的比特数：1=BPSK，偶数=QAM
static constexpr bool        kNormalizeExtrinsic            = false;      // 是否对 extrinsic 做归一化
static constexpr bool        kGenerateRandomBits            = true;       // true=随机信息比特，false=全 0 比特
static constexpr bool        kNormalizeKnownPrefixTail      = false;      // 是否对 known-prefix 之后的尾部做归一化
static constexpr float       kQuantClipRatio                = 0.5f;       // 动态 clip 比例，0 表示禁用
static constexpr std::size_t kLlrBits                       = 6;          // LLR 位宽：16=浮点，2~15=qfloat
static constexpr bool        kQuietConsole                  = false;      // true=减少控制台打印，false=保留详细日志

// 早停固定配置：决定“停不停”和“停了以后怎么办”
static constexpr bool        kEnableEarlyStop               = true;       // 是否启用 early-stop 总开关
static constexpr int         kEarlyStopConditionMode        = 1;          // 早停条件编号：1=v1，2=v2
static constexpr int         kEarlyStopActionMode           = 1;          // 早停命中后的动作编号：1=sign beta，2=residual only

// 条件 1（v1）参数
static constexpr bool        kEarlyStopCondV1RequireBch     = true;       // 条件1里是否要求 BCH syndrome 为 0
static constexpr bool        kEarlyStopCondV1RequireOverall = true;       // 条件1里是否要求 overall parity 一致

// 条件 2（v2）参数
static constexpr float       kEarlyStopV2LlrAbsThreshold    = 0.5f;       // 条件2中判不可靠 bit 的 |LLR| 阈值
static constexpr int         kEarlyStopV2MaxUnreliableBits  = 8;          // 条件2允许的不可靠 bit 数上限
static constexpr bool        kEarlyStopCondV2IncludeOverall = true;       // 条件2统计时是否把 overall bit 纳入

// 早停动作参数
static constexpr float       kEarlyStopActionResidualDivisor = 1.0f;      // 动作2里 residual 的除数
static constexpr float       kEarlyStopActionHardLlrMag     = 1.0f;       // 预留给硬输出类早停动作的 LLR 幅度

// MUX / 调度参数
static const std::vector<int> kSisoActiveList               = {32, 32, 32, 16}; // 每个 tile 的 SISO 行数预算
static constexpr int         kMuxGroupG                     = 1;          // MUX 分组粒度，1 表示全局池化
static constexpr bool        kMuxEnableReconfig             = false;      // 是否启用重配置版 MUX 调度
static constexpr int         kMuxBypassScheme               = 1;          // 旁路边集合方案编号：1=scheme1，2=scheme2

// Chase / 外信息扫描候选
static const std::vector<float> kAlphaStartCandidates                = utils::linspace(0.0f, 0.2f, 2); // alpha 起点候选
static const std::vector<float> kAlphaStepCandidates                 = utils::linspace(0.0f, 0.2f, 2); // alpha 步进候选
static const std::vector<float> kBetaStartCandidates                 = utils::linspace(0.0f, 0.2f, 2); // Chase beta 起点候选
static const std::vector<float> kBetaStepCandidates                  = utils::linspace(0.0f, 0.2f, 2); // Chase beta 步进候选
static const std::vector<int>   kChaseLCandidates                    = {6};                             // Chase L 候选列表

// 早停条件 / 动作模式扫描候选
static const std::vector<int>   kEarlyStopConditionCandidates        = {1, 2};                          // 早停条件候选列表
static const std::vector<int>   kEarlyStopActionCandidates           = {1, 2};                          // 早停动作候选列表

// 条件 1（v1）扫描候选
static const std::vector<bool>  kEarlyStopCondV1RequireBchCandidates = {};                              // 条件1里是否要求 BCH 的候选，空表示沿用固定值
static const std::vector<bool>  kEarlyStopCondV1RequireOverallCandidates = {};                          // 条件1里是否要求 overall 的候选，空表示沿用固定值

// 条件 2（v2）扫描候选
static const std::vector<float> kEarlyStopV2LlrAbsThresholdCandidates = {};                             // 条件2阈值候选，空表示沿用固定值
static const std::vector<int>   kEarlyStopV2MaxUnreliableBitsCandidates = {};                           // 条件2不可靠 bit 上限候选，空表示沿用固定值

// 早停动作扫描候选
static const std::vector<float> kEarlyStopActionBetaStartCandidates  = {};                              // early-stop 动作 beta 起点候选，空表示不单独扫描
static const std::vector<float> kEarlyStopActionBetaStepCandidates   = {};                              // early-stop 动作 beta 步进候选，空表示不单独扫描

// Monte Carlo 随机种子设置（为空则自动生成）
static constexpr int kBitgenSeedCount  = 1; // 自动生成的 bitgen seed 数量
static constexpr int kChannelSeedCount = 1; // 自动生成的 channel seed 数量

// Eb/N0 扫描区间
static constexpr float kEbN0Start  = 3.07f; // 扫描起始 Eb/N0
static constexpr float kEbN0End    = 3.07f; // 扫描结束 Eb/N0
static constexpr int   kEbN0Points = 1;     // Eb/N0 采样点数

// 显式列表方案（可选）：直接给定每个 tile 的 alpha / beta / early-stop beta
static const std::vector<ofec_sweep::ExplicitAlphaBetaPattern> kExplicitAlphaBetaSets = { // 显式给出 alpha/beta/early-stop beta 列表的方案集合
 {"custom_label",
  {0.342857,0.387439,0.435806,0.485714},
  {8.571428,10.037715,16.865997,31.428572},
  {}},
};

// Decoder 调试跟踪配置
static constexpr bool kDecoderTraceEnable      = false; // 是否启用 decoder trace
static constexpr long kDecoderTraceRow         = -1;    // 需要跟踪的全局行号，-1 表示不指定
static constexpr long kDecoderTraceCol         = -1;    // 需要跟踪的全局列号，-1 表示不指定
static constexpr bool kDecoderTraceLogRead     = false; // 是否打印 tile 读取映射
static constexpr bool kDecoderTraceLogWrite    = false; // 是否打印 tile 写回映射
static constexpr bool kDecoderTraceLogMismatch = false; // 是否打印写回冲突告警

// ==================================

int main() {
  const auto& selected_mux_bypass_edges =
      app_mux::bypass_edges_for_scheme(kMuxBypassScheme);
  ofec_sweep::SweepParameterConfig config;

  // 基础链路配置
  config.interleaver_name = kInterleaverName;
  config.decoder_name = kDecoderName;
  config.bits_per_symbol = kBitsPerSymbol;
  config.normalize_extrinsic = kNormalizeExtrinsic;
  config.generate_random_bits = kGenerateRandomBits;
  config.normalize_known_prefix_tail = kNormalizeKnownPrefixTail;
  config.quant_clip_ratio = kQuantClipRatio;
  config.quiet_pipeline = kQuietConsole;
  config.quiet_logs = kQuietConsole;

  // 早停固定配置
  config.enable_early_stop = kEnableEarlyStop;
  config.early_stop_condition_mode = kEarlyStopConditionMode;
  config.early_stop_action_mode = kEarlyStopActionMode;
  config.early_stop_cond_v1_require_bch = kEarlyStopCondV1RequireBch;
  config.early_stop_cond_v1_require_overall = kEarlyStopCondV1RequireOverall;
  config.early_stop_v2_llr_abs_threshold = kEarlyStopV2LlrAbsThreshold;
  config.early_stop_v2_max_unreliable_bits = kEarlyStopV2MaxUnreliableBits;
  config.early_stop_cond_v2_include_overall = kEarlyStopCondV2IncludeOverall;
  config.early_stop_action_residual_divisor = kEarlyStopActionResidualDivisor;
  config.early_stop_action_hard_llr_mag = kEarlyStopActionHardLlrMag;

  // 扫描候选：Chase / 外信息
  config.alpha_start_candidates = kAlphaStartCandidates;
  config.alpha_step_candidates = kAlphaStepCandidates;
  config.beta_start_candidates = kBetaStartCandidates;
  config.beta_step_candidates = kBetaStepCandidates;

  // 扫描候选：early-stop 动作 beta
  config.early_stop_action_beta_start_candidates = kEarlyStopActionBetaStartCandidates;
  config.early_stop_action_beta_step_candidates = kEarlyStopActionBetaStepCandidates;

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
  config.mux_enable_reconfig = kMuxEnableReconfig;
  config.mux_extra_bypass_edges = selected_mux_bypass_edges;
  config.chase_l_candidates = kChaseLCandidates;

  // Monte Carlo seed 配置
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
