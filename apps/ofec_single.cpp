#include <vector>

#include "mux_bypass_edges.hpp"
#include "newcode/ofec_single_runner.hpp"

// ======== 用户可改区域 ========
// 只需要改这里的常量/列表即可完成一次“单次调试运行”的配置

// 发射端参数
static constexpr const char* kLabel              = "debug_L6";    // 运行标签：日志名、输出文件名前缀都会带这个名字
static constexpr int         kChaseL_override    = 6;             // Chase L，-1 表示使用 Params 里的默认值
static constexpr int         kBitgenSeed         = 1521867291;    // 比特生成随机种子，固定后可复现实验
static constexpr bool        kGenerateRandomBits = true;          // true=发送随机信息比特，false=发送全 0 比特

// 信道参数
static constexpr float       kEbN0_db                      = 3.57f;   // 信道 Eb/N0，单位 dB
static constexpr int         kChannelSeed                  = 998258255; // 信道噪声随机种子，固定后可复现实验
static constexpr unsigned    kBitsPerSymbol                = 1;       // 每个调制符号携带的比特数：1=BPSK，偶数=QAM

//早停参数
static constexpr bool        kEnableEarlyStop              = false;   // true=启用早停，false=完全关闭早停路径
static constexpr int         kEarlyStopConditionMode       = 2;       // 早停条件编号：1=v1，2=v2
static constexpr int         kEarlyStopActionMode          = 1;       // 早停命中后的动作：1=sign beta，2=residual only

static constexpr bool        kEarlyStopCondV1RequireBch    = true;    // 条件1里是否要求 BCH syndrome 为 0
static constexpr bool        kEarlyStopCondV1RequireOverall = true;   // 条件1里是否要求 overall parity 一致

static constexpr float       kEarlyStopV2LlrAbsThreshold   = 26.0f;   // 条件2里把 bit 视为“不可靠”的 |LLR| 阈值
static constexpr int         kEarlyStopV2MaxUnreliableBits = 25;      // 条件2里允许的不可靠 bit 数上限
static constexpr bool        kEarlyStopCondV2IncludeOverall = true;   // 条件2统计不可靠 bit 时是否把 overall bit 算进去

// LLR 量化相关参数
static constexpr std::size_t kLlrBits        = 6;      // LLR 位宽：16=浮点，2~15=qfloat
static constexpr float       kQuantClipRatio = 0.5f;   // 动态裁剪比例，0 表示禁用自适应 clip

// 解码主参数
static constexpr const char* kInterleaverName          = "identity"; // 交织器名称，identity 表示不改变比特顺序
static constexpr const char* kDecoderName              = "plain";    // 解码器名称，当前常用 plain
static constexpr bool        kNormalizeExtrinsic       = false;      // 是否对 Chase 输出的 extrinsic 做归一化
static constexpr bool        kNormalizeKnownPrefixTail = false;      // 是否对 known-prefix 之后的尾部 LLR 做归一化

// 方式 A：统一填充值（长度自动取 Params::TILES_PER_WIN）
static constexpr float kAlpha_fill                     = 1.0f;   // 每个 tile 共用的 extrinsic 缩放系数 alpha
static constexpr float kBeta_fill                      = 0.40f;  // 每个 tile 共用的 Chase/fallback beta
static constexpr float kEarlyStopActionBeta_fill       = 0.40f;  // 每个 tile 共用的 early-stop 动作 beta
static constexpr float kEarlyStopActionResidualDivisor = 1.0f;   // 动作2里 residual 的除数
static constexpr float kEarlyStopActionHardLlrMag      = 1.0f;   // 预留给硬输出类早停动作的 LLR 幅度

// 方式 B：显式列表（若非空，将覆盖填充值；长度必须等于 TILES_PER_WIN）
static const std::vector<float> kAlpha_explicit = {      // 每个 tile 的 alpha 显式列表
   0.342857,0.387439,0.435806,0.485714
};
static const std::vector<float> kBeta_explicit = {       // 每个 tile 的 Chase/fallback beta 显式列表
  8.571428,10.037715,16.865997,31.428572
};
static const std::vector<float> kEarlyStopActionBeta_explicit = { // 每个 tile 的 early-stop 动作 beta 显式列表
  8.571428,10.037715,16.865997,31.428572
};
static const std::vector<int> kSisoActiveList = {32, 32, 32, 32}; // 每个 tile 允许参与 SISO 的行数预算
static constexpr int  kMuxGroupG          = 1;                     // MUX 分组粒度，1 表示全局池化
static constexpr bool kMuxEnableReconfig  = false;                 // true 表示启用重配置版 MUX 调度
static constexpr int  kMuxBypassScheme    = 1;                     // 旁路边集合方案编号：1=scheme1，2=scheme2

// LLR 导出相关
static constexpr bool        kDumpQuantizedLlr = true;                         // 是否导出量化后的信道 LLR
static constexpr const char* kQuantizedLlrPath = "data/llr/quantized_llr.txt"; // 量化 LLR 输出路径
static constexpr bool        kDumpWorkLlr      = true;                         // 是否保存最终累积得到的 work_llr
static constexpr const char* kWorkLlrPath      = "data/llr/work_llr.txt";      // work_llr 输出路径

namespace {
constexpr long BitIndexToRow(long bit_index) { return bit_index / 111 + 352; } // 把 info bit 编号映射回矩阵行号
constexpr long BitIndexToCol(long bit_index) { return bit_index % 111; }        // 把 info bit 编号映射回矩阵列号
struct TraceBitSpec {
  long bit_index;
  const char* label;
};
constexpr TraceBitSpec kTraceBitSpecs[] = { // 需要重点跟踪的目标比特列表
    {1143884, "bit1143884"},
    {1148431, "bit1148431"},
    {5544508, "bit5544508"},
    {5563798, "bit5563798"},
};
} 

static const std::vector<newcode::Params::DebugTraceConfig::TraceTarget>
    kDecoderTraceTargets = []() {
      std::vector<newcode::Params::DebugTraceConfig::TraceTarget> targets;
      for (const auto& spec : kTraceBitSpecs) {
        newcode::Params::DebugTraceConfig::TraceTarget entry;
        entry.row = BitIndexToRow(spec.bit_index);
        entry.col = BitIndexToCol(spec.bit_index);
        entry.bit_index = spec.bit_index;
        if (spec.label) entry.label = spec.label;
        targets.push_back(entry);
      }
      return targets;
    }();

// Decoder 调试跟踪配置
// 用于定位单个比特的映射/写回问题：row/col 为全局 work_llr 坐标（0-based）
// enable 与 row/col 同时满足才会输出日志
// log_read_mapping   输出 tile -> global 的读取坐标展开
// log_write_mapping  输出解码回写对应的坐标及 LLR
// log_mismatch       当窗口内多个 tile 写回同一坐标且值不同时报错提示
// log_chase_detail   追踪对应比特在 Chase 内部（plain/ebchPF）的输入/输出
// dump_chase_csv     每次进入 Chase 时导出 256 码字的 LLR/ω/ML/硬判决到 CSV
// chase_csv_dir      CSV 导出目录（可直接用 Excel 打开）

static constexpr bool        kDecoderTraceEnable      = true;             // 是否启用 decoder trace 总开关
static constexpr bool        kDecoderTraceLogRead     = true;             // 是否打印 tile 读取映射
static constexpr bool        kDecoderTraceLogWrite    = true;             // 是否打印 tile 写回映射
static constexpr bool        kDecoderTraceLogMismatch = true;             // 是否打印同坐标写回不一致告警
static constexpr bool        kDecoderTraceLogChase    = true;             // 是否打印 Chase 内部细节
static constexpr bool        kDecoderTraceDumpCsv     = true;             // 是否把 Chase 细节落成 CSV
static constexpr const char* kDecoderTraceCsvDir      = "data/chase_csv"; // Chase CSV 导出目录


// =============================

int main() {
  const auto& selected_mux_bypass_edges =
      app_mux::bypass_edges_for_scheme(kMuxBypassScheme);
  ofec_single::Config config{
    .label = kLabel,
    .ebn0_db = kEbN0_db,
    .chaseL_override = kChaseL_override,
    .normalize_extrinsic = kNormalizeExtrinsic,
    .bits_per_symbol = kBitsPerSymbol,
    .bitgen_seed = kBitgenSeed,
    .channel_seed = kChannelSeed,
    .enable_early_stop = kEnableEarlyStop,
    .early_stop_condition_mode = kEarlyStopConditionMode,
    .early_stop_action_mode = kEarlyStopActionMode,
    .early_stop_cond_v1_require_bch = kEarlyStopCondV1RequireBch,
    .early_stop_cond_v1_require_overall = kEarlyStopCondV1RequireOverall,
    .early_stop_v2_llr_abs_threshold = kEarlyStopV2LlrAbsThreshold,
    .early_stop_v2_max_unreliable_bits = kEarlyStopV2MaxUnreliableBits,
    .early_stop_cond_v2_include_overall = kEarlyStopCondV2IncludeOverall,
    .alpha_fill = kAlpha_fill,
    .beta_fill = kBeta_fill,
    .alpha_explicit = kAlpha_explicit,
    .beta_explicit = kBeta_explicit,
    .early_stop_action_sign_beta_fill = kEarlyStopActionBeta_fill,
    .early_stop_action_sign_beta_explicit = kEarlyStopActionBeta_explicit,
    .early_stop_action_residual_divisor = kEarlyStopActionResidualDivisor,
    .early_stop_action_hard_llr_mag = kEarlyStopActionHardLlrMag,
    .siso_active_list = kSisoActiveList,
    .mux_group_g = kMuxGroupG,
    .mux_enable_reconfig = kMuxEnableReconfig,
    .mux_extra_bypass_edges = selected_mux_bypass_edges,
    .interleaver_name = kInterleaverName,
    .decoder_name = kDecoderName,
    .generate_random_bits = kGenerateRandomBits,
    .normalize_known_prefix_tail = kNormalizeKnownPrefixTail,
    .quant_clip_ratio = kQuantClipRatio,
    .llr_bits = kLlrBits,
    .dump_quantized_llr = kDumpQuantizedLlr,
    .quantized_llr_output_path = kQuantizedLlrPath,
    .dump_work_llr = kDumpWorkLlr,
    .work_llr_output_path = kWorkLlrPath,
    .debug_trace = newcode::Params::DebugTraceConfig{
      .enable = kDecoderTraceEnable,
      .log_read_mapping = kDecoderTraceLogRead,
      .log_write_mapping = kDecoderTraceLogWrite,
      .log_mismatch = kDecoderTraceLogMismatch,
      .log_chase_detail = kDecoderTraceLogChase,
      .dump_chase_csv = kDecoderTraceDumpCsv,
      .row = -1,
      .col = -1,
      .chase_decoder_row = -1,
      .chase_decoder_col = -1,
      .chase_tile_index = -1,
      .chase_invocation = -1,
      .chase_csv_dir = kDecoderTraceCsvDir,
      .chase_expected_bits = {},
      .targets = kDecoderTraceTargets,
      .active_chase_entries = {},
    },
  };
  return ofec_single::run_ofec_single(config);
}
