#include <vector>

#include "newcode/ofec_single_runner.hpp"

// ======== 用户可改区域 ========
// 只需要改这里的常量/列表即可完成一次“单次调试运行”的配置

//发射端参数
static constexpr const char* kLabel             = "debug_L6";   // 运行标签：日志/输出文件标识
static constexpr int         kChaseL_override   = 6;            // Chase L，-1 表示使用默认
static constexpr int         kBitgenSeed        = 115948;    // 比特生成随机种子
static constexpr bool        kGenerateRandomBits = false;             // true=随机比特，false=全 0

//信道相关参数
static constexpr float       kEbN0_db           = 3.17f;        // 信道 Eb/N0 (dB)
static constexpr int         kChannelSeed       = 1481598;   // 信道噪声随机种子
static constexpr unsigned    kBitsPerSymbol     = 1;            // 每符号比特数：1=BPSK，偶数=QAM

//量化相关参数
static constexpr std::size_t kLlrBits =16;                           // LLR 位宽：16=浮点，2~15=qfloat
static constexpr float       kQuantClipRatio = 0.0f;                 // 动态裁剪比例，0=禁用

//解码相关参数
static constexpr const char* kInterleaverName = "identity";          // 交织器名称
static constexpr const char* kDecoderName     = "plain";             // 解码器名称
static constexpr bool        kNormalizeExtrinsic = false;       // 是否对外信息做归一化
static constexpr bool        kNormalizeKnownPrefixTail = false;      // known_prefix 外是否归一化

  // 方式 A：统一填充值（长度自动取 Params::TILES_PER_WIN）
  static constexpr float kAlpha_fill = 1.0f;   // extrinsic 缩放系数 α（统一填充）
  static constexpr float kBeta_fill  = 0.40f;  // fallback 可靠度 β（统一填充）

  // 方式 B：显式列表（若非空，将覆盖填充值；长度必须等于 TILES_PER_WIN）
  static const std::vector<float> kAlpha_explicit = {
    0.2,0.4,0.8,1
    //0.4f,1.0f
    //0.6f
    //0.606316
  };
  static const std::vector<float> kBeta_explicit = {
    0.2,0.4,0.6,0.8
    //1.4f
    //24.0f,1.0f
    //1.4f
    //1.794872*31
  };


//调试相关参数

//llr导出相关
  static constexpr bool        kDumpQuantizedLlr = true;               // 是否导出量化后 LLR
  static constexpr const char* kQuantizedLlrPath = "data/llr/quantized_llr.txt"; // 量化 LLR 输出路径
  static constexpr bool        kDumpWorkLlr = true;                    // 是否保存最后的的 work_llr
  static constexpr const char* kWorkLlrPath = "data/llr/work_llr.txt"; // work_llr 输出路径

namespace {
constexpr long BitIndexToRow(long bit_index) { return bit_index / 111 + 352; }
constexpr long BitIndexToCol(long bit_index) { return bit_index % 111; }
struct TraceBitSpec {
  long bit_index;
  const char* label;
};
constexpr TraceBitSpec kTraceBitSpecs[] = {
    {704369, "bit704369"},
    {727555, "bit727555"},
    {981571, "bit981571"},
    {1369219, "bit1369219"},
    {2403048, "bit2403048"},
    {2634830, "bit2634830"},
};
} // namespace

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
static constexpr bool kDecoderTraceEnable        = true;
static constexpr bool kDecoderTraceLogRead       = true;
static constexpr bool kDecoderTraceLogWrite      = true;
static constexpr bool kDecoderTraceLogMismatch   = true;
static constexpr bool kDecoderTraceLogChase      = true;
static constexpr bool kDecoderTraceDumpCsv       = true;
static constexpr const char* kDecoderTraceCsvDir = "data/chase_csv";


// =============================

int main() {
  ofec_single::Config config{
    .label = kLabel,
    .ebn0_db = kEbN0_db,
    .chaseL_override = kChaseL_override,
    .normalize_extrinsic = kNormalizeExtrinsic,
    .bits_per_symbol = kBitsPerSymbol,
    .bitgen_seed = kBitgenSeed,
    .channel_seed = kChannelSeed,
    .alpha_fill = kAlpha_fill,
    .beta_fill = kBeta_fill,
    .alpha_explicit = kAlpha_explicit,
    .beta_explicit = kBeta_explicit,
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
      .targets = kDecoderTraceTargets,
    },
  };
  return ofec_single::run_ofec_single(config);
}
