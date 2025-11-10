#include <vector>

#include "newcode/ofec_single_runner.hpp"

// ======== 用户可改区域 ========
// 只需要改这里的常量/列表即可完成一次“单次调试运行”的配置
static constexpr const char* kLabel             = "debug_L6";
static constexpr float       kEbN0_db           = 3.52f;
static constexpr int         kChaseL_override   = 6;    // 设为 -1 则沿用 Params 默认
static constexpr bool        kNormalizeExtrinsic = false;
static constexpr unsigned    kBitsPerSymbol     = 1;    // 设为 1 使用 BPSK，>=2 且偶数使用 QAM
static constexpr int         kBitgenSeed        = 885968131;
static constexpr int         kChannelSeed       = 1278801575;

// 方式 A：统一填充值（长度自动取 Params::TILES_PER_WIN）
static constexpr float kAlpha_fill = 1.0f;
static constexpr float kBeta_fill  = 0.40f;

// 方式 B：显式列表（若非空，将覆盖填充值；长度必须等于 TILES_PER_WIN）
static const std::vector<float> kAlpha_explicit = {
  //0.3f,0.45f,0.60f,0.9f,0.5f
  //0.3f,0.45f
  0.4f,0.5f
};
static const std::vector<float> kBeta_explicit = {
  //0.2f,0.225f,0.30f,0.4f,0.0f
  //0.2f,0.225f
  1.3f,2.3f
};

// static const std::vector<float> kAlpha_explicit = {
//   0.0f
// };
// static const std::vector<float> kBeta_explicit = {
//   0.0f
// };

static constexpr const char* kInterleaverName = "identity";
static constexpr const char* kDecoderName     = "plain";
static constexpr bool        kGenerateRandomBits = true;
static constexpr bool        kNormalizeKnownPrefixTail = false;

// Decoder 调试跟踪配置
static constexpr bool kDecoderTraceEnable        = false;
static constexpr long kDecoderTraceRow           = -1;
static constexpr long kDecoderTraceCol           = -1;
static constexpr bool kDecoderTraceLogRead       = false;
static constexpr bool kDecoderTraceLogWrite      = false;
static constexpr bool kDecoderTraceLogMismatch   = false;
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
    .debug_trace = newcode::Params::DebugTraceConfig{
      .enable = kDecoderTraceEnable,
      .log_read_mapping = kDecoderTraceLogRead,
      .log_write_mapping = kDecoderTraceLogWrite,
      .log_mismatch = kDecoderTraceLogMismatch,
      .row = kDecoderTraceRow,
      .col = kDecoderTraceCol,
    },
  };
  return ofec_single::run_ofec_single(config);
}
