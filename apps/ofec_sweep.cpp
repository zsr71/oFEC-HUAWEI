#include <vector>

#include "newcode/linspace.hpp"
#include "newcode/ofec_sweep_runner.hpp"

// ======== 用户可调参数区域 ========
static constexpr const char* kInterleaverName = "identity";
static constexpr const char* kDecoderName = "plain";
static constexpr unsigned    kBitsPerSymbol = 1; // 设为 1 使用 BPSK，>=2 且偶数使用 QAM
static constexpr bool        kNormalizeExtrinsic = false;
static constexpr bool kGenerateRandomBits        = true;
static constexpr bool kNormalizeKnownPrefixTail  = false;
static constexpr float kQuantClipRatio           = 0.5f; // 0 表示禁用动态 clip
static constexpr std::size_t kLlrBits            = 6;
static constexpr bool kQuietConsole              = false;

// Alpha/Beta 扫描候选
static const std::vector<float> kAlphaStartCandidates = newcode::linspace(0.0f, 0.2f, 2);
static const std::vector<float> kAlphaStepCandidates  = newcode::linspace(0.0f, 0.2f, 2);
static const std::vector<float> kBetaStartCandidates  = newcode::linspace(0.0f, 0.2f, 2);
static const std::vector<float> kBetaStepCandidates   = newcode::linspace(0.0f, 0.2f, 2);
static const std::vector<int>   kChaseLCandidates     = {6};

// 随机种子（为空则自动生成 bitgen/channel seeds 数量）
static constexpr int kBitgenSeedCount  = 1;
static constexpr int kChannelSeedCount = 1;

// Eb/N0 设置
static constexpr float kEbN0Start = 3.07f;
static constexpr float kEbN0End   = 3.37f;
static constexpr int   kEbN0Points = 14;

// 显式 alpha/beta 模式（可选）
static const std::vector<ofec_sweep::ExplicitAlphaBetaPattern> kExplicitAlphaBetaSets = {
  {"custom_label", {0.333333,0.416667,0.500000}, {4.444445,15.274477,31.111113}},
};

// Decoder 调试跟踪配置
static constexpr bool kDecoderTraceEnable        = false;
static constexpr long kDecoderTraceRow           = -1;
static constexpr long kDecoderTraceCol           = -1;
static constexpr bool kDecoderTraceLogRead       = false;
static constexpr bool kDecoderTraceLogWrite      = false;
static constexpr bool kDecoderTraceLogMismatch   = false;

// ==================================

int main() {
  ofec_sweep::SweepParameterConfig config;
  config.interleaver_name = kInterleaverName;
  config.decoder_name = kDecoderName;
  config.bits_per_symbol = kBitsPerSymbol;
  config.normalize_extrinsic = kNormalizeExtrinsic;
  config.quiet_pipeline = kQuietConsole;
  config.quiet_logs = kQuietConsole;

  config.alpha_start_candidates = kAlphaStartCandidates;
  config.alpha_step_candidates = kAlphaStepCandidates;
  config.beta_start_candidates = kBetaStartCandidates;
  config.beta_step_candidates = kBetaStepCandidates;
  config.chase_l_candidates = kChaseLCandidates;

  config.bitgen_seed_count = kBitgenSeedCount;
  config.channel_seed_count = kChannelSeedCount;

  config.ebn0_start = kEbN0Start;
  config.ebn0_end = kEbN0End;
  config.ebn0_points = kEbN0Points;

  config.explicit_patterns = kExplicitAlphaBetaSets;
  config.generate_random_bits = kGenerateRandomBits;
  config.normalize_known_prefix_tail = kNormalizeKnownPrefixTail;
  config.quant_clip_ratio = kQuantClipRatio;

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
  };
  config.base_params.BITGEN_RANDOM_BITS = kGenerateRandomBits;
  config.base_params.NORMALIZE_KNOWN_PREFIX_TAIL = kNormalizeKnownPrefixTail;
  config.base_params.LLR_CLIP_RATIO = kQuantClipRatio;
  config.base_params.LLR_BITS = kLlrBits;

  return ofec_sweep::run_sweep(config);
}
