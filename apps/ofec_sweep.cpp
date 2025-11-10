#include <vector>

#include "newcode/ofec_sweep_runner.hpp"

// ======== 用户可调参数区域 ========
static constexpr const char* kInterleaverName = "identity";
static constexpr const char* kDecoderName = "plain";
static constexpr unsigned    kBitsPerSymbol = 2; // 设为 1 使用 BPSK，>=2 且偶数使用 QAM
static constexpr bool        kNormalizeExtrinsic = true;

// Alpha/Beta 扫描候选
static const std::vector<float> kAlphaStartCandidates = {0.01f, 0.1f, 0.15f, 0.3f, 0.5f, 0.7f};
static const std::vector<float> kAlphaStepCandidates  = {0.0f, 0.025f, 0.05f, 0.1f, 0.15f, 0.2f};
static const std::vector<float> kBetaStartCandidates  = {0.1f, 0.2f, 0.3f, 0.4f, 0.5f, 0.7f, 0.9f, 1.2f, 1.7f};
static const std::vector<float> kBetaStepCandidates   = {0.0f, 0.025f, 0.05f, 0.1f, 0.2f};
static const std::vector<int>   kChaseLCandidates     = {6};

// 随机种子（为空则自动生成 bitgen/channel seeds 数量）
static constexpr int kBitgenSeedCount  = 8;
static constexpr int kChannelSeedCount = 8;

// Eb/N0 设置
static constexpr float kEbN0Start = 3.27f;
static constexpr float kEbN0End   = 3.27f;
static constexpr int   kEbN0Points = 1;

// 显式 alpha/beta 模式（可选）
static const std::vector<ofec_sweep::ExplicitAlphaBetaPattern> kExplicitAlphaBetaSets = {
  {"custom_label", {0.3f, 0.45f, 0.60f, 0.9f, 0.5f}, {0.2f, 0.225f, 0.250f, 0.275f, 0.8f}},
};

// Decoder 调试跟踪配置
static constexpr bool kDecoderTraceEnable        = false;
static constexpr long kDecoderTraceRow           = -1;
static constexpr long kDecoderTraceCol           = -1;
static constexpr bool kDecoderTraceLogRead       = false;
static constexpr bool kDecoderTraceLogWrite      = false;
static constexpr bool kDecoderTraceLogMismatch   = false;
static constexpr bool kGenerateRandomBits        = true;
static constexpr bool kNormalizeKnownPrefixTail  = true;
// ==================================

int main() {
  ofec_sweep::SweepParameterConfig config;
  config.interleaver_name = kInterleaverName;
  config.decoder_name = kDecoderName;
  config.bits_per_symbol = kBitsPerSymbol;
  config.normalize_extrinsic = kNormalizeExtrinsic;

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

  config.base_params.debug_trace = newcode::Params::DebugTraceConfig{
    .enable = kDecoderTraceEnable,
    .log_read_mapping = kDecoderTraceLogRead,
    .log_write_mapping = kDecoderTraceLogWrite,
    .log_mismatch = kDecoderTraceLogMismatch,
    .row = kDecoderTraceRow,
    .col = kDecoderTraceCol,
  };
  config.base_params.BITGEN_RANDOM_BITS = kGenerateRandomBits;
  config.base_params.NORMALIZE_KNOWN_PREFIX_TAIL = kNormalizeKnownPrefixTail;

  return ofec_sweep::run_sweep(config);
}
