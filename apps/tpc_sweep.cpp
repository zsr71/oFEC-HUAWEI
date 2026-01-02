#include <vector>

#include "newcode/tpc_sweep_runner.hpp"

// ======== 用户可调参数区域 ========
static constexpr unsigned    kBitsPerSymbol = 1; // 1=BPSK，偶数=QAM
static constexpr int         kMaxIters = 4;
static constexpr int         kNumBlocks = 10;
static constexpr bool        kGenerateRandomBits = true;
static constexpr bool        kQuietPipeline = true;
static constexpr bool        kQuietLogs = false;

// alpha/beta schedule（长度必须为 2 * kMaxIters）
static const std::vector<float> kAlphaSchedule = {
  0.0f, 0.2f, 0.3f, 0.5f, 0.7f, 0.9f, 1.0f, 1.0f
};
static const std::vector<float> kBetaSchedule = {
  0.2f, 0.4f, 0.6f, 0.8f, 1.0f, 1.0f, 1.0f, 1.0f
};

// 随机种子（为空则自动生成 bitgen/channel seeds 数量）
static constexpr int kBitgenSeedCount = 4;
static constexpr int kChannelSeedCount = 4;

// Eb/N0 设置
static constexpr float kEbN0Start = 3.9f - 0.59684f;
static constexpr float kEbN0End = 3.9f - 0.59684f;
static constexpr int kEbN0Points = 1;

int main() {
  tpc_sweep::SweepParameterConfig config;
  config.bits_per_symbol = kBitsPerSymbol;
  config.max_iters = kMaxIters;
  config.num_blocks = kNumBlocks;
  config.alpha_schedule = kAlphaSchedule;
  config.beta_schedule = kBetaSchedule;
  config.generate_random_bits = kGenerateRandomBits;
  config.quiet_pipeline = kQuietPipeline;
  config.quiet_logs = kQuietLogs;

  config.bitgen_seed_count = kBitgenSeedCount;
  config.channel_seed_count = kChannelSeedCount;
  config.ebn0_start = kEbN0Start;
  config.ebn0_end = kEbN0End;
  config.ebn0_points = kEbN0Points;

  return tpc_sweep::run_sweep(config);
}
