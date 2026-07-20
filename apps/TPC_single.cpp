#include <vector>

#include "newcode/tpc_single_runner.hpp"

// ======== 用户可改区域 ========
// 只需要改这里的常量即可完成一次“单次调试运行”的配置

// 发射端参数
static constexpr const char* kLabel              = "tpc_debug";
static constexpr int         kBitgenSeed         = 11598;
static constexpr bool        kGenerateRandomBits = true;

// 信道相关参数
static constexpr float       kEbN0_db        = 3.9f - 0.59684f;
static constexpr int         kChannelSeed    = 148898101;
static constexpr unsigned    kBitsPerSymbol  = 1; // 1=BPSK，偶数=QAM

// 解码相关参数
static constexpr int         kMaxIters = 4;
static constexpr int         kNumBlocks = 100;

// 每次行/列解码的系数表，长度需为 2 * kMaxIters（先行后列）
static const std::vector<float> kAlpha_schedule = {
  0.0f, 0.2f, 0.3f, 0.5f, 0.7f, 0.9f, 1.0f, 1.0f
  //0.0f, 0.2f, 0.3f, 0.5f, 0.7f, 0.9f
};
static const std::vector<float> kBeta_schedule = {
  0.2f, 0.4f, 0.6f, 0.8f, 1.0f, 1.0f, 1.0f, 1.0f
  //0.2f, 0.4f, 0.6f, 0.8f, 1.0f, 1.0f
};

// 量化相关参数
static constexpr std::size_t kLlrBits = 16; // 仅支持 16(浮点)

int main() {
  tpc_single::Config config{
    .label = kLabel,
    .ebn0_db = kEbN0_db,
    .max_iters = kMaxIters,
    .num_blocks = kNumBlocks,
    .bits_per_symbol = kBitsPerSymbol,
    .bitgen_seed = kBitgenSeed,
    .channel_seed = kChannelSeed,
    .generate_random_bits = kGenerateRandomBits,
    .llr_bits = kLlrBits,
    .alpha_schedule = kAlpha_schedule,
    .beta_schedule = kBeta_schedule,
  };
  return tpc_single::run_tpc_single(config);
}
