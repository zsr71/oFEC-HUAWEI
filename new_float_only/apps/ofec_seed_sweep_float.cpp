#include <iostream>

#include "new_float_only/seed_sweep_runner.hpp"

namespace {
constexpr const char* kLabel = "seed_sweep_float_demo";
constexpr float kEbN0Db = 3.13f;
constexpr std::size_t kTrialCount = 6;
constexpr int kBitgenSeedBase = 1521867291;
constexpr int kChannelSeedBase = 998258255;
constexpr unsigned kBitsPerSymbol = 1;
constexpr bool kGenerateRandomBits = true;
constexpr bool kNormalizeExtrinsic = true;
constexpr unsigned kMaxWorkersOverride = 6;
constexpr bool kQuietPipeline = true;
constexpr bool kQuietLogs = false;
constexpr bool kWriteSummaryCsv = true;
constexpr bool kWriteTrialCsv = true;

constexpr bool kNormalizeKnownPrefixTail = false;
constexpr int kChaseL = 6;
constexpr float kAlpha0 = 0.2f;
constexpr float kAlpha1 = 0.4f;
constexpr float kAlpha2 = 0.6f;
constexpr float kAlpha3 = 0.8f;
constexpr float kBeta0 = 0.2f;
constexpr float kBeta1 = 0.4f;
constexpr float kBeta2 = 0.6f;
constexpr float kBeta3 = 0.8f;
}  // namespace

/**
 * 多核 seed sweep 示例程序。
 * 使用方法是直接修改本文件顶部常量，然后重新构建运行；
 * 程序会在固定参数下并行跑多组 seed，并把汇总 BER 写到 data/ 下的 CSV 文件。
 */
int main() {
  new_float_only::SeedSweepConfig config;
  config.label = kLabel;
  config.ebn0_db = kEbN0Db;
  config.trial_count = kTrialCount;
  config.bitgen_seed_base = kBitgenSeedBase;
  config.channel_seed_base = kChannelSeedBase;
  config.bits_per_symbol = kBitsPerSymbol;
  config.generate_random_bits = kGenerateRandomBits;
  config.normalize_extrinsic = kNormalizeExtrinsic;
  config.max_workers_override = kMaxWorkersOverride;
  config.quiet_pipeline = kQuietPipeline;
  config.quiet_logs = kQuietLogs;
  config.write_summary_csv = kWriteSummaryCsv;
  config.write_trial_csv = kWriteTrialCsv;

  // 这里直接写死需要比较的一组 float plain 解码参数，便于重复扫描不同 seed。
  config.decoder.NORMALIZE_KNOWN_PREFIX_TAIL = kNormalizeKnownPrefixTail;
  config.decoder.CHASE_L = kChaseL;
  config.decoder.CHASE_NTEST = 1 << config.decoder.CHASE_L;
  config.decoder.ALPHA_LIST = {kAlpha0, kAlpha1, kAlpha2, kAlpha3};
  config.decoder.beta_list = {kBeta0, kBeta1, kBeta2, kBeta3};
  config.decoder.HARD_DECODE_DEFAULT = false;
  config.decoder.HARD_TILE_LIST = {0, 0, 0, 0};
  config.decoder.DUMP_WORK_LLR = false;
  config.decoder.debug_trace.enable = false;

  try {
    const auto result = new_float_only::run_seed_sweep(config);
    std::cout << "[APP] done: trials=" << result.trials_completed
              << ", post-FEC BER=" << result.post_fec.ber << "\n";
    return 0;
  } catch (const std::exception& ex) {
    std::cerr << "[APP] failed: " << ex.what() << "\n";
    return 2;
  }
}
