#include <iostream>

#include "new_float_only/seed_sweep_runner.hpp"

namespace {
// 这组常量控制 seed sweep 本身：跑多少次 trial、使用哪组基础 seed、
// 是否静默 pipeline，以及是否把汇总结果写到 CSV。
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

// 这组常量控制固定不变的 float plain 解码参数。
// kAlphaList / kBetaList 按 tile 顺序给出每个 tile 使用的 alpha/beta。
constexpr bool kNormalizeKnownPrefixTail = false;
constexpr int kChaseL = 6;
const std::vector<float> kAlphaList = {0.2f, 0.4f, 0.6f, 0.8f};
const std::vector<float> kBetaList = {0.2f, 0.4f, 0.6f, 0.8f};
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
  config.decoder.ALPHA_LIST = kAlphaList;
  config.decoder.beta_list = kBetaList;
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
