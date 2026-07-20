#include <iostream>

#include "new_float_only/ebn0_sweep_runner.hpp"

namespace {

constexpr const char* kLabel = "ofec_ebn0_sweep_float_demo"; // 本次 sweep 的标签，会出现在日志和 CSV 中
constexpr float kEbN0Start = 3.07f; // Eb/N0 扫描起点，单位 dB
constexpr float kEbN0End = 3.37f; // Eb/N0 扫描终点，单位 dB
constexpr int kEbN0Points = 30; // 在 [kEbN0Start, kEbN0End] 之间均匀采样的点数
constexpr std::size_t kTrialCount = 1; // 每个 Eb/N0 点重复跑多少个 seed trial
constexpr int kBitgenSeedBase = 1521867291; // 自动生成 bit 源 seed 时使用的起始值
constexpr int kChannelSeedBase = 998258255; // 自动生成信道噪声 seed 时使用的起始值
constexpr unsigned kBitsPerSymbol = 1; // 每个调制符号携带的 bit 数，1 表示 BPSK
constexpr bool kGenerateRandomBits = true; // true=随机信息比特，false=全 0 比特
constexpr bool kNormalizeExtrinsic = true; // 是否对 decoder 输出的 extrinsic 做归一化
constexpr unsigned kMaxWorkersOverride = 0; // 并行 worker 上限；0 表示按默认/环境变量自动决定
constexpr bool kQuietPipeline = true; // true=单个 task 内部少打印 pipeline 细节
constexpr bool kQuietLogs = false; // true=减少 sweep 总控日志与进度输出
constexpr bool kWriteSummaryCsv = true; // 是否输出按 Eb/N0 聚合后的 summary CSV
constexpr bool kWriteTrialCsv = true; // 是否输出每个 (Eb/N0, seed) task 的 trial CSV

constexpr bool kNormalizeKnownPrefixTail = false; // 是否对 known-prefix 之后的尾部 LLR 做归一化
constexpr int kChaseL = 6; // Chase decoder 选取的最不可靠 bit 数量
const std::vector<float> kAlphaList = {0.2f, 0.4f, 0.6f, 0.8f}; // 每个 tile 使用的 alpha 列表
const std::vector<float> kBetaList = {0.2f, 0.4f, 0.6f, 0.8f}; // 每个 tile 使用的 beta 列表

}  // namespace

int main() {
  new_float_only::Ebn0SweepConfig config;
  config.label = kLabel;
  config.ebn0_start = kEbN0Start;
  config.ebn0_end = kEbN0End;
  config.ebn0_points = kEbN0Points;
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
    const auto result = new_float_only::run_ebn0_sweep(config);
    std::cout << "[APP] done: Eb/N0 points=" << result.point_results.size()
              << ", trials=" << result.trial_results.size() << "\n";
    return 0;
  } catch (const std::exception& ex) {
    std::cerr << "[APP] failed: " << ex.what() << "\n";
    return 2;
  }
}
