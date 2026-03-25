#include <iostream>
#include <vector>

#include "new_float_only/single_runner.hpp"

namespace {
constexpr const char* kLabel = "float_plain_demo";
constexpr float kEbN0Db = 3.07f;
constexpr int kBitgenSeed = 1521867291;
constexpr int kChannelSeed = 998258255;
constexpr unsigned kBitsPerSymbol = 1;
constexpr bool kGenerateRandomBits = true;
constexpr bool kNormalizeExtrinsic = true;
constexpr bool kNormalizeKnownPrefixTail = false;
constexpr bool kDumpWorkLlr = false;
constexpr const char* kWorkLlrPath = "data/llr/work_llr_float.txt";
constexpr long kTraceBitIndex = 1048514;
constexpr const char* kTraceLabel = "bit1048514";
constexpr const char* kTraceCsvDir = "data/chase_csv/bit1048514";

/**
 * 老库 ofec_single 的 bit-index 调试口径：
 * bit_index 先映射到发送矩阵的信息区，再由 Tile 输入重排逻辑追到 Chase 的 (row,k)。
 */
constexpr long BitIndexToRow(long bit_index) {
  // 功能：把信息比特的一维 bit_index 映射到全局 work_llr 的目标行号。
  // 输入：bit_index 为按信息位顺序编号的索引。
  // 输出：返回追踪该 bit 时使用的全局矩阵行坐标。
  return bit_index / 111 + 352;
}

constexpr long BitIndexToCol(long bit_index) {
  // 功能：把信息比特的一维 bit_index 映射到全局 work_llr 的目标列号。
  // 输入：bit_index 为按信息位顺序编号的索引。
  // 输出：返回追踪该 bit 时使用的全局矩阵列坐标。
  return bit_index % 111;
}
}  // namespace

/**
 * 新库的最小示例程序。
 * 这里只保留 float plain 单次运行需要的配置，不再出现量化、MUX、early-stop 等选项。
 */
int main() {
  // 功能：执行一次单独的 float plain 链路实验，并打开 work_llr / chase trace 输出。
  // 输入：不接收命令行参数，所有实验配置都来自文件顶部常量。
  // 输出：成功时返回 0，并在 data/ 下写日志和调试文件；失败时返回 2。
  new_float_only::SingleRunConfig config;
  config.label = kLabel;
  config.ebn0_db = kEbN0Db;
  config.bits_per_symbol = kBitsPerSymbol;
  config.bitgen_seed = kBitgenSeed;
  config.channel_seed = kChannelSeed;
  config.generate_random_bits = kGenerateRandomBits;
  config.normalize_extrinsic = kNormalizeExtrinsic;

  // 这里直接配置 float plain 解码器参数。
  config.decoder.NORMALIZE_KNOWN_PREFIX_TAIL = kNormalizeKnownPrefixTail;
  config.decoder.CHASE_L = 6;
  config.decoder.CHASE_NTEST = 1 << config.decoder.CHASE_L;
  config.decoder.ALPHA_LIST = {0.2f, 0.4f, 0.6f, 0.8f};
  config.decoder.beta_list = {0.2f, 0.4f, 0.6f, 0.8f};
  config.decoder.DUMP_WORK_LLR = kDumpWorkLlr;
  config.decoder.WORK_LLR_OUTPUT_PATH = kWorkLlrPath;
  config.decoder.debug_trace.enable = true;
  config.decoder.debug_trace.dump_chase_csv = true;
  config.decoder.debug_trace.chase_csv_dir = kTraceCsvDir;
  config.decoder.debug_trace.targets.push_back(
      new_float_only::Params::DebugTraceConfig::TraceTarget{
          .row = BitIndexToRow(kTraceBitIndex),
          .col = BitIndexToCol(kTraceBitIndex),
          .bit_index = kTraceBitIndex,
          .label = kTraceLabel,
      });

  try {
    const auto result = new_float_only::run_single(config);
    std::cout << "[APP] done: post-FEC BER=" << result.post_fec.ber << "\n";
    return 0;
  } catch (const std::exception& ex) {
    std::cerr << "[APP] failed: " << ex.what() << "\n";
    return 2;
  }
}
