#pragma once

#include <cstddef>
#include <cstdint>
#include <memory>
#include <string>
#include <vector>

namespace new_float_only {

/**
 * 浮点 oFEC 解码器配置。
 * 这个结构同时承载编码布局常量、窗口/Tile 组织参数、Chase 参数、
 * 硬判决参数以及调试跟踪开关。它会被 TX/RX 单次链路共同使用。
 */
struct DecoderConfig {
  static constexpr std::size_t NUM_SUBBLOCK_COLS = 8;
  static constexpr std::size_t INFO_SUBROWS_PER_CODE = 16;
  static constexpr std::size_t BITS_PER_SUBBLOCK_DIM = 16;

  static constexpr std::size_t BCH_N = 256;
  static constexpr std::size_t BCH_K = 239;
  static constexpr std::size_t BCH_PARITY_BITS = BCH_N - BCH_K - 1;
  static constexpr std::size_t BCH_OVERALL_IDX = BCH_N - 1;

  std::size_t NUM_INFO_BITS = 6 * 88 * 16 * 111;
  bool NORMALIZE_KNOWN_PREFIX_TAIL = false;
  std::size_t NUM_GUARD_SUBROWS = 2;
  bool DUMP_WORK_LLR = false;
  std::string WORK_LLR_OUTPUT_PATH;

  std::size_t TILES_PER_WIN = 4;
  std::size_t TILE_OVERLAP_BR = 0;
  std::size_t TILE_HEIGHT_BR = 22;
  std::size_t WINDOW_POP_PUSH = 2;

  int CHASE_L = 6;
  int CHASE_NTEST = 64;
  int CHASE_SBR = 2;

  float beta = 0.35f;
  float ALPHA = 1.0f;
  std::vector<float> ALPHA_LIST = {0.2f, 0.4f, 0.6f, 0.8f};
  std::vector<float> beta_list = {0.2f, 0.4f, 0.6f, 0.8f};

  bool HARD_DECODE_DEFAULT = false;
  std::vector<int> HARD_TILE_LIST = {0, 0, 0, 0, 1};
  float HARD_LLR_MAG = 1.0f;

  /**
   * Chase 级别调试开关。
   * row/col 为全局 work_llr 坐标，targets 可用于一次跟踪多个比特。
   */
  struct DebugTraceConfig {
    bool enable = false;
    bool log_read_mapping = false;
    bool log_write_mapping = false;
    bool log_mismatch = false;
    bool log_chase_detail = false;
    bool dump_chase_csv = false;
    long row = -1;
    long col = -1;
    int chase_decoder_row = -1;
    int chase_decoder_col = -1;
    int chase_tile_index = -1;
    int chase_invocation = -1;
    std::string chase_csv_dir;
    std::shared_ptr<std::vector<std::vector<int8_t>>> chase_expected_bits;
    const std::vector<int8_t>* chase_expected_bits_row = nullptr;

    struct TraceTarget {
      long row = -1;
      long col = -1;
      long bit_index = -1;
      std::string label;
    };

    struct ChaseTraceEntry {
      int row_index = -1;
      int k = -1;
      long global_row = -1;
      long global_col = -1;
      long bit_index = -1;
      std::string label;
      int expected_bit = -1;
      std::vector<uint8_t> cplus_bits;
      std::vector<uint8_t> cminus_bits;
    };

    std::vector<TraceTarget> targets;
    std::vector<ChaseTraceEntry> active_chase_entries;
  };

  DebugTraceConfig debug_trace{};

  /// 返回单个 Tile 的比特行高。
  constexpr std::size_t tile_height_rows() const {
    return TILE_HEIGHT_BR * BITS_PER_SUBBLOCK_DIM;
  }

  /// 返回相邻 Tile 在窗口内滑动时的比特行步长。
  constexpr std::size_t tile_stride_rows() const {
    return (TILE_HEIGHT_BR - TILE_OVERLAP_BR) * BITS_PER_SUBBLOCK_DIM;
  }

  /// 返回单个窗口覆盖的 sub-block-row 数。
  constexpr std::size_t win_height_br() const {
    return TILES_PER_WIN * TILE_HEIGHT_BR - (TILES_PER_WIN - 1) * TILE_OVERLAP_BR;
  }

  /// 返回单个窗口覆盖的比特行数。
  constexpr std::size_t win_height_rows() const {
    return win_height_br() * BITS_PER_SUBBLOCK_DIM;
  }

  /// 返回窗口前移一次时的比特行数。
  constexpr std::size_t pop_push_rows() const {
    return WINDOW_POP_PUSH * BITS_PER_SUBBLOCK_DIM;
  }

  /// 当前实现从矩阵顶部开始滑窗，因此初始起点固定为 0。
  constexpr std::size_t initial_win_start_rows() const {
    return 0;
  }

  /// 检查窗口/Tile/Chase 参数是否满足基本约束。
  constexpr bool valid() const {
    return (TILES_PER_WIN >= 1) &&
           (TILE_OVERLAP_BR < TILE_HEIGHT_BR) &&
           (CHASE_L >= 1) &&
           (CHASE_NTEST >= 1) &&
           (CHASE_SBR == 1 || CHASE_SBR == 2);
  }
};

using Params = DecoderConfig;

/**
 * 单次完整链路的运行配置。
 * 这层负责比特生成、信道、调制解调等实验级参数，
 * 真正的 oFEC 解码细节由 decoder 子配置承载。
 */
struct SingleRunConfig {
  std::string label = "float_only";
  float ebn0_db = 3.24f;
  unsigned bits_per_symbol = 2;
  int bitgen_seed = 56456;
  int channel_seed = 57112;
  bool generate_random_bits = true;
  bool normalize_extrinsic = true;
  DecoderConfig decoder{};
};

}  // namespace new_float_only
