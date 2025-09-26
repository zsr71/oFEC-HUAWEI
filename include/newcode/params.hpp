#pragma once
#include <cstddef>

namespace newcode {

struct Params {
  // ===== 编码结构常量（与 oFEC 布局一致） =====
  static constexpr size_t NUM_SUBBLOCK_COLS     = 8;   // 每行的子块数
  static constexpr size_t INFO_SUBROWS_PER_CODE = 16;  // BCH 前半信息对应的子块行数
  static constexpr size_t BITS_PER_SUBBLOCK_DIM = 16;  // 子块维度(子块大小=16x16 bits)

  // ===== BCH 码参数 =====
  static constexpr size_t BCH_N = 256; // 码长
  static constexpr size_t BCH_K = 239; // 信息长

  // ===== 运行/仿真参数 =====
  size_t NUM_INFO_BITS     = 362304; // 信息比特总数
  int    BITGEN_SEED       = 42;     // 随机种子
  size_t NUM_GUARD_SUBROWS = 2;      // 保护块子行数 G

  // ===== 解码组织参数（单位：sub-block rows）=====
  // 与现有实现保持一致的默认值
  size_t TILES_PER_WIN   = 5;  // 每个 window 含多少个 tile（自下而上处理）
  size_t TILE_OVERLAP_BR = 2;  // 相邻 tile 在 sub-block-row 维度上的重叠
  size_t TILE_HEIGHT_BR  = 22; // 单个 tile 的高度（sub-block-row）
  size_t WINDOW_POP_PUSH = 2;  // window 每次滑动的 sub-block-row 数（pop/push）

  size_t WINDOW_ITERS    = 4;  // 每个 window 内重复迭代的轮数（>=1）
  size_t LLR_BITS = 5;   // 4 或 5
  float  LLR_CLIP = 8.0f;
  // ===== 便捷派生（统一换算到“比特行 rows”）=====
  // tile 高度（比特行）
  constexpr size_t tile_height_rows() const {
    return TILE_HEIGHT_BR * BITS_PER_SUBBLOCK_DIM;
  }
  // tile 步幅（比特行）
  constexpr size_t tile_stride_rows() const {
    return (TILE_HEIGHT_BR - TILE_OVERLAP_BR) * BITS_PER_SUBBLOCK_DIM;
  }
  // window 高度（以 sub-block-row 计）
  constexpr size_t win_height_br() const {
    return TILES_PER_WIN * TILE_HEIGHT_BR - (TILES_PER_WIN - 1) * TILE_OVERLAP_BR;
  }
  // window 高度（比特行）
  constexpr size_t win_height_rows() const {
    return win_height_br() * BITS_PER_SUBBLOCK_DIM;
  }
  // window 每次滑动的行数（比特行）
  constexpr size_t pop_push_rows() const {
    return WINDOW_POP_PUSH * BITS_PER_SUBBLOCK_DIM;
  }
  // 初始 window 起始行（比特行）：跳过 warmup 的 B 行 + 2G 行保护
  constexpr size_t initial_win_start_rows() const {
    return (BITS_PER_SUBBLOCK_DIM + 2 * NUM_GUARD_SUBROWS) * BITS_PER_SUBBLOCK_DIM;
  }

  // 基本有效性检查（可选）
  constexpr bool valid() const {
    return (TILES_PER_WIN >= 1) &&
           (TILE_OVERLAP_BR >= 1) &&
           (TILE_OVERLAP_BR < TILE_HEIGHT_BR) &&
           (WINDOW_ITERS >= 1);
  }
};

} // namespace newcode
