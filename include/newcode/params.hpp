#pragma once
#include <cstddef>

namespace newcode {

struct Params {
  // ===== 编码结构常量（与 oFEC 布局一致） =====
  static constexpr size_t NUM_SUBBLOCK_COLS     = 8;   // 每行的子块数
  static constexpr size_t INFO_SUBROWS_PER_CODE = 16;  // BCH 前半信息对应的子块行数
  static constexpr size_t BITS_PER_SUBBLOCK_DIM = 16;  // 子块维度(子块大小=16x16 bits)

  // ===== BCH 码参数 =====
  static constexpr size_t BCH_N = 256; // 码长（255 + 1 个整体奇偶位）
  static constexpr size_t BCH_K = 239; // 信息长

  // ===== 运行/仿真参数 =====
  size_t NUM_INFO_BITS     = 362304; // 信息比特总数
  int    BITGEN_SEED       = 42;     // 随机种子
  size_t NUM_GUARD_SUBROWS = 2;      // 保护块子行数 G

  // ===== 解码组织参数（单位：sub-block rows）=====
  size_t TILES_PER_WIN   = 5;  // 每个 window 含多少个 tile（自下而上处理）
  size_t TILE_OVERLAP_BR = 2;  // 相邻 tile 在 sub-block-row 维度上的重叠
  size_t TILE_HEIGHT_BR  = 22; // 单个 tile 的高度（sub-block-row）
  size_t WINDOW_POP_PUSH = 2;  // window 每次滑动的 sub-block-row 数（pop/push）

  size_t WINDOW_ITERS    = 4;  // 每个 window 内重复迭代的轮数（>=1）

  // ===== LLR 量化（与 qfloat/int8_t 打包相关）=====
  size_t LLR_BITS = 5;   // 4 或 5 常见；与硬件位宽一致
  float  LLR_CLIP = 8.0f; // 量化/反量化剪裁幅度（与 qfloat 默认一致）

  // ===== Chase–Pyndiah（chase256 用到的参数）=====
  // L：选择的“最不可靠位置”个数；NTEST：候选数；NCOMP：参与外信息的竞争者数（<= NTEST）
  int   CHASE_L     = 4;
  int   CHASE_NTEST = 16;
  int   CHASE_NCOMP = 8;

  // Pyndiah 系数：Y2 = reliability - a*Y；竞争者相对度量缩放：b；
  // beta = sum(|最不可靠|[前 e 个]) - c * metric(DW)；若该位所有竞争者与 DW 同判：reliability = beta + d*|Y|
  float CP_A = 0.7f;  // 去通道权重 a
  float CP_B = 1.0f;  // 相对度量缩放 b
  float CP_C = 0.0f;  // beta 中的 -c*metric(DW)
  float CP_D = 0.0f;  // 同判兜底项中的 d*|Y|
  int   CP_E = 0;     // 参与 beta 求和的最不可靠位置个数上限（0 表示用 L）

  // ===== 便捷派生（统一换算到“比特行 rows”）=====
  constexpr size_t tile_height_rows() const {
    return TILE_HEIGHT_BR * BITS_PER_SUBBLOCK_DIM;
  }
  constexpr size_t tile_stride_rows() const {
    return (TILE_HEIGHT_BR - TILE_OVERLAP_BR) * BITS_PER_SUBBLOCK_DIM;
  }
  constexpr size_t win_height_br() const {
    return TILES_PER_WIN * TILE_HEIGHT_BR - (TILES_PER_WIN - 1) * TILE_OVERLAP_BR;
  }
  constexpr size_t win_height_rows() const {
    return win_height_br() * BITS_PER_SUBBLOCK_DIM;
  }
  constexpr size_t pop_push_rows() const {
    return WINDOW_POP_PUSH * BITS_PER_SUBBLOCK_DIM;
  }
  // 初始 window 起始行（跳过 warmup 的 B 行 + 2G 行保护）
  constexpr size_t initial_win_start_rows() const {
    return (BITS_PER_SUBBLOCK_DIM + 2 * NUM_GUARD_SUBROWS) * BITS_PER_SUBBLOCK_DIM;
  }

  // ===== 基本有效性检查（可选）=====
  constexpr bool valid() const {
    const bool tiling_ok =
      (TILES_PER_WIN >= 1) &&
      (TILE_OVERLAP_BR >= 1) &&
      (TILE_OVERLAP_BR < TILE_HEIGHT_BR) &&
      (WINDOW_ITERS >= 1);

    const bool chase_ok =
      (CHASE_L >= 1) && (CHASE_L <= static_cast<int>(BCH_N) - 1) && // 不含整体奇偶位
      (CHASE_NTEST >= 1) &&
      (CHASE_NCOMP >= 1) && (CHASE_NCOMP <= CHASE_NTEST);

    return tiling_ok && chase_ok;
  }
};

} // namespace newcode
