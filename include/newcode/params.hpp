#pragma once
#include <cstddef> 

namespace newcode {
struct Params {
  //常数
  static constexpr size_t NUM_SUBBLOCK_COLS     = 8;   // oFEC无限行矩阵中每行包含的子块数
  static constexpr size_t INFO_SUBROWS_PER_CODE = 16;  // BCH 前半信息部分对应的子块行数
  static constexpr size_t BITS_PER_SUBBLOCK_DIM = 16;  // 子块维度 (子块大小 = 16x16 bits)

  // 仿真
  size_t NUM_INFO_BITS     = 65536;  // 信息比特总数
  int    BITGEN_SEED       = 42;     // 随机比特生成器种子
  int    NUM_GUARD_SUBROWS = 2;      // 保护块子行数

};
}
