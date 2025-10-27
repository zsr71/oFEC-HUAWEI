#include "newcode/interleaver.hpp"
#include <random>
#include <iostream>
#include <cassert>

int main() {
  using namespace ofec;

  // OpenROADM oFEC 常用参数
  const int R = 84, C = 8, H = 16, W = 16;
  auto itv = Interleaver::build_from_spec(R, C, H, W);

  const size_t N = itv.size();
  std::vector<uint8_t> x(N), y, x_hat;

  // 随机比特
  std::mt19937 rng(0);
  std::uniform_int_distribution<int> dist(0,1);
  for (auto &b : x) b = static_cast<uint8_t>(dist(rng));

  // 交织 + 解交织
  y     = itv.interleave(x);
  x_hat = itv.deinterleave(y);

  // 校验
  if (x_hat != x) {
    std::cerr << "[CHECK] FAILED: round-trip mismatch\n";
    return 1;
  }
  std::cout << "[CHECK] PASSED: TX interleave + RX deinterleave\n";

  // （可选）打印一些索引检查
  std::cout << "[INFO] N=" << N << "  first 8 idx_in: ";
  for (int i=0;i<8;++i) std::cout << itv.idx_in[i] << (i==7?'\n':' ');

  return 0;
}
