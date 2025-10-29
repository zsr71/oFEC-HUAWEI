#include <iostream>
#include "newcode/frontend/interleaver.hpp"
#include "newcode/decoder_api.hpp"
using namespace newcode;
int main() {
  auto itv = make_interleaver("identity");
  auto dec = make_decoder("plain");   // 先占位，稍后合并另一个分支时补实现
  Matrix<float> lin,lch,lout;         // 你的前处理会填充它们
  (void)itv; (void)dec; (void)lin; (void)lch; (void)lout;
  std::cout << "plain app skeleton\n";
  return 0;
}
