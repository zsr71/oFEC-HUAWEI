#include <iostream>
#include "newcode/frontend/interleaver.hpp"
#include "newcode/decoder_api.hpp"
using namespace newcode;
int main() {
  auto itv = make_interleaver("ofec");
  auto dec = make_decoder("ebchPF");  // 占位，稍后补实现
  Matrix<float> lin0,lch,lin,lout;    // 你的前处理会填充它们
  if (itv) itv->interleave(lin0, lin);
  (void)dec; (void)lch; (void)lout;
  std::cout << "interleaved app skeleton\n";
  return 0;
}
