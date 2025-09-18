#pragma once
#include <cstddef> 

namespace newcode {
struct Params {
  // 仿真
  size_t Global_Information_bits = 65536;  // 信息比特数
  int Bit_Gen_seed = 42;
};
}
