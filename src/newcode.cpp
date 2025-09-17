#include "newcode/newcode.hpp"
#include <bit>      // std::popcount
#include <string>

namespace newcode {

Version version() { return {0, 1, 0}; }

int parity(unsigned x) {
  return std::popcount(x) & 1u; // 1 表示奇校验
}

std::string hello(const std::string& name) {
  return "hello, " + name;
}

} // namespace newcode
