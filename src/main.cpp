#include "newcode/newcode.hpp"
#include <iostream>

int main() {
  using namespace newcode;

  std::cout << hello("oFEC") << "\n";

  const auto v = version();
  std::cout << "version " << v.major << "." << v.minor << "." << v.patch << "\n";

  unsigned x = 0xA5;
  std::cout << "parity(" << x << ") = " << parity(x) << "\n";
  return 0;
}
