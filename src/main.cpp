#include "newcode/params.hpp"
#include "newcode/bitgen.hpp"
#include <iostream>

int main() {
    using namespace newcode;

    Params p;
    p.Global_Information_bits = 16;   // 例如生成16个比特
    auto bits = generate_bits(p);
    std::cout << "Generated bits: ";
    for (auto b : bits) std::cout << int(b);
    std::cout << "\n";

    return 0;
}
