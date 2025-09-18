#include "newcode/params.hpp"
#include "newcode/bitgen.hpp"
#include "newcode/ofec_encoder.hpp"
#include <iostream>

int main() {
    using namespace newcode;

    Params p;
    auto bits = generate_bits(p);

    auto code_matrix = ofec_encode(bits, p);

    std::cout << "oFEC encoded matrix: " << code_matrix.rows()
              << " x " << code_matrix.cols() << std::endl;

    return 0;
}