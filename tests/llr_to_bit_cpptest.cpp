#include "newcode/llr_to_bit.hpp"
#include "newcode/matrix.hpp"   // 需要 Matrix 的定义
#include <cassert>
#include <iostream>
#include <cstdint>

using namespace newcode;

int main() {
    Matrix<float> llr(2, 3);
    llr[0][0] = +0.1f;  llr[0][1] = -0.2f; llr[0][2] =  0.0f;
    llr[1][0] = -1.5f;  llr[1][1] = +2.3f; llr[1][2] = -0.0f; // -0.0f 也会判为 0

    auto bits = llr_to_bit(llr);
    assert(bits.rows() == 2 && bits.cols() == 3);

    const uint8_t expect[2][3] = {{0,1,0},{1,0,0}};
    for (size_t r = 0; r < 2; ++r) {
        for (size_t c = 0; c < 3; ++c) {
            if (bits[r][c] != expect[r][c]) {
                std::cerr << "Mismatch at (" << r << "," << c << "): got "
                          << unsigned(bits[r][c]) << " expected "
                          << unsigned(expect[r][c]) << "\n";
                return 1;
            }
        }
    }

    std::cout << "llr_to_bit cpptest: PASS\n";
    return 0;
}
