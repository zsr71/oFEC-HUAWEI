#include "newcode/matrix.hpp"
#include <cassert>
#include <iostream>

int main() {
    newcode::Matrix<int> M = newcode::Matrix<int>::zero(2, 3);
    assert(M.rows() == 2);
    assert(M.cols() == 3);

    M.add_row();
    M.add_col();
    assert(M.rows() == 3);
    assert(M.cols() == 4);

    M.erase_row(0);
    M.erase_col(0);
    assert(M.rows() == 2);
    assert(M.cols() == 3);

    std::cout << "All Matrix tests passed!\n";
    return 0;
}
