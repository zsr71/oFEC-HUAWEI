#pragma once
#include <vector>
#include <iostream>
#include <stdexcept>

namespace newcode {

template<typename T>
class Matrix {
public:
    Matrix() = default;
    Matrix(size_t n_rows, size_t n_cols)
        : data_(n_rows, std::vector<T>(n_cols, T{})) {}

    static Matrix<T> zero(size_t n_rows, size_t n_cols) {
        return Matrix<T>(n_rows, n_cols);
    }

    void add_row() {
        size_t n_cols = cols();
        data_.emplace_back(n_cols, T{});
    }

    void add_col() {
        for (auto &row : data_) row.push_back(T{});
    }

    void erase_row(size_t row_index, size_t n_rows = 1) {
        if (row_index + n_rows > rows())
            throw std::out_of_range("Row index out of range");
        data_.erase(data_.begin() + row_index, data_.begin() + row_index + n_rows);
    }

    void erase_col(size_t col_index, size_t n_cols = 1) {
        if (col_index + n_cols > cols())
            throw std::out_of_range("Column index out of range");
        for (auto &row : data_)
            row.erase(row.begin() + col_index, row.begin() + col_index + n_cols);
    }

    size_t rows() const { return data_.size(); }
    size_t cols() const { return data_.empty() ? 0 : data_[0].size(); }

    std::vector<T>& operator[](size_t r) { return data_[r]; }
    const std::vector<T>& operator[](size_t r) const { return data_[r]; }

    void print(std::ostream& os = std::cout) const {
        for (auto &row : data_) {
            for (auto &val : row) os << val << " ";
            os << "\n";
        }
    }

private:
    std::vector<std::vector<T>> data_;
};

} // namespace newcode
