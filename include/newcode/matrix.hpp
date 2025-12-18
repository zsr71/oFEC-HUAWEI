#pragma once
#include <vector>
#include <iostream>
#include <stdexcept>
#include <type_traits>
#include <utility>

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

// Flatten a Matrix in row-major order into a contiguous vector.
template <typename T>
std::vector<T> flatten_row_major(const Matrix<T>& matrix)
{
    std::vector<T> out;
    out.reserve(matrix.rows() * matrix.cols());
    for (size_t r = 0; r < matrix.rows(); ++r)
        for (size_t c = 0; c < matrix.cols(); ++c)
            out.push_back(matrix[r][c]);
    return out;
}

// Flatten a Matrix in row-major order and apply a transform to each element.
// The transform should be callable as: U f(const T&).
template <typename T, typename Transform>
auto flatten_row_major(const Matrix<T>& matrix, Transform&& transform)
    -> std::vector<std::decay_t<decltype(std::declval<Transform&>()(std::declval<const T&>()))>>
{
    using U = std::decay_t<decltype(std::declval<Transform&>()(std::declval<const T&>()))>;
    std::vector<U> out;
    out.reserve(matrix.rows() * matrix.cols());
    for (size_t r = 0; r < matrix.rows(); ++r) {
        for (size_t c = 0; c < matrix.cols(); ++c) {
            out.push_back(transform(matrix[r][c]));
        }
    }
    return out;
}

} // namespace newcode
