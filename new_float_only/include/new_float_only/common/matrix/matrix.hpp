#pragma once
#include <vector>
#include <iostream>
#include <stdexcept>
#include <type_traits>
#include <utility>

namespace matrix {

template<typename T>
class Matrix {
public:
    class Slice {
    public:
        Slice(Matrix& matrix, size_t r0, size_t r1, size_t c0, size_t c1)
            : matrix_(&matrix), r0_(r0), r1_(r1), c0_(c0), c1_(c1) {
            const size_t max_r = matrix_->rows();
            const size_t max_c = matrix_->cols();
            if (r0_ > r1_ || c0_ > c1_)
                throw std::out_of_range("Matrix slice bounds are inverted");
            if (r1_ > max_r || c1_ > max_c)
                throw std::out_of_range("Matrix slice out of range");
        }

        size_t rows() const { return (r1_ >= r0_) ? (r1_ - r0_) : 0; }
        size_t cols() const { return (c1_ >= c0_) ? (c1_ - c0_) : 0; }

        T& at(size_t r, size_t c) { return (*matrix_)[r0_ + r][c0_ + c]; }
        const T& at(size_t r, size_t c) const { return (*matrix_)[r0_ + r][c0_ + c]; }

    private:
        Matrix* matrix_;
        size_t r0_;
        size_t r1_;
        size_t c0_;
        size_t c1_;
    };

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

    Slice slice(size_t r0, size_t r1, size_t c0, size_t c1) {
        return Slice(*this, r0, r1, c0, c1);
    }

    void print(std::ostream& os = std::cout) const {
        for (auto &row : data_) {
            for (auto &val : row) os << val << " ";
            os << "\n";
        }
    }

private:
    std::vector<std::vector<T>> data_;
};

template <typename T>
void scale(typename Matrix<T>::Slice slice, T factor)
{
    const size_t rN = slice.rows();
    const size_t cN = slice.cols();
    for (size_t r = 0; r < rN; ++r)
        for (size_t c = 0; c < cN; ++c)
            slice.at(r, c) *= factor;
}

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

} // namespace new_float_only
