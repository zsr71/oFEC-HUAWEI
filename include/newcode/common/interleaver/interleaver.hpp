#pragma once

#include <cstddef>
#include <memory>
#include <stdexcept>
#include <string>
#include <vector>

#include "newcode/ofec_llr_matrix.hpp"

namespace interleaver {

// 统一接口：矩阵版交织
struct IInterleaver {
  virtual ~IInterleaver() = default;
  virtual void interleave(const newcode::Matrix<float>& in, newcode::Matrix<float>& out) const = 0;
  virtual void deinterleave(const newcode::Matrix<float>& in, newcode::Matrix<float>& out) const = 0;
  virtual const std::vector<int>& forward_mapping() const = 0; // y[i] = x[forward_mapping()[i]]
  virtual const std::vector<int>& inverse_mapping() const = 0; // inverse permutation
};

// 工厂：按名字创建
std::unique_ptr<IInterleaver> make_interleaver(const std::string& name);

// 工厂（带形状信息）
std::unique_ptr<IInterleaver> make_interleaver(const std::string& name,
                                               std::size_t R, std::size_t C,
                                               std::size_t H, std::size_t W);

// 兼容旧 API 的薄封装，把 unique_ptr<IInterleaver> 适配成“老接口”
struct Interleaver {
  struct Handle {
    std::unique_ptr<IInterleaver> impl;
    std::size_t N = 0;                // 总元素数（例如 R*C*H*W）
    std::vector<int> idx_in;          // y[pos] = x[idx_in[pos]]
    std::vector<int> idx_out;         // x[pos] = y[idx_out[pos]]

    template <typename T>
    std::vector<T> interleave(const std::vector<T>& x) const;

    template <typename T>
    std::vector<T> deinterleave(const std::vector<T>& y) const;

    template <typename T>
    std::vector<T> interleave_chunks(const std::vector<T>& v) const;
    template <typename T>
    std::vector<T> deinterleave_chunks(const std::vector<T>& v) const;

    // 矩阵版交织：直接委托给实现；若 impl 为空则直通
    void interleave(const newcode::Matrix<float>& in, newcode::Matrix<float>& out) const;
    void deinterleave(const newcode::Matrix<float>& in, newcode::Matrix<float>& out) const;

    std::size_t size() const;
    bool has_mapping() const;
  };

  // 旧 API：按形状构建（默认 kind="ofec"）
  static Handle build_from_shape(std::size_t rows, std::size_t cols,
                                 std::size_t H, std::size_t W,
                                 const std::string& kind = "ofec");

  // 旧 API：兼容的另一个入口名
  static Handle build_from_spec(std::size_t R, std::size_t C,
                                std::size_t H, std::size_t W,
                                const std::string& kind = "ofec");
};

} // namespace interleaver

#include "interleaver/interleaver.ipp"

