#pragma once

#include <cstddef>
#include <memory>
#include <stdexcept>
#include <string>
#include <vector>

#include "newcode/ofec_llr_matrix.hpp"

namespace newcode {

// 统一接口：矩阵版交织
struct IInterleaver {
  virtual ~IInterleaver() = default;
  virtual void interleave(const Matrix<float>& in, Matrix<float>& out) const = 0;
  virtual void deinterleave(const Matrix<float>& in, Matrix<float>& out) const = 0;
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
    std::vector<T> interleave(const std::vector<T>& x) const {
      if (idx_in.empty()) return x;
      if (x.size() != idx_in.size())
        throw std::invalid_argument("interleave: input size mismatch");
      std::vector<T> y(idx_in.size());
      for (std::size_t i = 0; i < idx_in.size(); ++i)
        y[i] = x[static_cast<std::size_t>(idx_in[i])];
      return y;
    }

    template <typename T>
    std::vector<T> deinterleave(const std::vector<T>& y) const {
      if (idx_out.empty()) return y;
      if (y.size() != idx_out.size())
        throw std::invalid_argument("deinterleave: input size mismatch");
      std::vector<T> x(idx_out.size());
      for (std::size_t i = 0; i < idx_out.size(); ++i)
        x[i] = y[static_cast<std::size_t>(idx_out[i])];
      return x;
    }

    template <typename T>
    std::vector<T> interleave_chunks(const std::vector<T>& v) const {
      return interleave(v);
    }
    template <typename T>
    std::vector<T> deinterleave_chunks(const std::vector<T>& v) const {
      return deinterleave(v);
    }

    // 矩阵版交织：直接委托给实现；若 impl 为空则直通
    void interleave(const Matrix<float>& in, Matrix<float>& out) const {
      if (impl) impl->interleave(in, out); else out = in;
    }
    void deinterleave(const Matrix<float>& in, Matrix<float>& out) const {
      if (impl) impl->deinterleave(in, out); else out = in;
    }

    std::size_t size() const { return idx_in.size(); }
    bool has_mapping() const { return !idx_in.empty() && idx_in.size() == idx_out.size(); }
  };

  // 旧 API：按形状构建（默认 kind="ofec"）
  static Handle build_from_shape(std::size_t rows, std::size_t cols,
                                 std::size_t H, std::size_t W,
                                 const std::string& kind = "ofec")
  {
    Handle h;
    if (rows % H != 0 || cols % W != 0)
      throw std::invalid_argument("build_from_shape: rows/cols not divisible by H/W");

    const std::size_t blocks_R = rows / H;
    const std::size_t blocks_C = cols / W;

    h.N = rows * cols;
    h.impl = make_interleaver(kind, blocks_R, blocks_C, H, W);
    if (h.impl) {
      h.idx_in  = h.impl->forward_mapping();
      h.idx_out = h.impl->inverse_mapping();
    }
    if (!h.has_mapping()) {
      h.idx_in.resize(h.N);
      h.idx_out.resize(h.N);
      for (std::size_t i = 0; i < h.N; ++i) {
        h.idx_in[i] = static_cast<int>(i);
        h.idx_out[i] = static_cast<int>(i);
      }
    }
    return h;
  }

  // 旧 API：兼容的另一个入口名
  static Handle build_from_spec(std::size_t R, std::size_t C,
                                std::size_t H, std::size_t W,
                                const std::string& kind = "ofec")
  {
    return build_from_shape(R * H, C * W, H, W, kind);
  }
};

} // namespace newcode
