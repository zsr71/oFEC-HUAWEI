#pragma once
#include <memory>
#include <string>
#include <vector>
#include <cstddef>
#include "newcode/ofec_llr_matrix.hpp"

namespace newcode {

// 统一接口：矩阵版交织
struct IInterleaver {
  virtual ~IInterleaver() = default;
  virtual void interleave  (const Matrix<float>& in, Matrix<float>& out) = 0;
  virtual void deinterleave(const Matrix<float>& in, Matrix<float>& out) = 0;
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
    std::size_t N = 0;             // 总元素数（例如 R*C*H*W）
    std::vector<int> idx_in;       // 旧测试会打印它；先给 identity

    // 旧代码常用：一维向量交织/去交织（此处先直通，后续你可替换为真实映射）
    template <typename T>
    std::vector<T> interleave(const std::vector<T>& x) const {
      return x; // TODO: 用 idx_in 做真实映射
    }
    template <typename T>
    std::vector<T> deinterleave(const std::vector<T>& y) const {
      return y; // TODO: 用 idx_in 的逆映射
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
    void interleave  (const Matrix<float>& in, Matrix<float>& out) const {
      if (impl) impl->interleave(in, out); else out = in;
    }
    void deinterleave(const Matrix<float>& in, Matrix<float>& out) const {
      if (impl) impl->deinterleave(in, out); else out = in;
    }

    std::size_t size() const { return N; }
  };

  // 旧 API：按形状构建（默认 kind="ofec"）
  static Handle build_from_shape(std::size_t R, std::size_t C,
                                 std::size_t H, std::size_t W,
                                 const std::string& kind = "ofec")
  {
    Handle h;
    h.impl = make_interleaver(kind, R, C, H, W);
    h.N = R * C * H * W;
    h.idx_in.resize(h.N);
    for (std::size_t i = 0; i < h.N; ++i) h.idx_in[i] = static_cast<int>(i); // identity
    return h;
  }

  // 旧 API：兼容的另一个入口名
  static Handle build_from_spec(std::size_t R, std::size_t C,
                                std::size_t H, std::size_t W,
                                const std::string& kind = "ofec")
  {
    return build_from_shape(R, C, H, W, kind);
  }
};

} // namespace newcode
