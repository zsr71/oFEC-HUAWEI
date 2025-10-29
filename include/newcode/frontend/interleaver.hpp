#pragma once
#include <memory>
#include <string>
#include <cstddef>
#include "newcode/ofec_llr_matrix.hpp"

namespace newcode {

// 统一的交织接口
struct IInterleaver {
  virtual ~IInterleaver() = default;
  virtual void interleave  (const Matrix<float>& in, Matrix<float>& out) = 0;
  virtual void deinterleave(const Matrix<float>& in, Matrix<float>& out) = 0;
};

// 轻量工厂：按名字创建
std::unique_ptr<IInterleaver> make_interleaver(const std::string& name);

// 带形状信息的重载（老代码常见：R,C,H,W）
std::unique_ptr<IInterleaver> make_interleaver(const std::string& name,
                                               std::size_t R, std::size_t C,
                                               std::size_t H, std::size_t W);

// 兼容旧 API 的薄封装（静态构建）
struct Interleaver {
  // 旧代码里常见的两个入口名，内部转到工厂
  static std::unique_ptr<IInterleaver> build_from_shape(std::size_t R, std::size_t C,
                                                        std::size_t H, std::size_t W,
                                                        const std::string& kind = "ofec")
  {
    return make_interleaver(kind, R, C, H, W);
  }
  static std::unique_ptr<IInterleaver> build_from_spec(std::size_t R, std::size_t C,
                                                       std::size_t H, std::size_t W,
                                                       const std::string& kind = "ofec")
  {
    return make_interleaver(kind, R, C, H, W);
  }
};

} // namespace newcode
