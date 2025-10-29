#include "newcode/frontend/interleaver.hpp"

namespace newcode {

struct IdentityInterleaver : IInterleaver {
  void interleave  (const Matrix<float>& in, Matrix<float>& out) override { out = in; }
  void deinterleave(const Matrix<float>& in, Matrix<float>& out) override { out = in; }
};

// 无形状版本：直接返回
std::unique_ptr<IInterleaver> make_interleaver(const std::string& name) {
  if (name == "identity") return std::make_unique<IdentityInterleaver>();
  return nullptr; // 其它名字由 ofec.cpp 处理
}

// 带形状版本：identity 忽略形状
std::unique_ptr<IInterleaver> make_interleaver(const std::string& name,
                                               std::size_t /*R*/, std::size_t /*C*/,
                                               std::size_t /*H*/, std::size_t /*W*/) {
  if (name == "identity") return std::make_unique<IdentityInterleaver>();
  return nullptr;
}

} // namespace newcode
