#include "newcode/frontend/interleaver.hpp"

namespace newcode {

struct OfecInterleaver : IInterleaver {
  void interleave  (const Matrix<float>& in, Matrix<float>& out) override {
    // TODO: ofec 交织映射（根据 R,C,H,W 的规则；临时先直接拷贝以便编过）
    out = in;
  }
  void deinterleave(const Matrix<float>& in, Matrix<float>& out) override {
    // TODO: ofec 逆映射
    out = in;
  }
};

std::unique_ptr<IInterleaver> make_interleaver(const std::string& name) {
  if (name == "ofec") return std::make_unique<OfecInterleaver>();
  return nullptr;
}

std::unique_ptr<IInterleaver> make_interleaver(const std::string& name,
                                               std::size_t /*R*/, std::size_t /*C*/,
                                               std::size_t /*H*/, std::size_t /*W*/) {
  if (name == "ofec") return std::make_unique<OfecInterleaver>();
  return nullptr;
}

} // namespace newcode
