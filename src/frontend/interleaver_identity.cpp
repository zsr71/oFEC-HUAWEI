#include "newcode/frontend/interleaver.hpp"
namespace newcode {
struct IdentityInterleaver : IInterleaver {
  void interleave  (const Matrix<float>& in, Matrix<float>& out) override { out = in; }
  void deinterleave(const Matrix<float>& in, Matrix<float>& out) override { out = in; }
};
std::unique_ptr<IInterleaver> make_interleaver(const std::string& name) {
  if (name == "identity") return std::make_unique<IdentityInterleaver>();
  return nullptr; // 让 ofec 的实现去接住
}
} // namespace newcode
