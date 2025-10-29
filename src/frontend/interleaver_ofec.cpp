#include "newcode/frontend/interleaver.hpp"
namespace newcode {
struct OfecInterleaver : IInterleaver {
  void interleave  (const Matrix<float>& in, Matrix<float>& out) override { /* TODO: ofec 映射 */ }
  void deinterleave(const Matrix<float>& in, Matrix<float>& out) override { /* TODO: 逆映射 */ }
};
// 与 identity.cpp 共同“拼接”工厂
std::unique_ptr<IInterleaver> make_interleaver(const std::string& name);
std::unique_ptr<IInterleaver> make_interleaver(const std::string& name) {
  if (name == "ofec")     return std::make_unique<OfecInterleaver>();
  return nullptr; // 其它名字交给 identity.cpp
}
} // namespace newcode
