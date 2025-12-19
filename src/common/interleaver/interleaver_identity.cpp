#include "newcode/common/interleaver/interleaver.hpp"

namespace interleaver {
namespace {

class IdentityInterleaver final : public IInterleaver {
public:
  IdentityInterleaver() = default;
  explicit IdentityInterleaver(std::size_t size) { set_size(size); }

  void set_size(std::size_t size) {
    mapping_.resize(size);
    for (std::size_t i = 0; i < size; ++i)
      mapping_[i] = static_cast<int>(i);
  }

  void interleave(const newcode::Matrix<float>& in, newcode::Matrix<float>& out) const override { out = in; }
  void deinterleave(const newcode::Matrix<float>& in, newcode::Matrix<float>& out) const override { out = in; }

  const std::vector<int>& forward_mapping() const override { return mapping_; }
  const std::vector<int>& inverse_mapping() const override { return mapping_; }

private:
  std::vector<int> mapping_;
};

} // namespace

std::unique_ptr<IInterleaver> make_interleaver(const std::string& name) {
  if (name == "identity")
    return std::make_unique<IdentityInterleaver>();
  return nullptr;
}

std::unique_ptr<IInterleaver> make_interleaver(const std::string& name,
                                               std::size_t R, std::size_t C,
                                               std::size_t H, std::size_t W) {
  if (name != "identity") return nullptr;
  auto ptr = std::make_unique<IdentityInterleaver>();
  ptr->set_size(R * C * H * W);
  return ptr;
}

} // namespace newcode
