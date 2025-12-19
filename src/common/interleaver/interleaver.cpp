#include "newcode/common/interleaver/interleaver.hpp"

namespace interleaver {

void Interleaver::Handle::interleave(const newcode::Matrix<float>& in, newcode::Matrix<float>& out) const {
  if (impl) {
    impl->interleave(in, out);
  } else {
    out = in; 
  }
}

void Interleaver::Handle::deinterleave(const newcode::Matrix<float>& in, newcode::Matrix<float>& out) const {
  if (impl) {
    impl->deinterleave(in, out);
  } else {
    out = in;
  }
}

std::size_t Interleaver::Handle::size() const {
  return idx_in.size();
}

bool Interleaver::Handle::has_mapping() const {
  return !idx_in.empty() && idx_in.size() == idx_out.size();
}

Interleaver::Handle Interleaver::build_from_shape(std::size_t rows, std::size_t cols,
                                                  std::size_t H, std::size_t W,
                                                  const std::string& kind)
{
  Handle h;
  if (rows % H != 0 || cols % W != 0) {
    throw std::invalid_argument("build_from_shape: rows/cols not divisible by H/W");
  }

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

Interleaver::Handle Interleaver::build_from_spec(std::size_t R, std::size_t C,
                                                 std::size_t H, std::size_t W,
                                                 const std::string& kind)
{
  return build_from_shape(R * H, C * W, H, W, kind);
}

} // namespace newcode
