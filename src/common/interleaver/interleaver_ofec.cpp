#include "newcode/common/interleaver/interleaver.hpp"

#include <algorithm>
#include <array>
#include <stdexcept>
#include <vector>

namespace interleaver {
namespace {

using PermTable = std::vector<int>;

static std::size_t buf_lin(int R, int C, int r, int c, int C_tot, int H, int W) {
  const std::size_t block_index = static_cast<std::size_t>(R) * static_cast<std::size_t>(C_tot)
                                + static_cast<std::size_t>(C);
  const std::size_t inner_index = static_cast<std::size_t>(r) * static_cast<std::size_t>(W)
                                + static_cast<std::size_t>(c);
  return block_index * static_cast<std::size_t>(H * W) + inner_index;
}

const std::vector<uint16_t>& perm16_dst_to_src() {
  static std::vector<uint16_t> P;
  if (!P.empty()) return P;

  static const int sr[16][16] = {
      { 0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15 },
      { 14,15,0,1,2,3,4,5,6,7,8,9,10,11,12,13 },
      { 12,13,14,15,0,1,2,3,4,5,6,7,8,9,10,11 },
      { 10,11,12,13,14,15,0,1,2,3,4,5,6,7,8,9 },
      { 8,9,10,11,12,13,14,15,0,1,2,3,4,5,6,7 },
      { 6,7,8,9,10,11,12,13,14,15,0,1,2,3,4,5 },
      { 4,5,6,7,8,9,10,11,12,13,14,15,0,1,2,3 },
      { 2,3,4,5,6,7,8,9,10,11,12,13,14,15,0,1 },
      { 15,0,1,2,3,4,5,6,7,8,9,10,11,12,13,14 },
      { 13,14,15,0,1,2,3,4,5,6,7,8,9,10,11,12 },
      { 11,12,13,14,15,0,1,2,3,4,5,6,7,8,9,10 },
      { 9,10,11,12,13,14,15,0,1,2,3,4,5,6,7,8 },
      { 7,8,9,10,11,12,13,14,15,0,1,2,3,4,5,6 },
      { 5,6,7,8,9,10,11,12,13,14,15,0,1,2,3,4 },
      { 3,4,5,6,7,8,9,10,11,12,13,14,15,0,1,2 },
      { 1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,0 }
  };
  static const int sc[16][16] = {
      { 0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15 },
      { 15,0,1,2,3,4,5,6,7,8,9,10,11,12,13,14 },
      { 14,15,0,1,2,3,4,5,6,7,8,9,10,11,12,13 },
      { 13,14,15,0,1,2,3,4,5,6,7,8,9,10,11,12 },
      { 12,13,14,15,0,1,2,3,4,5,6,7,8,9,10,11 },
      { 11,12,13,14,15,0,1,2,3,4,5,6,7,8,9,10 },
      { 10,11,12,13,14,15,0,1,2,3,4,5,6,7,8,9 },
      { 9,10,11,12,13,14,15,0,1,2,3,4,5,6,7,8 },
      { 7,8,9,10,11,12,13,14,15,0,1,2,3,4,5,6 },
      { 6,7,8,9,10,11,12,13,14,15,0,1,2,3,4,5 },
      { 5,6,7,8,9,10,11,12,13,14,15,0,1,2,3,4 },
      { 4,5,6,7,8,9,10,11,12,13,14,15,0,1,2,3 },
      { 3,4,5,6,7,8,9,10,11,12,13,14,15,0,1,2 },
      { 2,3,4,5,6,7,8,9,10,11,12,13,14,15,0,1 },
      { 1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,0 },
      { 0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15 }
  };

  P.resize(16 * 16);
  for (int r = 0; r < 16; ++r) {
    for (int c = 0; c < 16; ++c) {
      const int dst = r * 16 + c;
      const int src = sr[r][c] * 16 + sc[r][c];
      P[dst] = static_cast<uint16_t>(src);
    }
  }
  return P;
}

template <typename T>
void apply_permutation(const std::vector<T>& src, std::vector<T>& dst, const PermTable& perm) {
  dst.resize(perm.size());
  for (std::size_t i = 0; i < perm.size(); ++i) {
    dst[i] = src[static_cast<std::size_t>(perm[i])];
  }
}

void matrix_from_flat(const std::vector<float>& flat, std::size_t rows, std::size_t cols,
                      newcode::Matrix<float>& out) {
  out = newcode::Matrix<float>(rows, cols);
  for (std::size_t r = 0; r < rows; ++r)
    for (std::size_t c = 0; c < cols; ++c)
      out[r][c] = flat[r * cols + c];
}

void matrix_to_flat(const newcode::Matrix<float>& in, std::vector<float>& flat) {
  flat.resize(in.rows() * in.cols());
  for (std::size_t r = 0; r < in.rows(); ++r)
    for (std::size_t c = 0; c < in.cols(); ++c)
      flat[r * in.cols() + c] = in[r][c];
}

class OfecInterleaver final : public IInterleaver {
public:
  OfecInterleaver(std::size_t R, std::size_t C, std::size_t H, std::size_t W) {
    if (R == 0 || C == 0 || H == 0 || W == 0)
      throw std::invalid_argument("OfecInterleaver: invalid dimension");
    build_mapping(R, C, H, W);
    rows_ = R * H;
    cols_ = C * W;
  }

  void interleave(const newcode::Matrix<float>& in, newcode::Matrix<float>& out) const override {
    if (forward_.empty()) { out = in; return; }
    std::vector<float> flat;
    matrix_to_flat(in, flat);
    std::vector<float> perm;
    apply_permutation(flat, perm, forward_);
    matrix_from_flat(perm, in.rows(), in.cols(), out);
  }

  void deinterleave(const newcode::Matrix<float>& in, newcode::Matrix<float>& out) const override {
    if (inverse_.empty()) { out = in; return; }
    std::vector<float> flat;
    matrix_to_flat(in, flat);
    std::vector<float> perm;
    apply_permutation(flat, perm, inverse_);
    matrix_from_flat(perm, in.rows(), in.cols(), out);
  }

  const std::vector<int>& forward_mapping() const override { return forward_; }
  const std::vector<int>& inverse_mapping() const override { return inverse_; }

private:
  void build_mapping(std::size_t R, std::size_t C, std::size_t H, std::size_t W) {
    const std::size_t N = R * C * H * W;
    forward_.resize(N);
    inverse_.resize(N);

    const auto& dst2src = perm16_dst_to_src();
    std::array<uint16_t, 16 * 16> src2dst{};
    for (int dst = 0; dst < 256; ++dst)
      src2dst[dst2src[dst]] = static_cast<uint16_t>(dst);

    std::vector<std::size_t> idx_buf(N);
    std::size_t in_pos = 0;
    for (std::size_t r = 0; r < R; ++r) {
      for (std::size_t c = 0; c < C; ++c) {
        for (std::size_t sr = 0; sr < H; ++sr) {
          for (std::size_t sc = 0; sc < W; ++sc) {
            const int src_lin = static_cast<int>(sr * W + sc);
            const int dst_lin = src2dst[static_cast<std::size_t>(src_lin)];
            const int dr = dst_lin / static_cast<int>(W);
            const int dc = dst_lin % static_cast<int>(W);
            idx_buf[in_pos++] = buf_lin(static_cast<int>(r), static_cast<int>(c),
                                        dr, dc, static_cast<int>(C), static_cast<int>(H), static_cast<int>(W));
          }
        }
      }
    }

    std::vector<int> S0, S1, S2, S3;
    S0.reserve(R / 2); S1.reserve(R / 2); S2.reserve(R / 2); S3.reserve(R / 2);
    for (std::size_t r = 0; r < R / 2; ++r)    ((r % 2) == 0 ? S0 : S1).push_back(static_cast<int>(r));
    for (std::size_t r = R / 2; r < R; ++r)    (((r % 2) == 0) ? S2 : S3).push_back(static_cast<int>(r));

    std::vector<int> order_rows; order_rows.reserve(R);
    std::size_t t = 0;
    while (order_rows.size() < R) {
      if (t < S0.size()) order_rows.push_back(S0[t]);
      if (t < S1.size()) order_rows.push_back(S1[t]);
      if (t < S2.size()) order_rows.push_back(S2[t]);
      if (t < S3.size()) order_rows.push_back(S3[t]);
      ++t;
    }

    std::vector<std::size_t> idx_read(N);
    std::size_t out_pos = 0;
    for (std::size_t k = 0; k < C * W; ++k) {
      const int sc = static_cast<int>(k / W);
      const int ic = static_cast<int>(k % W);
      for (int rr : order_rows) {
        for (std::size_t i = 0; i < H; ++i) {
          idx_read[out_pos++] = buf_lin(rr, sc, static_cast<int>(i), ic,
                                        static_cast<int>(C), static_cast<int>(H), static_cast<int>(W));
        }
      }
    }

    std::vector<std::size_t> buf_to_in(N);
    for (std::size_t p = 0; p < N; ++p)
      buf_to_in[idx_buf[p]] = p;
    for (std::size_t p = 0; p < N; ++p) {
      forward_[p] = static_cast<int>(buf_to_in[idx_read[p]]);
    }
    for (std::size_t p = 0; p < N; ++p) {
      const std::size_t src_idx = static_cast<std::size_t>(forward_[p]);
      inverse_[src_idx] = static_cast<int>(p);
    }
  }

  std::size_t rows_ = 0;
  std::size_t cols_ = 0;
  PermTable forward_;
  PermTable inverse_;
};

} // namespace

std::unique_ptr<IInterleaver> make_interleaver(const std::string& name) {
  if (name == "ofec")
    throw std::invalid_argument("make_interleaver: ofec requires shape (R,C,H,W)");
  return nullptr;
}

std::unique_ptr<IInterleaver> make_interleaver(const std::string& name,
                                               std::size_t R, std::size_t C,
                                               std::size_t H, std::size_t W) {
  if (name != "ofec") return nullptr;
  return std::make_unique<OfecInterleaver>(R, C, H, W);
}

} // namespace newcode
