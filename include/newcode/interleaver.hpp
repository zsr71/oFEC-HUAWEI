#pragma once
#include <cstdint>
#include <vector>
#include <stdexcept>

namespace newcode {

// 交织器：按 R×C 个 H×W 小块组织（OpenROADM oFEC: R=84, C=8, H=W=16）
struct Interleaver {
  int R;   // 小块行数
  int C;   // 小块列数
  int H;   // 小块内行
  int W;   // 小块内列

  // y[pos] = x[idx_in[pos]]; idx_out 作为对称性保留（恒等 0..N-1）
  std::vector<uint32_t> idx_in;   // 输入位置映射
  std::vector<uint32_t> idx_out;  // 输出位置（恒等）

  // 从规范的 16×16 表（dst <- src）构建（显式指定 R,C,H,W）
  static Interleaver build_from_spec(int R, int C, int H, int W);

  // 从编码矩阵形状构建（rows = R*H, cols = C*W；H/W 默认为 16）
  static Interleaver build_from_shape(int rows, int cols, int H = 16, int W = 16);

  // 单块 bit 数
  inline size_t size() const { return static_cast<size_t>(R) * C * H * W; }
  inline size_t block_size() const { return size(); }

  // —— 单块交织/解交织 —— //
  template <typename T>
  std::vector<T> interleave(const std::vector<T>& x) const {
    if (x.size() != size()) throw std::runtime_error("interleave: size mismatch");
    std::vector<T> y(x.size());
    for (size_t i = 0; i < x.size(); ++i) y[i] = x[ static_cast<size_t>(idx_in[i]) ];
    return y;
  }

  template <typename T>
  std::vector<T> deinterleave(const std::vector<T>& y) const {
    if (y.size() != size()) throw std::runtime_error("deinterleave: size mismatch");
    std::vector<T> x(y.size());
    for (size_t i = 0; i < y.size(); ++i) x[ static_cast<size_t>(idx_in[i]) ] = y[i];
    return x;
  }

  // —— 多块交织/解交织（输入长度必须是 block_size 的整数倍）—— //
  template <typename T>
  std::vector<T> interleave_chunks(const std::vector<T>& x) const {
    const size_t B = block_size();
    if (x.size() % B != 0) throw std::runtime_error("interleave_chunks: size not multiple of block");
    std::vector<T> y(x.size());
    const size_t nblk = x.size() / B;
    for (size_t b = 0; b < nblk; ++b) {
      const size_t off = b * B;
      for (size_t i = 0; i < B; ++i)
        y[off + i] = x[off + static_cast<size_t>(idx_in[i])];
    }
    return y;
  }

  template <typename T>
  std::vector<T> deinterleave_chunks(const std::vector<T>& y) const {
    const size_t B = block_size();
    if (y.size() % B != 0) throw std::runtime_error("deinterleave_chunks: size not multiple of block");
    std::vector<T> x(y.size());
    const size_t nblk = y.size() / B;
    for (size_t b = 0; b < nblk; ++b) {
      const size_t off = b * B;
      for (size_t i = 0; i < B; ++i)
        x[off + static_cast<size_t>(idx_in[i])] = y[off + i];
    }
    return x;
  }
};

// 返回规范 16×16 的“dst←src”表（线性 0..255，行主序）
// P[dst_linear] = src_linear
const std::vector<uint16_t>& perm16_dst_to_src();

// 工具：把 4D 坐标映射成 0-based 线性下标
inline uint32_t buf_lin(int r, int c, int i, int j, int C, int H, int W) {
  return static_cast<uint32_t>(((r * C) + c) * H * W + i * W + j);
}

} // namespace newcode
