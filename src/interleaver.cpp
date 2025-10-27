#include "newcode/interleaver.hpp"
#include <array>
#include <algorithm>
#include <cassert>

namespace newcode {

// ---------------- 规范表（Table 7）：dst(i,j) <- src(sr,sc) ----------------
// P[dst] = src （线性 0..255，行主序）
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

// ---------- 构造可逆映射（显式 R,C,H,W） ----------
Interleaver Interleaver::build_from_spec(int R, int C, int H, int W) {
  if (R <= 0 || C <= 0 || H <= 0 || W <= 0)
    throw std::invalid_argument("build_from_spec: invalid dims");

  const size_t N = static_cast<size_t>(R) * C * H * W;
  Interleaver itv{R, C, H, W};
  itv.idx_out.resize(N);
  itv.idx_in.resize(N);

  // 1) 块内置换：准备 src->dst 的逆表（线性 0..255）
  const auto& dst2src = perm16_dst_to_src(); // dst <- src
  std::array<uint16_t, 16 * 16> src2dst{};
  for (int dst = 0; dst < 256; ++dst) {
    src2dst[ dst2src[dst] ] = static_cast<uint16_t>(dst);
  }

  // 2) 写入映射：输入序列位置 -> 缓冲线性位置
  std::vector<uint32_t> idx_buf(N);        // in_pos -> buffer_pos
  size_t in_pos = 0;
  for (int r = 0; r < R; ++r) {
    for (int c = 0; c < C; ++c) {
      for (int sr = 0; sr < H; ++sr) {
        for (int sc = 0; sc < W; ++sc) {
          const int src_lin = sr * W + sc;              // 0..255
          const int dst_lin = src2dst[src_lin];         // 0..255
          const int dr = dst_lin / W;
          const int dc = dst_lin % W;
          idx_buf[in_pos++] = buf_lin(r, c, dr, dc, C, H, W);
        }
      }
    }
  }
  assert(in_pos == N);

  // 3) 读出顺序（列优先 + 四子集轮询）：buffer_pos 序列
  std::vector<int> S0, S1, S2, S3;
  S0.reserve(R/2); S1.reserve(R/2); S2.reserve(R/2); S3.reserve(R/2);
  for (int r = 0; r < R/2; ++r)    ((r % 2) == 0 ? S0 : S1).push_back(r);       // 上半
  for (int r = R/2; r < R; ++r)    (((r % 2) == 0) ? S2 : S3).push_back(r);     // 下半

  std::vector<int> order_rows; order_rows.reserve(R);
  size_t t = 0;
  while (order_rows.size() < static_cast<size_t>(R)) {
    if (t < S0.size()) order_rows.push_back(S0[t]);
    if (t < S1.size()) order_rows.push_back(S1[t]);
    if (t < S2.size()) order_rows.push_back(S2[t]);
    if (t < S3.size()) order_rows.push_back(S3[t]);
    ++t;
  }

  std::vector<uint32_t> idx_read(N);       // out_pos -> buffer_pos
  size_t out_pos = 0;
  for (int k = 0; k < C * W; ++k) {
    const int sc = k / W;       // 块列
    const int ic = k % W;       // 小块内列
    for (int rr : order_rows) {
      for (int i = 0; i < H; ++i) {
        idx_read[out_pos++] = buf_lin(rr, sc, i, ic, C, H, W);
      }
    }
  }
  assert(out_pos == N);

  // 4) 合成最终排列：y[pos] = x[idx_in[pos]]
  //    即：idx_in[pos] = which input wrote the same buffer_pos as idx_read[pos]
  std::vector<uint32_t> buf_to_in(N);  // buffer_pos -> in_pos
  for (size_t p = 0; p < N; ++p) buf_to_in[ idx_buf[p] ] = static_cast<uint32_t>(p);
  for (size_t p = 0; p < N; ++p) {
    itv.idx_out[p] = static_cast<uint32_t>(p);               // 恒等
    itv.idx_in [p] = buf_to_in[ idx_read[p] ];
  }

  // 自检：应为排列
  std::vector<uint32_t> check(N);
  for (size_t i = 0; i < N; ++i) check[ itv.idx_in[i] ] = static_cast<uint32_t>(i);
  for (size_t i = 0; i < N; ++i) assert(check[i] < N);

  return itv;
}

// ---------- 由矩阵形状构建 ----------
Interleaver Interleaver::build_from_shape(int rows, int cols, int H, int W) {
  if (rows % H != 0 || cols % W != 0)
    throw std::invalid_argument("build_from_shape: rows/cols not divisible by H/W");
  const int R = rows / H;
  const int C = cols / W;
  return Interleaver::build_from_spec(R, C, H, W);
}

} // namespace newcode
