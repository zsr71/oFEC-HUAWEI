#include "newcode/chase255.hpp"

#include "newcode/bch_255_239.hpp"
#include "newcode/params.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <cstdint>
#include <limits>
#include <vector>

namespace newcode {
namespace detail255 {

constexpr int BCH_N = 255;

template <typename LLR>
static void find_least_reliable(const LLR* llrs,
                                int L,
                                std::vector<int>& pos,
                                std::vector<float>& absval)
{
  struct Node { float a; int i; };
  std::vector<Node> nodes;
  nodes.reserve(BCH_N);
  for (int i = 0; i < BCH_N; ++i) {
    nodes.push_back({ std::fabs(static_cast<float>(llrs[i])), i });
  }

  const int take = std::min(L, static_cast<int>(nodes.size()));
  if (take < static_cast<int>(nodes.size())) {
    std::nth_element(nodes.begin(), nodes.begin() + take, nodes.end(),
                     [](const Node& x, const Node& y){ return x.a < y.a; });
    std::sort(nodes.begin(), nodes.begin() + take,
              [](const Node& x, const Node& y){ return x.a < y.a; });
  } else {
    std::sort(nodes.begin(), nodes.end(),
              [](const Node& x, const Node& y){ return x.a < y.a; });
  }

  pos.resize(take);
  absval.resize(take);
  for (int k = 0; k < take; ++k) {
    pos[k] = nodes[k].i;
    absval[k] = nodes[k].a;
  }
}

static void gen_test_patterns(int L, int n_test, std::vector<std::vector<bool>>& patt)
{
  if (n_test <= 0) n_test = 1;
  patt.assign(n_test, std::vector<bool>(std::max(L,0), false));

  if (L >= 0 && L < 31) {
    const int full = 1 << L;
    if (n_test <= full) {
      for (int c = 0; c < n_test; ++c)
        for (int j = 0; j < L; ++j)
          patt[c][j] = ((c >> j) & 1) != 0;
      return;
    }
  }

  int c = 0;
  if (c < n_test) c++;
  for (int layer = 1; c < n_test && layer <= L; ++layer) {
    for (int j = layer - 1; j < L && c < n_test; ++j) {
      std::fill(patt[c].begin(), patt[c].end(), false);
      for (int k = 0; k < layer - 1 && k < L; ++k) patt[c][k] = true;
      if (j < L) patt[c][j] = true;
      ++c;
    }
  }
}

} // namespace detail255

void chase_decode_255(const float* Lin255, float* extr255, const Params& p)
{
  using namespace detail255;

  const int L = std::max(1, p.CHASE_L);
  const int NTEST = std::max(1, p.CHASE_NTEST);
  const float beta = p.beta;
  const float alpha = p.ALPHA;

  std::array<float, BCH_N> y{};
  std::array<uint8_t, BCH_N> hard_ch{};
  for (int i = 0; i < BCH_N; ++i) {
    const float v = Lin255[i];
    y[i] = v;
    hard_ch[i] = (v >= 0.0f) ? 0u : 1u;
  }

  std::vector<int> lrp_pos;
  std::vector<float> lrp_abs;
  find_least_reliable(Lin255, L, lrp_pos, lrp_abs);
  const int L_eff = static_cast<int>(lrp_pos.size());

  std::vector<std::vector<bool>> patterns;
  gen_test_patterns(L_eff, NTEST, patterns);

  std::vector<std::vector<uint8_t>> candidates(NTEST, std::vector<uint8_t>(BCH_N, 0));
  struct Comp { float score; int idx; bool good; };
  std::vector<Comp> comps;
  comps.reserve(NTEST);

  std::array<float, BCH_N> abs_y{};
  for (int i = 0; i < BCH_N; ++i) abs_y[i] = std::fabs(y[i]);

  std::vector<uint8_t> tmp(BCH_N);
  std::array<uint8_t, BCH_N> decoded{};

  for (int c = 0; c < NTEST; ++c) {
    std::copy(hard_ch.begin(), hard_ch.end(), tmp.begin());
    for (int j = 0; j < L_eff; ++j) {
      if (patterns[c][j]) {
        const int pos = lrp_pos[j];
        if (pos >= 0 && pos < BCH_N) tmp[pos] ^= 1u;
      }
    }

    const bool ok = bch_255_239_decode_hiho_cw_255(tmp.data(), decoded.data());

    auto& CW = candidates[c];
    std::copy(decoded.begin(), decoded.end(), CW.begin());

    float dist = 0.f;
    for (int k = 0; k < BCH_N; ++k) {
      const uint8_t diff = (hard_ch[k] ^ CW[k]);
      dist += abs_y[k] * (diff ? 1.f : 0.f);
    }
    const float score = -dist;
    comps.push_back({score, c, ok});
  }

  int ml_idx = -1;
  float ml_score = -std::numeric_limits<float>::infinity();
  for (const auto& cp : comps) {
    if (!cp.good) continue;
    if (cp.score > ml_score) {
      ml_score = cp.score;
      ml_idx = cp.idx;
    }
  }

  std::vector<uint8_t> ML(BCH_N, 0);
  if (ml_idx >= 0) {
    ML = candidates[ml_idx];
  } else {
    ML.assign(hard_ch.begin(), hard_ch.end());
    ml_score = 0.f;
    for (int k = 0; k < BCH_N; ++k) {
      ml_score += y[k] * (ML[k] ? -1.f : +1.f);
    }
  }

  std::vector<float> omega(BCH_N, std::numeric_limits<float>::quiet_NaN());
  for (int j = 0; j < BCH_N; ++j) {
    int best_idx_plus = -1;
    float best_S_plus = -std::numeric_limits<float>::infinity();
    int best_idx_minus = -1;
    float best_S_minus = -std::numeric_limits<float>::infinity();

    for (const auto& cp : comps) {
      if (!cp.good) continue;
      const auto& C = candidates[cp.idx];
      if (C[j] == 0) {
        if (cp.score > best_S_plus) {
          best_S_plus = cp.score;
          best_idx_plus = cp.idx;
        }
      } else {
        if (cp.score > best_S_minus) {
          best_S_minus = cp.score;
          best_idx_minus = cp.idx;
        }
      }
    }

    if (best_idx_plus >= 0 && best_idx_minus >= 0) {
      const auto& Cplus = candidates[best_idx_plus];
      const auto& Cminus = candidates[best_idx_minus];
      float wj = 0.f;
      for (int l = 0; l < BCH_N; ++l) {
        if (l == j) continue;
        const float r_l = y[l];
        const float c_plus_l  = Cplus[l] ? -1.f : +1.f;
        const float c_minus_l = Cminus[l] ? -1.f : +1.f;
        const int p_l = (c_plus_l != c_minus_l) ? 1 : 0;
        wj += r_l * c_plus_l * static_cast<float>(p_l);
      }
      omega[j] = wj;
    } else {
      omega[j] = std::numeric_limits<float>::quiet_NaN();
    }

    if (std::isnan(omega[j])) {
      const float sgn = ML[j] ? -1.f : +1.f;
      omega[j] = beta * sgn-Lin255[j];
    }
  }

  for (int j = 0; j < BCH_N; ++j) {
    extr255[j] = alpha * omega[j];
  }
}

} // namespace newcode
