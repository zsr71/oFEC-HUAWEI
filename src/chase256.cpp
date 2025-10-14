// chase256.cpp
// ------------------------------------------------------------
// - nth_element 边界修复
// - 分层翻转模式（接近 AFF3CT）
// - 预计算 |Lch| 做 metric
// - 写回时对整型 LLR 饱和
// - 固定 Chase-Pyndiah 系数 (A,B,C,D,E)
// - ★ 新增三参 API：chase_decode_256(Lin, Lch, Y2, p)
//   * metric/LRP 用 Lch（纯通道）
//   * 减项用 -A*Lin（本轮先验 = Lch + B*E_other + C）
//   * metric 仅累计 0..254，overall 位(255)的外信息恒设为 0
// ------------------------------------------------------------

#include "newcode/chase256.hpp"
#include "newcode/bch_255_239.hpp"
#include "newcode/params.hpp"
#include "newcode/qfloat.hpp"

#include <vector>
#include <algorithm>
#include <cmath>
#include <limits>
#include <cstdint>
#include <type_traits>

namespace newcode {
namespace detail {

constexpr int BCH_N_TOTAL = 256; // 255 + overall parity
constexpr int BCH_N_CORE  = 255;
constexpr int PAR_IDX     = 255;

// -------- LLR 类型适配：to-float / from-float --------
template<typename LLR>
inline float llr_to_float(LLR x) { return static_cast<float>(x); }

// 浮点：直转
template<typename LLR, typename std::enable_if<std::is_floating_point<LLR>::value, int>::type = 0>
inline LLR llr_from_float(float x) { return static_cast<LLR>(x); }

// 有符号整型：饱和 + 四舍五入
template<typename LLR, typename std::enable_if<std::is_integral<LLR>::value && std::is_signed<LLR>::value, int>::type = 0>
inline LLR llr_from_float(float x)
{
    const float lo = (float)std::numeric_limits<LLR>::min();
    const float hi = (float)std::numeric_limits<LLR>::max();
    if (x < lo) x = lo;
    if (x > hi) x = hi;
    return (LLR)std::lrintf(x);
}

// 非算术（如 qfloat）：走自定义 from_float()
template<typename LLR, typename std::enable_if<!std::is_arithmetic<LLR>::value, int>::type = 0>
inline LLR llr_from_float(float x) { return LLR::from_float(x); }

// ---- 选取（A,B,C,D,E）系数 ----
static inline void pick_cp(const Params& p, float& A, float& B, float& C, float& D, int& E)
{
    A = p.CP_A;
    B = p.CP_B;
    C = p.CP_C;
    D = p.CP_D;
    E = p.CP_E;
}

// ---- LRP：在 0..254 里选 L 个最小 |Lch| ----
template<typename LLR>
static void find_least_reliable(const LLR* LLR_list, int L,
                                std::vector<int>& pos, std::vector<float>& absval)
{
    struct Node { float a; int i; };
    std::vector<Node> v; v.reserve(BCH_N_CORE);
    for (int i = 0; i < BCH_N_CORE; ++i)
        v.push_back({ std::fabs(llr_to_float(LLR_list[i])), i });

    const int take = std::min(L, (int)v.size());
    if (take < (int)v.size())
        std::nth_element(v.begin(), v.begin() + take, v.end(),
                         [](const Node& x, const Node& y){ return x.a < y.a; });
    std::sort(v.begin(), v.begin() + take, [](const Node& x, const Node& y){ return x.a < y.a; });

    pos.resize(take);
    absval.resize(take);
    for (int k = 0; k < take; ++k) { pos[k] = v[k].i; absval[k] = v[k].a; }
}

// ---- 生成翻转模式（分层推进，覆盖更接近 AFF3CT）----
static void gen_test_patterns(int L, int n_test, std::vector<std::vector<bool>>& patt)
{
    if (n_test <= 0) n_test = 1;
    patt.assign(n_test, std::vector<bool>(L, false));

    const bool can_full = (L >= 0 && L < 31);
    if (can_full) {
        const int full = 1 << L;
        if (n_test <= full) {
            for (int c = 0; c < n_test; ++c)
                for (int j = 0; j < L; ++j)
                    patt[c][j] = ((c >> j) & 1) != 0;
            return;
        }
    }

    int c = 0;
    if (c < n_test) c++; // 0 翻：全 false

    for (int layer = 1; c < n_test && layer <= L; ++layer)
        for (int j = layer - 1; j < L && c < n_test; ++j) {
            std::fill(patt[c].begin(), patt[c].end(), false);
            for (int k = 0; k < layer - 1; ++k) patt[c][k] = true;
            patt[c][j] = true;
            c++;
        }
}

// ---- overall 偶校验位（让 256 位异或和为 0）----
inline uint8_t parity256_from255(const uint8_t* cw255) {
    uint8_t acc = 0;
    for (int i = 0; i < BCH_N_CORE; ++i) acc ^= cw255[i];
    return acc;
}

// ---- metric（仅 0..254；用 |Lch|）----
static inline float chase_metric_core(const uint8_t* hard_ch_256, const uint8_t* DW_256, const float* absLc)
{
    float m = 0.f;
    for (int i = 0; i < BCH_N_CORE; ++i)
        if ((hard_ch_256[i] ^ DW_256[i]) & 1u) m += absLc[i];
    return m;
}

} // namespace detail

// ======================== 主函数（3参：Lin + Lch） ========================
template<typename LLR>
void chase_decode_256(const LLR* Lin256,   // 本轮先验 = Lch + B*E_other (+C)
                      const LLR* Lch256,   // 纯通道
                      LLR* Y2_256,         // 输出纯外信息
                      const Params& p)
{
    using namespace detail;

    float A,B,C,D; int E;
    pick_cp(p, A,B,C,D,E);

    const int L      = std::max(1, p.CHASE_L);
    const int NTEST  = std::max(1, p.CHASE_NTEST);
    const int NCOMP  = std::max(1, std::min(p.CHASE_NCOMP, NTEST));

    uint8_t hard_ch[BCH_N_TOTAL];
    std::vector<float> absLc(BCH_N_TOTAL);
    for (int i = 0; i < BCH_N_TOTAL; ++i) {
        const float lc = llr_to_float(Lin256[i]);
        hard_ch[i] = (lc >= 0.f) ? 0u : 1u;
        absLc[i]   = std::fabs(lc);
    }

    // LRP（基于 |Lch|）
    std::vector<int>   lrp_pos;  lrp_pos.reserve(L);
    std::vector<float> lrp_abs;  lrp_abs.reserve(L);
    find_least_reliable(Lin256, L, lrp_pos, lrp_abs);
    const int L_eff = (int)lrp_pos.size();

    // 测试模式
    std::vector<std::vector<bool>> patt;
    gen_test_patterns(L_eff, NTEST, patt);

    // 生成候选 + metric
    std::vector<std::vector<uint8_t>> DW_all(NTEST, std::vector<uint8_t>(BCH_N_TOTAL, 0));
    struct Comp { float metric; int idx; };
    std::vector<Comp> comps; comps.reserve(NTEST);

    std::vector<uint8_t> hcand_256(BCH_N_TOTAL);
    std::vector<uint8_t> out255(BCH_N_CORE);

    for (int c = 0; c < NTEST; ++c)
    {
        std::copy(hard_ch, hard_ch + BCH_N_TOTAL, hcand_256.begin());
        for (int j = 0; j < L_eff; ++j)
            if (patt[c][j]) hcand_256[ lrp_pos[j] ] ^= 1u;

        bool ok = bch_255_239_decode_hiho_cw_255(hcand_256.data(), out255.data());
        if (!ok) continue;

        auto& DW = DW_all[c];
        std::copy(out255.begin(), out255.end(), DW.begin());
        DW[PAR_IDX] = parity256_from255(DW.data());

        float m = chase_metric_core(hard_ch, DW.data(), absLc.data());
        comps.push_back({m, c});
    }

    if (comps.empty()) {
        // 回退：用 D*|Lch| 的符号外信息，再减 A*Lin；overall 外信息设 0
        for (int i = 0; i < BCH_N_CORE; ++i) {
            float rel = D * absLc[i];
            if (hard_ch[i]) rel = -rel;
            const float apr = llr_to_float(Lin256[i]) - llr_to_float(Lch256[i]); // 只取先验部分
            const float y2  = rel - A * apr;            Y2_256[i] = llr_from_float<LLR>(y2);
        }
        Y2_256[PAR_IDX] = llr_from_float<LLR>(0.0f);
        return;
    }

    std::sort(comps.begin(), comps.end(),
              [](const Comp& a, const Comp& b){ return a.metric < b.metric; });
    const int NG = std::min((int)comps.size(), NCOMP);
    const auto& DW0 = DW_all[ comps[0].idx ];
    const float    M0  = comps[0].metric;

    std::vector<float> comp_delta(NG, 0.f);
    for (int j = 1; j < NG; ++j) comp_delta[j] = (comps[j].metric - M0) * B;

    const int e_cap = (E > 0) ? std::min(E, L_eff) : L_eff;
    float beta = 0.f;
    for (int i = 0; i < e_cap; ++i) beta += lrp_abs[i];
    beta -= C * M0;
    beta = std::max(0.f, beta - C * M0);

    for (int i = 0; i < BCH_N_CORE; ++i)
    {
        const uint8_t db = DW0[i];
        int j = 1;
        for (; j < NG; ++j) {
            const auto& DWj = DW_all[ comps[j].idx ];
            if (DWj[i] != db) break;
        }
        float reliability = 0.f;
        if (j < NG) {
            reliability = comp_delta[j];
        } else {
            reliability = beta + D * absLc[i];
            if (reliability < 0.f) reliability = 0.f;
        }
        if (db) reliability = -reliability;

        const float apr = llr_to_float(Lin256[i]) - llr_to_float(Lch256[i]);
        const float y2  = reliability - A * apr;
        Y2_256[i] = llr_from_float<LLR>(y2);
    }
    Y2_256[PAR_IDX] = llr_from_float<LLR>(0.0f);
}

// ======================== 兼容旧接口（2参） ========================
template<typename LLR>
void chase_decode_256(const LLR* Y256, LLR* Y2_256, const Params& p)
{
    chase_decode_256<LLR>(Y256, Y256, Y2_256, p);
}

// ===== 显式实例化 =====
template void chase_decode_256<float >(const float*,  const float*,  float*,  const Params&);
template void chase_decode_256<int8_t>(const int8_t*, const int8_t*, int8_t*, const Params&);
template void chase_decode_256<newcode::qfloat<4>>(const newcode::qfloat<4>*,
                                                   const newcode::qfloat<4>*,
                                                   newcode::qfloat<4>*,
                                                   const Params&);
template void chase_decode_256<newcode::qfloat<5>>(const newcode::qfloat<5>*,
                                                   const newcode::qfloat<5>*,
                                                   newcode::qfloat<5>*,
                                                   const Params&);

// 旧接口实例化
template void chase_decode_256<float >(const float*,  float*,  const Params&);
template void chase_decode_256<int8_t>(const int8_t*, int8_t*, const Params&);
template void chase_decode_256<newcode::qfloat<4>>(const newcode::qfloat<4>*,
                                                   newcode::qfloat<4>*,
                                                   const Params&);
template void chase_decode_256<newcode::qfloat<5>>(const newcode::qfloat<5>*,
                                                   newcode::qfloat<5>*,
                                                   const Params&);

} // namespace newcode
