// chase256.cpp
// ------------------------------------------------------------
// - Fix: nth_element boundary when L == v.size()
// - Improve: gen_test_patterns coverage (layered growth like AFF3CT)
// - Improve: precompute |LLR| for metric
// - Improve: integral LLR saturation on write-back
// - NEW: iteration schedule (A,B,C,D,E) via Params::CP_*_SCHED and CP_ITER
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
#include <numeric>

namespace newcode {
namespace detail {

constexpr int BCH_N_TOTAL = 256; // 255 + overall parity
constexpr int BCH_N_CORE  = 255;
constexpr int PAR_IDX     = 255;

// -------- LLR 类型适配：to-float / from-float --------

// to-float：要求 LLR 支持 explicit operator float()
template<typename LLR>
inline float llr_to_float(LLR x) { return static_cast<float>(x); }

// from-float：浮点 → 直接转换
template<typename LLR, typename std::enable_if<std::is_floating_point<LLR>::value, int>::type = 0>
inline LLR llr_from_float(float x) { return static_cast<LLR>(x); }

// from-float：有符号整型 → 饱和 + 四舍五入
template<typename LLR, typename std::enable_if<std::is_integral<LLR>::value && std::is_signed<LLR>::value, int>::type = 0>
inline LLR llr_from_float(float x)
{
    const float lo = (float)std::numeric_limits<LLR>::min();
    const float hi = (float)std::numeric_limits<LLR>::max();
    if (x < lo) x = lo;
    if (x > hi) x = hi;
    return (LLR)std::lrintf(x);
}

// from-float：非算术类型（如 qfloat）→ 走自定义 from_float()
template<typename LLR, typename std::enable_if<!std::is_arithmetic<LLR>::value, int>::type = 0>
inline LLR llr_from_float(float x) { return LLR::from_float(x); }

// ---- 选取（A,B,C,D,E）系数（支持迭代日程）----
static inline void pick_cp(const Params& p, float& A, float& B, float& C, float& D, int& E)
{
    if (p.CP_USE_SCHED && p.CP_SCHED_LEN > 0)
    {
        const int it = std::min(std::max(0, p.CP_ITER), p.CP_SCHED_LEN - 1);
        A = p.CP_A_SCHED[it];
        B = p.CP_B_SCHED[it];
        C = p.CP_C_SCHED[it];
        D = p.CP_D_SCHED[it];
        E = p.CP_E_SCHED[it];
    }
    else
    {
        A = p.CP_A;
        B = p.CP_B;
        C = p.CP_C;
        D = p.CP_D;
        E = p.CP_E;
    }
}

// ---- 最不可靠位置（在 0..254 中选择 L 个）----
template<typename LLR>
static void find_least_reliable(const LLR* Y256, int L,
                                std::vector<int>& pos, std::vector<float>& absval)
{
    struct Node { float a; int i; };
    std::vector<Node> v; v.reserve(BCH_N_CORE);
    for (int i = 0; i < BCH_N_CORE; ++i)
        v.push_back({ std::fabs(llr_to_float(Y256[i])), i });

    const int take = std::min(L, (int)v.size());
    if (take < (int)v.size()) // 只有当 nth 在 [begin, end) 内时才调用
        std::nth_element(v.begin(), v.begin() + take, v.end(),
                         [](const Node& x, const Node& y){ return x.a < y.a; });
    std::sort(v.begin(), v.begin() + take, [](const Node& x, const Node& y){ return x.a < y.a; });

    pos.resize(take);
    absval.resize(take);
    for (int k = 0; k < take; ++k) { pos[k] = v[k].i; absval[k] = v[k].a; }
}

// ---- 生成翻转模式（分层推进，覆盖更接近 AFF3CT） ----
static void gen_test_patterns(int L, int n_test, std::vector<std::vector<bool>>& patt)
{
    if (n_test <= 0) n_test = 1;
    patt.assign(n_test, std::vector<bool>(L, false));

    // 如果要求完整枚举且可行，则直接枚举（只生成前 n_test 个）
    const bool can_full = (L >= 0 && L < 31); // 防止 (1<<L) 溢出
    if (can_full)
    {
        const int full = 1 << L;
        if (n_test <= full)
        {
            for (int c = 0; c < n_test; ++c)
                for (int j = 0; j < L; ++j)
                    patt[c][j] = ((c >> j) & 1) != 0;
            return;
        }
    }

    // 分层策略：0翻 → 1翻 → “固定前k-1位为真 + 在其余位上单翻” → … 直到凑够 n_test
    int c = 0;

    // (0) 不翻
    if (c < n_test) c++; // patt[0] 默认全 false

    // layer 从 1 开始（layer 表示此层总共翻转的位数）
    for (int layer = 1; c < n_test && layer <= L; ++layer)
    {
        for (int j = layer - 1; j < L && c < n_test; ++j)
        {
            std::fill(patt[c].begin(), patt[c].end(), false);
            for (int k = 0; k < layer - 1; ++k) patt[c][k] = true; // 固定 0..(layer-2)
            patt[c][j] = true; // 再在第 j 位上单翻
            c++;
        }
    }
}

// ---- 计算整体偶校验位（让 256 位异或和为 0）----
inline uint8_t parity256_from255(const uint8_t* cw255) {
    uint8_t acc = 0;
    for (int i = 0; i < BCH_N_CORE; ++i) acc ^= cw255[i];
    return acc;
}

// ---- Chase 度量（预计算 absY 版本）----
static inline float chase_metric_fast(const uint8_t* hard0_256, const uint8_t* DW_256, const float* absY)
{
    float m = 0.f;
    for (int i = 0; i < BCH_N_TOTAL; ++i)
        if ((hard0_256[i] ^ DW_256[i]) & 1u)
            m += absY[i];
    return m;
}

} // namespace detail

// ======================== 主函数（流程 2~6） ========================
template<typename LLR>
void chase_decode_256(const LLR* Y256, LLR* Y2_256, const Params& p)
{
    using namespace detail;

    // —— 依据迭代日程选 a..e —— //
    float A,B,C,D; int E;
    pick_cp(p, A,B,C,D,E);

    const int L      = std::max(1, p.CHASE_L);
    const int NTEST  = std::max(1, p.CHASE_NTEST);
    const int NCOMP  = std::max(1, std::min(p.CHASE_NCOMP, NTEST));

    // 硬判决基线（含整体奇偶位）
    uint8_t hard0[BCH_N_TOTAL];
    for (int i = 0; i < BCH_N_TOTAL; ++i)
        hard0[i] = (llr_to_float(Y256[i]) >= 0.f) ? 0u : 1u;

    // 预计算 |LLR|
    std::vector<float> absY(BCH_N_TOTAL);
    for (int i = 0; i < BCH_N_TOTAL; ++i)
        absY[i] = std::fabs(llr_to_float(Y256[i]));

    // (2) LRP
    std::vector<int>   lrp_pos;  lrp_pos.reserve(L);
    std::vector<float> lrp_abs;  lrp_abs.reserve(L);
    find_least_reliable(Y256, L, lrp_pos, lrp_abs);
    const int L_eff = (int)lrp_pos.size();

    // (3) patt
    std::vector<std::vector<bool>> patt;
    gen_test_patterns(L_eff, NTEST, patt);

    // (4) 候选：翻转 -> BCH(255) 硬译码 -> 组装 256 -> metric
    std::vector<uint8_t> DW_all(NTEST * BCH_N_TOTAL, 0);
    struct Comp { float metric; int idx; };
    std::vector<Comp> comps; comps.reserve(NTEST);

    std::vector<uint8_t> hcand_256(BCH_N_TOTAL);
    std::vector<uint8_t> out255(BCH_N_CORE);

    for (int c = 0; c < NTEST; ++c)
    {
        // 4a 硬判决翻转（仅 0..254）
        std::copy(hard0, hard0 + BCH_N_TOTAL, hcand_256.begin());
        for (int j = 0; j < L_eff; ++j)
            if (patt[c][j]) hcand_256[ lrp_pos[j] ] ^= 1u;

        // 4b BCH 硬译码（255→255）
        bool ok = bch_255_239_decode_hiho_cw_255(hcand_256.data(), out255.data());
        if (!ok) continue; // 非法候选：跳过或视为大度量

        // 4c 组装 256（重算第 256 位）
        uint8_t* DW = &DW_all[c * BCH_N_TOTAL];
        std::copy(out255.begin(), out255.end(), DW);
        DW[PAR_IDX] = parity256_from255(DW);

        // 4d metric
        float m = chase_metric_fast(hard0, DW, absY.data());
        comps.push_back({m, c});
    }

    if (comps.empty()) {
        // 所有候选都失败：Y2 = -a*Y
        for (int i = 0; i < BCH_N_TOTAL; ++i)
            Y2_256[i] = llr_from_float<LLR>(-A * llr_to_float(Y256[i]));
        return;
    }

    // (5) 排序 + 竞争者
    std::sort(comps.begin(), comps.end(),
              [](const Comp& a, const Comp& b){ return a.metric < b.metric; });
    const int NG = std::min((int)comps.size(), NCOMP);
    const uint8_t* DW0 = &DW_all[ comps[0].idx * BCH_N_TOTAL ];
    const float    M0  = comps[0].metric;

    std::vector<float> comp_delta(NG, 0.f);
    for (int j = 1; j < NG; ++j) comp_delta[j] = (comps[j].metric - M0) * B;

    // beta = sum(|最不可靠|[前 e 个]) - c*M0
    const int e_cap = (E > 0) ? std::min(E, L_eff) : L_eff;
    float beta = 0.f;
    for (int i = 0; i < e_cap; ++i) beta += lrp_abs[i];
    beta -= C * M0;

    // (6) 外信息
    for (int i = 0; i < BCH_N_TOTAL; ++i)
    {
        const uint8_t db = DW0[i];
        int j = 1;
        for (; j < NG; ++j) {
            const uint8_t* DWj = &DW_all[ comps[j].idx * BCH_N_TOTAL ];
            if (DWj[i] != db) break;
        }
        float reliability = 0.f;
        if (j < NG) {
            reliability = comp_delta[j];
        } else {
            reliability = beta + D * absY[i];
            if (reliability < 0.f) reliability = 0.f;
        }
        if (db) reliability = -reliability;

        const float y2 = reliability - A * llr_to_float(Y256[i]);
        Y2_256[i] = llr_from_float<LLR>(y2);
    }
}

// ===== 显式实例化（与你工程常用类型一致）=====
template void chase_decode_256<float >(const float*,  float*,  const Params&);
template void chase_decode_256<int8_t>(const int8_t*, int8_t*, const Params&);
template void chase_decode_256<newcode::qfloat<4>>(const newcode::qfloat<4>*,
                                                   newcode::qfloat<4>*,
                                                   const Params&);
template void chase_decode_256<newcode::qfloat<5>>(const newcode::qfloat<5>*,
                                                   newcode::qfloat<5>*,
                                                   const Params&);

} // namespace newcode
