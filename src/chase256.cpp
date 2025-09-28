#include "newcode/chase256.hpp"
#include "newcode/bch_255_239.hpp"
#include <vector>
#include <algorithm>
#include <cmath>
#include <limits>
#include <cstdint>
#include <type_traits>
// chase256.cpp 顶部
#include "newcode/chase256.hpp"
#include "newcode/params.hpp"
#include "newcode/qfloat.hpp"   // ← 新增，保证能看到 qfloat 模板
// 其它 include ...

namespace newcode {
namespace detail {

constexpr int BCH_N_TOTAL = 256; // 255 + overall parity
constexpr int BCH_N_CORE  = 255;
constexpr int PAR_IDX     = 255;

// -------- LLR 类型适配：to-float / from-float --------
// to-float：要求 LLR 支持 explicit operator float()（你的 qfloat 已支持）
template<typename LLR>
inline float llr_to_float(LLR x) { return static_cast<float>(x); }

// from-float：算术类型直接 cast；非算术类型要求有 static from_float(float)
template<typename LLR, typename std::enable_if<std::is_arithmetic<LLR>::value, int>::type = 0>
inline LLR llr_from_float(float x) { return static_cast<LLR>(x); }

template<typename LLR, typename std::enable_if<!std::is_arithmetic<LLR>::value, int>::type = 0>
inline LLR llr_from_float(float x) { return LLR::from_float(x); }

// ---- 最不可靠位置（在 0..254 中选择 L 个）----
template<typename LLR>
static void find_least_reliable(const LLR* Y256, int L,
                                std::vector<int>& pos, std::vector<float>& absval)
{
    struct Node { float a; int i; };
    std::vector<Node> v; v.reserve(BCH_N_CORE);
    for (int i = 0; i < BCH_N_CORE; ++i)
        v.push_back({ std::fabs(llr_to_float(Y256[i])), i });

    std::nth_element(v.begin(), v.begin() + std::min(L, (int)v.size()), v.end(),
                     [](const Node& x, const Node& y){ return x.a < y.a; });
    const int take = std::min(L, (int)v.size());
    std::sort(v.begin(), v.begin() + take, [](auto& x, auto& y){ return x.a < y.a; });

    pos.resize(take);
    absval.resize(take);
    for (int k = 0; k < take; ++k) { pos[k] = v[k].i; absval[k] = v[k].a; }
}

// ---- 生成翻转模式 ----
static void gen_test_patterns(int L, int n_test, std::vector<std::vector<bool>>& patt)
{
    const int full = (L <= 20) ? (1 << L) : 0;
    if (n_test <= 0) n_test = 1;

    patt.assign(n_test, std::vector<bool>(L, false));

    if (L <= 15 && n_test <= full) {
        for (int c = 0; c < n_test; ++c)
            for (int j = 0; j < L; ++j)
                patt[c][j] = ((c >> j) & 1) != 0;
        return;
    }

    int c = 0;
    // (0) 不翻
    if (c < n_test) c++;
    // (1) 单翻
    for (int j = 0; j < L && c < n_test; ++j) { patt[c][j] = true; c++; }
    // (2+) 固定锚点 + 单翻
    for (int a = 0; a + 1 < L && c < n_test; ++a)
        for (int j = a + 1; j < L && c < n_test; ++j) {
            patt[c].assign(L, false);
            patt[c][a] = true; patt[c][j] = true; c++;
        }
}

// ---- 计算整体偶校验位（让 256 位异或和为 0）----
inline uint8_t parity256_from255(const uint8_t* cw255) {
    uint8_t acc = 0;
    for (int i = 0; i < BCH_N_CORE; ++i) acc ^= cw255[i];
    return acc;
}

// ---- Chase 度量：sum [hard0^DW] * |Y|（含第 256 位）----
template<typename LLR>
static float chase_metric(const uint8_t* hard0_256, const uint8_t* DW_256, const LLR* Y256)
{
    float m = 0.f;
    for (int i = 0; i < BCH_N_TOTAL; ++i)
        if ((hard0_256[i] ^ DW_256[i]) & 1u)
            m += std::fabs(llr_to_float(Y256[i]));
    return m;
}

} // namespace detail

// ======================== 主函数（流程 2~6） ========================
template<typename LLR>
void chase_decode_256(const LLR* Y256, LLR* Y2_256, const Params& p)
{
    using namespace detail;

    const int L      = std::max(1, p.CHASE_L);
    const int NTEST  = std::max(1, p.CHASE_NTEST);
    const int NCOMP  = std::max(1, std::min(p.CHASE_NCOMP, NTEST));

    // 硬判决基线（含整体奇偶位）
    uint8_t hard0[BCH_N_TOTAL];
    for (int i = 0; i < BCH_N_TOTAL; ++i)
        hard0[i] = (llr_to_float(Y256[i]) >= 0.f) ? 0u : 1u;

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
        float m = chase_metric(hard0, DW, Y256);
        comps.push_back({m, c});
    }

    if (comps.empty()) {
        // 所有候选都失败：Y2 = -a*Y
        for (int i = 0; i < BCH_N_TOTAL; ++i)
            Y2_256[i] = llr_from_float<LLR>(-p.CP_A * llr_to_float(Y256[i]));
        return;
    }

    // (5) 排序 + 竞争者
    std::sort(comps.begin(), comps.end(),
              [](const Comp& a, const Comp& b){ return a.metric < b.metric; });
    const int NG = std::min((int)comps.size(), NCOMP);
    const uint8_t* DW0 = &DW_all[ comps[0].idx * BCH_N_TOTAL ];
    const float    M0  = comps[0].metric;

    std::vector<float> relM(NG, 0.f);
    for (int j = 1; j < NG; ++j) relM[j] = (comps[j].metric - M0) * p.CP_B;

    // beta = sum(|最不可靠|[前 e 个]) - c*M0
    int e_cap = (p.CP_E > 0) ? std::min(p.CP_E, L_eff) : L_eff;
    float beta = 0.f;
    for (int i = 0; i < e_cap; ++i) beta += lrp_abs[i];
    beta -= p.CP_C * M0;

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
        if (j < NG) reliability = relM[j];
        else {
            reliability = beta + p.CP_D * std::fabs(llr_to_float(Y256[i]));
            if (reliability < 0.f) reliability = 0.f;
        }
        if (db) reliability = -reliability;

        const float y2 = reliability - p.CP_A * llr_to_float(Y256[i]);
        Y2_256[i] = llr_from_float<LLR>(y2); // 关键：算术类型直接 cast；qfloat 走 from_float()
    }
}

// ===== 显式实例化（与你工程常用类型一致）=====
template void chase_decode_256<float >(const float*,  float*,  const Params&);
template void chase_decode_256<int8_t>(const int8_t*, int8_t*, const Params&);
// chase256.cpp 末尾（在命名空间 newcode 内）
// 注意：类型要精确匹配错误里显示的类型（你的默认 Store=short）
template void chase_decode_256<newcode::qfloat<4>>(const newcode::qfloat<4>*,
                                                   newcode::qfloat<4>*,
                                                   const Params&);
template void chase_decode_256<newcode::qfloat<5>>(const newcode::qfloat<5>*,
                                                   newcode::qfloat<5>*,
                                                   const Params&);

} // namespace newcode
