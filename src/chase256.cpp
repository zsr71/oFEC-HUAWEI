// chase256.cpp — Pyndiah '98 SISO (eqs. (14)–(17), (20), (21))
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
#include <array> // <<< 新增

namespace newcode {
namespace detail {

constexpr int BCH_N_TOTAL = 256; // 255 core + 1 overall parity (extended code)
constexpr int BCH_N_CORE  = 255;
constexpr int PAR_IDX     = 255;

// ----- LLR adapters -----
template<typename LLR>
inline float llr_to_float(LLR x) { return static_cast<float>(x); }

template<typename LLR, typename std::enable_if<std::is_floating_point<LLR>::value, int>::type = 0>
inline LLR llr_from_float(float x) { return static_cast<LLR>(x); }

template<typename LLR, typename std::enable_if<std::is_integral<LLR>::value && std::is_signed<LLR>::value, int>::type = 0>
inline LLR llr_from_float(float x)
{
    const float lo = (float)std::numeric_limits<LLR>::min();
    const float hi = (float)std::numeric_limits<LLR>::max();
    if (x < lo) x = lo;
    if (x > hi) x = hi;
    return (LLR)std::lrintf(x);
}

template<typename LLR, typename std::enable_if<!std::is_arithmetic<LLR>::value, int>::type = 0>
inline LLR llr_from_float(float x) { return LLR::from_float(x); }

// ----- select coefficients from Params -----
static inline void pick_cp(const Params& p, float& beta, float& alpha)
{
    beta = p.beta;   // fallback magnitude
    alpha = p.ALPHA; // scaling of extrinsic
}

// ----- find L least reliable core positions (indices 0..254) by |LLR| -----
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

// ----- generate Chase test patterns over the L unreliable positions -----
static void gen_test_patterns(int L, int n_test, std::vector<std::vector<bool>>& patt)
{
    if (n_test <= 0) n_test = 1;
    patt.assign(n_test, std::vector<bool>(L, false));

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
    if (c < n_test) c++; // all-zero pattern

    for (int layer = 1; c < n_test && layer <= L; ++layer)
        for (int j = layer - 1; j < L && c < n_test; ++j) {
            std::fill(patt[c].begin(), patt[c].end(), false);
            for (int k = 0; k < layer - 1; ++k) patt[c][k] = true;
            patt[c][j] = true;
            c++;
        }
}

// ----- overall (even) parity from 255-bit core (extend to 256 bits) -----
inline uint8_t parity256_from255(const uint8_t* cw255) {
    uint8_t acc = 0;
    for (int i = 0; i < BCH_N_CORE; ++i) acc ^= cw255[i];
    return acc;
}

// ===== 新增：0/1 -> ±1；构建/去重 codeset（仅 BCH 成功） =====
static inline float bit01_to_pm1(uint8_t b) { return b ? -1.f : +1.f; }

static void build_unique_codeset_pm1(
    const std::vector<std::vector<uint8_t>>& CW_all,   // NTEST × 256 (0/1)
    const std::vector<bool>& ok_mask,                  // NTEST
    std::vector<std::array<float,256>>& codeset_pm1    // 输出：M × 256 (±1)
){
    std::vector<std::array<float,256>> tmp;
    tmp.reserve(CW_all.size());
    for (size_t i = 0; i < CW_all.size(); ++i) {
        if (!ok_mask[i]) continue;
        std::array<float,256> v{};
        for (int j = 0; j < 256; ++j) v[(size_t)j] = bit01_to_pm1(CW_all[i][(size_t)j]);
        tmp.push_back(v);
    }
    if (tmp.empty()) { codeset_pm1.clear(); return; }

    auto less_arr = [](const std::array<float,256>& a, const std::array<float,256>& b){
        for (int i = 0; i < 256; ++i)
            if (a[(size_t)i] != b[(size_t)i]) return a[(size_t)i] < b[(size_t)i];
        return false;
    };
    std::sort(tmp.begin(), tmp.end(), less_arr);
    tmp.erase(std::unique(tmp.begin(), tmp.end(),
                          [](const auto& x, const auto& y){
                              for (int i = 0; i < 256; ++i)
                                  if (x[(size_t)i] != y[(size_t)i]) return false;
                              return true;
                          }), tmp.end());
    codeset_pm1.swap(tmp);
}

// ===== 新增：与 MATLAB 一致的欧氏距离：|softin - code|^2 =====
static inline float eucldist_sq_256(const float* softin, const std::array<float,256>& code_pm1)
{
    float acc = 0.f;
    for (int i = 0; i < 256; ++i) {
        const float d = softin[(size_t)i] - code_pm1[(size_t)i];
        acc += d * d;
    }
    return acc;
}

} // namespace detail

template<typename LLR>
void chase_decode_256(const LLR* Lin256,
                      const LLR* /*Lch256*/,
                      LLR* Y2_256,
                      const Params& p)
{
    using namespace detail;

    float beta, alpha;
    pick_cp(p, beta, alpha);

    const int L     = std::max(1, p.CHASE_L);
    const int NTEST = std::max(1, p.CHASE_NTEST);

    // 当前行输入软值（本轮的“总 LLR”）与硬判
    float y[BCH_N_TOTAL];
    uint8_t hard_ch[BCH_N_TOTAL];
    for (int i = 0; i < BCH_N_TOTAL; ++i) {
        const float v = llr_to_float(Lin256[i]);
        y[i]       = v;
        hard_ch[i] = (v >= 0.f) ? 0u : 1u;
    }

    // LRP 仅在 0..254 上（不包含整体位）
    std::vector<int>   lrp_pos; lrp_pos.reserve(L);
    std::vector<float> lrp_abs; lrp_abs.reserve(L);
    find_least_reliable(Lin256, L, lrp_pos, lrp_abs);
    const int L_eff = (int)lrp_pos.size();

    // 生成测试翻转模式
    std::vector<std::vector<bool>> patt;
    gen_test_patterns(L_eff, NTEST, patt);

    // 生成候选：对 255 位做 BCH 硬判译码，补整体位 -> 256 位候选（0/1）
    std::vector<std::vector<uint8_t>> CW_all(NTEST, std::vector<uint8_t>(BCH_N_TOTAL, 0));
    std::vector<bool> ok_mask(NTEST, false);
    std::vector<uint8_t> tmp_in(BCH_N_TOTAL), cw255(BCH_N_CORE);

    for (int c = 0; c < NTEST; ++c)
    {
        // apply flips on unreliable set (only 0..254)
        std::copy(hard_ch, hard_ch + BCH_N_TOTAL, tmp_in.begin());
        for (int j = 0; j < L_eff; ++j)
            if (patt[c][j]) tmp_in[ lrp_pos[j] ] ^= 1u;

        // BCH over first 255 bits
        bool ok = bch_255_239_decode_hiho_cw_255(tmp_in.data(), cw255.data());

        // build full 256-bit candidate
        auto& CW = CW_all[c];
        std::copy(cw255.begin(), cw255.end(), CW.begin());
        CW[PAR_IDX] = parity256_from255(CW.data());

        ok_mask[(size_t)c] = ok; // 标记合法候选
    }

    // 构建 ±1 候选集（仅保留合法，去重）
    std::vector<std::array<float,256>> codeset; // M × 256, ±1
    build_unique_codeset_pm1(CW_all, ok_mask, codeset);

    // softin = y（本轮输入），做行级均值绝对值归一化
    float softin[256];
    for (int i = 0; i < 256; ++i) softin[i] = y[i];
    {
        float m = 0.f;
        for (int i = 0; i < 256; ++i) m += std::fabs(softin[i]);
        m /= 256.f;
        if (m > 0.f) for (int i = 0; i < 256; ++i) softin[i] /= m;
    }

    // 若没有任何合法候选：按 MATLAB chasedec fallback，取“0翻转”硬判并补整体位
    if (codeset.empty()) {
        std::array<float,256> v{};
        for (int j = 0; j < 255; ++j) v[(size_t)j] = (hard_ch[j] ? -1.f : +1.f);
        v[255] = (parity256_from255(hard_ch) ? -1.f : +1.f);
        codeset.push_back(v);
    }

    const int M = (int)codeset.size();

    // 计算所有候选的欧氏距离 |softin - code|^2
    std::vector<float> dists(M);
    for (int m_i = 0; m_i < M; ++m_i)
        dists[(size_t)m_i] = eucldist_sq_256(softin, codeset[(size_t)m_i]);

    // 选最小距离的判决码字 d
    int d_idx = 0;
    for (int m_i = 1; m_i < M; ++m_i)
        if (dists[(size_t)m_i] < dists[(size_t)d_idx]) d_idx = m_i;

    const auto& d = codeset[(size_t)d_idx];
    const float distD = dists[(size_t)d_idx];

    // 计算外信息：tpcrel
    std::vector<float> omega(256, 0.f);

    if (M == 1) {
        // 只有一行：Eq.(19) 回退 β·d
        for (int j = 0; j < 256; ++j)
            omega[(size_t)j] = beta * d[(size_t)j];
    } else {
        for (int j = 0; j < 256; ++j) {
            // 在“j位与 d 不同”的候选里找最小距离
            float bestComp = std::numeric_limits<float>::infinity();
            for (int m_i = 0; m_i < M; ++m_i) {
                if (m_i == d_idx) continue;
                if (codeset[(size_t)m_i][(size_t)j] != d[(size_t)j])
                    if (dists[(size_t)m_i] < bestComp) bestComp = dists[(size_t)m_i];
            }

            if (std::isfinite(bestComp)) {
                const float delta = bestComp - distD;          // |R-C|^2 - |R-D|^2
                omega[(size_t)j] = (delta * 0.25f) * d[(size_t)j] - softin[(size_t)j]; // Eq.(18) + Eq.(16)
            } else {
                omega[(size_t)j] = beta * d[(size_t)j];        // 无竞争者：Eq.(19)
            }
        }
    }

    // 外信息归一化：排除 |w|==β 的点，使 mean(|w|)=1
    {
        float acc = 0.f; int cnt = 0;
        for (int j = 0; j < 256; ++j) {
            if (std::fabs(omega[(size_t)j]) != beta) { acc += std::fabs(omega[(size_t)j]); cnt++; }
        }
        if (cnt > 0) {
            const float g = acc / cnt;
            if (g > 0.f) for (int j = 0; j < 256; ++j) omega[(size_t)j] /= g;
        }
    }

    // 输出 α·ω（上层用“加等于”累加到总 LLR：L_next = L_cur + α·ω）
    for (int j = 0; j < 256; ++j)
        Y2_256[(size_t)j] = llr_from_float<LLR>(alpha * omega[(size_t)j]);
}

// ======================== 2-arg wrapper (kept for API parity) ========================
template<typename LLR>
void chase_decode_256(const LLR* Y256, LLR* Y2_256, const Params& p)
{
    chase_decode_256<LLR>(Y256, Y256, Y2_256, p);
}

// ======================== explicit instantiations ========================
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

template void chase_decode_256<float >(const float*,  float*,  const Params&);
template void chase_decode_256<int8_t>(const int8_t*, int8_t*, const Params&);
template void chase_decode_256<newcode::qfloat<4>>(const newcode::qfloat<4>*,
                                                   newcode::qfloat<4>*,
                                                   const Params&);
template void chase_decode_256<newcode::qfloat<5>>(const newcode::qfloat<5>*,
                                                   newcode::qfloat<5>*,
                                                   const Params&);

} // namespace newcode
