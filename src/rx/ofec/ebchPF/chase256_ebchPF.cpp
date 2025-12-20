// chase256.cpp — Pyndiah '98 SISO (eqs. (14)–(17), (20), (21))
// -------------------------------------------------------------------------------------
// This file implements the soft-input/soft-output (SISO) Chase component decoder
// exactly following Pyndiah’s derivation:
//   • (14)–(17): soft output Λ_j = y_j + ω_j, where ω_j is computed from the
//     “best ML codeword” and the “best competing codeword with flipped bit j”
//     using correlation (inner-product) metrics equivalent to Euclidean criteria.
//   • (20): when no competing codeword exists for bit j, use a constant-reliability
//     fallback L0 whose magnitude reflects the average reliability; sign is that of
//     the ML decision at bit j.
//   • (21): the decoder must output ONLY EXTRINSIC information ω_j; the caller
//     shall form the next input as y(next) = y(channel) + α·ω, with α being a schedule.
// -------------------------------------------------------------------------------------
#include "newcode/params.hpp"
#include "newcode/rx/ofec/chase/chase256.hpp"
#include "newcode/common/bch/bch_255_239.hpp"
#include "newcode/common/qfloat/qfloat.hpp"

#include <filesystem>
#include <fstream>
#include <string>
#include <cctype>
#include <vector>
#include <algorithm>
#include <cmath>
#include <limits>
#include <cstdint>
#include <type_traits>
#include <array>

namespace chase {
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

// ----- select coefficients from newcode::Params -----
static inline void pick_cp(const newcode::Params& p, float& beta, float& alpha)
{
    // p.beta is reused as the fallback magnitude scale for L0 in (20)
    beta = p.beta;
    alpha = p.ALPHA;
}

// ----- find L least reliable core positions (indices 0..255) by |LLR| -----
template<typename LLR>
static void find_least_reliable(const LLR* LLR_list, int L,
                                std::vector<int>& pos, std::vector<float>& absval)
{
    struct Node { float a; int i; };
    std::vector<Node> v; v.reserve(BCH_N_TOTAL);
    for (int i = 0; i < BCH_N_TOTAL; ++i)
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

    // complete enumeration if small enough
    if (L >= 0 && L < 31) {
        const int full = 1 << L;
        if (n_test <= full) {
            for (int c = 0; c < n_test; ++c)
                for (int j = 0; j < L; ++j)
                    patt[c][j] = ((c >> j) & 1) != 0;
            return;
        }
    }

    // layered growth: 0-flip, then all single-flips, then a few double-flips, ...
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

static void dump_chase_csv(const newcode::Params::DebugTraceConfig& trace,
                           const float* y,
                           const uint8_t* hard_ch,
                           const uint8_t* ML,
                           const float* omega)
{
    if (!trace.enable || !trace.dump_chase_csv) return;
    if (trace.chase_tile_index < 0 || trace.chase_invocation < 0) return;

    std::vector<newcode::Params::DebugTraceConfig::ChaseTraceEntry> entries =
        trace.active_chase_entries;
    if (entries.empty()) {
        if (trace.chase_decoder_col < 0) return;
        newcode::Params::DebugTraceConfig::ChaseTraceEntry fallback;
        fallback.row_index = trace.chase_decoder_row;
        fallback.k = trace.chase_decoder_col;
        fallback.global_row = trace.row;
        fallback.global_col = trace.col;
        fallback.bit_index = -1;
        entries.push_back(fallback);
    }

    const int tile_idx = trace.chase_tile_index;
    const int inv_id   = trace.chase_invocation;

    namespace fs = std::filesystem;
    const fs::path dir = trace.chase_csv_dir.empty()
                           ? fs::path("data/chase_csv")
                           : fs::path(trace.chase_csv_dir);
    std::error_code ec;
    fs::create_directories(dir, ec);

    for (const auto& entry : entries) {
        if (entry.row_index < 0 || entry.k < 0) continue;
        std::string label = entry.label;
        if (label.empty()) {
            if (entry.bit_index >= 0) {
                label = "bit" + std::to_string(entry.bit_index);
            } else {
                label = "row" + std::to_string(entry.global_row) +
                        "_col" + std::to_string(entry.global_col);
            }
        }
        std::string safe_label = label;
        for (char& ch : safe_label) {
            if (!std::isalnum(static_cast<unsigned char>(ch)) && ch != '-') {
                ch = '_';
            }
        }
        const fs::path file = dir / ("tile" + std::to_string(tile_idx)
                                     + "_call" + std::to_string(inv_id)
                                     + "_" + safe_label + ".csv");
        std::ofstream out(file);
        if (!out) continue;
        out << "target_label," << label << '\n';
        out << "target_bit_index," << entry.bit_index << '\n';
        out << "target_global_row," << entry.global_row << '\n';
        out << "target_global_col," << entry.global_col << '\n';
        out << "target_expected_bit," << entry.expected_bit << '\n';
        out << "lin_matrix_row_index," << entry.row_index << '\n';
        out << "lin_matrix_k," << entry.k << '\n';
        const std::vector<int8_t>* expected_bits_row = trace.chase_expected_bits_row;
        const bool has_competing =
            entry.cplus_bits.size() == detail::BCH_N_TOTAL &&
            entry.cminus_bits.size() == detail::BCH_N_TOTAL;
        out << "Index,LLR,Omega,ML,HardCh,ExpectedBit";
        if (has_competing) {
            out << ",CplusBit,CminusBit";
        }
        out << '\n';
        for (int k = 0; k < BCH_N_TOTAL; ++k) {
            int expected_bit = -1;
            if (expected_bits_row && static_cast<size_t>(k) < expected_bits_row->size()) {
                expected_bit = (*expected_bits_row)[static_cast<size_t>(k)];
            }
            out << k << ',' << y[k] << ',' << omega[k] << ','
                << int(ML[k]) << ',' << int(hard_ch[k]) << ','
                << expected_bit;
            if (has_competing) {
                out << ',' << int(entry.cplus_bits[static_cast<size_t>(k)])
                    << ',' << int(entry.cminus_bits[static_cast<size_t>(k)]);
            }
            out << '\n';
        }
    }
}

} // namespace detail

template<typename LLR>
void chase_decode_256_ebchPF(const LLR* Lin256,
                            const LLR* /*Lch256*/,
                            float* Y2_256,
                            const newcode::Params& p)
{
    using namespace detail;

    float beta;  // (20) fallback magnitude scale
    float alpha; // (21) extrinsic scaling (applied by caller)
    pick_cp(p, beta, alpha);
    (void)alpha; // scaling moved to decoder wrapper

    const int L      = std::max(1, p.CHASE_L);
    const int NTEST  = std::max(1, p.CHASE_NTEST);

    // y_k = LLR inputs for correlation metric in (14))
    float y[BCH_N_TOTAL];
    uint8_t hard_ch[BCH_N_TOTAL];
    for (int i = 0; i < BCH_N_TOTAL; ++i) {
        const float v = llr_to_float(Lin256[i]);
        y[i]       = v;
        hard_ch[i] = (v >= 0.f) ? 0u : 1u;
    }
    std::array<float, BCH_N_TOTAL> abs_y{};
    for (int k = 0; k < BCH_N_TOTAL; ++k) abs_y[k] = std::fabs(y[k]);
    // unreliable set over core (0..255)
    std::vector<int>   lrp_pos; lrp_pos.reserve(L);
    std::vector<float> lrp_abs; lrp_abs.reserve(L);
    find_least_reliable(Lin256, L, lrp_pos, lrp_abs);
    const int L_eff = (int)lrp_pos.size();

    // test patterns
    std::vector<std::vector<bool>> patt;
    gen_test_patterns(L_eff, NTEST, patt);

    // generate candidates, BCH hard-decode, extend to 256, and compute S(c)
    std::vector<std::vector<uint8_t>> CW_all(NTEST, std::vector<uint8_t>(BCH_N_TOTAL, 0));
    struct Comp {
        float score;
        int idx;
        bool good;
        int corrected_errors;
    };
    std::vector<Comp> comps; comps.reserve(NTEST);

    std::vector<uint8_t> tmp_in(BCH_N_TOTAL), cw255(BCH_N_CORE);

    for (int c = 0; c < NTEST; ++c)
    {
        // apply flips on unreliable set
        std::copy(hard_ch, hard_ch + BCH_N_TOTAL, tmp_in.begin());
        for (int j = 0; j < L_eff; ++j)
            if (patt[c][j]) tmp_in[ lrp_pos[j] ] ^= 1u;

        // BCH decode over 255 (hard-input, hard-output)
        int corrected_errors = 0;
        bool ok =bch::bch_255_239_decode_hiho_cw_255(tmp_in.data(),
                                                 cw255.data(),
                                                 &corrected_errors);

        // build full 256-bit codeword
        auto& CW = CW_all[c];
        std::copy(cw255.begin(), cw255.end(), CW.begin());
        uint8_t parity = parity256_from255(CW.data());

        if (parity != tmp_in[PAR_IDX] &&(corrected_errors==2) ) 
        { ok = false; corrected_errors = -1; }

        CW[PAR_IDX] = parity256_from255(CW.data());

        float dist = 0.f;
        for (int k = 0; k < BCH_N_TOTAL; ++k) {
            const uint8_t diff = (hard_ch[k] ^ CW[k]); // 1=不一致，0=一致
            dist += abs_y[k] * (diff ? 1.f : 0.f);
        }
        float score = -dist;  // 越大越好（等价于最小化 dist）

        comps.push_back({score, c, ok, corrected_errors});
    }

    // pick ML among valid decodes; if none valid, fall back to channel hard word
    int ml_idx = -1;
    float ml_S = -std::numeric_limits<float>::infinity();
    for (auto &cp : comps) {
        if (!cp.good) continue;
        if (cp.score > ml_S) { ml_S = cp.score; ml_idx = cp.idx; }
    }

    std::vector<uint8_t> ML(BCH_N_TOTAL, 0);
    if (ml_idx >= 0) {
        ML = CW_all[ml_idx];
    } else {
        std::fprintf(stderr,
                 "[chase256] Warning: no valid BCH candidate; "
                 "fallback to channel hard decisions (extend parity).\n");
        // no valid codeword: take channel hard decisions (extend parity)
        std::copy(hard_ch, hard_ch + BCH_N_CORE, ML.begin());
        ML[PAR_IDX] = parity256_from255(ML.data());
        
        ml_S = 0.f;
        for (int k = 0; k < BCH_N_TOTAL; ++k)
            ml_S += y[k] * (ML[k] ? -1.f : +1.f);
    }

    auto trace_copy = p.debug_trace;
    if (!trace_copy.active_chase_entries.empty()) {
        for (auto& entry : trace_copy.active_chase_entries) {
            entry.cplus_bits.clear();
            entry.cminus_bits.clear();
        }
    }
    std::vector<float> omega(BCH_N_TOTAL, std::numeric_limits<float>::quiet_NaN());
    for (int j = 0; j < BCH_N_TOTAL; ++j) {
        // ----- find best competing codeword for bit j -----
        int   best_idx_plus  = -1;                    
        float best_S_plus    = -std::numeric_limits<float>::infinity();
        int   best_idx_minus = -1;                     
        float best_S_minus   = -std::numeric_limits<float>::infinity();

        for (const auto &cp : comps) {
            if (!cp.good) continue; // 只考虑 BCH 成功的候选

            const auto &C = CW_all[cp.idx];
            if (C[j] == 0) { // j 位为 0 -> BPSK(+1)
                if (cp.score > best_S_plus) {
                    best_S_plus = cp.score;
                    best_idx_plus = cp.idx;
                }
            } else { // j 位为 1 -> BPSK(-1)
                if (cp.score > best_S_minus) {
                    best_S_minus = cp.score;
                    best_idx_minus = cp.idx;
                }
            }
        }

        if (best_idx_plus >= 0 && best_idx_minus >= 0) {
            const auto &Cplus  = CW_all[best_idx_plus];  // c^{+1(j)}
            const auto &Cminus = CW_all[best_idx_minus]; // c^{-1(j)}
            if (!trace_copy.active_chase_entries.empty()) {
                for (auto& entry : trace_copy.active_chase_entries) {
                    if (entry.k == j) {
                        entry.cplus_bits.assign(Cplus.begin(), Cplus.end());
                        entry.cminus_bits.assign(Cminus.begin(), Cminus.end());
                    }
                }
            }

            float wj = 0.f;
            for (int l = 0; l < BCH_N_TOTAL; ++l) {
                if (l == j) continue;

                const float r_l = y[l];
                const float c_plus_l  = Cplus [l] ? -1.f : +1.f;
                const float c_minus_l = Cminus[l] ? -1.f : +1.f;
                const int   p_l = (c_plus_l != c_minus_l) ? 1 : 0;

                wj += r_l * c_plus_l * static_cast<float>(p_l);
            }

            omega[j] = wj; // 只输出外信息（不含 r_j）
        } else {
            omega[j] = std::numeric_limits<float>::quiet_NaN(); // 交由后面 L0 回退填充
        }
        if (std::isnan(omega[j])) {
            const float sgn = ML[j] ? -1.f : +1.f;
            //omega[j] = beta * sgn-llr_to_float(Lin256[j]);
            omega[j] = beta * sgn;
        }
    }

    dump_chase_csv(trace_copy, y, hard_ch, ML.data(), omega.data());

    // Output EXTRINSIC  Per (21), caller shall form y(next) = y(ch) + α·ω.
    for (int j = 0; j < BCH_N_TOTAL; ++j)
        Y2_256[j] = omega[j];
}

// ======================== 2-arg wrapper (kept for API parity) ========================
template<typename LLR>
void chase_decode_256_ebchPF(const LLR* Y256, float* Y2_256, const newcode::Params& p)
{
    chase_decode_256_ebchPF<LLR>(Y256, Y256, Y2_256, p);
}

// ======================== explicit instantiations ========================
template void chase_decode_256_ebchPF<float >(const float*,  const float*,  float*,  const newcode::Params&);
template void chase_decode_256_ebchPF<int8_t>(const int8_t*, const int8_t*, float*, const newcode::Params&);
template void chase_decode_256_ebchPF<float >(const float*,  float*,  const newcode::Params&);
template void chase_decode_256_ebchPF<int8_t>(const int8_t*, float*, const newcode::Params&);

#define INSTANTIATE_CHASE256_ebchPF_QFLOAT(N) \
template void chase_decode_256_ebchPF<qfloat::qfloat<N>>( \
    const qfloat::qfloat<N>*, const qfloat::qfloat<N>*, float*, const newcode::Params&); \
template void chase_decode_256_ebchPF<qfloat::qfloat<N>>( \
    const qfloat::qfloat<N>*, float*, const newcode::Params&);

INSTANTIATE_CHASE256_ebchPF_QFLOAT(2)
INSTANTIATE_CHASE256_ebchPF_QFLOAT(3)
INSTANTIATE_CHASE256_ebchPF_QFLOAT(4)
INSTANTIATE_CHASE256_ebchPF_QFLOAT(5)
INSTANTIATE_CHASE256_ebchPF_QFLOAT(6)
INSTANTIATE_CHASE256_ebchPF_QFLOAT(7)
INSTANTIATE_CHASE256_ebchPF_QFLOAT(8)
INSTANTIATE_CHASE256_ebchPF_QFLOAT(9)
INSTANTIATE_CHASE256_ebchPF_QFLOAT(10)
INSTANTIATE_CHASE256_ebchPF_QFLOAT(11)
INSTANTIATE_CHASE256_ebchPF_QFLOAT(12)
INSTANTIATE_CHASE256_ebchPF_QFLOAT(13)
INSTANTIATE_CHASE256_ebchPF_QFLOAT(14)
INSTANTIATE_CHASE256_ebchPF_QFLOAT(15)

#undef INSTANTIATE_CHASE256_ebchPF_QFLOAT

} // namespace chase
