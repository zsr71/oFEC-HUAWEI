#pragma once

namespace chase {

/**
 * Plain Chase(256) float-only 实现。
 * Lin256 是当前输入的行级 LLR，Y2_256 输出 extrinsic。
 */
void chase_decode_256_plain(const float* lin256,
                            const float* /*lch256*/,
                            float* Y2_256,
                            const new_float_only::Params& p)
{
    using namespace new_float_only::detail;

    float beta;  // (20) fallback magnitude scale
    float alpha; // (21) extrinsic scaling (applied by caller)
    pick_cp(p, beta, alpha);
    (void)alpha; // scaling moved to decoder wrapper

    const int L      = std::max(1, p.CHASE_L);
    const int NTEST  = std::max(1, p.CHASE_NTEST);

    // y_k = LLR inputs for correlation metric in (14)–(17)
    float y[BCH_N_TOTAL];
    uint8_t hard_ch[BCH_N_TOTAL];
    for (int i = 0; i < BCH_N_TOTAL; ++i) {
        const float v = lin256[i];
        y[i]       = v;
        hard_ch[i] = (v >= 0.f) ? 0u : 1u;
    }
    std::array<float, BCH_N_TOTAL> abs_y{};
    for (int k = 0; k < BCH_N_TOTAL; ++k) abs_y[k] = std::fabs(y[k]);
    // unreliable set over core (0..254)
    std::vector<int>   lrp_pos; lrp_pos.reserve(L);
    std::vector<float> lrp_abs; lrp_abs.reserve(L);
    find_least_reliable(lin256, L, lrp_pos, lrp_abs);
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
        bool ok = bch::bch_255_239_decode_hiho_cw_255(tmp_in.data(),
                                                 cw255.data(),
                                                 &corrected_errors);

        // build full 256-bit codeword
        auto& CW = CW_all[c];
        std::copy(cw255.begin(), cw255.end(), CW.begin());
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
            omega[j] = beta * sgn;
        }
    }

    dump_chase_csv(trace_copy, y, hard_ch, ML.data(), omega.data());

    // Output EXTRINSIC  Per (21), caller shall form y(next) = y(ch) + α·ω.
    for (int j = 0; j < BCH_N_TOTAL; ++j)
        Y2_256[j] = omega[j];
}

// ======================== 2-arg wrapper (kept for API parity) ========================
void chase_decode_256_plain(const float* Y256, float* Y2_256, const new_float_only::Params& p)
{
    chase_decode_256_plain(Y256, Y256, Y2_256, p);
}

} // namespace chase
