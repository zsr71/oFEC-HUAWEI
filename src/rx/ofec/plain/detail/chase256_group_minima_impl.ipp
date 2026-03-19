#pragma once

namespace chase {

template<typename LLR>
void chase_decode_256_group_minima(const LLR* Lin256,
                                   const LLR* /*Lch256*/,
                                   float* Y2_256,
                                   const newcode::Params& p)
{
    using namespace newcode::detail;

    float beta;
    float alpha;
    pick_cp(p, beta, alpha);
    (void)alpha;

    const int L = std::max(1, p.CHASE_L);
    const int NTEST = std::max(1, p.CHASE_NTEST);

    float y[BCH_N_TOTAL];
    uint8_t hard_ch[BCH_N_TOTAL];
    for (int i = 0; i < BCH_N_TOTAL; ++i) {
        const float v = llr_to_float(Lin256[i]);
        y[i] = v;
        hard_ch[i] = (v >= 0.f) ? 0u : 1u;
    }

    std::array<float, BCH_N_TOTAL> abs_y{};
    for (int k = 0; k < BCH_N_TOTAL; ++k) {
        abs_y[k] = std::fabs(y[k]);
    }

    std::vector<int> lrp_pos;
    std::vector<float> lrp_abs;
    lrp_pos.reserve(L);
    lrp_abs.reserve(L);
    find_least_reliable(Lin256, L, lrp_pos, lrp_abs);
    const int L_eff = static_cast<int>(lrp_pos.size());

    std::vector<std::vector<bool>> patt;
    gen_test_patterns(L_eff, NTEST, patt);

    std::vector<std::vector<uint8_t>> CW_all(
        static_cast<std::size_t>(NTEST),
        std::vector<uint8_t>(BCH_N_TOTAL, 0));
    struct Comp {
        float score;
        int idx;
        bool good;
        int corrected_errors;
    };
    std::vector<Comp> comps;
    comps.reserve(static_cast<std::size_t>(NTEST));

    std::vector<uint8_t> tmp_in(BCH_N_TOTAL), cw255(BCH_N_CORE);
    for (int c = 0; c < NTEST; ++c)
    {
        std::copy(hard_ch, hard_ch + BCH_N_TOTAL, tmp_in.begin());
        for (int j = 0; j < L_eff; ++j) {
            if (patt[static_cast<std::size_t>(c)][j]) {
                tmp_in[static_cast<std::size_t>(lrp_pos[j])] ^= 1u;
            }
        }

        int corrected_errors = 0;
        const bool ok = bch::bch_255_239_decode_hiho_cw_255(tmp_in.data(),
                                                            cw255.data(),
                                                            &corrected_errors);

        auto& CW = CW_all[static_cast<std::size_t>(c)];
        std::copy(cw255.begin(), cw255.end(), CW.begin());
        CW[PAR_IDX] = parity256_from255(CW.data());

        float dist = 0.f;
        for (int k = 0; k < BCH_N_TOTAL; ++k) {
            const uint8_t diff = (hard_ch[k] ^ CW[k]);
            dist += abs_y[k] * (diff ? 1.f : 0.f);
        }
        comps.push_back({-dist, c, ok, corrected_errors});
    }

    const int group_bits = std::min(std::max(0, p.CHASE_GROUP_MINIMA_BITS), L_eff);
    const int group_count = 1 << group_bits;
    std::vector<bool> group_has_best(static_cast<std::size_t>(group_count), false);
    std::vector<Comp> group_best(static_cast<std::size_t>(group_count),
                                 Comp{-std::numeric_limits<float>::infinity(), -1, false, 0});

    for (const auto& cp : comps) {
        if (!cp.good) {
            continue;
        }

        int group_key = 0;
        for (int j = 0; j < group_bits; ++j) {
            if (patt[static_cast<std::size_t>(cp.idx)][j]) {
                group_key |= (1 << j);
            }
        }

        auto& best = group_best[static_cast<std::size_t>(group_key)];
        const bool should_replace =
            !group_has_best[static_cast<std::size_t>(group_key)] ||
            cp.score > best.score ||
            (cp.score == best.score && cp.idx < best.idx);
        if (should_replace) {
            best = cp;
            group_has_best[static_cast<std::size_t>(group_key)] = true;
        }
    }

    std::vector<Comp> kept_comps;
    kept_comps.reserve(static_cast<std::size_t>(group_count));
    for (int group = 0; group < group_count; ++group) {
        if (group_has_best[static_cast<std::size_t>(group)]) {
            kept_comps.push_back(group_best[static_cast<std::size_t>(group)]);
        }
    }
    std::sort(kept_comps.begin(), kept_comps.end(),
              [](const Comp& lhs, const Comp& rhs) {
                  if (lhs.score != rhs.score) {
                      return lhs.score > rhs.score;
                  }
                  return lhs.idx < rhs.idx;
              });

    int ml_idx = -1;
    float ml_S = -std::numeric_limits<float>::infinity();
    for (const auto& cp : kept_comps) {
        if (cp.score > ml_S) {
            ml_S = cp.score;
            ml_idx = cp.idx;
        }
    }

    std::vector<uint8_t> ML(BCH_N_TOTAL, 0);
    if (ml_idx >= 0) {
        ML = CW_all[static_cast<std::size_t>(ml_idx)];
    } else {
        std::fprintf(stderr,
                     "[chase_group_minima] Warning: no valid BCH candidate after grouped good-only pruning; "
                     "fallback to channel hard decisions (extend parity).\n");
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
        int best_idx_plus = -1;
        float best_S_plus = -std::numeric_limits<float>::infinity();
        int best_idx_minus = -1;
        float best_S_minus = -std::numeric_limits<float>::infinity();

        for (const auto& cp : kept_comps) {
            const auto& C = CW_all[static_cast<std::size_t>(cp.idx)];
            if (C[j] == 0u) {
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
            const auto& Cplus = CW_all[static_cast<std::size_t>(best_idx_plus)];
            const auto& Cminus = CW_all[static_cast<std::size_t>(best_idx_minus)];
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
                const float c_plus_l = Cplus[l] ? -1.f : +1.f;
                const float c_minus_l = Cminus[l] ? -1.f : +1.f;
                const int p_l = (c_plus_l != c_minus_l) ? 1 : 0;

                wj += r_l * c_plus_l * static_cast<float>(p_l);
            }
            omega[j] = wj;
        }

        if (std::isnan(omega[j])) {
            const float sgn = ML[j] ? -1.f : +1.f;
            omega[j] = beta * sgn;
        }
    }

    dump_chase_csv(trace_copy, y, hard_ch, ML.data(), omega.data());

    for (int j = 0; j < BCH_N_TOTAL; ++j) {
        Y2_256[j] = omega[j];
    }
}

template<typename LLR>
void chase_decode_256_group_minima(const LLR* Y256,
                                   float* Y2_256,
                                   const newcode::Params& p)
{
    chase_decode_256_group_minima<LLR>(Y256, Y256, Y2_256, p);
}

} // namespace chase
