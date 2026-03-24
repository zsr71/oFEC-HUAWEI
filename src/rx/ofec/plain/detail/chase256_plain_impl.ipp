#pragma once

namespace chase {

template<typename LLR>
void chase_decode_256_plain(const LLR* Lin256,
                            const LLR* /*Lch256*/,
                            float* Y2_256,
                            const newcode::Params& p)
{
    // 输入:
    // - Lin256: 当前 256 维输入 LLR。
    // - Lch256: 信道 LLR；本实现未使用，保留接口兼容性。
    // - Y2_256: 输出外信息数组。
    // - p: Chase/Pyndiah 参数。
    // 输出:
    // - 在 Y2_256 中写入 256 维外信息 omega。
    // 用途:
    // - 这是最基础的 Plain Chase-256 实现：
    //   选最不可靠位 -> 枚举翻转模式 -> BCH 硬译码 -> 选 ML 码字 -> 计算每位外信息。
    using namespace newcode::detail;

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
        // 输入同时保留 float 形式 y 和硬判形式 hard_ch。
        const float v = llr_to_float(Lin256[i]);
        y[i]       = v;
        hard_ch[i] = (v >= 0.f) ? 0u : 1u;
    }
    std::array<float, BCH_N_TOTAL> abs_y{};
    for (int k = 0; k < BCH_N_TOTAL; ++k) abs_y[k] = std::fabs(y[k]);
    // unreliable set over core (0..254)
    std::vector<int>   lrp_pos; lrp_pos.reserve(L);
    std::vector<float> lrp_abs; lrp_abs.reserve(L);
    find_least_reliable(Lin256, L, lrp_pos, lrp_abs);
    const int L_eff = (int)lrp_pos.size();

    // test patterns
    std::vector<std::vector<bool>> patt;
    gen_test_patterns(L_eff, NTEST, patt);

    // 生成候选、执行 BCH 硬译码、补 overall parity，并对每个候选计算度量分数。
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
        // 对最不可靠位施加当前测试模式的翻转。
        std::copy(hard_ch, hard_ch + BCH_N_TOTAL, tmp_in.begin());
        for (int j = 0; j < L_eff; ++j)
            if (patt[c][j]) tmp_in[ lrp_pos[j] ] ^= 1u;

        // 在 255 位 BCH 核心码上做硬输入/硬输出译码。
        int corrected_errors = 0;
        bool ok = bch::bch_255_239_decode_hiho_cw_255(tmp_in.data(),
                                                 cw255.data(),
                                                 &corrected_errors);

        // 补成完整 256 位码字。
        auto& CW = CW_all[c];
        std::copy(cw255.begin(), cw255.end(), CW.begin());
        CW[PAR_IDX] = parity256_from255(CW.data());

        // 度量等于与原始硬判不同位置上的 |LLR| 之和；分数越大表示距离越近。
        float dist = 0.f;
        for (int k = 0; k < BCH_N_TOTAL; ++k) {
            const uint8_t diff = (hard_ch[k] ^ CW[k]); // 1=不一致，0=一致
            dist += abs_y[k] * (diff ? 1.f : 0.f);
        }
        float score = -dist;  // 越大越好（等价于最小化 dist）

        comps.push_back({score, c, ok, corrected_errors});
    }

    // 从所有 BCH 成功的候选中选分数最高的 ML 码字；若没有合法候选则回退到输入硬判。
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
        // 对每个 bit j，分别寻找“j 位为 +1/-1”时分数最好的两个竞争码字。
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

            // 用竞争码字差异构造第 j 位的外信息，不包含本位的 channel 项。
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
            // 如果找不到有效竞争者，则交给 fallback 规则处理。
            omega[j] = std::numeric_limits<float>::quiet_NaN(); // 交由后面 L0 回退填充
        }
        if (std::isnan(omega[j])) {
            const float sgn = ML[j] ? -1.f : +1.f;
            // fallback: 只保留 ML 判决符号和固定幅度 beta。
            omega[j] = beta * sgn;
        }
    }

    dump_chase_csv(trace_copy, y, hard_ch, ML.data(), omega.data());

    // 这里只输出纯外信息；上层负责再按 y(next) = y(ch) + alpha * omega 合成下一轮输入。
    for (int j = 0; j < BCH_N_TOTAL; ++j)
        Y2_256[j] = omega[j];
}

// ======================== 2-arg wrapper (kept for API parity) ========================
template<typename LLR>
void chase_decode_256_plain(const LLR* Y256, float* Y2_256, const newcode::Params& p)
{
    // 输入:
    // - Y256: 输入 LLR。
    // - Y2_256: 输出外信息。
    // - p: 参数集合。
    // 输出:
    // - 调用三参版本完成写回。
    // 用途:
    // - 提供与其他 decoder 接口一致的简化包装。
    chase_decode_256_plain<LLR>(Y256, Y256, Y2_256, p);
}

} // namespace chase
