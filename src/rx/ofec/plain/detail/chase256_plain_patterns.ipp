#pragma once

namespace newcode {
namespace detail {

// ----- generate Chase test patterns over the L unreliable positions -----
static void gen_test_patterns(int L, int n_test, std::vector<std::vector<bool>>& patt)
{
    // 输入:
    // - L: 最不可靠位数量。
    // - n_test: 需要生成的测试模式数。
    // - patt: 输出容器引用。
    // 输出:
    // - patt[c][j] 表示第 c 个模式是否翻转第 j 个最不可靠位。
    // 用途:
    // - 为 Chase 生成 bit-flipping 模式，后续每个模式都会送入 BCH 硬译码形成候选码字。
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

} // namespace detail
} // namespace newcode
