#include <iostream>
#include <vector>
#include <cassert>
#include "newcode/info_extract.hpp"
#include "newcode/matrix.hpp"
#include "newcode/params.hpp"

using namespace newcode;

static Params make_params_skip_demo(size_t& warmup_rows_out) {
    Params p;
    // 典型 oFEC 配置
    p.BITS_PER_SUBBLOCK_DIM = 16;  // B
    p.NUM_SUBBLOCK_COLS     = 8;   // => N = 128
    // 设置一个可控的 warmup：warmup_rows = (2G + INFO)*B = (0 + 1)*16 = 16
    p.NUM_GUARD_SUBROWS     = 0;   // G
    p.INFO_SUBROWS_PER_CODE = 1;   // INFO
    warmup_rows_out = (2 * p.NUM_GUARD_SUBROWS + p.INFO_SUBROWS_PER_CODE) * p.BITS_PER_SUBBLOCK_DIM;
    return p;
}

int main() {
    size_t warmup_rows = 0;
    const Params p = make_params_skip_demo(warmup_rows);

    const int N = static_cast<int>(p.NUM_SUBBLOCK_COLS * p.BITS_PER_SUBBLOCK_DIM); // 128
    const int K = 239;
    const int TAKE_BITS = K - N; // 111

    // 构造：总行数 = warmup_rows(被跳过) + 2(真实信息行)
    const size_t total_rows = warmup_rows + 2;
    Matrix<float> llr(total_rows, static_cast<size_t>(N));

    // warmup 区（前 warmup_rows 行）：随便放点值（即使判为 1/0 也应被忽略）
    for (size_t r = 0; r < warmup_rows; ++r) {
        for (int c = 0; c < TAKE_BITS; ++c) {
            llr[r][static_cast<size_t>(c)] = (c % 3 == 0) ? -1.0f : +0.5f; // 无所谓
        }
        for (int c = TAKE_BITS; c < N; ++c) llr[r][static_cast<size_t>(c)] = 0.0f;
    }

    // 第 1 个真实信息行：交替 0/1（LLR 的符号控制）
    const size_t r0 = warmup_rows + 0;
    for (int c = 0; c < TAKE_BITS; ++c) {
        llr[r0][static_cast<size_t>(c)] = (c % 2 == 0) ? +2.5f : -1.3f;
    }
    for (int c = TAKE_BITS; c < N; ++c) llr[r0][static_cast<size_t>(c)] = 0.0f;

    // 第 2 个真实信息行：全 1（LLR<0）
    const size_t r1 = warmup_rows + 1;
    for (int c = 0; c < TAKE_BITS; ++c) llr[r1][static_cast<size_t>(c)] = -0.7f;
    for (int c = TAKE_BITS; c < N; ++c) llr[r1][static_cast<size_t>(c)] = 0.0f;

    // 抽取
    std::vector<uint8_t> rx = rx_info_from_bit_llr(llr, p);

    // 期望长度：只包含真实信息行
    assert(rx.size() == 2u * static_cast<size_t>(TAKE_BITS));

    // 校验真实第 1 行：交替 0/1 → 0,1,0,1,...
    for (int c = 0; c < TAKE_BITS; ++c) {
        uint8_t expect = (c % 2 == 0) ? 0u : 1u;
        if (rx[static_cast<size_t>(c)] != expect) {
            std::cerr << "Row0 bit mismatch at " << c << ": got "
                      << int(rx[static_cast<size_t>(c)]) << " expect " << int(expect) << "\n";
            return 1;
        }
    }
    // 校验真实第 2 行：全 1
    for (int c = 0; c < TAKE_BITS; ++c) {
        uint8_t expect = 1u;
        size_t idx = static_cast<size_t>(TAKE_BITS) + static_cast<size_t>(c);
        if (rx[idx] != expect) {
            std::cerr << "Row1 bit mismatch at " << c << ": got "
                      << int(rx[idx]) << " expect " << int(expect) << "\n";
            return 1;
        }
    }

    // 额外小测：如果没有真实行（全是 warmup），应返回空
    {
        Matrix<float> only_warm(warmup_rows, static_cast<size_t>(N));
        auto rx2 = rx_info_from_bit_llr(only_warm, p);
        assert(rx2.empty());
    }

    std::cout << "[OK] rx_info_from_bit_llr with warmup-skip test passed.\n";
    return 0;
}
