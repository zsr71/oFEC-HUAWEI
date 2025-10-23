#include "newcode/ofec_encoder.hpp"
#include "newcode/bch_255_239.hpp"

#include <vector>
#include <cstddef>
#include <cstdint>

namespace newcode {

Matrix<uint8_t> ofec_encode(const std::vector<uint8_t>& bits, const Params& p)
{
    // 行宽 N = 子块列数 × 子块维度
    const int B = static_cast<int>(p.BITS_PER_SUBBLOCK_DIM);
    const int N = static_cast<int>(p.NUM_SUBBLOCK_COLS * p.BITS_PER_SUBBLOCK_DIM);
    const int G = static_cast<int>(p.NUM_GUARD_SUBROWS);

    // BCH(255,239) 参数（与实现保持一致）
    const int K        = static_cast<int>(Params::BCH_K);       // 239
    const int PAR_LEN  = 16;                                    // 16 parity bits per row
    const int TAKE_BITS = K - N;                                // 每行新增的信息位数 (239-128=111)

    // 预留 (2*G + INFO_SUBROWS_PER_CODE) 行空间，可在后续填充
    Matrix<uint8_t> mat = Matrix<uint8_t>::zero(p.win_height_rows(), N);

    size_t bit_pos    = 0;                                // bits 游标
    size_t global_row = static_cast<size_t>(p.win_height_rows());

    while (bit_pos < bits.size()) {
        const long R = static_cast<long>(global_row / B); // block row
        const int  r = static_cast<int>(global_row % B);  // bit-row in block

        std::vector<uint8_t> row(static_cast<size_t>(N), 0);
        std::vector<uint8_t> info239;
        info239.reserve(static_cast<size_t>(K));

        // (1) 搬运旧信息 (先验) 128 位
        for (int k = 0; k < N; ++k) {
            const long br = (R ^ 1L) - 2 * G - 2 * (N / B) + 2 * (k / B);
            const long bc = (k / B);
            const long bit_row_in_block = (k % B) ^ r;
            const long bit_col_in_block = r;

            const size_t rr = static_cast<size_t>(br * B + bit_row_in_block);
            const size_t cc = static_cast<size_t>(bc * B + bit_col_in_block);

            info239.push_back(mat[rr][cc]);
        }

        // (2) 附加新的 111 位信息
        for (int i = 0; i < TAKE_BITS; ++i) {
            uint8_t v = (bit_pos < bits.size()) ? bits[bit_pos++] : 0;
            info239.push_back(v);
            row[static_cast<size_t>(i)] = v;
        }

        // (3) 计算 16 位 BCH parity
        auto parity16 = bch_255_239_parity(info239);
        for (int i = 0; i < PAR_LEN; ++i)
            row[static_cast<size_t>(TAKE_BITS + i)] = parity16[static_cast<size_t>(i)];

        // (4) 整体偶校验位
        uint8_t overall = 0;
        for (uint8_t b : info239)    overall ^= (b & 1u);
        for (uint8_t pbit : parity16) overall ^= (pbit & 1u);
        row[static_cast<size_t>(TAKE_BITS + PAR_LEN)] = overall;

        mat.add_row();
        mat[global_row] = row;
        ++global_row;
    }

    return mat;
}

} // namespace newcode
