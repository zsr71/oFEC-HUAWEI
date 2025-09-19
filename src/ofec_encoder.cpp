#include "newcode/ofec_encoder.hpp"
#include "newcode/bch_255_239.hpp"  // ← 新增：使用(255,239) BCH 实现
#include <vector>
#include <cstddef>

namespace newcode {

Matrix<uint8_t> ofec_encode(const std::vector<uint8_t>& bits, const Params& p)
{
    // 行宽 N = 子块列数 × 子块维度
    const int B = static_cast<int>(p.BITS_PER_SUBBLOCK_DIM);
    const int N = static_cast<int>(p.NUM_SUBBLOCK_COLS * p.BITS_PER_SUBBLOCK_DIM); // 典型为 8*16=128
    const int G = static_cast<int>(p.NUM_GUARD_SUBROWS);

    // BCH(255,239) 参数（与实现保持一致）
    const int K          = 239;     // 信息位数
    const int PAR_LEN    = 16;      // BCH 校验位数
    const int OVERALLLEN = 1;       // 整体偶校验位数
    const int TAKE_BITS  = K - N;   // 每行从 bits 取的新信息比特数（典型 239-128=111）

    // 预放入“启动”所需的行：(2G + INFO_SUBROWS_PER_CODE) 个子块行 × 每子块 B 行
    Matrix<uint8_t> mat = Matrix<uint8_t>::zero((2 * G + p.INFO_SUBROWS_PER_CODE) * B, N);

    size_t bit_pos    = 0;                                               // bits 游标
    size_t global_row = static_cast<size_t>((2 * G + p.INFO_SUBROWS_PER_CODE) * B); // 追加行起始

    while (bit_pos < bits.size())
    {
        // 当前行对应的 (R, r)
        const long R = static_cast<long>(global_row / B); // block row
        const int  r = static_cast<int>(global_row % B);  // bit-row in block

        // 准备一行缓冲（长度 N）
        std::vector<uint8_t> row(N, 0);

        // ---------------- 先抽取 N 个“旧信息位” → BCH 信息的前半部分 ----------------
        std::vector<uint8_t> info239;
        info239.reserve(K); // K = N + (K - N)

        for (int k = 0; k < N; ++k)
        {
            // 公式 (1):
            // { (R ^ 1) – 2G − 2 N/B + 2 ⌊k/B⌋ , ⌊k/B⌋ , (k % B) ^ r , r }
            const long br = (R ^ 1L) - 2 * G - 2 * (N / B) + 2 * (k / B); // block row
            const long bc = (k / B);                                      // block col
            const long bit_row_in_block = (k % B) ^ r;                    // 子块内行
            const long bit_col_in_block = r;                              // 子块内列

            const size_t rr = static_cast<size_t>(br * B + bit_row_in_block); // 全局行
            const size_t cc = static_cast<size_t>(bc * B + bit_col_in_block); // 全局列

            // 不做边界保护，直接访问（按你的要求）
            info239.push_back(mat[rr][cc]);
        }

        // ---------------- 再取 (K - N) = 111 个“新信息位” ----------------
        for (int i = 0; i < TAKE_BITS; ++i)
        {
            uint8_t v = (bit_pos < bits.size()) ? bits[bit_pos++] : 0;
            info239.push_back(v); // 进入 BCH 信息 239
            row[i] = v;           // 行前 (K - N) 位 = 新信息
        }

        // ---------------- 计算 16 位 BCH 校验 + 整体偶校验位 ----------------
        auto parity16 = bch_255_239_parity(info239);

        // 16 位 BCH parity 放在 row 的紧随其后位置
        for (int i = 0; i < PAR_LEN; ++i)
            row[TAKE_BITS + i] = parity16[i];

        // 整体偶校验位：对 239 信息 + 16 parity 取异或，使 256 位异或和为 0
        uint8_t overall = 0;
        for (auto b : info239) overall ^= (b & 1u);
        for (auto pbit : parity16) overall ^= (pbit & 1u);

        row[TAKE_BITS + PAR_LEN] = overall; // 补齐到 N（典型 111 + 16 + 1 = 128）

        // ---------------- 追加到矩阵 ----------------
        mat.add_row();
        mat[global_row] = row;
        ++global_row;
    }

    return mat;
}

} // namespace newcode
