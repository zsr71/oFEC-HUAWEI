#include "newcode/ofec_decoder.hpp"
#include <stdexcept>
#include <cassert>

namespace newcode {

static inline uint8_t hard_decide(float L) noexcept { return L >= 0.f ? 0u : 1u; }

Matrix<uint8_t> ofec_decode_llr(const Matrix<float>& llr_mat, const Params& p)
{
    // 与编码保持一致的参数
    const int B = static_cast<int>(p.BITS_PER_SUBBLOCK_DIM);
    const int N = static_cast<int>(p.NUM_SUBBLOCK_COLS * p.BITS_PER_SUBBLOCK_DIM); // 典型 128
    const int G = static_cast<int>(p.NUM_GUARD_SUBROWS);

    // BCH(255,239) + overall(1) 的布局假设
    const int K          = 239;     // 信息位数
    const int PAR_LEN    = 16;      // BCH 校验位
    const int OVERALLLEN = 1;       // 整体偶校验
    const int TAKE_BITS  = K - N;   // 每行“新信息位”的数量（典型 111）

    // 关键一致性检查
    if (2 * N != (K + PAR_LEN + OVERALLLEN))
        throw std::invalid_argument("ofec_decode_llr: layout mismatch, require 2N = K+PAR_LEN+OVERALLLEN (i.e., N=128).");
    if (B <= 0 || (N % B) != 0)
        throw std::invalid_argument("ofec_decode_llr: N must be multiple of B.");

    const size_t RROWS = llr_mat.rows();
    const size_t CCOLS = llr_mat.cols();
    if (CCOLS != static_cast<size_t>(N))
        throw std::invalid_argument("ofec_decode_llr: llr_mat cols != N.");

    Matrix<uint8_t> out = Matrix<uint8_t>::zero(RROWS, CCOLS); // 默认 0，仅填每行前 111 列

    // 从最后一行往前解码（按你的要求）
    for (size_t global_row = (RROWS == 0 ? 0 : RROWS - 1); global_row  >(p.INFO_SUBROWS_PER_CODE+2*G)*B-1; )
    {
        const long R = static_cast<long>(global_row / B); // block row
        const int  r = static_cast<int>(global_row % B);  // bit-row in block

        // 组装 239 信息位 + 16 parity + 1 overall 的 LLR（按编码布局）
        std::vector<float> info239_llr; info239_llr.reserve(K);
        // 1) 先取 N 个“旧信息位”的 LLR：按与编码相同的式(1)索引
        for (int k = 0; k < N; ++k)
        {
            const long br = (R ^ 1L) - 2 * G - 2 * (N / B) + 2 * (k / B); // block row
            const long bc = (k / B);                                      // block col
            const long bit_row_in_block = (k % B) ^ r;                    // 子块内行
            const long bit_col_in_block = r;                              // 子块内列

            if (br < 0)
                throw std::runtime_error("ofec_decode_llr: negative back-reference (br<0). Increase INFO_SUBROWS_PER_CODE.");
            const size_t rr = static_cast<size_t>(br * B + bit_row_in_block);
            const size_t cc = static_cast<size_t>(bc * B + bit_col_in_block);

            if (rr >= RROWS || cc >= CCOLS)
                throw std::runtime_error("ofec_decode_llr: back-reference out of range.");
            info239_llr.push_back(llr_mat[rr][cc]);
        }

        // 2) 再取本行前 111 列（新信息位）的 LLR
        for (int i = 0; i < TAKE_BITS; ++i)
            info239_llr.push_back(llr_mat[global_row][i]);

        // 3) 取 parity 与 overall（若将来做真 BCH 可用；此占位解码实际未使用）
        std::vector<float> p16_llr; p16_llr.reserve(PAR_LEN);
        for (int i = 0; i < PAR_LEN; ++i)
            p16_llr.push_back(llr_mat[global_row][TAKE_BITS + i]);
        const float overall_llr = llr_mat[global_row][TAKE_BITS + PAR_LEN];

        // ---- 占位版“BCH 硬解码”：仅做硬判决并截取前 239 位 ----
        std::vector<uint8_t> info239_bits; info239_bits.reserve(K);
        for (int i = 0; i < K; ++i)
            info239_bits.push_back(hard_decide(info239_llr[i]));

        // 取“后 111 位”作为本行的新信息位（与编码一致）
        for (int i = 0; i < TAKE_BITS; ++i)
            out[global_row][i] = info239_bits[N + i];

        // 递减到上一行（size_t 下行写法避免 underflow）
        if (global_row == 0) break;
        --global_row;
    }

    return out;
}

} // namespace newcode
