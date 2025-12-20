#include <stdexcept>
#include <algorithm>
#include "newcode/tx/ofecencoder/ofec_encoder.hpp"
#include "newcode/common/bch/bch_255_239.hpp"
#include <vector>
#include <cstddef>
#include <cstdint>

namespace ofecencoder {

matrix::Matrix<uint8_t> ofec_encode(const std::vector<uint8_t>& bits, const newcode::Params& p)
{
    // 基本尺寸
    const int B  = static_cast<int>(p.BITS_PER_SUBBLOCK_DIM);                         // 16
    const int N  = static_cast<int>(p.NUM_SUBBLOCK_COLS * p.BITS_PER_SUBBLOCK_DIM);   // 128
    const int G  = static_cast<int>(p.NUM_GUARD_SUBROWS);                             // e.g. 2
    const int NB = static_cast<int>(p.NUM_SUBBLOCK_COLS);                             // 8 (= N/B)

    // BCH(255,239)
    const int K         = static_cast<int>(newcode::Params::BCH_K); // 239
    const int PAR_LEN   = 16;                              // parity bits
    const int TAKE_BITS = K - N;                           // 111

    // 一个输入矩形：32x111 = 3552
    const int ROWS_RECT = 2 * B;                           // 32
    const int RECT_BITS = ROWS_RECT * TAKE_BITS;           // 3552

    // (3) 入口检查 & 形成全局 u(i) = [ zero_prefix , bits ]
    if (bits.size() % static_cast<size_t>(RECT_BITS) != 0)
        throw std::runtime_error("ofec_encode: input 'bits' length must be a multiple of 3552.");

    const int TILE_HEIGHT_BR = static_cast<int>(p.tile_height_rows()) / B;   // e.g. 22
    const int PAD_RECTS      = TILE_HEIGHT_BR / 2;                           // e.g. 11
    const size_t ZERO_PREFIX = static_cast<size_t>(PAD_RECTS) * RECT_BITS;

    std::vector<uint8_t> u;
    u.resize(ZERO_PREFIX + bits.size(), 0);
    std::copy(bits.begin(), bits.end(), u.begin() + ZERO_PREFIX);

    // 将 V(R,C,r,c) 展开为二维：行 = R*B + r，列 = C*B + c
    matrix::Matrix<uint8_t> mat = matrix::Matrix<uint8_t>::zero(p.tile_height_rows(), N);

    // —— 左半历史位读取（按你给的式子，带 -2*(N/B) 项）——
    auto read_hist_bit = [&](long R, int r, int k) -> uint8_t {
        // { (R^1) − 2G − 2*(N/B) + 2*floor(k/B) , floor(k/B) , (k%B)^r , r } for k < N
        const long br = (R ^ 1L) - 2L * G - 2L * NB + 2L * (k / B);
        const long bc = (k / B);
        const long rr_in_blk = (k % B) ^ r;
        const long cc_in_blk = r;

        if (br < 0) return 0; // 启动期/负索引置 0
        const size_t rr = static_cast<size_t>(br) * static_cast<size_t>(B) + static_cast<size_t>(rr_in_blk);
        const size_t cc = static_cast<size_t>(bc) * static_cast<size_t>(B) + static_cast<size_t>(cc_in_blk);
        if (rr >= mat.rows() || cc >= mat.cols()) {
            throw std::out_of_range("Error: Index out of bounds. rr = " + std::to_string(rr) + ", cc = " + std::to_string(cc) + ", matrix size: (" + std::to_string(mat.rows()) + ", " + std::to_string(mat.cols()) + ")");
        }        
        return mat[rr][cc];
    };

    // —— 右半 111 的全局 u(i) 索引（k=0..110 对应 WR(128+k)）——
    // u( floor(R/2)*32*111 + ((R%2)*16 + r)*(16 - floor(k/96)) + floor(k/16)*512 + (k%16) )
    const int K96 = (NB - 2) * B; // 96
    auto u_index = [&](long R, int r, int k) -> size_t {
        const size_t P        = static_cast<size_t>(R / 2);
        const size_t row_in_P = static_cast<size_t>((R % 2) * B + r);         // 0..31
        const size_t shrink   = static_cast<size_t>( (k >= K96) ? 1 : 0 );    // floor(k/96)
        const size_t tile_col = static_cast<size_t>(k / B);                   // 0..6
        const size_t col_in_t = static_cast<size_t>(k % B);
        const size_t stride512= static_cast<size_t>(ROWS_RECT * B);           // 32*16 = 512

        return P * static_cast<size_t>(RECT_BITS)
             + row_in_P * static_cast<size_t>(B - shrink)
             + tile_col * stride512
             + col_in_t;
    };

    // —— 将 WR(128..255) 写回到 V 展开矩阵：与 MATLAB 一致 —— 
    auto write_right_to_mat = [&](long R, int r, int k128, uint8_t bit) {
        // MATLAB:
        //   Ct = floor((k-128)/16);
        //   rt = r;
        //   ct = bitxor(mod(k,16), r);
        //   V(R+1, Ct+1, rt+1, ct+1) = WR(k+1);
        const int k = k128; // 128..255
        const size_t Ct = static_cast<size_t>((k - N) / B);
        const size_t ct = static_cast<size_t>((k % B) ^ r);

        const size_t rr = static_cast<size_t>(R) * static_cast<size_t>(B) + static_cast<size_t>(r); // 行 = R*B+r
        const size_t cc = Ct * static_cast<size_t>(B) + ct;                                         // 列 = C*B+c

        // 确保行够：如果不够就逐行扩
        while (rr >= mat.rows()) mat.add_row();
        mat[rr][cc] = bit;
    };

    // 逐“全局行”（按 r 上增）编码（每行消耗 111 个系统位；系统位来自 u(i)）
    size_t produced_rows = 0;
    size_t global_row    = static_cast<size_t>(p.tile_height_rows()); // 保持你的起点
    const size_t rows_to_make = bits.size() / static_cast<size_t>(TAKE_BITS);

    while (produced_rows < rows_to_make)
    {
        const long R = static_cast<long>(global_row / static_cast<size_t>(B));
        const int  r = static_cast<int>(global_row % static_cast<size_t>(B));

        // (1) 左半 128 历史位
        std::vector<uint8_t> left128; left128.reserve(static_cast<size_t>(N));
        for (int k = 0; k < N; ++k) left128.push_back(read_hist_bit(R, r, k));

        // (2) 右半 111 系统位：严格从 u(i) 抽
        std::vector<uint8_t> right111; right111.reserve(static_cast<size_t>(TAKE_BITS));
        for (int k = 0; k < TAKE_BITS; ++k) {
            const size_t iu = u_index(R, r, k);
        if (iu < ZERO_PREFIX) {
        throw std::runtime_error("ofec_encode: u_index points into the zero prefix (iu < ZERO_PREFIX).");
        }
        if (iu >= u.size()) {
            throw std::runtime_error("ofec_encode: u_index out of range.");
        }
            right111.push_back(u[iu]);
        }

        // (3) textbook BCH(255,239) + 扩展偶校验
        std::vector<uint8_t> msg239; msg239.reserve(static_cast<size_t>(K));
        msg239.insert(msg239.end(), left128.begin(),  left128.end());
        msg239.insert(msg239.end(), right111.begin(), right111.end());

        auto parity16 = bch::bch_255_239_parity(msg239);
        uint8_t overall = 0;
        for (uint8_t b : msg239)     overall ^= (b & 1u);
        for (uint8_t pbit : parity16) overall ^= (pbit & 1u);

        // (4) 写回 V
        // 128..238 : 111 系统位
        for (int kk = 0; kk < TAKE_BITS; ++kk)
            write_right_to_mat(R, r, N + kk, right111[static_cast<size_t>(kk)]);
        // 239..254 : 16 parity
        for (int i = 0; i < PAR_LEN; ++i)
            write_right_to_mat(R, r, N + TAKE_BITS + i, parity16[static_cast<size_t>(i)]);
        // 255 : overall
        write_right_to_mat(R, r, N + TAKE_BITS + PAR_LEN, overall);

        ++global_row;
        ++produced_rows;
    }

    return mat;
}

} // namespace ofecencoder
