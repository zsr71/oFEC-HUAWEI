#include <iostream>
#include <algorithm>
#include "newcode/params.hpp"
#include "newcode/bitgen.hpp"
#include "newcode/ofec_encoder.hpp"
#include "newcode/qam.hpp"
#include "newcode/awgn.hpp"
#include "newcode/qam_llr.hpp"
#include "newcode/ofec_llr_matrix.hpp"
#include "newcode/ofec_decoder.hpp"

using namespace newcode;

static std::vector<uint8_t> flatten_row_major(const Matrix<uint8_t>& M)
{
    std::vector<uint8_t> out;
    out.reserve(M.rows() * M.cols());
    for (size_t r = 0; r < M.rows(); ++r)
        for (size_t c = 0; c < M.cols(); ++c)
            out.push_back(M[r][c] & 1u);
    return out;
}

int main()
{
    Params p;

    // 1) 随机信息比特
    auto info_bits = generate_bits(p);
    std::cout << "[INFO] Generated bits: " << info_bits.size() << "\n";

    // 2) oFEC 编码
    auto code_matrix = ofec_encode(info_bits, p);
    std::cout << "[INFO] oFEC matrix: " << code_matrix.rows()
              << " x " << code_matrix.cols() << "\n";

    // 3) 展平 → QAM 调制
    auto coded_bits = flatten_row_major(code_matrix);
    const unsigned n_bps = 2; // QPSK
    auto tx_syms = qam_modulate(coded_bits, n_bps);

    // 4) AWGN
    const float ebn0_dB = 8.0f;
    const float code_rate = 1.0f; // 仅作为 LLR 缩放；如需考虑码率请设置实际 R
    const uint32_t awgn_seed = static_cast<uint32_t>(p.BITGEN_SEED + 100);
    auto rx_syms = add_awgn(tx_syms, ebn0_dB, n_bps, code_rate, awgn_seed);

    // 5) 软解调 LLR → 还原为矩阵
    auto llr = qam_llr_from_ebn0(rx_syms, n_bps, ebn0_dB, code_rate);
    Matrix<float> llr_mat = llr_to_matrix_row_major(llr,
                                                    code_matrix.rows(),
                                                    code_matrix.cols());
    std::cout << "[INFO] LLR matrix: " << llr_mat.rows()
              << " x " << llr_mat.cols() << "\n";

    // 6) oFEC 简易“解码”（占位版 BCH：硬判决 + 取后 111 位）
    auto bit_llr_mat = ofec_decode_llr(llr_mat, p);
    std::cout << "[INFO] bit_llr_mat: " << bit_llr_mat.rows()
              << " x " << bit_llr_mat.cols() << "\n";

    // 打印末尾 2 行的前 16 个新信息位，目测检查
    const size_t R = bit_llr_mat.rows();
    const size_t show_rows = std::min<size_t>(2, R);
    for (size_t i = 0; i < show_rows; ++i)
    {
        size_t gr = R - 1 - i;
        std::cout << "[ROW " << gr << "] first 16 new-info bits: ";
        for (int c = 0; c < 16 && c < (int)bit_llr_mat.cols(); ++c)
            std::cout << (int)bit_llr_mat[gr][c] << (c+1<16? ' ' : '\n');
    }

    std::cout << "[DONE] bits -> oFEC -> QAM -> AWGN -> LLR -> (placeholder) decode\n";
    return 0;
}
