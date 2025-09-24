#include <iostream>
#include <vector>
#include <complex>
#include <algorithm>

#include "newcode/params.hpp"
#include "newcode/bitgen.hpp"
#include "newcode/ofec_encoder.hpp"
#include "newcode/qam.hpp"
#include "newcode/awgn.hpp"
#include "newcode/qam_llr.hpp"

#include "newcode/ofec_llr_matrix.hpp"   // << 新增：LLR 还原

using namespace newcode;

// 将编码后的矩阵按行展平为比特流（0/1）——与你之前一致
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

    // 2) oFEC 编码 → 得到编码矩阵
    auto code_matrix = ofec_encode(info_bits, p);
    std::cout << "[INFO] oFEC matrix: " << code_matrix.rows()
              << " x " << code_matrix.cols() << "\n";

    // 3) 将编码矩阵展平成比特流（逐行）
    auto coded_bits = flatten_row_major(code_matrix);
    std::cout << "[INFO] Coded bits (flattened): " << coded_bits.size() << "\n";

    // 4) QAM 调制
    const unsigned n_bps = 2; // QPSK
    auto tx_syms = qam_modulate(coded_bits, n_bps);
    std::cout << "[INFO] Modulated symbols: " << tx_syms.size() << " (Es≈1)\n";

    // 5) AWGN
    const float ebn0_dB   = 8.0f;
    const float code_rate = 1.0f; // 如需考虑码率，改成实际 R
    const uint32_t awgn_seed = static_cast<uint32_t>(p.BITGEN_SEED + 100);
    auto rx_syms = add_awgn(tx_syms, ebn0_dB, n_bps, code_rate, awgn_seed);

    // 6) QAM 解调：输出比特 LLR（>0 表示更像 0）
    auto llr = qam_llr_from_ebn0(rx_syms, n_bps, ebn0_dB, code_rate);
    std::cout << "[INFO] LLR count: " << llr.size()
              << " (should be tx_syms.size()*n_bps)\n";

    // 7) 还原成与编码矩阵相同形状（行优先）
    Matrix<float> llr_mat = llr_to_matrix_row_major(llr,
                                                    code_matrix.rows(),
                                                    code_matrix.cols());
    std::cout << "[INFO] LLR matrix: " << llr_mat.rows()
              << " x " << llr_mat.cols() << "\n";

    // 8) 简单 sanity check：取前 8 个位置，比较符号和原始比特的一致性
    size_t mism = 0;
    for (size_t i = 0; i < std::min<size_t>(8, llr.size()); ++i) {
        int hard = (llr[i] < 0.f) ? 1 : 0;   // LLR>0 认为是 0，LLR<0 认为是 1
        if (hard != int(coded_bits[i])) ++mism;
    }
    std::cout << "[INFO] First-8 hard decisions mismatches vs. coded bits: "
              << mism << "\n";

    // 9) 打印 LLR 矩阵前 2 行、前 8 列
    std::cout << "[INFO] LLR matrix [0:2, 0:8):\n";
    for (size_t r = 0; r < std::min<size_t>(2, llr_mat.rows()); ++r) {
        std::cout << "  r=" << r << ": ";
        for (size_t c = 0; c < std::min<size_t>(8, llr_mat.cols()); ++c) {
            std::cout << llr_mat[r][c] << (c+1<8? ", ":"\n");
        }
    }

    std::cout << "[DONE] LLR reshape OK.\n";
    return 0;
}
