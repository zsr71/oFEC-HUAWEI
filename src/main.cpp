#include <iostream>
#include <vector>
#include <complex>
#include <algorithm>  // for std::min

#include "newcode/params.hpp"
#include "newcode/bitgen.hpp"
#include "newcode/ofec_encoder.hpp"
#include "newcode/qam.hpp"
#include "newcode/awgn.hpp"
#include "newcode/qam_llr.hpp"   // << 新增：QAM 解调（LLR）

using namespace newcode;

// 将编码后的矩阵按行展平为比特流（0/1）
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

    // 4) QAM 调制（示例：QPSK/4-QAM => n_bps = 2；16-QAM 可改为 4，依此类推）
    const unsigned n_bps   = 2;      // 每符号比特数（QPSK=2, 16QAM=4, 64QAM=6, ...）
    auto tx_syms = qam_modulate(coded_bits, n_bps);
    std::cout << "[INFO] Modulated symbols: " << tx_syms.size() << " (Es≈1)\n";

    // 5) 通过 AWGN 信道（按 Eb/N0 设置噪声）
    const float    ebn0_dB   = 8.0f; // 这里设置目标 Eb/N0(dB)
    const float    code_rate = 1.0f; // 若考虑编码增益，请改为实际码率 R（例如 0.9、0.5 等）
    const uint32_t awgn_seed = static_cast<uint32_t>(p.BITGEN_SEED + 100); // 随机数种子（可改）
    auto rx_syms = add_awgn(tx_syms, ebn0_dB, n_bps, code_rate, awgn_seed);

    // 简单打印前 3 个符号对比
    std::cout << "[INFO] Example symbols (TX -> RX):\n";
    for (size_t i = 0; i < std::min<size_t>(3, tx_syms.size()); ++i)
    {
        std::cout << "  " << i
                  << ": (" << tx_syms[i].real() << ", " << tx_syms[i].imag() << ")"
                  << " -> (" << rx_syms[i].real() << ", " << rx_syms[i].imag() << ")\n";
    }

    // 6) QAM 解调：输出比特 LLR（>0 表示更像 0）
    //    这里用 Eb/N0 和码率换算 sigma，再做 Max-Log LLR
    auto llr = qam_llr_from_ebn0(rx_syms, n_bps, ebn0_dB, code_rate);
    std::cout << "[INFO] LLR count: " << llr.size() << " (should be tx_syms.size()*n_bps)\n";
    std::cout << "[INFO] First few LLRs: ";
    for (size_t i = 0; i < std::min<size_t>(8, llr.size()); ++i)
        std::cout << llr[i] << (i + 1 < std::min<size_t>(8, llr.size()) ? ", " : "\n");

    std::cout << "[DONE] Pipeline: bits -> oFEC -> QAM -> AWGN -> QAM LLR\n";
    return 0;
}
