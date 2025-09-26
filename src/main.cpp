#include <iostream>
#include <vector>
#include <complex>
#include <algorithm>  // for std::min
#include <cstdint>

#include "newcode/params.hpp"
#include "newcode/bitgen.hpp"
#include "newcode/ofec_encoder.hpp"
#include "newcode/qam.hpp"
#include "newcode/awgn.hpp"
#include "newcode/qam_llr.hpp"        // QAM 解调（LLR）
#include "newcode/ofec_llr_matrix.hpp"// LLR 一维->矩阵
#include "newcode/info_extract.hpp"   // 从 bit_llr_mat 提取信息比特（跳过 warmup）
#include "newcode/ber.hpp"            // BER 统计
#include "newcode/ofec_decoder.hpp"

// ★ 新增：定点量化打包（int8_t），以及回到 float 的适配
#include "newcode/llr_qpack.hpp"

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

    // 4) QAM 调制（示例：QPSK/4-QAM => n_bps = 2；16-QAM 可改为 4）
    const unsigned n_bps = 2; // QPSK
    auto tx_syms = qam_modulate(coded_bits, n_bps);
    std::cout << "[INFO] Modulated symbols: " << tx_syms.size() << " (Es≈1)\n";

    // 5) 通过 AWGN 信道（按 Eb/N0 设置噪声）
    const float ebn0_dB = 20.0f; // 目标 Eb/N0(dB)
    // 真实码率：每行 111 个信息位 / 128 个发送位
    const int   N        = static_cast<int>(p.NUM_SUBBLOCK_COLS * p.BITS_PER_SUBBLOCK_DIM); // 典型 128
    const int   K        = 239;
    const int   TAKEBITS = K - N; // 典型 111
    const float code_rate = static_cast<float>(TAKEBITS) / static_cast<float>(N);
    const uint32_t awgn_seed = static_cast<uint32_t>(p.BITGEN_SEED + 100);
    auto rx_syms = add_awgn(tx_syms, ebn0_dB, n_bps, code_rate, awgn_seed);

    // 简单打印前 3 个符号对比
    std::cout << "[INFO] Example symbols (TX -> RX):\n";
    for (size_t i = 0; i < std::min<size_t>(3, tx_syms.size()); ++i) {
        std::cout << "  " << i
                  << ": (" << tx_syms[i].real() << ", " << tx_syms[i].imag() << ")"
                  << " -> (" << rx_syms[i].real() << ", " << rx_syms[i].imag() << ")\n";
    }

    // 6) QAM 解调：输出比特 LLR（>0 表示更像 0）
    auto llr = qam_llr_from_ebn0(rx_syms, n_bps, ebn0_dB, code_rate);
    std::cout << "[INFO] LLR count: " << llr.size() << " (should be tx_syms.size()*n_bps)\n";
    std::cout << "[INFO] First few LLRs: ";
    for (size_t i = 0; i < std::min<size_t>(8, llr.size()); ++i)
        std::cout << llr[i] << (i + 1 < std::min<size_t>(8, llr.size()) ? ", " : "\n");

    // 7) 展成行×列 LLR 矩阵（float）
    Matrix<float> llr_mat = llr_to_matrix_row_major(llr, code_matrix.rows(), code_matrix.cols());

    // 7.1) ★ 量化为定点 LLR（int8_t），以模拟硬件 4/5 位
    Matrix<int8_t> qllr_mat = quantize_llr_to_int8(llr_mat, p);
    std::cout << "[INFO] Quantized LLR: bits=" << p.LLR_BITS
              << " clip=" << p.LLR_CLIP << "  (int8_t in [-Q,+Q])\n";

    // 8) oFEC 解码（定点 LLR 版本，模板会推导出 int8_t）
    Matrix<int8_t> dec_qllr = ofec_decode_llr(qllr_mat, p);

    // 9) 提取信息位（Pre-FEC / Post-FEC）
    //    rx_info_from_bit_llr 需要 float LLR；这里做轻量转换以复用接口。
    auto llr_mat_pre_f  = cast_qllr_to_float(qllr_mat);   // 或 dequantize_llr_to_float(...)
    auto llr_mat_post_f = cast_qllr_to_float(dec_qllr);

    auto rx_info_bits_pre  = rx_info_from_bit_llr(llr_mat_pre_f,  p);
    auto rx_info_bits_post = rx_info_from_bit_llr(llr_mat_post_f, p);

    std::cout << "[INFO] rx_info_bits: " << rx_info_bits_pre.size()
              << " (flattened, warmup skipped)\n";

    // === 10) BER 统计 ===
    auto ber_pre  = compute_and_print_ber(info_bits, rx_info_bits_pre,  "Pre-FEC");
    auto ber_post = compute_and_print_ber(info_bits, rx_info_bits_post, "Post-FEC");

    std::cout << "[DONE] Pipeline: bits -> oFEC -> QAM -> AWGN -> QAM LLR"
                 " -> quant(int8) -> decode(int8) -> info extract & BER\n";
    return 0;
}
