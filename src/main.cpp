#include <iostream>
#include <vector>
#include <complex>
#include <algorithm>
#include <cstdint>

#include "newcode/params.hpp"
#include "newcode/bitgen.hpp"
#include "newcode/ofec_encoder.hpp"
#include "newcode/qam.hpp"
#include "newcode/awgn.hpp"
#include "newcode/qam_llr.hpp"
#include "newcode/ofec_llr_matrix.hpp"
#include "newcode/info_extract.hpp"
#include "newcode/ber.hpp"
#include "newcode/ofec_decoder.hpp"

// A) qfloat 路径（推荐）：像浮点一样用 4/5 位“类浮点”
#include "newcode/qfloat.hpp"
// B) 备选 int8_t 路径（回退用）
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

// —— 用 qfloat<NBITS> 跑完整解码 + BER —— //
template<int NBITS>
static void run_with_qfloat(const Matrix<float>& llr_mat_f,
                            const Params& p,
                            const std::vector<uint8_t>& info_bits)
{
    // 量化为 qfloat<NBITS>
    auto q_mat = quantize_matrix_to_qfloat<NBITS>(llr_mat_f, /*clip=*/qfloat<NBITS>::DEFAULT_CLIP);
    std::cout << "[INFO] Quantized qfloat<" << NBITS << ">: clip=" << qfloat<NBITS>::DEFAULT_CLIP
              << "  code range [-" << qfloat<NBITS>::Q() << ", +" << qfloat<NBITS>::Q() << "]\n";

    // oFEC 解码（模板推导 qfloat<NBITS>）
    auto q_dec = ofec_decode_llr(q_mat, p);

    // 为复用 info_extract/BER（需要 float LLR），做轻量类型转换（阈值仍为 0）
    auto pre_f  = cast_matrix_from_qfloat(q_mat);
    auto post_f = cast_matrix_from_qfloat(q_dec);

    auto rx_info_bits_pre  = rx_info_from_bit_llr(pre_f,  p);
    auto rx_info_bits_post = rx_info_from_bit_llr(post_f, p);

    std::cout << "[INFO] rx_info_bits: " << rx_info_bits_pre.size()
              << " (flattened, warmup skipped)\n";

    compute_and_print_ber(info_bits, rx_info_bits_pre,  "Pre-FEC",  p);
    compute_and_print_ber(info_bits, rx_info_bits_post, "Post-FEC", p);

    std::cout << "[DONE] Pipeline(qfloat<" << NBITS
              << ">) bits -> oFEC -> QAM -> AWGN -> QAM LLR -> quant(qfloat) -> decode -> info extract & BER\n";
}

// —— 用 int8_t（线性量化）跑完整解码 + BER（通用回退） —— //
static void run_with_int8(const Matrix<float>& llr_mat_f,
                          const Params& p,
                          const std::vector<uint8_t>& info_bits)
{
    // 量化为 int8（使用 p.LLR_BITS 和 p.LLR_CLIP，Q = 2^(bits-1)-1）
    auto q8_mat = quantize_llr_to_int8(llr_mat_f, p);
    std::cout << "[INFO] Quantized int8: bits=" << p.LLR_BITS << " clip=" << p.LLR_CLIP << "\n";

    // int8_t 版本解码
    auto q8_dec = ofec_decode_llr(q8_mat, p);
    // 转回 float 以复用 info_extract/BER
    auto pre_f  = cast_qllr_to_float(q8_mat);
    auto post_f = cast_qllr_to_float(q8_dec);

    auto rx_info_bits_pre  = rx_info_from_bit_llr(pre_f,  p);
    auto rx_info_bits_post = rx_info_from_bit_llr(post_f, p);

    std::cout << "[INFO] rx_info_bits: " << rx_info_bits_pre.size()
              << " (flattened, warmup skipped)\n";

    compute_and_print_ber(info_bits, rx_info_bits_pre,  "Pre-FEC",  p);
    compute_and_print_ber(info_bits, rx_info_bits_post, "Post-FEC", p);

    std::cout << "[DONE] Pipeline(int8) bits -> oFEC -> QAM -> AWGN -> QAM LLR -> quant(int8) -> decode -> info extract & BER\n";
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

    // 4) QAM 调制（示例：QPSK/4-QAM => n_bps = 2）
    const unsigned n_bps = 2; // QPSK
    auto tx_syms = qam_modulate(coded_bits, n_bps);
    std::cout << "[INFO] Modulated symbols: " << tx_syms.size() << " (Es≈1)\n";

    // 5) AWGN（按 Eb/N0 设置噪声）
    const float ebn0_dB = 4.0f;
    const int   N        = static_cast<int>(p.NUM_SUBBLOCK_COLS * p.BITS_PER_SUBBLOCK_DIM); // 128
    const int   K        = 239;
    const int   TAKEBITS = K - N; // 111
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

    // 8) 分发：优先走 qfloat<4>/qfloat<5>，否则回退 int8 路径
    if (p.LLR_BITS == 5) {
        run_with_qfloat<5>(llr_mat, p, info_bits);
    } else if (p.LLR_BITS == 4) {
        run_with_qfloat<4>(llr_mat, p, info_bits);
    } else {
        run_with_int8(llr_mat, p, info_bits);
    }

    return 0;
}
