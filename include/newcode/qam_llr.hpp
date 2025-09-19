#pragma once
#include <vector>
#include <complex>
#include <cstdint>

namespace newcode {

// 基础接口：已知噪声标准差 sigma（每一维），计算 QAM LLR（Max-Log）
// - y: 接收复符号（经过 AWGN 后）
// - n_bps: 每符号比特数（偶数，QPSK=2, 16QAM=4, 64QAM=6, ...）
// - sigma: AWGN 每一维标准差（与 add_awgn 用的一致）
std::vector<float>
qam_llr_maxlog(const std::vector<std::complex<float>>& y,
               unsigned n_bps,
               float sigma);

// 便捷包装：用 Eb/N0 与码率换算 sigma 后再解调（等价于先调用 ebn0_to_sigma 再 qam_llr_maxlog）
std::vector<float>
qam_llr_from_ebn0(const std::vector<std::complex<float>>& y,
                  unsigned n_bps,
                  float ebn0_dB,
                  float code_rate);

} // namespace newcode
