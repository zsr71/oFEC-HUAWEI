#pragma once
#include <vector>
#include <complex>
#include <cstdint>

namespace mod {

// 将 0/1 比特流按 n_bps 分组做调制：
// - n_bps = 1 时执行 BPSK（实轴 ±1）；
// - n_bps 为偶数时执行方形 Gray QAM。
// 输出单位平均能量的复符号。
// - bits: 0/1 比特（其他值会按最低位取 &1）
// - n_bps: 每个符号的比特数（1 或正偶数，如 2=QPSK, 4=16-QAM, 6=64-QAM, ...）
// - 返回：长度 = ceil(bits.size() / n_bps) 的 std::vector<std::complex<float>>
std::vector<std::complex<float>>
qam_modulate(const std::vector<uint8_t>& bits, unsigned n_bps);

} // namespace newcode
