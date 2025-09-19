#pragma once
#include <vector>
#include <complex>
#include <cstdint>

namespace newcode {

// 将 0/1 比特流按 n_bps 分组做方形 QAM（格雷映射），输出单位平均能量的调制符号。
// - bits: 0/1 比特（其他值会按最低位取 &1）
// - n_bps: 每个符号的比特数，必须为偶数（如 2=QPSK/4-QAM, 4=16-QAM, 6=64-QAM, ...）
// - 返回：长度 = ceil(bits.size() / n_bps) 的 std::vector<std::complex<float>>
std::vector<std::complex<float>>
qam_modulate(const std::vector<uint8_t>& bits, unsigned n_bps);

} // namespace newcode
