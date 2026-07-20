// path: newcode/rx/demod/qam_llr.hpp
#pragma once
#include <vector>
#include <complex>
#include <cstdint>

namespace demod {



// 精确版：BPSK 走解析 LLR，其余调制走 log-sum-exp（数值稳定、准确）
std::vector<float>
qam_llr_logsumexp(const std::vector<std::complex<float>>& y,
                  unsigned n_bps,
                  float sigma);

// 便捷包装：用 Eb/N0 与码率换算 sigma 后再解调
// 【注意】默认采用精确（log-sum-exp）实现
std::vector<float>
qam_llr_from_ebn0(const std::vector<std::complex<float>>& y,
                  unsigned n_bps,
                  float ebn0_dB,
                  float code_rate);

} // namespace newcode
