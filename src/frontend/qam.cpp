#include "newcode/qam.hpp"
#include <stdexcept>
#include <cmath>

namespace newcode {

static inline std::complex<float>
bits_to_symbol_gray_qam(const uint8_t* b, unsigned n_bps)
{
    // 参考 AFF3CT：前一半比特映 I，后一半比特映 Q；逐级展开得到 {±1, ±3, ...} 幅度
    const unsigned m = n_bps / 2;

    float I = 1.0f - 2.0f * float(b[0] & 1u);          // +1 or -1
    float Q = 1.0f - 2.0f * float(b[m] & 1u);          // +1 or -1

    for (unsigned j = 1; j < m; ++j)
    {
        I = (1.0f - 2.0f * float(b[j]   & 1u)) * (float(1u << j) - I);
        Q = (1.0f - 2.0f * float(b[m+j] & 1u)) * (float(1u << j) - Q);
    }

    // 能量归一化：Es = 2*(M-1)/3, sqrt(Es) = sqrt(2*(M-1)/3)
    const unsigned M = 1u << n_bps;
    const float sqrt_es = std::sqrt(2.0f * float(M - 1) / 3.0f);
    return std::complex<float>(I / sqrt_es, Q / sqrt_es);
}

std::vector<std::complex<float>>
qam_modulate(const std::vector<uint8_t>& bits, unsigned n_bps)
{
    if (n_bps == 0 || (n_bps & 1u))
        throw std::invalid_argument("qam_modulate: n_bps must be a positive even number.");

    const size_t n_sym = (bits.size() + n_bps - 1) / n_bps; // ceil
    std::vector<std::complex<float>> syms;
    syms.reserve(n_sym);

    std::vector<uint8_t> buf(n_bps, 0); // 每符号的局部比特缓存（不足补 0）
    size_t p = 0;
    for (size_t s = 0; s < n_sym; ++s)
    {
        // 拷贝 n_bps 个比特，不足补 0（只取最低位）
        for (unsigned k = 0; k < n_bps; ++k)
            buf[k] = (p < bits.size()) ? (bits[p++] & 1u) : 0u;

        syms.emplace_back(bits_to_symbol_gray_qam(buf.data(), n_bps));
    }

    return syms;
}

} // namespace newcode
