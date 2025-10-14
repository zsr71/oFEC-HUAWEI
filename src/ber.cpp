#include "newcode/ber.hpp"
#include "newcode/params.hpp"
#include <algorithm>
#include <iostream>

namespace newcode {

static inline std::size_t saturating_mul(std::size_t a, std::size_t b, std::size_t cap)
{
    if (a == 0 || b == 0) return 0;
    if (a > cap / b) return cap; // 防溢出：上饱和到 cap
    return a * b;
}

BerStats compute_ber(const std::vector<uint8_t>& ref_bits,
                     const std::vector<uint8_t>& rx_bits,
                     const Params& p)
{
    const std::size_t L = std::min(ref_bits.size(), rx_bits.size());

    // 每行比特数 = 8 * 16 = 128（由 Params 得出，不再依赖矩阵形状）
    const std::size_t row_bits = p.NUM_SUBBLOCK_COLS * p.BITS_PER_SUBBLOCK_DIM;

    // 窗口高度（比特行） * 每行比特数 = 一个 window 覆盖的比特数
    const std::size_t win_rows  = p.win_height_rows(); // 已是“比特行”数量
    const std::size_t win_bits  = saturating_mul(win_rows, row_bits, L);

    // 去掉首尾各一个 window 覆盖的比特
    const std::size_t skip_prefix = std::min(L, 2*win_bits);
    const std::size_t skip_suffix = std::min(L - skip_prefix, win_bits);

    const std::size_t start = skip_prefix;
    const std::size_t stop  = L - skip_suffix;

    std::size_t err = 0;
    for (std::size_t i = start; i < stop; ++i) {
        err += (ref_bits[i] ^ rx_bits[i]) & 1u;
    }

    BerStats s;
    s.errors = err;
    s.total  = (stop > start) ? (stop - start) : 0;
    s.ber    = (s.total == 0) ? 0.0 : static_cast<double>(err) / static_cast<double>(s.total);
    return s;
}

BerStats compute_and_print_ber(const std::vector<uint8_t>& ref_bits,
                               const std::vector<uint8_t>& rx_bits,
                               const char* label,
                               const Params& p)
{
    BerStats s = compute_ber(ref_bits, rx_bits, p);

    // 为了可见性，把被剔除的前后窗口比特数也打印出来
    const std::size_t L = std::min(ref_bits.size(), rx_bits.size());
    const std::size_t row_bits = 111;
    const std::size_t win_bits = std::min(L, p.win_height_rows() * row_bits);
    const std::size_t cut_total = std::min(L, win_bits) + std::min(L > win_bits ? (L - win_bits) : 0, win_bits);

    std::cout << "[RESULT] " << (label ? label : "BER")
              << " BER=" << s.ber
              << "  (errs=" << s.errors << " / " << s.total << " compared"
              << ", cut=" << cut_total << " of " << L << ")\n";
    return s;
}

} // namespace newcode
