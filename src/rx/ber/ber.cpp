#include "newcode/rx/ber/ber.hpp"
#include "newcode/params.hpp"
#include <algorithm>
#include <iostream>

namespace newcode {
namespace {

constexpr std::size_t kInfoBitsPerRow = 111;
// Keep a single fixed comparison interval for every buffered-FIFO size.
// Two trailing windows cover the observed undrained FIFO tail for R_buf <= 32.
constexpr std::size_t kBerSkipPrefixWindows = 4;
constexpr std::size_t kBerSkipSuffixWindows = 2;

std::vector<WindowBerStats> compute_ber_per_segment_bits(
    const std::vector<uint8_t>& ref_bits,
    const std::vector<uint8_t>& rx_bits,
    std::size_t segment_bits)
{
    const std::size_t L = std::min(ref_bits.size(), rx_bits.size());
    std::vector<WindowBerStats> segments;
    if (L == 0 || segment_bits == 0) {
        return segments;
    }

    const std::size_t num_segments = (L + segment_bits - 1) / segment_bits;
    segments.reserve(num_segments);
    for (std::size_t idx = 0; idx < num_segments; ++idx) {
        const std::size_t start = idx * segment_bits;
        const std::size_t stop = std::min(L, start + segment_bits);
        std::size_t err = 0;
        const std::size_t total = stop - start;
        for (std::size_t i = start; i < stop; ++i) {
            if ((ref_bits[i] ^ rx_bits[i]) & 1u) {
                ++err;
            }
        }

        WindowBerStats s;
        s.window_idx = idx;
        s.errors = err;
        s.total = total;
        s.ber = (s.total == 0) ? 0.0
                               : static_cast<double>(s.errors) / static_cast<double>(s.total);
        segments.push_back(s);
    }
    return segments;
}

} // namespace

BerStats compute_ber(const std::vector<uint8_t>& ref_bits,
                     const std::vector<uint8_t>& rx_bits,
                     const Params& p,
                     std::vector<std::size_t>* error_positions)
{
    const std::size_t L = std::min(ref_bits.size(), rx_bits.size());

    // 每行比特数  = 111（由 Params 得出，不再依赖矩阵形状）
    const std::size_t row_bits = kInfoBitsPerRow;

    // 窗口高度（比特行） * 每行比特数 = 一个 window 覆盖的比特数
    const std::size_t win_rows  = p.win_height_rows(); // 已是“比特行”数量
    const std::size_t win_bits  = win_rows * row_bits;

    // 固定裁剪区间：前四个 warm-up window，以及两个帧尾 window。
    // 所有 R_buf 使用同一范围，不能随 FIFO 深度动态改变 BER 分母。
    const std::size_t skip_prefix =
        std::min(L, kBerSkipPrefixWindows * win_bits);
    const std::size_t skip_suffix = std::min(
        L - skip_prefix, kBerSkipSuffixWindows * win_bits);

    const std::size_t start = skip_prefix;
    const std::size_t stop  = L - skip_suffix;

    if (error_positions) error_positions->clear();
    std::size_t err = 0;
    const std::size_t total = stop - start;
    for (std::size_t i = start; i < stop; ++i) {
        if ((ref_bits[i] ^ rx_bits[i]) & 1u) {
            ++err;
            if (error_positions) error_positions->push_back(i);
        }
    }

    BerStats s;
    s.errors = err;
    s.total  = total;
    s.ber    = (s.total == 0) ? 0.0 : static_cast<double>(err) / static_cast<double>(s.total);
    return s;
}

BerStats compute_and_print_ber(const std::vector<uint8_t>& ref_bits,
                               const std::vector<uint8_t>& rx_bits,
                               const char* label,
                               const Params& p,
                               std::vector<std::size_t>* error_positions,
                               bool quiet)
{
    BerStats s = compute_ber(ref_bits, rx_bits, p, error_positions);

    // 为了可见性，按 compute_ber() 的实际规则打印固定边界裁剪量。
    const std::size_t L = std::min(ref_bits.size(), rx_bits.size());
    const std::size_t row_bits = kInfoBitsPerRow;
    const std::size_t win_bits = p.win_height_rows() * row_bits;
    const std::size_t skip_prefix =
        std::min(L, kBerSkipPrefixWindows * win_bits);
    const std::size_t skip_suffix = std::min(
        L - skip_prefix, kBerSkipSuffixWindows * win_bits);
    const std::size_t cut_total = skip_prefix + skip_suffix;

    if (!quiet) {
        std::cout << "[RESULT] " << (label ? label : "BER")
                  << " BER=" << s.ber
                  << "  (errs=" << s.errors << " / " << s.total << " compared"
                  << ", cut=" << cut_total << " of " << L << ")\n";
    }
    return s;
}

std::vector<WindowBerStats> compute_ber_per_window(
    const std::vector<uint8_t>& ref_bits,
    const std::vector<uint8_t>& rx_bits,
    const Params& p)
{
    const std::size_t win_bits = p.win_height_rows() * kInfoBitsPerRow;
    return compute_ber_per_segment_bits(ref_bits, rx_bits, win_bits);
}

std::vector<WindowBerStats> compute_ber_per_tile_window(
    const std::vector<uint8_t>& ref_bits,
    const std::vector<uint8_t>& rx_bits,
    const Params& p)
{
    const std::size_t tile_bits = p.tile_height_rows() * kInfoBitsPerRow;
    return compute_ber_per_segment_bits(ref_bits, rx_bits, tile_bits);
}

} // namespace newcode
