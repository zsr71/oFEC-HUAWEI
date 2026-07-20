#pragma once
#include <vector>
#include <cstdint>
#include "newcode/params.hpp"

namespace newcode {

struct BerStats {
    std::size_t errors{0};
    std::size_t total{0};
    double      ber{0.0};
};

struct WindowBerStats {
    std::size_t window_idx{0};
    std::size_t errors{0};
    std::size_t total{0};
    double      ber{0.0};
};

/**
 * 计算 BER（仅比较中间有效区间）：
 * 根据 Params 的窗口高度，丢弃首尾各一个 window 覆盖的比特，
 * 在剩余比特上计算误码率。若长度不足将自动做边界保护。
 *
 * 不打印，仅返回统计结果。
 */
BerStats compute_ber(const std::vector<uint8_t>& ref_bits,
                     const std::vector<uint8_t>& rx_bits,
                     const Params& p,
                     std::vector<std::size_t>* error_positions = nullptr);

/**
 * 计算并打印 BER（同上规则丢弃首尾 window）。
 * 打印格式：
 *   "[RESULT] <label> BER=<ber>  (errs=<errors> / <total> compared, cut=<cut> of <L>)"
 * 返回统计结果。
 */
BerStats compute_and_print_ber(const std::vector<uint8_t>& ref_bits,
                               const std::vector<uint8_t>& rx_bits,
                               const char* label,
                               const Params& p,
                               std::vector<std::size_t>* error_positions = nullptr,
                               bool quiet = false);

/**
 * 按 window 输出 BER 统计，不做现有 compute_ber() 的首尾 window 裁剪。
 *
 * 这里的输入应当已经是“完成 info_extract 之后”的一维信息比特流，
 * 即已保持当前 warm-up 机制，但尚未做 BER 的前后 window 截断。
 */
std::vector<WindowBerStats> compute_ber_per_window(
    const std::vector<uint8_t>& ref_bits,
    const std::vector<uint8_t>& rx_bits,
    const Params& p);

/**
 * 按 tile 高度输出 BER 统计。
 *
 * 与 compute_ber_per_window() 类似，但窗口长度改为一个 tile 覆盖的比特数。
 * 这里仍然基于一维信息比特流做局部切片统计，适合做 tile 尺度的时域分析。
 */
std::vector<WindowBerStats> compute_ber_per_tile_window(
    const std::vector<uint8_t>& ref_bits,
    const std::vector<uint8_t>& rx_bits,
    const Params& p);

} // namespace newcode
