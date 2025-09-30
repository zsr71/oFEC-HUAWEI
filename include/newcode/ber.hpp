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

/**
 * 计算 BER（仅比较中间有效区间）：
 * 根据 Params 的窗口高度，丢弃首尾各一个 window 覆盖的比特，
 * 在剩余比特上计算误码率。若长度不足将自动做边界保护。
 *
 * 不打印，仅返回统计结果。
 */
BerStats compute_ber(const std::vector<uint8_t>& ref_bits,
                     const std::vector<uint8_t>& rx_bits,
                     const Params& p);

/**
 * 计算并打印 BER（同上规则丢弃首尾 window）。
 * 打印格式：
 *   "[RESULT] <label> BER=<ber>  (errs=<errors> / <total> compared, cut=<cut> of <L>)"
 * 返回统计结果。
 */
BerStats compute_and_print_ber(const std::vector<uint8_t>& ref_bits,
                               const std::vector<uint8_t>& rx_bits,
                               const char* label,
                               const Params& p);

} // namespace newcode
