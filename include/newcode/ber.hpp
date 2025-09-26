#pragma once
#include <vector>
#include <cstdint>

namespace newcode {

struct BerStats {
    std::size_t errors{0};
    std::size_t total{0};
    double      ber{0.0};
};

/**
 * 计算 BER（按两个比特序列的重叠最短长度比较）。
 * 不打印，仅返回统计结果。
 */
BerStats compute_ber(const std::vector<uint8_t>& ref_bits,
                     const std::vector<uint8_t>& rx_bits);

/**
 * 计算并打印 BER。
 * 打印格式： "[RESULT] <label> BER=<ber>  (errs=<errors> / <total> compared)"
 * 返回统计结果。
 */
BerStats compute_and_print_ber(const std::vector<uint8_t>& ref_bits,
                               const std::vector<uint8_t>& rx_bits,
                               const char* label = "BER");

} // namespace newcode
