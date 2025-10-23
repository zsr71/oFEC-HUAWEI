#include "newcode/info_extract.hpp"
#include <cassert>
#include <algorithm>

namespace newcode {

static inline uint8_t hard_decide(float L) noexcept { return (L >= 0.f) ? 0u : 1u; }

std::vector<uint8_t> rx_info_from_bit_llr(const Matrix<float>& bit_llr_mat, const Params& p)
{
    // 基本参数
    const int B = static_cast<int>(p.BITS_PER_SUBBLOCK_DIM);
    const int N = static_cast<int>(p.NUM_SUBBLOCK_COLS * p.BITS_PER_SUBBLOCK_DIM); // 典型 128
    const int K          = 239;   // BCH 信息位

    const int TAKE_BITS  = K - N; // 每行新信息位列数（典型 111）


    if (B <= 0 || (N % B) != 0) {
        throw std::invalid_argument("rx_info_from_bit_llr: N must be a multiple of B.");
    }
    if (bit_llr_mat.cols() != static_cast<size_t>(N)) {
        throw std::invalid_argument("rx_info_from_bit_llr: llr_mat.cols != N.");
    }
    if (TAKE_BITS <= 0 || TAKE_BITS > N) {
        throw std::invalid_argument("rx_info_from_bit_llr: invalid TAKE_BITS.");
    }

    const size_t RROWS = bit_llr_mat.rows();

    // 跳过前置初始化行（不是“真实信息行”）
    size_t warmup_rows = 0;
    {
        const long tmp = static_cast<long>(p.win_height_rows());
        if (tmp > 0) warmup_rows = static_cast<size_t>(tmp);
        warmup_rows = std::min(warmup_rows, RROWS); // 防溢出
    }

    const size_t useful_rows = RROWS - warmup_rows;

    std::vector<uint8_t> rx_info;
    rx_info.reserve(useful_rows * static_cast<size_t>(TAKE_BITS));

    // 逐行提取前 TAKE_BITS 列（信息位），并硬判决，展平为一维数组
    for (size_t r = warmup_rows; r < RROWS; ++r) {
        for (int c = 0; c < TAKE_BITS; ++c) {
            rx_info.push_back(hard_decide(bit_llr_mat[r][static_cast<size_t>(c)]));
        }
    }

    return rx_info;
}

} // namespace newcode
