#include "new_float_only/common/matrix/info_extract.hpp"
#include <cassert>
#include <algorithm>

namespace matrix {

// 把单个比特位的 LLR 硬判成 0/1。
// 当前约定是：LLR >= 0 判 0，LLR < 0 判 1。
static inline uint8_t hard_decide(float L) noexcept { return (L >= 0.f) ? 0u : 1u; }

// 从按 bit 排布的 LLR 矩阵中提取最终的信息比特序列。
//
// 这一步假设输入矩阵已经和发送端的 code_matrix 对齐：
// - 每一行对应一个 oFEC 行码字在接收端恢复后的 bit LLR
// - 每一行只取前 TAKE_BITS 列作为“新信息位”
// - 前面的 warmup_rows 行属于滑窗初始化区域，不参与真实信息输出
//
// 返回值是一维比特数组，按“逐行、每行从左到右”的顺序展平。
std::vector<uint8_t> rx_info_from_bit_llr(const matrix::Matrix<float>& bit_llr_mat, const new_float_only::Params& p)
{
    // 基本布局参数。N 是每行的比特列数，K 是 BCH 信息位数。
    const int B = static_cast<int>(p.BITS_PER_SUBBLOCK_DIM);
    const int N = static_cast<int>(p.NUM_SUBBLOCK_COLS * p.BITS_PER_SUBBLOCK_DIM); // 典型 128
    const int K          = 239;   // BCH 信息位

    // 每行真正承载“新信息”的列数。当前布局下等于 239 - 128 = 111。
    const int TAKE_BITS  = K - N; // 每行新信息位列数（典型 111）


    // 先做尺寸合法性检查，避免后面按固定布局取列时越界。
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

    // 跳过前置初始化行。发送端最前面有一块 warmup 区域，用来给滑窗解码提供
    // 初始已知历史；这些行不对应真正需要输出的信息比特。
    size_t warmup_rows = 0;
    {
        const long tmp = static_cast<long>(p.tile_height_rows());
        if (tmp > 0) warmup_rows = static_cast<size_t>(tmp);
        warmup_rows = std::min(warmup_rows, RROWS); // 防溢出
    }

    const size_t useful_rows = RROWS - warmup_rows;

    std::vector<uint8_t> rx_info;
    rx_info.reserve(useful_rows * static_cast<size_t>(TAKE_BITS));

    // 逐行提取前 TAKE_BITS 列，并把每个 LLR 硬判成 0/1 后展平成一维数组。
    // 这里不看 BCH parity、history 区或 overall parity，只取发送端定义的信息位部分。
    for (size_t r = warmup_rows; r < RROWS; ++r) {
        for (int c = 0; c < TAKE_BITS; ++c) {
            rx_info.push_back(hard_decide(bit_llr_mat[r][static_cast<size_t>(c)]));
        }
    }

    return rx_info;
}

} // namespace new_float_only
