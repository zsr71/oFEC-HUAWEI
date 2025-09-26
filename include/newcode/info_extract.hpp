#pragma once
#include <vector>
#include <cstdint>
#include <stdexcept>

#include "newcode/matrix.hpp"
#include "newcode/params.hpp"

namespace newcode {

/**
 * 从 bit_llr_mat 提取信息比特（按 oFEC 行布局：前 TAKE_BITS 列为新信息位），
 * 展平成一维数组（行主序），并且会跳过前置初始化的行：
 *   warmup_rows = (2*G + INFO_SUBROWS_PER_CODE) * B
 *
 * 约定与检查：
 *   - B = p.BITS_PER_SUBBLOCK_DIM
 *   - N = p.NUM_SUBBLOCK_COLS * p.BITS_PER_SUBBLOCK_DIM
 *   - K=239, PAR_LEN=16, OVERALL=1，需满足 2N = K+PAR_LEN+OVERALL（典型 N=128）
 *   - 每行有效信息列数 TAKE_BITS = K - N（典型 111）
 *   - bit_llr_mat.cols() 必须等于 N
 *
 * 硬判决：LLR >= 0 → 0；LLR < 0 → 1。
 * 返回：按行主序拼接的“真实信息行”的信息比特。
 */
std::vector<uint8_t> rx_info_from_bit_llr(const Matrix<float>& bit_llr_mat, const Params& p);

} // namespace newcode
