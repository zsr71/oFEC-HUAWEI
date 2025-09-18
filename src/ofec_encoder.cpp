#include "newcode/ofec_encoder.hpp"

namespace newcode {

Matrix<uint8_t> ofec_encode(const std::vector<uint8_t>& bits, const Params& p)
{
    // TODO: 实现 oFEC 编码逻辑
    // 当前仅返回一个占位矩阵，尺寸由系统参数决定：
    // 行数 = (信息子行 + 上下保护子行) × 每子块大小
    // 列数 = 子块列数 × 每子块大小

    Matrix<uint8_t> mat = Matrix<uint8_t>::zero(
        (p.INFO_SUBROWS_PER_CODE + 2 * p.NUM_GUARD_SUBROWS) * p.BITS_PER_SUBBLOCK_DIM,
        p.NUM_SUBBLOCK_COLS * p.BITS_PER_SUBBLOCK_DIM
    );

    // TODO: 使用 bits 填充 mat

    return mat;
}

} // namespace newcode
