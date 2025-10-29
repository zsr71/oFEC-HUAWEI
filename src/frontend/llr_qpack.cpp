#include "newcode/llr_qpack.hpp"
#include "newcode/matrix.hpp"
#include "newcode/params.hpp"
#include <algorithm>
#include <cmath>

namespace newcode {

static inline float clampf(float x, float lo, float hi) {
    return std::max(lo, std::min(hi, x));
}

Matrix<float> dequantize_llr_to_float(const Matrix<int8_t>& in, const Params& p) {
    Matrix<float> out(in.rows(), in.cols());
    const std::size_t llr_bits = p.LLR_BITS;
    const float clip = p.LLR_CLIP;

    if (llr_bits < 2 || clip <= 0.0f) {
        for (size_t r = 0; r < in.rows(); ++r)
            for (size_t c = 0; c < in.cols(); ++c)
                out[r][c] = static_cast<float>(in[r][c]);
        return out;
    }

    const int Q = (1 << (llr_bits - 1)) - 1;
    const float inv_scale = clip / static_cast<float>(Q);

    for (size_t r = 0; r < in.rows(); ++r)
        for (size_t c = 0; c < in.cols(); ++c)
            out[r][c] = static_cast<float>(in[r][c]) * inv_scale;

    return out;
}

Matrix<float> cast_qllr_to_float(const Matrix<int8_t>& in) {
    Matrix<float> out(in.rows(), in.cols());
    for (size_t r = 0; r < in.rows(); ++r)
        for (size_t c = 0; c < in.cols(); ++c)
            out[r][c] = static_cast<float>(in[r][c]);
    return out;
}

} // namespace newcode
