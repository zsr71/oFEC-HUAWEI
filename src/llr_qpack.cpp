#include "newcode/llr_qpack.hpp"
#include "newcode/matrix.hpp"
#include "newcode/params.hpp"
#include <algorithm>
#include <cmath>

namespace newcode {

static inline float clampf(float x, float lo, float hi) {
    return std::max(lo, std::min(hi, x));
}

Matrix<int8_t> quantize_llr_to_int8(const Matrix<float>& in, const Params& p) {
    Matrix<int8_t> out(in.rows(), in.cols());
    const std::size_t llr_bits = p.LLR_BITS;
    const float clip = p.LLR_CLIP;

    if (llr_bits < 2 || clip <= 0.0f) {
        // 退化：不量化，直接截断到 int8 范围
        for (size_t r = 0; r < in.rows(); ++r)
            for (size_t c = 0; c < in.cols(); ++c) {
                float x = in[r][c];
                if (x > 127.f) x = 127.f;
                if (x < -128.f) x = -128.f;
                out[r][c] = static_cast<int8_t>(std::nearbyint(x));
            }
        return out;
    }

    const int Q = (1 << (llr_bits - 1)) - 1;   // 5bit -> 15, 4bit -> 7
    const float scale = static_cast<float>(Q) / clip;

    for (size_t r = 0; r < in.rows(); ++r) {
        for (size_t c = 0; c < in.cols(); ++c) {
            const float x  = clampf(in[r][c], -clip, +clip);
            int qi = static_cast<int>(std::nearbyint(x * scale));
            if (qi >  Q) qi =  Q;
            if (qi < -Q) qi = -Q;
            out[r][c] = static_cast<int8_t>(qi);  // 存在 int8_t 中
        }
    }
    return out;
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
