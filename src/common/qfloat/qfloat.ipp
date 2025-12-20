#pragma once

#include "newcode/common/qfloat/qfloat.hpp"

#include <algorithm>
#include <cmath>
#include <functional>
#include <limits>
#include <utility>
#include <vector>

namespace qfloat {

template<int NBITS, typename Store>
constexpr int qfloat<NBITS, Store>::Q() { return (1 << (NBITS - 1)) - 1; }

template<int NBITS, typename Store>
constexpr int qfloat<NBITS, Store>::LO() { return -Q(); }

template<int NBITS, typename Store>
constexpr int qfloat<NBITS, Store>::HI() { return +Q(); }

template<int NBITS, typename Store>
float qfloat<NBITS, Store>::current_clip() { return current_clip_ref(); }

template<int NBITS, typename Store>
void qfloat<NBITS, Store>::set_clip(float clip) { current_clip_ref() = clip; }

template<int NBITS, typename Store>
qfloat<NBITS, Store>::qfloat() : code_(0) {}

template<int NBITS, typename Store>
qfloat<NBITS, Store>::qfloat(int code, float clip) : code_(sat(code)) { (void)clip; }

template<int NBITS, typename Store>
qfloat<NBITS, Store>::qfloat(float x, float clip) : code_(0) { code_ = quantize(x, clip); }

template<int NBITS, typename Store>
qfloat<NBITS, Store> qfloat<NBITS, Store>::from_float(float x, float clip)
{
    return qfloat<NBITS, Store>(x, clip);
}

template<int NBITS, typename Store>
float qfloat<NBITS, Store>::to_float() const
{
    const float clip = current_clip();
    return static_cast<float>(code_) * (clip / static_cast<float>(Q()));
}

template<int NBITS, typename Store>
float qfloat<NBITS, Store>::clip() { return current_clip(); }

template<int NBITS, typename Store>
qfloat<NBITS, Store> qfloat<NBITS, Store>::operator-() const
{
    return qfloat<NBITS, Store>(sat(-code_));
}

template<int NBITS, typename Store>
qfloat<NBITS, Store>& qfloat<NBITS, Store>::operator+=(qfloat rhs)
{
    code_ = sat(int(code_) + int(rhs.code_));
    return *this;
}

template<int NBITS, typename Store>
qfloat<NBITS, Store>& qfloat<NBITS, Store>::operator-=(qfloat rhs)
{
    code_ = sat(int(code_) - int(rhs.code_));
    return *this;
}

template<int NBITS, typename Store>
qfloat<NBITS, Store>& qfloat<NBITS, Store>::operator*=(float k)
{
    code_ = sat(static_cast<int>(std::lrint(static_cast<float>(code_) * k)));
    return *this;
}

template<int NBITS, typename Store>
qfloat<NBITS, Store>& qfloat<NBITS, Store>::operator/=(float k)
{
    if (k != 0.0f)
        code_ = sat(static_cast<int>(std::lrint(static_cast<float>(code_) / k)));
    return *this;
}

template<int NBITS, typename Store>
qfloat<NBITS, Store> operator+(qfloat<NBITS, Store> a, qfloat<NBITS, Store> b)
{
    a += b;
    return a;
}

template<int NBITS, typename Store>
qfloat<NBITS, Store> operator-(qfloat<NBITS, Store> a, qfloat<NBITS, Store> b)
{
    a -= b;
    return a;
}

template<int NBITS, typename Store>
qfloat<NBITS, Store> operator*(qfloat<NBITS, Store> a, float k)
{
    a *= k;
    return a;
}

template<int NBITS, typename Store>
qfloat<NBITS, Store> operator*(float k, qfloat<NBITS, Store> a)
{
    a *= k;
    return a;
}

template<int NBITS, typename Store>
qfloat<NBITS, Store> operator/(qfloat<NBITS, Store> a, float k)
{
    a /= k;
    return a;
}

template<int NBITS, typename Store>
qfloat<NBITS, Store> operator*(qfloat<NBITS, Store> a, qfloat<NBITS, Store> b)
{
    float y = a.to_float() * b.to_float();
    return qfloat<NBITS, Store>::from_float(y);
}

template<int NBITS, typename Store>
qfloat<NBITS, Store> operator/(qfloat<NBITS, Store> a, qfloat<NBITS, Store> b)
{
    float den = b.to_float();
    if (den == 0.0f) return a;
    float y = a.to_float() / den;
    return qfloat<NBITS, Store>::from_float(y);
}

template<int NBITS, typename Store>
bool operator<(qfloat<NBITS, Store> a, qfloat<NBITS, Store> b) { return a.code() <  b.code(); }

template<int NBITS, typename Store>
bool operator>(qfloat<NBITS, Store> a, qfloat<NBITS, Store> b) { return a.code() >  b.code(); }

template<int NBITS, typename Store>
bool operator<=(qfloat<NBITS, Store> a, qfloat<NBITS, Store> b) { return a.code() <= b.code(); }

template<int NBITS, typename Store>
bool operator>=(qfloat<NBITS, Store> a, qfloat<NBITS, Store> b) { return a.code() >= b.code(); }

template<int NBITS, typename Store>
bool operator==(qfloat<NBITS, Store> a, qfloat<NBITS, Store> b) { return a.code() == b.code(); }

template<int NBITS, typename Store>
bool operator!=(qfloat<NBITS, Store> a, qfloat<NBITS, Store> b) { return a.code() != b.code(); }

template<int NBITS, typename Store>
qfloat<NBITS, Store>::operator float() const { return to_float(); }

template<int NBITS, typename Store>
qfloat<NBITS, Store>::operator int() const { return static_cast<int>(code_); }

template<int NBITS, typename Store>
int qfloat<NBITS, Store>::code() const { return static_cast<int>(code_); }

template<int NBITS, typename Store>
void qfloat<NBITS, Store>::set_code(int c) { code_ = sat(c); }

template<int NBITS, typename Store>
float& qfloat<NBITS, Store>::current_clip_ref()
{
    thread_local float clip = 0.0f;
    return clip;
}

template<int NBITS, typename Store>
int qfloat<NBITS, Store>::sat(int x)
{
    if (x > HI()) return HI();
    if (x < LO()) return LO();
    return x;
}

template<int NBITS, typename Store>
Store qfloat<NBITS, Store>::quantize(float x, float clip)
{
    if (clip <= 0.f) return static_cast<Store>(0);
    float xc = std::max(-clip, std::min(+clip, x));
    int   c  = static_cast<int>(std::lrint(xc * (static_cast<float>(Q()) / clip)));
    return static_cast<Store>(sat(c));
}

// Batch helpers (Matrix interop)
template<int NBITS>
matrix::Matrix< qfloat<NBITS> > quantize_matrix_to_qfloat(const matrix::Matrix<float>& in,
                                                          float clip)
{
    matrix::Matrix< qfloat<NBITS> > out(in.rows(), in.cols());
    for (size_t r = 0; r < in.rows(); ++r)
        for (size_t c = 0; c < in.cols(); ++c)
            out[r][c] = qfloat<NBITS>::from_float(in[r][c], clip);
    return out;
}

template<int NBITS>
matrix::Matrix<float> dequantize_matrix_from_qfloat(const matrix::Matrix< qfloat<NBITS> >& in,
                                                     float clip_unused)
{
    (void)clip_unused;
    matrix::Matrix<float> out(in.rows(), in.cols());
    for (size_t r = 0; r < in.rows(); ++r)
        for (size_t c = 0; c < in.cols(); ++c)
            out[r][c] = in[r][c].to_float();
    return out;
}

template<int NBITS>
matrix::Matrix<float> cast_matrix_from_qfloat(const matrix::Matrix< qfloat<NBITS> >& in)
{
    matrix::Matrix<float> out(in.rows(), in.cols());
    for (size_t r = 0; r < in.rows(); ++r)
        for (size_t c = 0; c < in.cols(); ++c)
            out[r][c] = static_cast<float>(static_cast<int>(in[r][c]));
    return out;
}

template<typename InputIt>
float compute_clip_from_ratio(InputIt first, InputIt last, float ratio)
{
    std::vector<float> mags;
    for (auto it = first; it != last; ++it) {
        mags.push_back(std::fabs(static_cast<float>(*it)));
    }

    if (mags.empty()) return 0.0f;

    ratio = std::clamp(ratio, 0.0f, 1.0f);
    const std::size_t count = mags.size();
    std::size_t target = static_cast<std::size_t>(std::ceil(ratio * count));
    if (target == 0) target = 1;
    const std::size_t index = target - 1;

    std::nth_element(mags.begin(), mags.begin() + index, mags.end(), std::greater<float>());
    return mags[index];
}

} // namespace qfloat
