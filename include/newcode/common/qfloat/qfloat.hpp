#pragma once
#include <cstdint>
#include <cmath>
#include <algorithm>
#include <type_traits>
#include <vector>
#include <limits>
#include "newcode/common/matrix/matrix.hpp"
namespace qfloat {

// 线性量化的“类浮点”类型：NBITS ∈ [2,15]，mid-tread 对称量化到 [-Q,+Q]
// 内部存码值（整型），对外提供 + - * /、比较、与 float 互转。
template<int NBITS, typename Store=int16_t>
class qfloat {
    static_assert(NBITS >= 2 && NBITS <= 15, "qfloat: NBITS must be in [2,15].");
    static_assert(std::is_signed<Store>::value && sizeof(Store) >= 2,
                  "qfloat: Store must be a signed integer type of at least 16 bits.");
public:
    static constexpr int   Q()  { return (1 << (NBITS - 1)) - 1; }    //最大正整数编码整个定点数用 1 位符号 + (NBITS-1) 位数值 来表示 小数部分，整数部分则由剩余的 Store 位数决定
    static constexpr int   LO() { return -Q(); }
    static constexpr int   HI() { return +Q(); }
    static float  current_clip() { return current_clip_ref(); }
    static void   set_clip(float clip) { current_clip_ref() = clip; }

    qfloat() : code_(0) {}
    explicit qfloat(int code, float clip = current_clip()) : code_(sat(code)) { (void)clip; }
    explicit qfloat(float x, float clip = current_clip()) : code_(0) { code_ = quantize(x, clip); }

    static qfloat from_float(float x, float clip = current_clip()) { return qfloat(x, clip); }
    float to_float() const {
        const float clip = current_clip();
        return static_cast<float>(code_) * (clip / static_cast<float>(Q()));
    }
    static float clip() { return current_clip(); }

    // 算术（码值域饱和）
    qfloat  operator-() const { return qfloat(sat(-code_)); }
    qfloat& operator+=(qfloat rhs) { code_ = sat(int(code_) + int(rhs.code_)); return *this; }
    qfloat& operator-=(qfloat rhs) { code_ = sat(int(code_) - int(rhs.code_)); return *this; }
    friend qfloat operator+(qfloat a, qfloat b) { a += b; return a; }
    friend qfloat operator-(qfloat a, qfloat b) { a -= b; return a; }

    qfloat& operator*=(float k) {
        code_ = sat(static_cast<int>(std::lrint(static_cast<float>(code_) * k)));
        return *this;
    }
    qfloat& operator/=(float k) {
        if (k != 0.0f)
            code_ = sat(static_cast<int>(std::lrint(static_cast<float>(code_) / k)));
        return *this;
    }
    friend qfloat operator*(qfloat a, float k) { a *= k; return a; }
    friend qfloat operator*(float k, qfloat a) { a *= k; return a; }
    friend qfloat operator/(qfloat a, float k) { a /= k; return a; }

    // 同类型相乘/除：按真实幅值相乘/除，再重量化
    friend qfloat operator*(qfloat a, qfloat b) {
        float y = a.to_float() * b.to_float();
        return qfloat::from_float(y);
    }
    friend qfloat operator/(qfloat a, qfloat b) {
        float den = b.to_float();
        if (den == 0.0f) return a; // 或者返回饱和值，按需改
        float y = a.to_float() / den;
        return qfloat::from_float(y);
    }

    // 比较（码值比较）
    friend bool operator<(qfloat a, qfloat b)  { return a.code_ <  b.code_; }
    friend bool operator>(qfloat a, qfloat b)  { return a.code_ >  b.code_; }
    friend bool operator<=(qfloat a, qfloat b) { return a.code_ <= b.code_; }
    friend bool operator>=(qfloat a, qfloat b) { return a.code_ >= b.code_; }
    friend bool operator==(qfloat a, qfloat b) { return a.code_ == b.code_; }
    friend bool operator!=(qfloat a, qfloat b) { return a.code_ != b.code_; }

    // 转换/访问
    explicit operator float() const { return to_float(); }
    explicit operator int()   const { return static_cast<int>(code_); }
    int   code() const { return static_cast<int>(code_); }
    void  set_code(int c) { code_ = sat(c); }

private:
    Store code_ = 0;

    static float& current_clip_ref() {
        thread_local float clip = 0.0f;
        return clip;
    }

    static int   sat(int x) {
        if (x > HI()) return HI();
        if (x < LO()) return LO();
        return x;
    }
    static Store quantize(float x, float clip) {
        if (clip <= 0.f) return static_cast<Store>(0);
        float xc = std::max(-clip, std::min(+clip, x));
        int   c  = static_cast<int>(std::lrint(xc * (static_cast<float>(Q()) / clip)));
        return static_cast<Store>(sat(c));
    }
};

// 便捷别名
using float_3  = qfloat<3>;
using float_4  = qfloat<4>;
using float_5  = qfloat<5>;
using float_6  = qfloat<6>;
using float_7  = qfloat<7>;
using float_8  = qfloat<8>;
using float_9  = qfloat<9>;
using float_10 = qfloat<10>;

// 批量工具（与 Matrix 协作）
template<typename T> class Matrix;

template<int NBITS>
matrix::Matrix< qfloat<NBITS> > quantize_matrix_to_qfloat(const matrix::Matrix<float>& in,
                                                  float clip = qfloat<NBITS>::current_clip())
{
    matrix::Matrix< qfloat<NBITS> > out(in.rows(), in.cols());
    for (size_t r=0; r<in.rows(); ++r)
        for (size_t c=0; c<in.cols(); ++c)
            out[r][c] = qfloat<NBITS>::from_float(in[r][c], clip);
    return out;
}

template<int NBITS>
matrix::Matrix<float> dequantize_matrix_from_qfloat(const matrix::Matrix< qfloat<NBITS> >& in,
                                            float clip_unused = qfloat<NBITS>::current_clip())
{
    (void)clip_unused;
    matrix::Matrix<float> out(in.rows(), in.cols());
    for (size_t r=0; r<in.rows(); ++r)
        for (size_t c=0; c<in.cols(); ++c)
            out[r][c] = in[r][c].to_float();
    return out;
}

// 仅类型转换为 float（直接用码值，不做幅度缩放；阈值仍是 0）
template<int NBITS>
matrix::Matrix<float> cast_matrix_from_qfloat(const matrix::Matrix< qfloat<NBITS> >& in)
{
    matrix::Matrix<float> out(in.rows(), in.cols());
    for (size_t r=0; r<in.rows(); ++r)
        for (size_t c=0; c<in.cols(); ++c)
            out[r][c] = static_cast<float>(static_cast<int>(in[r][c]));
    return out;
}

} // namespace newcode

namespace newcode {

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

} // namespace newcode
