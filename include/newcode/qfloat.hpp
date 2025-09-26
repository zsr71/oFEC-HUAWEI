#pragma once
#include <cstdint>
#include <cmath>
#include <algorithm>
#include <type_traits>

namespace newcode {

// 线性量化的“类浮点”类型：NBITS ∈ [3,10]，mid-tread 对称量化到 [-Q,+Q]
// 内部存码值（整型），对外提供 + - * /、比较、与 float 互转。
template<int NBITS, typename Store=int16_t>
class qfloat {
    static_assert(NBITS >= 3 && NBITS <= 10, "qfloat: NBITS must be in [3,10].");
    static_assert(std::is_signed<Store>::value && sizeof(Store) >= 2,
                  "qfloat: Store must be a signed integer type of at least 16 bits.");
public:
    static constexpr float DEFAULT_CLIP = 8.0f;
    static constexpr int   Q()  { return (1 << (NBITS - 1)) - 1; }
    static constexpr int   LO() { return -Q(); }
    static constexpr int   HI() { return +Q(); }

    qfloat() = default;
    explicit qfloat(int code) : code_(sat(code)) {}
    explicit qfloat(float x, float clip = DEFAULT_CLIP) { code_ = quantize(x, clip); }

    static qfloat from_float(float x, float clip = DEFAULT_CLIP) { return qfloat(x, clip); }
    float to_float(float clip = DEFAULT_CLIP) const {
        return static_cast<float>(code_) * (clip / static_cast<float>(Q()));
    }

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
Matrix< qfloat<NBITS> > quantize_matrix_to_qfloat(const Matrix<float>& in,
                                                  float clip = qfloat<NBITS>::DEFAULT_CLIP)
{
    Matrix< qfloat<NBITS> > out(in.rows(), in.cols());
    for (size_t r=0; r<in.rows(); ++r)
        for (size_t c=0; c<in.cols(); ++c)
            out[r][c] = qfloat<NBITS>::from_float(in[r][c], clip);
    return out;
}

template<int NBITS>
Matrix<float> dequantize_matrix_from_qfloat(const Matrix< qfloat<NBITS> >& in,
                                            float clip = qfloat<NBITS>::DEFAULT_CLIP)
{
    Matrix<float> out(in.rows(), in.cols());
    for (size_t r=0; r<in.rows(); ++r)
        for (size_t c=0; c<in.cols(); ++c)
            out[r][c] = in[r][c].to_float(clip);
    return out;
}

// 仅类型转换为 float（直接用码值，不做幅度缩放；阈值仍是 0）
template<int NBITS>
Matrix<float> cast_matrix_from_qfloat(const Matrix< qfloat<NBITS> >& in)
{
    Matrix<float> out(in.rows(), in.cols());
    for (size_t r=0; r<in.rows(); ++r)
        for (size_t c=0; c<in.cols(); ++c)
            out[r][c] = static_cast<float>(static_cast<int>(in[r][c]));
    return out;
}

} // namespace newcode
