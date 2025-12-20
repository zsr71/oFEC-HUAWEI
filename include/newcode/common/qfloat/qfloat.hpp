#pragma once
#include <cstdint>
#include <cmath>
#include <algorithm>
#include <type_traits>
#include <vector>
#include <limits>
#include "newcode/common/matrix/matrix.hpp"

namespace qfloat {

// Linear quantized "float-like" type, NBITS in [2, 15].
template<int NBITS, typename Store = int16_t>
class qfloat {
    static_assert(NBITS >= 2 && NBITS <= 15, "qfloat: NBITS must be in [2,15].");
    static_assert(std::is_signed<Store>::value && sizeof(Store) >= 2,
                  "qfloat: Store must be a signed integer type of at least 16 bits.");
public:
    static constexpr int   Q();
    static constexpr int   LO();
    static constexpr int   HI();
    static float  current_clip();
    static void   set_clip(float clip);

    qfloat();
    explicit qfloat(int code, float clip = current_clip());
    explicit qfloat(float x, float clip = current_clip());

    static qfloat from_float(float x, float clip = current_clip());
    float to_float() const;
    static float clip();

    qfloat  operator-() const;
    qfloat& operator+=(qfloat rhs);
    qfloat& operator-=(qfloat rhs);

    qfloat& operator*=(float k);
    qfloat& operator/=(float k);

    explicit operator float() const;
    explicit operator int()   const;
    int   code() const;
    void  set_code(int c);

private:
    Store code_ = 0;

    static float& current_clip_ref();
    static int   sat(int x);
    static Store quantize(float x, float clip);
};

template<int NBITS, typename Store>
qfloat<NBITS, Store> operator+(qfloat<NBITS, Store> a, qfloat<NBITS, Store> b);

template<int NBITS, typename Store>
qfloat<NBITS, Store> operator-(qfloat<NBITS, Store> a, qfloat<NBITS, Store> b);

template<int NBITS, typename Store>
qfloat<NBITS, Store> operator*(qfloat<NBITS, Store> a, float k);

template<int NBITS, typename Store>
qfloat<NBITS, Store> operator*(float k, qfloat<NBITS, Store> a);

template<int NBITS, typename Store>
qfloat<NBITS, Store> operator/(qfloat<NBITS, Store> a, float k);

template<int NBITS, typename Store>
qfloat<NBITS, Store> operator*(qfloat<NBITS, Store> a, qfloat<NBITS, Store> b);

template<int NBITS, typename Store>
qfloat<NBITS, Store> operator/(qfloat<NBITS, Store> a, qfloat<NBITS, Store> b);

template<int NBITS, typename Store>
bool operator<(qfloat<NBITS, Store> a, qfloat<NBITS, Store> b);

template<int NBITS, typename Store>
bool operator>(qfloat<NBITS, Store> a, qfloat<NBITS, Store> b);

template<int NBITS, typename Store>
bool operator<=(qfloat<NBITS, Store> a, qfloat<NBITS, Store> b);

template<int NBITS, typename Store>
bool operator>=(qfloat<NBITS, Store> a, qfloat<NBITS, Store> b);

template<int NBITS, typename Store>
bool operator==(qfloat<NBITS, Store> a, qfloat<NBITS, Store> b);

template<int NBITS, typename Store>
bool operator!=(qfloat<NBITS, Store> a, qfloat<NBITS, Store> b);

using float_3  = qfloat<3>;
using float_4  = qfloat<4>;
using float_5  = qfloat<5>;
using float_6  = qfloat<6>;
using float_7  = qfloat<7>;
using float_8  = qfloat<8>;
using float_9  = qfloat<9>;
using float_10 = qfloat<10>;

// Batch helpers (Matrix interop)
template<int NBITS>
matrix::Matrix< qfloat<NBITS> > quantize_matrix_to_qfloat(const matrix::Matrix<float>& in,
                                                          float clip = qfloat<NBITS>::current_clip());

template<int NBITS>
matrix::Matrix<float> dequantize_matrix_from_qfloat(const matrix::Matrix< qfloat<NBITS> >& in,
                                                     float clip_unused = qfloat<NBITS>::current_clip());

template<int NBITS>
matrix::Matrix<float> cast_matrix_from_qfloat(const matrix::Matrix< qfloat<NBITS> >& in);

// Explicit instantiations (class template)
extern template class qfloat<2,  int16_t>;
extern template class qfloat<3,  int16_t>;
extern template class qfloat<4,  int16_t>;
extern template class qfloat<5,  int16_t>;
extern template class qfloat<6,  int16_t>;
extern template class qfloat<7,  int16_t>;
extern template class qfloat<8,  int16_t>;
extern template class qfloat<9,  int16_t>;
extern template class qfloat<10, int16_t>;
extern template class qfloat<11, int16_t>;
extern template class qfloat<12, int16_t>;
extern template class qfloat<13, int16_t>;
extern template class qfloat<14, int16_t>;
extern template class qfloat<15, int16_t>;

// Explicit instantiations (function templates)
extern template matrix::Matrix< qfloat<2> > quantize_matrix_to_qfloat<2>(const matrix::Matrix<float>&, float);
extern template matrix::Matrix< qfloat<3> > quantize_matrix_to_qfloat<3>(const matrix::Matrix<float>&, float);
extern template matrix::Matrix< qfloat<4> > quantize_matrix_to_qfloat<4>(const matrix::Matrix<float>&, float);
extern template matrix::Matrix< qfloat<5> > quantize_matrix_to_qfloat<5>(const matrix::Matrix<float>&, float);
extern template matrix::Matrix< qfloat<6> > quantize_matrix_to_qfloat<6>(const matrix::Matrix<float>&, float);
extern template matrix::Matrix< qfloat<7> > quantize_matrix_to_qfloat<7>(const matrix::Matrix<float>&, float);
extern template matrix::Matrix< qfloat<8> > quantize_matrix_to_qfloat<8>(const matrix::Matrix<float>&, float);
extern template matrix::Matrix< qfloat<9> > quantize_matrix_to_qfloat<9>(const matrix::Matrix<float>&, float);
extern template matrix::Matrix< qfloat<10> > quantize_matrix_to_qfloat<10>(const matrix::Matrix<float>&, float);
extern template matrix::Matrix< qfloat<11> > quantize_matrix_to_qfloat<11>(const matrix::Matrix<float>&, float);
extern template matrix::Matrix< qfloat<12> > quantize_matrix_to_qfloat<12>(const matrix::Matrix<float>&, float);
extern template matrix::Matrix< qfloat<13> > quantize_matrix_to_qfloat<13>(const matrix::Matrix<float>&, float);
extern template matrix::Matrix< qfloat<14> > quantize_matrix_to_qfloat<14>(const matrix::Matrix<float>&, float);
extern template matrix::Matrix< qfloat<15> > quantize_matrix_to_qfloat<15>(const matrix::Matrix<float>&, float);

extern template matrix::Matrix<float> dequantize_matrix_from_qfloat<2>(const matrix::Matrix< qfloat<2> >&,
                                                                       float);
extern template matrix::Matrix<float> dequantize_matrix_from_qfloat<3>(const matrix::Matrix< qfloat<3> >&,
                                                                       float);
extern template matrix::Matrix<float> dequantize_matrix_from_qfloat<4>(const matrix::Matrix< qfloat<4> >&,
                                                                       float);
extern template matrix::Matrix<float> dequantize_matrix_from_qfloat<5>(const matrix::Matrix< qfloat<5> >&,
                                                                       float);
extern template matrix::Matrix<float> dequantize_matrix_from_qfloat<6>(const matrix::Matrix< qfloat<6> >&,
                                                                       float);
extern template matrix::Matrix<float> dequantize_matrix_from_qfloat<7>(const matrix::Matrix< qfloat<7> >&,
                                                                       float);
extern template matrix::Matrix<float> dequantize_matrix_from_qfloat<8>(const matrix::Matrix< qfloat<8> >&,
                                                                       float);
extern template matrix::Matrix<float> dequantize_matrix_from_qfloat<9>(const matrix::Matrix< qfloat<9> >&,
                                                                       float);
extern template matrix::Matrix<float> dequantize_matrix_from_qfloat<10>(const matrix::Matrix< qfloat<10> >&,
                                                                        float);
extern template matrix::Matrix<float> dequantize_matrix_from_qfloat<11>(const matrix::Matrix< qfloat<11> >&,
                                                                        float);
extern template matrix::Matrix<float> dequantize_matrix_from_qfloat<12>(const matrix::Matrix< qfloat<12> >&,
                                                                        float);
extern template matrix::Matrix<float> dequantize_matrix_from_qfloat<13>(const matrix::Matrix< qfloat<13> >&,
                                                                        float);
extern template matrix::Matrix<float> dequantize_matrix_from_qfloat<14>(const matrix::Matrix< qfloat<14> >&,
                                                                        float);
extern template matrix::Matrix<float> dequantize_matrix_from_qfloat<15>(const matrix::Matrix< qfloat<15> >&,
                                                                        float);

extern template matrix::Matrix<float> cast_matrix_from_qfloat<2>(const matrix::Matrix< qfloat<2> >&);
extern template matrix::Matrix<float> cast_matrix_from_qfloat<3>(const matrix::Matrix< qfloat<3> >&);
extern template matrix::Matrix<float> cast_matrix_from_qfloat<4>(const matrix::Matrix< qfloat<4> >&);
extern template matrix::Matrix<float> cast_matrix_from_qfloat<5>(const matrix::Matrix< qfloat<5> >&);
extern template matrix::Matrix<float> cast_matrix_from_qfloat<6>(const matrix::Matrix< qfloat<6> >&);
extern template matrix::Matrix<float> cast_matrix_from_qfloat<7>(const matrix::Matrix< qfloat<7> >&);
extern template matrix::Matrix<float> cast_matrix_from_qfloat<8>(const matrix::Matrix< qfloat<8> >&);
extern template matrix::Matrix<float> cast_matrix_from_qfloat<9>(const matrix::Matrix< qfloat<9> >&);
extern template matrix::Matrix<float> cast_matrix_from_qfloat<10>(const matrix::Matrix< qfloat<10> >&);
extern template matrix::Matrix<float> cast_matrix_from_qfloat<11>(const matrix::Matrix< qfloat<11> >&);
extern template matrix::Matrix<float> cast_matrix_from_qfloat<12>(const matrix::Matrix< qfloat<12> >&);
extern template matrix::Matrix<float> cast_matrix_from_qfloat<13>(const matrix::Matrix< qfloat<13> >&);
extern template matrix::Matrix<float> cast_matrix_from_qfloat<14>(const matrix::Matrix< qfloat<14> >&);
extern template matrix::Matrix<float> cast_matrix_from_qfloat<15>(const matrix::Matrix< qfloat<15> >&);

#define DECLARE_QFLOAT_OPS(N) \
extern template qfloat<N> operator+(qfloat<N>, qfloat<N>); \
extern template qfloat<N> operator-(qfloat<N>, qfloat<N>); \
extern template qfloat<N> operator*(qfloat<N>, float); \
extern template qfloat<N> operator*(float, qfloat<N>); \
extern template qfloat<N> operator/(qfloat<N>, float); \
extern template qfloat<N> operator*(qfloat<N>, qfloat<N>); \
extern template qfloat<N> operator/(qfloat<N>, qfloat<N>); \
extern template bool operator<(qfloat<N>, qfloat<N>); \
extern template bool operator>(qfloat<N>, qfloat<N>); \
extern template bool operator<=(qfloat<N>, qfloat<N>); \
extern template bool operator>=(qfloat<N>, qfloat<N>); \
extern template bool operator==(qfloat<N>, qfloat<N>); \
extern template bool operator!=(qfloat<N>, qfloat<N>);

DECLARE_QFLOAT_OPS(2)
DECLARE_QFLOAT_OPS(3)
DECLARE_QFLOAT_OPS(4)
DECLARE_QFLOAT_OPS(5)
DECLARE_QFLOAT_OPS(6)
DECLARE_QFLOAT_OPS(7)
DECLARE_QFLOAT_OPS(8)
DECLARE_QFLOAT_OPS(9)
DECLARE_QFLOAT_OPS(10)
DECLARE_QFLOAT_OPS(11)
DECLARE_QFLOAT_OPS(12)
DECLARE_QFLOAT_OPS(13)
DECLARE_QFLOAT_OPS(14)
DECLARE_QFLOAT_OPS(15)

#undef DECLARE_QFLOAT_OPS



template<typename InputIt>
float compute_clip_from_ratio(InputIt first, InputIt last, float ratio);

extern template float compute_clip_from_ratio<std::vector<float>::iterator>(
    std::vector<float>::iterator, std::vector<float>::iterator, float);
extern template float compute_clip_from_ratio<std::vector<float>::const_iterator>(
    std::vector<float>::const_iterator, std::vector<float>::const_iterator, float);

} // namespace qfloat
