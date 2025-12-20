#include "newcode/common/qfloat/qfloat.hpp"

#include "qfloat.ipp"

namespace qfloat {

// Explicit instantiation for NBITS 2..15 (Store = int16_t)
template class qfloat<2,  int16_t>;
template class qfloat<3,  int16_t>;
template class qfloat<4,  int16_t>;
template class qfloat<5,  int16_t>;
template class qfloat<6,  int16_t>;
template class qfloat<7,  int16_t>;
template class qfloat<8,  int16_t>;
template class qfloat<9,  int16_t>;
template class qfloat<10, int16_t>;
template class qfloat<11, int16_t>;
template class qfloat<12, int16_t>;
template class qfloat<13, int16_t>;
template class qfloat<14, int16_t>;
template class qfloat<15, int16_t>;

// Explicit instantiation for matrix helpers
template matrix::Matrix< qfloat<2> > quantize_matrix_to_qfloat<2>(const matrix::Matrix<float>&, float);
template matrix::Matrix< qfloat<3> > quantize_matrix_to_qfloat<3>(const matrix::Matrix<float>&, float);
template matrix::Matrix< qfloat<4> > quantize_matrix_to_qfloat<4>(const matrix::Matrix<float>&, float);
template matrix::Matrix< qfloat<5> > quantize_matrix_to_qfloat<5>(const matrix::Matrix<float>&, float);
template matrix::Matrix< qfloat<6> > quantize_matrix_to_qfloat<6>(const matrix::Matrix<float>&, float);
template matrix::Matrix< qfloat<7> > quantize_matrix_to_qfloat<7>(const matrix::Matrix<float>&, float);
template matrix::Matrix< qfloat<8> > quantize_matrix_to_qfloat<8>(const matrix::Matrix<float>&, float);
template matrix::Matrix< qfloat<9> > quantize_matrix_to_qfloat<9>(const matrix::Matrix<float>&, float);
template matrix::Matrix< qfloat<10> > quantize_matrix_to_qfloat<10>(const matrix::Matrix<float>&, float);
template matrix::Matrix< qfloat<11> > quantize_matrix_to_qfloat<11>(const matrix::Matrix<float>&, float);
template matrix::Matrix< qfloat<12> > quantize_matrix_to_qfloat<12>(const matrix::Matrix<float>&, float);
template matrix::Matrix< qfloat<13> > quantize_matrix_to_qfloat<13>(const matrix::Matrix<float>&, float);
template matrix::Matrix< qfloat<14> > quantize_matrix_to_qfloat<14>(const matrix::Matrix<float>&, float);
template matrix::Matrix< qfloat<15> > quantize_matrix_to_qfloat<15>(const matrix::Matrix<float>&, float);

template matrix::Matrix<float> dequantize_matrix_from_qfloat<2>(const matrix::Matrix< qfloat<2> >&, float);
template matrix::Matrix<float> dequantize_matrix_from_qfloat<3>(const matrix::Matrix< qfloat<3> >&, float);
template matrix::Matrix<float> dequantize_matrix_from_qfloat<4>(const matrix::Matrix< qfloat<4> >&, float);
template matrix::Matrix<float> dequantize_matrix_from_qfloat<5>(const matrix::Matrix< qfloat<5> >&, float);
template matrix::Matrix<float> dequantize_matrix_from_qfloat<6>(const matrix::Matrix< qfloat<6> >&, float);
template matrix::Matrix<float> dequantize_matrix_from_qfloat<7>(const matrix::Matrix< qfloat<7> >&, float);
template matrix::Matrix<float> dequantize_matrix_from_qfloat<8>(const matrix::Matrix< qfloat<8> >&, float);
template matrix::Matrix<float> dequantize_matrix_from_qfloat<9>(const matrix::Matrix< qfloat<9> >&, float);
template matrix::Matrix<float> dequantize_matrix_from_qfloat<10>(const matrix::Matrix< qfloat<10> >&, float);
template matrix::Matrix<float> dequantize_matrix_from_qfloat<11>(const matrix::Matrix< qfloat<11> >&, float);
template matrix::Matrix<float> dequantize_matrix_from_qfloat<12>(const matrix::Matrix< qfloat<12> >&, float);
template matrix::Matrix<float> dequantize_matrix_from_qfloat<13>(const matrix::Matrix< qfloat<13> >&, float);
template matrix::Matrix<float> dequantize_matrix_from_qfloat<14>(const matrix::Matrix< qfloat<14> >&, float);
template matrix::Matrix<float> dequantize_matrix_from_qfloat<15>(const matrix::Matrix< qfloat<15> >&, float);

template matrix::Matrix<float> cast_matrix_from_qfloat<2>(const matrix::Matrix< qfloat<2> >&);
template matrix::Matrix<float> cast_matrix_from_qfloat<3>(const matrix::Matrix< qfloat<3> >&);
template matrix::Matrix<float> cast_matrix_from_qfloat<4>(const matrix::Matrix< qfloat<4> >&);
template matrix::Matrix<float> cast_matrix_from_qfloat<5>(const matrix::Matrix< qfloat<5> >&);
template matrix::Matrix<float> cast_matrix_from_qfloat<6>(const matrix::Matrix< qfloat<6> >&);
template matrix::Matrix<float> cast_matrix_from_qfloat<7>(const matrix::Matrix< qfloat<7> >&);
template matrix::Matrix<float> cast_matrix_from_qfloat<8>(const matrix::Matrix< qfloat<8> >&);
template matrix::Matrix<float> cast_matrix_from_qfloat<9>(const matrix::Matrix< qfloat<9> >&);
template matrix::Matrix<float> cast_matrix_from_qfloat<10>(const matrix::Matrix< qfloat<10> >&);
template matrix::Matrix<float> cast_matrix_from_qfloat<11>(const matrix::Matrix< qfloat<11> >&);
template matrix::Matrix<float> cast_matrix_from_qfloat<12>(const matrix::Matrix< qfloat<12> >&);
template matrix::Matrix<float> cast_matrix_from_qfloat<13>(const matrix::Matrix< qfloat<13> >&);
template matrix::Matrix<float> cast_matrix_from_qfloat<14>(const matrix::Matrix< qfloat<14> >&);
template matrix::Matrix<float> cast_matrix_from_qfloat<15>(const matrix::Matrix< qfloat<15> >&);

#define INSTANTIATE_QFLOAT_OPS(N) \
template qfloat<N> operator+(qfloat<N>, qfloat<N>); \
template qfloat<N> operator-(qfloat<N>, qfloat<N>); \
template qfloat<N> operator*(qfloat<N>, float); \
template qfloat<N> operator*(float, qfloat<N>); \
template qfloat<N> operator/(qfloat<N>, float); \
template qfloat<N> operator*(qfloat<N>, qfloat<N>); \
template qfloat<N> operator/(qfloat<N>, qfloat<N>); \
template bool operator<(qfloat<N>, qfloat<N>); \
template bool operator>(qfloat<N>, qfloat<N>); \
template bool operator<=(qfloat<N>, qfloat<N>); \
template bool operator>=(qfloat<N>, qfloat<N>); \
template bool operator==(qfloat<N>, qfloat<N>); \
template bool operator!=(qfloat<N>, qfloat<N>);
INSTANTIATE_QFLOAT_OPS(2)
INSTANTIATE_QFLOAT_OPS(3)
INSTANTIATE_QFLOAT_OPS(4)
INSTANTIATE_QFLOAT_OPS(5)
INSTANTIATE_QFLOAT_OPS(6)
INSTANTIATE_QFLOAT_OPS(7)
INSTANTIATE_QFLOAT_OPS(8)
INSTANTIATE_QFLOAT_OPS(9)
INSTANTIATE_QFLOAT_OPS(10)
INSTANTIATE_QFLOAT_OPS(11)
INSTANTIATE_QFLOAT_OPS(12)
INSTANTIATE_QFLOAT_OPS(13)
INSTANTIATE_QFLOAT_OPS(14)
INSTANTIATE_QFLOAT_OPS(15)

#undef INSTANTIATE_QFLOAT_OPS



template float compute_clip_from_ratio<std::vector<float>::iterator>(
    std::vector<float>::iterator, std::vector<float>::iterator, float);
template float compute_clip_from_ratio<std::vector<float>::const_iterator>(
    std::vector<float>::const_iterator, std::vector<float>::const_iterator, float);

} // namespace newcode
