#include "newcode/qfloat.hpp"

namespace newcode {

// 为 2..15 位做显式实例化（Store=int16_t），可按需删减
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

} // namespace newcode
