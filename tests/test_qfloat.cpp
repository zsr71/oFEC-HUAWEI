#include "newcode/qfloat.hpp"
#include "newcode/matrix.hpp"

#include <algorithm>
#include <cassert>
#include <cmath>
#include <iostream>
#include <iomanip>

namespace {

template <typename QFloat>
int saturate_code(int value) {
  if (value > QFloat::HI()) return QFloat::HI();
  if (value < QFloat::LO()) return QFloat::LO();
  return value;
}

template <int NBITS>
void check_round_trip(float clip, float value) {
  using Q = newcode::qfloat<NBITS>;
  Q q = Q::from_float(value, clip);
  const float clamped = std::clamp(value, -clip, clip);
  const float step = clip / static_cast<float>(Q::Q());
  const float restored = q.to_float(clip);
  std::cout << "[round-trip qfloat<" << NBITS << ">] input=" << value
            << " clamped=" << clamped
            << " code=" << q.code()
            << " restored=" << restored
            << " step=" << step << "\n";
  assert(std::fabs(restored - clamped) <= step);
}

void test_scalar_quantization() {
  constexpr float clip = newcode::qfloat<5>::DEFAULT_CLIP;
  check_round_trip<5>(clip, 0.0f);
  check_round_trip<5>(clip, 2.5f);
  check_round_trip<5>(clip, -3.75f);
  auto hi = newcode::qfloat<5>::from_float(clip * 10.0f);
  auto lo = newcode::qfloat<5>::from_float(-clip * 10.0f);
  std::cout << "[saturation] hi input=" << clip * 10.0f
            << " -> code=" << hi.code() << " (HI=" << newcode::qfloat<5>::HI() << ")\n";
  std::cout << "[saturation] lo input=" << -clip * 10.0f
            << " -> code=" << lo.code() << " (LO=" << newcode::qfloat<5>::LO() << ")\n";
  assert(hi.code() == newcode::qfloat<5>::HI());
  assert(lo.code() == newcode::qfloat<5>::LO());
}

void test_unary_and_sign() {
  using Q = newcode::qfloat<5>;
  auto pos = Q::from_float(3.2f);
  auto neg = -pos;
  auto zero = Q::from_float(0.0f);
  std::cout << "[unary] pos=" << pos.code()
            << " neg=" << neg.code()
            << " zero=" << zero.code() << "\n";
  assert(neg.code() == -pos.code());
  assert(pos > zero);
  assert(neg < zero);
}

void test_add_sub_ops() {
  using Q = newcode::qfloat<4>;
  auto a = Q::from_float(2.0f);
  auto b = Q::from_float(-0.5f);
  auto sum = a + b;
  auto diff = a - b;
  std::cout << "[add/sub] a=" << a.code()
            << " b=" << b.code()
            << " sum=" << sum.code()
            << " diff=" << diff.code() << "\n";
  assert(sum.code() == saturate_code<Q>(a.code() + b.code()));
  assert(diff.code() == saturate_code<Q>(a.code() - b.code()));

  auto sat = Q::from_float(Q::DEFAULT_CLIP);
  auto sat_sum = sat + sat;
  std::cout << "[add/sub saturation] sat=" << sat.code()
            << " sat+sat=" << sat_sum.code() << "\n";
  assert(sat_sum.code() == Q::HI());
}

void test_float_scale_mul_div() {
  using Q = newcode::qfloat<6>;
  auto base = Q::from_float(1.5f);
  auto scaled = base * 2.5f;
  int expected = saturate_code<Q>(
      static_cast<int>(std::lrint(static_cast<float>(base.code()) * 2.5f)));
  std::cout << "[scale] base_code=" << base.code()
            << " *2.5 -> code=" << scaled.code()
            << " expected=" << expected << "\n";
  assert(scaled.code() == expected);

  auto rhs = Q::from_float(0.75f);
  auto prod = base * rhs;
  auto quot = base / rhs;
  const float tolerance = Q::DEFAULT_CLIP / static_cast<float>(Q::Q());
  const float prod_expected = base.to_float() * rhs.to_float();
  const float quot_expected = base.to_float() / rhs.to_float();
  std::cout << std::fixed << std::setprecision(6);
  std::cout << "[mul] " << base.to_float() << " * " << rhs.to_float()
            << " -> " << prod.to_float()
            << " (expected " << prod_expected << ", tol " << tolerance << ")\n";
  std::cout << "[div] " << base.to_float() << " / " << rhs.to_float()
            << " -> " << quot.to_float()
            << " (expected " << quot_expected << ", tol " << tolerance << ")\n";
  std::cout.unsetf(std::ios::floatfield);
  assert(std::fabs(prod.to_float() - prod_expected) <= tolerance);
  assert(std::fabs(quot.to_float() - quot_expected) <= tolerance);
}

void test_matrix_helpers() {
  using Q = newcode::qfloat<4>;
  const float clip = 4.0f;
  newcode::Matrix<float> mat(2, 3);
  mat[0][0] = -6.0f;
  mat[0][1] = -1.0f;
  mat[0][2] = 0.5f;
  mat[1][0] = 1.0f;
  mat[1][1] = 2.0f;
  mat[1][2] = 6.0f;

  auto qmat = newcode::quantize_matrix_to_qfloat<4>(mat, clip);
  std::cout << "[matrix quant]\n";
  for (size_t r = 0; r < mat.rows(); ++r) {
    for (size_t c = 0; c < mat.cols(); ++c) {
      std::cout << "  (" << r << "," << c << ") value=" << mat[r][c]
                << " code=" << qmat[r][c].code() << "\n";
    }
  }
  assert(qmat[0][0].code() == Q::LO());
  assert(qmat[1][2].code() == Q::HI());

  auto deq = newcode::dequantize_matrix_from_qfloat<4>(qmat, clip);
  auto casted = newcode::cast_matrix_from_qfloat(qmat);
  const float step = clip / static_cast<float>(Q::Q());

  std::cout << "[matrix dequant/cast]\n";
  for (size_t r = 0; r < mat.rows(); ++r) {
    for (size_t c = 0; c < mat.cols(); ++c) {
      const float clamped = std::clamp(mat[r][c], -clip, clip);
      std::cout << "  (" << r << "," << c << ") dequant=" << deq[r][c]
                << " cast=" << casted[r][c]
                << " clamped=" << clamped
                << " step=" << step << "\n";
      assert(std::fabs(deq[r][c] - clamped) <= step);
      assert(std::fabs(casted[r][c] - static_cast<float>(qmat[r][c].code())) <
             1e-6f);
    }
  }
}

} // namespace

int main() {
  test_scalar_quantization();
  test_unary_and_sign();
  test_add_sub_ops();
  test_float_scale_mul_div();
  test_matrix_helpers();
  std::cout << "All qfloat quantization tests passed!\n";
  return 0;
}
