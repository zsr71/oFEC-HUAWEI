#pragma once
#include <memory>
#include <string>
#include "newcode/ofec_llr_matrix.hpp"

namespace newcode {

struct IInterleaver {
  virtual ~IInterleaver() = default;
  virtual void interleave  (const Matrix<float>& in, Matrix<float>& out) = 0;
  virtual void deinterleave(const Matrix<float>& in, Matrix<float>& out) = 0;
};

std::unique_ptr<IInterleaver> make_interleaver(const std::string& name);

} // namespace newcode
