#pragma once
#include <memory>
#include <string>
#include "newcode/ofec_llr_matrix.hpp"

namespace newcode {

struct DecodeStats { int iters = 0; bool success = false; };

struct IDecoder {
  virtual ~IDecoder() = default;
  virtual DecodeStats decode(const Matrix<float>& lin,
                             const Matrix<float>& lch,
                             Matrix<float>& lout) = 0;
};

// 由各变体各自“部分实现”的工厂；某实现不认识的 name 返回 nullptr
std::unique_ptr<IDecoder> make_decoder(const std::string& name);

} // namespace newcode
