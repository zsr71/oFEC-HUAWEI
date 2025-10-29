#include "newcode/decoder_api.hpp"
namespace newcode {
struct EbchPFDecoder : IDecoder {
  DecodeStats decode(const Matrix<float>& lin,const Matrix<float>& lch,Matrix<float>& lout) override {
    // TODO: 调用 ebchPF 变体实际译码流程
    return {0, true};
  }
};
std::unique_ptr<IDecoder> make_decoder(const std::string& name);
std::unique_ptr<IDecoder> make_decoder(const std::string& name) {
  if (name == "ebchPF") return std::make_unique<EbchPFDecoder>();
  return nullptr;
}
} // namespace newcode
