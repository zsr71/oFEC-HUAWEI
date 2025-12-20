#pragma once

#include <string>

#include "newcode/decoder_api.hpp"

namespace newcode {

// 将矩阵展开写入文本文件（便于 MATLAB readmatrix + histogram）
// 返回保存的路径，若未保存则返回空字符串。
// enable 控制是否执行，path_override 为空时使用默认命名（带 suffix）。
std::string dump_quantized_llr(const matrix::Matrix<float>& llr_mat,
                               const DecodeRequest& request,
                               bool enable,
                               const std::string& path_override,
                               const std::string& default_suffix);

} // namespace newcode
