#include "newcode/params.hpp"
#include "newcode/rx/ofec/chase/decoder_core.hpp"
#include "newcode/rx/ofec/chase/chase256.hpp"
#include "newcode/ofec/earlystop/row_early_stop_process_1.hpp"
#include "newcode/ofec/earlystop/row_early_stop_process_2.hpp"
#include "newcode/ofec_decoder_hard.hpp"

#include <array>
#include <stdexcept>

namespace chase {
namespace {

// Chase解码函数指针类型定义
// LLR*: 线性矩阵输入
// LLR*: 信道矩阵输入
// float*: 输出结果
// const newcode::Params&: 解码参数
template<typename LLR>
using ChaseFn = void (*)(const LLR*, const LLR*, float*, const newcode::Params&);

// Chase解码核心实现函数
// 模板参数LLR: 对数似然比类型，可以是float、int8_t或量化类型qfloat
// 参数：
//   lin_matrix: 线性矩阵输入
//   lch_matrix: 信道矩阵输入
//   use_hard_decode: 是否使用硬解码模式
//   p: 解码参数
//   chase_fn: 具体的Chase解码函数指针
//   early_stop_row_flags: 提前终止行标志
// 返回值：包含解码结果的DecoderCoreResult对象
template<typename LLR>
DecoderCoreResult<LLR> Decoder_Core_impl(const matrix::Matrix<LLR>& lin_matrix,
                                         const matrix::Matrix<LLR>& lch_matrix,
                                         bool use_hard_decode,
                                         const newcode::Params& p,
                                         ChaseFn<LLR> chase_fn,
                                         const std::vector<bool>* early_stop_row_flags,
                                         const std::vector<uint8_t>* mux_state)
{
  const size_t rows = lin_matrix.rows();
  const size_t cols = lin_matrix.cols();
  
  // 计算预期的列数：2 * 子块列数 * 子块维度位数
  const size_t expected_cols =
      2 * newcode::Params::NUM_SUBBLOCK_COLS * newcode::Params::BITS_PER_SUBBLOCK_DIM;
  
  // 验证输入矩阵列数是否符合预期
  if (cols != expected_cols) {
    throw std::invalid_argument("Decoder_Core: unexpected column count");
  }

  // 初始化解码结果：
  // - 输出矩阵大小为rows x cols，元素类型为float
  // - produced_rows向量标记哪些行成功解码
  DecoderCoreResult<LLR> result{
      matrix::Matrix<float>(rows, cols),
      std::vector<bool>(rows, false)};

  // 逐行处理解码
  for (size_t row = 0; row < rows; ++row) {
    // 创建当前行的线性矩阵和信道矩阵向量
    std::array<LLR, newcode::Params::BCH_N> LinVec{};
    std::array<LLR, newcode::Params::BCH_N> LchVec{};

    // 复制当前行的数据到局部数组
    for (size_t col = 0; col < cols; ++col) {
      LinVec[col] = lin_matrix[row][col];
      LchVec[col] = lch_matrix[row][col];
    }

    // 创建输出结果数组
    std::array<float, newcode::Params::BCH_N> Y2{};
    // 标记当前行是否成功解码
    bool produced = false;

    // 配置当前行的解码参数，特别是调试跟踪信息
    newcode::Params row_params = p;
    const auto& trace_cfg = p.debug_trace;
    
    // 清除活动的Chase条目，只保留当前行相关的
    row_params.debug_trace.active_chase_entries.clear();
    row_params.debug_trace.chase_expected_bits = trace_cfg.chase_expected_bits;
    
    // 设置当前行的期望比特信息（如果有）
    if (trace_cfg.chase_expected_bits &&
        static_cast<size_t>(row) < trace_cfg.chase_expected_bits->size()) {
      row_params.debug_trace.chase_expected_bits_row =
          &(*trace_cfg.chase_expected_bits)[row];
    } else {
      row_params.debug_trace.chase_expected_bits_row = nullptr;
    }
    
    // 收集当前行相关的活动Chase条目
    for (const auto& entry : trace_cfg.active_chase_entries) {
      if (entry.row_index == static_cast<int>(row)) {
        row_params.debug_trace.active_chase_entries.push_back(entry);
      }
    }
    
    // 设置Chase解码器的行和列跟踪信息
    if (!row_params.debug_trace.active_chase_entries.empty()) {
      row_params.debug_trace.chase_decoder_row = static_cast<int>(row);
      row_params.debug_trace.chase_decoder_col =
          row_params.debug_trace.active_chase_entries.front().k;
    } else if (trace_cfg.chase_decoder_row >= 0 &&
               static_cast<int>(row) == trace_cfg.chase_decoder_row) {
      row_params.debug_trace.chase_decoder_row = static_cast<int>(row);
      row_params.debug_trace.chase_decoder_col = trace_cfg.chase_decoder_col;
    } else {
      // 禁用当前行的跟踪
      row_params.debug_trace.chase_decoder_row = -1;
      row_params.debug_trace.chase_decoder_col = -1;
    }

    const bool has_mux_state = mux_state &&
                               row < mux_state->size();
    const uint8_t mux_tag =
        has_mux_state ? (*mux_state)[row] : static_cast<uint8_t>(0xFF);

    // 根据配置选择解码模式
    if (has_mux_state && mux_tag == 2u) {
      // 本轮未分配到 SISO：保持未产生输出（不更新）
      produced = false;
    } else if (use_hard_decode) {
      // 使用硬解码
      produced = newcode::perform_hard_decode<LLR>(LinVec, LchVec, Y2, p);
    } else {
      // 使用软解码
      if (has_mux_state) {
        if (mux_tag == 1u) {
          //newcode::row_early_stop_process_2(LinVec.data(),LchVec.data(),Y2.data(),row_params);
          newcode::row_early_stop_process_1(LinVec.data(),
                                            LchVec.data(),
                                            Y2.data(),
                                            row_params);
          produced = true;
        } else {
          chase_fn(LinVec.data(), LchVec.data(), Y2.data(), row_params);
          produced = true;
        }
      } else {
        // 调用具体的Chase解码函数
        chase_fn(LinVec.data(), LchVec.data(), Y2.data(), row_params);
        produced = true;
      }
    }

    // 如果当前行成功解码，将结果保存到输出矩阵
    if (produced) {
      result.produced_rows[row] = true;
      for (size_t col = 0; col < cols; ++col) {
        result.lout[row][col] = Y2[col];
      }
    }
  }

  // 返回解码结果
  return result;
}

} // namespace

// 普通Chase解码器实现
// 使用chase_decode_256_plain作为具体的解码函数
// 模板参数LLR: 对数似然比类型
template<typename LLR>
DecoderCoreResult<LLR> Decoder_Core_plain(const matrix::Matrix<LLR>& lin_matrix,
                                          const matrix::Matrix<LLR>& lch_matrix,
                                          bool use_hard_decode,
                                          const newcode::Params& p,
                                          const std::vector<bool>* early_stop_row_flags,
                                          const std::vector<uint8_t>* mux_state)
{
  return Decoder_Core_impl(lin_matrix, lch_matrix, use_hard_decode, p,
                           &chase_decode_256_plain<LLR>,
                           early_stop_row_flags,
                           mux_state);
}

// ebchPF版本的Chase解码器实现
// 使用chase_decode_256_ebchPF作为具体的解码函数
// 模板参数LLR: 对数似然比类型
template<typename LLR>
DecoderCoreResult<LLR> Decoder_Core_ebchPF(const matrix::Matrix<LLR>& lin_matrix,
                                           const matrix::Matrix<LLR>& lch_matrix,
                                           bool use_hard_decode,
                                           const newcode::Params& p,
                                           const std::vector<bool>* early_stop_row_flags,
                                           const std::vector<uint8_t>* mux_state)
{
  return Decoder_Core_impl(lin_matrix, lch_matrix, use_hard_decode, p,
                           &chase_decode_256_ebchPF<LLR>,
                           early_stop_row_flags,
                           mux_state);
}

// 模板实例化：针对float类型的普通Chase解码器
template DecoderCoreResult<float> Decoder_Core_plain<float>(const matrix::Matrix<float>&,
                                                            const matrix::Matrix<float>&,
                                                            bool,
                                                            const newcode::Params&,
                                                            const std::vector<bool>*,
                                                            const std::vector<uint8_t>*);

// 模板实例化：针对int8_t类型的普通Chase解码器
template DecoderCoreResult<int8_t> Decoder_Core_plain<int8_t>(const matrix::Matrix<int8_t>&,
                                                              const matrix::Matrix<int8_t>&,
                                                              bool,
                                                              const newcode::Params&,
                                                              const std::vector<bool>*,
                                                              const std::vector<uint8_t>*);

// 模板实例化：针对float类型的ebchPF Chase解码器
template DecoderCoreResult<float> Decoder_Core_ebchPF<float>(const matrix::Matrix<float>&,
                                                             const matrix::Matrix<float>&,
                                                             bool,
                                                             const newcode::Params&,
                                                             const std::vector<bool>*,
                                                             const std::vector<uint8_t>*);

// 模板实例化：针对int8_t类型的ebchPF Chase解码器
template DecoderCoreResult<int8_t> Decoder_Core_ebchPF<int8_t>(const matrix::Matrix<int8_t>&,
                                                               const matrix::Matrix<int8_t>&,
                                                               bool,
                                                               const newcode::Params&,
                                                               const std::vector<bool>*,
                                                               const std::vector<uint8_t>*);

// 定义量化浮点数类型的模板实例化宏
// N: 量化位数
#define INSTANTIATE_DECODER_CORE_QFLOAT(N) \
template DecoderCoreResult<qfloat::qfloat<N>> Decoder_Core_plain<qfloat::qfloat<N>>( \
    const matrix::Matrix<qfloat::qfloat<N>>& lin_matrix, const matrix::Matrix<qfloat::qfloat<N>>& lch_matrix, bool, const newcode::Params&, \
    const std::vector<bool>*, const std::vector<uint8_t>*); \
template DecoderCoreResult<qfloat::qfloat<N>> Decoder_Core_ebchPF<qfloat::qfloat<N>>( \
    const matrix::Matrix<qfloat::qfloat<N>>& lin_matrix, const matrix::Matrix<qfloat::qfloat<N>>& lch_matrix, bool, const newcode::Params&, \
    const std::vector<bool>*, const std::vector<uint8_t>*);

// 为不同量化位数（2-15位）的qfloat类型实例化解码器
// 这些实例化允许解码器处理各种精度的量化输入
// 其中6位量化（INSTANTIATE_DECODER_CORE_QFLOAT(6)）已被证明在性能上优于浮点实现
INSTANTIATE_DECODER_CORE_QFLOAT(2)
INSTANTIATE_DECODER_CORE_QFLOAT(3)
INSTANTIATE_DECODER_CORE_QFLOAT(4)
INSTANTIATE_DECODER_CORE_QFLOAT(5)
INSTANTIATE_DECODER_CORE_QFLOAT(6)  // 6位量化，性能优于浮点
INSTANTIATE_DECODER_CORE_QFLOAT(7)
INSTANTIATE_DECODER_CORE_QFLOAT(8)
INSTANTIATE_DECODER_CORE_QFLOAT(9)
INSTANTIATE_DECODER_CORE_QFLOAT(10)
INSTANTIATE_DECODER_CORE_QFLOAT(11)
INSTANTIATE_DECODER_CORE_QFLOAT(12)
INSTANTIATE_DECODER_CORE_QFLOAT(13)
INSTANTIATE_DECODER_CORE_QFLOAT(14)
INSTANTIATE_DECODER_CORE_QFLOAT(15)

// 取消宏定义
#undef INSTANTIATE_DECODER_CORE_QFLOAT

} // namespace chase
