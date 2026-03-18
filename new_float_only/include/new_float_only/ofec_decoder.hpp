#pragma once

#include <cstddef>

#include "new_float_only/common/matrix/matrix.hpp"
#include "new_float_only/params.hpp"

namespace new_float_only {

/**
 * 单个 Tile 处理后的返回结果。
 * tile_out 是当前 Tile 范围内的新外信息矩阵，尺寸与输入 Tile 完全一致。
 */
struct TileProcessResult {
  matrix::Matrix<float> tile_out;
};

/**
 * 顶层浮点 plain 解码入口。
 * 输入为去交织后的信道 LLR 矩阵，输出为与输入同尺寸的译码后 LLR 矩阵。
 * 若提供 tx_llr_ref，则会在调试模式下为 Chase 跟踪补充期望比特信息。
 */
matrix::Matrix<float> decode_plain_llr(const matrix::Matrix<float>& channel_llr,
                                       const Params& config,
                                       bool normalize_extrinsic = true,
                                       const matrix::Matrix<float>* tx_llr_ref = nullptr);

/**
 * 处理一个窗口内的全部 Tile。
 * 该函数会原地更新 work_llr，并可选地记录最后一个 soft Tile 的 history。
 */
void process_window_plain(matrix::Matrix<float>& work_llr,
                          const matrix::Matrix<float>& channel_llr,
                          std::size_t win_start,
                          std::size_t win_end,
                          const Params& config,
                          std::size_t tile_height_rows,
                          std::size_t tile_stride_rows,
                          std::size_t tiles_per_window,
                          bool normalize_extrinsic = true,
                          const matrix::Matrix<float>* tx_llr_ref = nullptr,
                          matrix::Matrix<float>* last_tile_history_accum = nullptr);

/**
 * 处理单个 Tile。
 * 输入 tile_in/ch_tile 为当前 Tile 的历史外信息与信道 LLR，输出为新的 Tile 外信息。
 */
TileProcessResult process_tile_plain(const matrix::Matrix<float>& tile_in,
                                     const matrix::Matrix<float>& ch_tile,
                                     const Params& config,
                                     std::size_t tile_top_row_global,
                                     bool use_hard_decode,
                                     bool normalize_extrinsic,
                                     const matrix::Matrix<float>* tx_llr_ref = nullptr,
                                     matrix::Matrix<float>* last_tile_history_accum = nullptr,
                                     bool capture_last_tile_history = false);

}  // namespace new_float_only
