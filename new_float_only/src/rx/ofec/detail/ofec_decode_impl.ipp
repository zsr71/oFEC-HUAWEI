#pragma once

#include "ofec_window_impl.ipp"

#include <fstream>
#include <string>
#include <vector>

namespace new_float_only {
namespace detail {

inline void dump_matrix_text(const matrix::Matrix<float>& matrix_in,
                             const std::string& path) {
  if (path.empty()) {
    return;
  }
  std::ofstream out(path);
  for (std::size_t r = 0; r < matrix_in.rows(); ++r) {
    for (std::size_t c = 0; c < matrix_in.cols(); ++c) {
      if (c) {
        out << ' ';
      }
      out << matrix_in[r][c];
    }
    out << '\n';
  }
}

/**
 * 顶层滑窗解码实现。
 * 关键语句：
 * 1. work_llr 作为窗口内可被反复覆盖的工作矩阵，初值全 0；
 * 2. last_tile_history_llr 单独累计“最终用于 history 的值”；
 * 3. 最终输出继续保留旧口径：post = channel_llr + last_tile_history_llr。
 */
matrix::Matrix<float> decode_plain_llr_impl(const matrix::Matrix<float>& llr_mat,
                                            const new_float_only::Params& p,
                                            bool normalize_extrinsic,
                                            const matrix::Matrix<float>* tx_llr_ref,
                                            CoreFn core_fn)
{
  const size_t N = new_float_only::Params::NUM_SUBBLOCK_COLS * new_float_only::Params::BITS_PER_SUBBLOCK_DIM;

  const size_t RROWS = llr_mat.rows();
  const size_t CCOLS = llr_mat.cols();
  if (CCOLS != N)
      throw std::invalid_argument("decode_plain_llr: llr_mat cols != N.");

  assert(p.valid());
  const size_t TILE_HEIGHT_ROWS = p.tile_height_rows();
  const size_t TILE_STRIDE_ROWS = p.tile_stride_rows();
  const size_t WIN_HEIGHT_ROWS  = p.win_height_rows();
  const size_t POP_PUSH_ROWS    = p.pop_push_rows();
  const size_t TILES_PER_WIN    = p.TILES_PER_WIN;


  matrix::Matrix<float> channel_llr = llr_mat;
  matrix::Matrix<float> work_llr(RROWS, N);
  for (size_t r = 0; r < RROWS; ++r)
      for (size_t c = 0; c < N; ++c)
          work_llr[r][c] = 0.0f;
  matrix::Matrix<float> last_tile_history_llr(RROWS, N);

  if (RROWS < WIN_HEIGHT_ROWS) {
    return channel_llr;
  }

  size_t win_start     = p.initial_win_start_rows();
  const size_t last_ws = RROWS - WIN_HEIGHT_ROWS;

  while (win_start <= last_ws) {
    const size_t win_end = win_start + WIN_HEIGHT_ROWS - 1;
    process_window_impl(work_llr, channel_llr,
                        win_start, win_end, p,
                        TILE_HEIGHT_ROWS, TILE_STRIDE_ROWS, TILES_PER_WIN,
                        normalize_extrinsic,
                        tx_llr_ref,
                        core_fn,
                        &last_tile_history_llr);

    win_start += POP_PUSH_ROWS;
  }

  matrix::Matrix<float> out(RROWS, N);
  for (size_t r = 0; r < RROWS; ++r) {
    for (size_t c = 0; c < N; ++c) {
      out[r][c] = channel_llr[r][c] + last_tile_history_llr[r][c];
    }
  }

  if (p.DUMP_WORK_LLR) {
    dump_matrix_text(work_llr, p.WORK_LLR_OUTPUT_PATH);
  }

  return out;
}

} // namespace detail
} // namespace new_float_only
