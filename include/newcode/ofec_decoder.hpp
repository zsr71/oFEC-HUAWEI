#pragma once
#include <cstddef>
#include <cstdint>
#include <vector>
#include "newcode/matrix.hpp"
#include "newcode/params.hpp"
#include "newcode/qfloat.hpp"

namespace newcode {

struct TileEarlyStopCounter {
  std::size_t triggered = 0;
  std::size_t total = 0;
  std::size_t row_triggered = 0;
  std::size_t row_total = 0;
};

template <typename LLR>
struct TileProcessResult {
  Matrix<LLR> tile_out;
  bool early_stop_triggered = false;
  std::size_t rows_early_stop = 0;
  std::size_t rows_total = 0;
};

// 顶层解码（plain / ebchPF 两个变体分别导出）
template <typename LLR>
Matrix<LLR> ofec_decode_llr_plain(const Matrix<LLR>& llr_mat, const Params& p,
                                  std::vector<TileEarlyStopCounter>* tile_stats = nullptr,
                                  bool normalize_extrinsic = true,
                                  const Matrix<float>* tx_llr_ref = nullptr);

template <typename LLR>
Matrix<LLR> ofec_decode_llr_ebchPF(const Matrix<LLR>& llr_mat, const Params& p,
                                   std::vector<TileEarlyStopCounter>* tile_stats = nullptr,
                                   bool normalize_extrinsic = true,
                                   const Matrix<float>* tx_llr_ref = nullptr);

// 窗口处理：对 work_llr 的 [win_start, win_end] 做一轮（或多轮）tile 扫描（就地修改）
template <typename LLR>
void process_window_plain(Matrix<LLR>& work_llr,
                          const Matrix<LLR>& channel_llr,
                          std::size_t win_start, std::size_t win_end, const Params& p,
                          std::size_t tile_height_rows, std::size_t tile_stride_rows, std::size_t TILES_PER_WIN,
                          std::vector<TileEarlyStopCounter>* tile_stats = nullptr,
                          bool normalize_extrinsic = true,
                          const Matrix<float>* tx_llr_ref = nullptr);

template <typename LLR>
void process_window_ebchPF(Matrix<LLR>& work_llr,
                           const Matrix<LLR>& channel_llr,
                           std::size_t win_start, std::size_t win_end, const Params& p,
                           std::size_t tile_height_rows, std::size_t tile_stride_rows, std::size_t TILES_PER_WIN,
                           std::vector<TileEarlyStopCounter>* tile_stats = nullptr,
                           bool normalize_extrinsic = true,
                           const Matrix<float>* tx_llr_ref = nullptr);

template <typename LLR>
inline void process_window(Matrix<LLR>& work_llr,
                           const Matrix<LLR>& channel_llr,
                           std::size_t win_start, std::size_t win_end, const Params& p,
                           std::size_t tile_height_rows, std::size_t tile_stride_rows, std::size_t TILES_PER_WIN,
                           std::vector<TileEarlyStopCounter>* tile_stats = nullptr,
                           bool normalize_extrinsic = true,
                           const Matrix<float>* tx_llr_ref = nullptr)
{
  process_window_plain(work_llr, channel_llr, win_start, win_end, p,
                       tile_height_rows, tile_stride_rows, TILES_PER_WIN,
                       tile_stats, normalize_extrinsic, tx_llr_ref);
}

// Tile 处理：对单个 tile 的 LLR 做一次“Chase 解码/更新”
template <typename LLR>
TileProcessResult<LLR> process_tile_plain(const Matrix<LLR>& tile_in,
                                          const Matrix<LLR>& ch_tile,
                                          const Params& p,
                                          std::size_t tile_top_row_global,
                                          bool use_hard_decode,
                                          bool normalize_extrinsic,
                                          const Matrix<float>* tx_llr_ref = nullptr);

template <typename LLR>
TileProcessResult<LLR> process_tile_ebchPF(const Matrix<LLR>& tile_in,
                                          const Matrix<LLR>& ch_tile,
                                          const Params& p,
                                          std::size_t tile_top_row_global,
                                          bool use_hard_decode,
                                          bool normalize_extrinsic,
                                          const Matrix<float>* tx_llr_ref = nullptr);

template <typename LLR>
inline TileProcessResult<LLR> process_tile(const Matrix<LLR>& tile_in,
                                           const Matrix<LLR>& ch_tile,
                                           const Params& p,
                                           std::size_t tile_top_row_global,
                                           bool use_hard_decode,
                                           bool normalize_extrinsic,
                                           const Matrix<float>* tx_llr_ref = nullptr)
{
  return process_tile_plain(tile_in, ch_tile, p, tile_top_row_global,
                            use_hard_decode, normalize_extrinsic, tx_llr_ref);
}

// ===== extern template（减少重复实例化）=====
// 基础类型
extern template Matrix<float>  ofec_decode_llr_plain<float >(const Matrix<float>&,  const Params&, std::vector<TileEarlyStopCounter>*, bool, const Matrix<float>*);
extern template Matrix<int8_t> ofec_decode_llr_plain<int8_t>(const Matrix<int8_t>&, const Params&, std::vector<TileEarlyStopCounter>*, bool, const Matrix<float>*);

extern template Matrix<float>  ofec_decode_llr_ebchPF<float >(const Matrix<float>&,  const Params&, std::vector<TileEarlyStopCounter>*, bool, const Matrix<float>*);
extern template Matrix<int8_t> ofec_decode_llr_ebchPF<int8_t>(const Matrix<int8_t>&, const Params&, std::vector<TileEarlyStopCounter>*, bool, const Matrix<float>*);

extern template void process_window_plain<float >(Matrix<float>&,  const Matrix<float>&,
                                                  std::size_t, std::size_t, const Params&,
                                                  std::size_t, std::size_t, std::size_t,
                                                  std::vector<TileEarlyStopCounter>*, bool,
                                                  const Matrix<float>*);
extern template void process_window_plain<int8_t>(Matrix<int8_t>&, const Matrix<int8_t>&,
                                                  std::size_t, std::size_t, const Params&,
                                                  std::size_t, std::size_t, std::size_t,
                                                  std::vector<TileEarlyStopCounter>*, bool,
                                                  const Matrix<float>*);
extern template void process_window_ebchPF<float >(Matrix<float>&,  const Matrix<float>&,
                                                  std::size_t, std::size_t, const Params&,
                                                  std::size_t, std::size_t, std::size_t,
                                                  std::vector<TileEarlyStopCounter>*, bool,
                                                  const Matrix<float>*);
extern template void process_window_ebchPF<int8_t>(Matrix<int8_t>&, const Matrix<int8_t>&,
                                                  std::size_t, std::size_t, const Params&,
                                                  std::size_t, std::size_t, std::size_t,
                                                  std::vector<TileEarlyStopCounter>*, bool,
                                                  const Matrix<float>*);

extern template TileProcessResult<float>  process_tile_plain<float >(const Matrix<float>&,  const Matrix<float>&,
                                                                     const Params&, std::size_t, bool, bool,
                                                                     const Matrix<float>*);
extern template TileProcessResult<int8_t> process_tile_plain<int8_t>(const Matrix<int8_t>&, const Matrix<int8_t>&,
                                                                     const Params&, std::size_t, bool, bool,
                                                                     const Matrix<float>*);
extern template TileProcessResult<float>  process_tile_ebchPF<float >(const Matrix<float>&,  const Matrix<float>&,
                                                                      const Params&, std::size_t, bool, bool,
                                                                      const Matrix<float>*);
extern template TileProcessResult<int8_t> process_tile_ebchPF<int8_t>(const Matrix<int8_t>&, const Matrix<int8_t>&,
                                                                      const Params&, std::size_t, bool, bool,
                                                                      const Matrix<float>*);

// qfloat 量化类型（按需开启）
extern template Matrix<newcode::qfloat<4>> ofec_decode_llr_plain<newcode::qfloat<4>>(const Matrix<newcode::qfloat<4>>&,
                                                                                    const Params&, std::vector<TileEarlyStopCounter>*, bool,
                                                                                    const Matrix<float>*);
extern template Matrix<newcode::qfloat<5>> ofec_decode_llr_plain<newcode::qfloat<5>>(const Matrix<newcode::qfloat<5>>&, 
                                                                                    const Params&, std::vector<TileEarlyStopCounter>*, bool,
                                                                                    const Matrix<float>*);

extern template Matrix<newcode::qfloat<4>> ofec_decode_llr_ebchPF<newcode::qfloat<4>>(const Matrix<newcode::qfloat<4>>&,
                                                                                      const Params&, std::vector<TileEarlyStopCounter>*, bool,
                                                                                      const Matrix<float>*);
extern template Matrix<newcode::qfloat<5>> ofec_decode_llr_ebchPF<newcode::qfloat<5>>(const Matrix<newcode::qfloat<5>>&, 
                                                                                      const Params&, std::vector<TileEarlyStopCounter>*, bool,
                                                                                      const Matrix<float>*);

extern template void process_window<newcode::qfloat<4>>(Matrix<newcode::qfloat<4>>&,
                                                        const Matrix<newcode::qfloat<4>>&,
                                                        std::size_t, std::size_t, const Params&,
                                                        std::size_t, std::size_t, std::size_t,
                                                        std::vector<TileEarlyStopCounter>*, bool,
                                                        const Matrix<float>*);
extern template void process_window<newcode::qfloat<5>>(Matrix<newcode::qfloat<5>>&,
                                                        const Matrix<newcode::qfloat<5>>&,
                                                        std::size_t, std::size_t, const Params&,
                                                        std::size_t, std::size_t, std::size_t,
                                                        std::vector<TileEarlyStopCounter>*, bool,
                                                        const Matrix<float>*);
extern template void process_window_plain<newcode::qfloat<4>>(Matrix<newcode::qfloat<4>>&,
                                                              const Matrix<newcode::qfloat<4>>&,
                                                              std::size_t, std::size_t, const Params&,
                                                              std::size_t, std::size_t, std::size_t,
                                                              std::vector<TileEarlyStopCounter>*, bool,
                                                              const Matrix<float>*);
extern template void process_window_plain<newcode::qfloat<5>>(Matrix<newcode::qfloat<5>>&,
                                                              const Matrix<newcode::qfloat<5>>&,
                                                              std::size_t, std::size_t, const Params&,
                                                              std::size_t, std::size_t, std::size_t,
                                                              std::vector<TileEarlyStopCounter>*, bool,
                                                              const Matrix<float>*);
extern template void process_window_ebchPF<newcode::qfloat<4>>(Matrix<newcode::qfloat<4>>&,
                                                              const Matrix<newcode::qfloat<4>>&,
                                                              std::size_t, std::size_t, const Params&,
                                                              std::size_t, std::size_t, std::size_t,
                                                              std::vector<TileEarlyStopCounter>*, bool,
                                                              const Matrix<float>*);
extern template void process_window_ebchPF<newcode::qfloat<5>>(Matrix<newcode::qfloat<5>>&,
                                                              const Matrix<newcode::qfloat<5>>&,
                                                              std::size_t, std::size_t, const Params&,
                                                              std::size_t, std::size_t, std::size_t,
                                                              std::vector<TileEarlyStopCounter>*, bool,
                                                              const Matrix<float>*);

extern template TileProcessResult<newcode::qfloat<4>> process_tile_plain<newcode::qfloat<4>>(const Matrix<newcode::qfloat<4>>&,
                                                                                             const Matrix<newcode::qfloat<4>>&,
                                                                                             const Params&, std::size_t, bool, bool,
                                                                                             const Matrix<float>*);
extern template TileProcessResult<newcode::qfloat<5>> process_tile_plain<newcode::qfloat<5>>(const Matrix<newcode::qfloat<5>>&,
                                                                                             const Matrix<newcode::qfloat<5>>&,
                                                                                             const Params&, std::size_t, bool, bool,
                                                                                             const Matrix<float>*);
extern template TileProcessResult<newcode::qfloat<4>> process_tile_ebchPF<newcode::qfloat<4>>(const Matrix<newcode::qfloat<4>>&,
                                                                                              const Matrix<newcode::qfloat<4>>&,
                                                                                              const Params&, std::size_t, bool, bool,
                                                                                              const Matrix<float>*);
extern template TileProcessResult<newcode::qfloat<5>> process_tile_ebchPF<newcode::qfloat<5>>(const Matrix<newcode::qfloat<5>>&,
                                                                                              const Matrix<newcode::qfloat<5>>&,
                                                                                              const Params&, std::size_t, bool, bool,
                                                                                              const Matrix<float>*);
} // namespace newcode
