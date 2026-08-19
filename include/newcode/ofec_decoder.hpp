#pragma once
#include <array>
#include <cstddef>
#include <cstdint>
#include <string>
#include <vector>
#include "newcode/common/matrix/matrix.hpp"
#include "newcode/params.hpp"
#include "newcode/common/qfloat/qfloat.hpp"
#include "newcode/ofec/hybrid/hybrid_classifier.hpp"

namespace newcode {

enum class Level56ScheduleBranch : uint8_t {
  K0 = 0,
  KLessThan8,
  KEqual8,
  KGreaterThan8
};

enum class Level56TemporalBranch : uint8_t {
  Disabled = 0,
  NoHistory,
  SupplementHistory,
  CurrentFirst,
  NoFuture
};

enum class Level56TemporalInfoType : uint8_t {
  DecodeInfo = 0,
  EarlyStopInfo
};

enum class Level56DecodeStatus : uint8_t {
  NotDecoded = 0,
  Produced,
  ActionFailed,
  CompletedNoWriteback,
  ForcedEvicted
};

enum class Level56BufferedBatchRetirement : uint8_t {
  Normal = 0,
  FullEarlyStop,
  ForcedEvicted
};

struct Level56ScheduleRoundSample {
  std::size_t round_index = 0;
  std::size_t time_index = 1;
  std::array<std::size_t, 16> initial_counts{};
  std::size_t initial_nonzero_groups = 0;
  Level56ScheduleBranch branch = Level56ScheduleBranch::K0;
  std::array<std::size_t, 16> remaining_before{};
  std::vector<std::size_t> selected_groups;
  std::array<std::size_t, 16> remaining_after{};
  std::size_t used_entries_before = 0;
  std::size_t used_entries_after = 0;
};

struct Level56ScheduleCodeSample {
  std::size_t code_index = 0;
  std::size_t time_index = 1;
  int time_offset = 0;
  Level56TemporalInfoType info_type =
      Level56TemporalInfoType::DecodeInfo;
  Level56DecodeStatus decode_status = Level56DecodeStatus::NotDecoded;
  std::size_t source_level = 0;
  std::size_t source_local_row = 0;
  std::size_t source_global_row = 0;
  std::size_t group_index = 0;
  std::size_t position_in_group = 0;
  bool early_stop_hit = false;
  uint8_t hybrid_class = 0;
  uint8_t resource_eligibility = 0;
  bool planned_hiso = false;
  bool planned_siso = false;
  bool remaining_for_schedule = false;
  uint8_t final_action = 0;
  int assigned_entry_slot = -1;
  int assigned_core = -1;
  bool produced = false;
  bool already_decoded = false;
  bool pending = false;
  bool forced_evicted = false;
  bool writeback_complete = false;
};

struct Level56BufferedRetirementSample {
  std::size_t batch_id = 0;
  std::size_t arrival_time = 0;
  Level56BufferedBatchRetirement reason =
      Level56BufferedBatchRetirement::Normal;
};

struct Level56BufferedTimeSample {
  std::size_t service_time = 0;
  std::size_t arrived_batch_id = 0;
  std::size_t fifo_depth_before = 0;
  std::size_t fifo_depth_after = 0;
  std::size_t head_batch_id_before = 0;
  bool had_head_before = false;
  std::size_t pending_before = 0;
  std::size_t pending_after = 0;
  std::size_t completed_batches = 0;
  std::size_t full_early_stop_batches = 0;
  std::size_t forced_evicted_batches = 0;
  std::size_t window_start_before = 0;
  std::size_t window_start_after = 0;
  bool ordinary_service_used = false;
  std::size_t ordinary_batch_id = 0;
  std::size_t ordinary_schedule_invocation = 0;
  std::vector<Level56BufferedRetirementSample> retirements;
  std::vector<std::size_t> forced_evicted_global_rows;
};

struct Level56ScheduleSample {
  std::size_t invocation = 0;
  bool buffered_fifo_enabled = false;
  std::size_t buffered_service_time = 0;
  std::size_t buffered_batch_id = 0;
  bool temporal_lookahead_enabled = false;
  bool temporal_has_history = false;
  bool temporal_has_future = false;
  std::size_t temporal_x = 0;
  std::size_t temporal_k1 = 0;
  std::size_t temporal_k2 = 0;
  Level56TemporalBranch temporal_branch = Level56TemporalBranch::Disabled;
  std::size_t temporal_t0_group_entries = 0;
  std::size_t temporal_t1_group_entries = 0;
  std::size_t temporal_t0_pending_before = 0;
  std::size_t temporal_t0_pending_after = 0;
  std::size_t temporal_t1_pending_before = 0;
  std::size_t temporal_t1_pending_after = 0;
  std::size_t temporal_t0_new_produced = 0;
  // Temporal samples keep per-time group-entry counts. The legacy
  // group_entry_counts field remains the sum for compatibility.
  std::array<std::size_t, 16> temporal_t0_group_entry_counts{};
  std::array<std::size_t, 16> temporal_t1_group_entry_counts{};
  std::array<std::size_t, 16> initial_counts{};
  std::size_t initial_nonzero_groups = 0;
  Level56ScheduleBranch branch = Level56ScheduleBranch::K0;
  std::vector<Level56ScheduleRoundSample> rounds;
  std::array<std::size_t, 16> group_entry_counts{};
  std::size_t total_group_entries = 0;
  std::size_t planned_hiso_count = 0;
  std::size_t planned_siso_count = 0;
  std::vector<Level56ScheduleCodeSample> codes;
};

struct HybridClassCount {
  std::size_t invocation = 0;
  std::size_t tile_index = 0;
  std::size_t rows_seen_by_hybrid = 0;
  std::size_t class_none_count = 0;
  std::size_t class_bch_hard_decoded_count = 0;
  std::size_t class_clean_count = 0;
  std::size_t class_parity_only_count = 0;
  std::size_t class_one_main_count = 0;
  std::size_t class_one_main_plus_parity_count = 0;
  std::size_t class_two_main_count = 0;
  std::size_t class_suspicious_count = 0;
  std::size_t class_hard_fail_count = 0;
  std::size_t deferred_candidate_count = 0;
  std::size_t deferred_priority_0_count = 0;
  std::size_t deferred_priority_1_count = 0;
  std::size_t deferred_priority_2_count = 0;
  std::size_t deferred_priority_3_count = 0;
  std::size_t deferred_reclaimed_to_hard_finish_count = 0;
};

struct TileEarlyStopSample {
  std::size_t invocation = 0;
  std::size_t tile_index = 0;
  std::size_t rows_total = 0;
  std::size_t rows_passed = 0;
  std::size_t rows_hard_finish = 0;
  std::size_t rows_need_siso_before_mux = 0;
  std::size_t rows_unscheduled = 0;
};

struct TileEarlyStopGroupBindDebugSample {
  std::size_t invocation = 0;
  std::size_t tile_index = 0;
  int condition_mode = 0;
  int bind_group_size = 1;
  std::string raw_early_stop_flags;
  std::string bound_early_stop_flags;
};

struct TileEarlyStopCounter {
  std::size_t triggered = 0;
  std::size_t total = 0;
  std::size_t row_triggered = 0;
  std::size_t row_total = 0;
  std::size_t row_hard_finish = 0;
  std::size_t row_need_siso_before_mux = 0;
  std::size_t row_unscheduled = 0;
  std::vector<TileEarlyStopSample> samples;
  std::vector<TileEarlyStopGroupBindDebugSample> group_bind_debug_samples;
  std::vector<HybridClassCount> hybrid_class_counts;
  std::vector<Level56ScheduleSample> level56_schedule_samples;
  std::vector<Level56BufferedTimeSample> level56_buffered_time_samples;
};

template <typename LLR>
struct TileProcessResult {
  matrix::Matrix<LLR> tile_out;
  bool early_stop_triggered = false;
  std::size_t rows_early_stop = 0;
  std::size_t rows_total = 0;
  std::size_t rows_hard_finish = 0;
  std::size_t rows_need_siso_before_mux = 0;
  std::size_t rows_unscheduled = 0;
  HybridClassCount hybrid_class_count{};
  bool has_group_bind_debug_sample = false;
  TileEarlyStopGroupBindDebugSample group_bind_debug_sample{};
};

// 顶层解码（不同 Chase 变体分别导出）
template <typename LLR>
matrix::Matrix<LLR> ofec_decode_llr_plain(const matrix::Matrix<LLR>& llr_mat, const Params& p,
                                  std::vector<TileEarlyStopCounter>* tile_stats = nullptr,
                                  bool normalize_extrinsic = true,
                                  const matrix::Matrix<float>* tx_llr_ref = nullptr);

template <typename LLR>
matrix::Matrix<LLR> ofec_decode_llr_ebchPF(const matrix::Matrix<LLR>& llr_mat, const Params& p,
                                   std::vector<TileEarlyStopCounter>* tile_stats = nullptr,
                                   bool normalize_extrinsic = true,
                                   const matrix::Matrix<float>* tx_llr_ref = nullptr);

template <typename LLR>
matrix::Matrix<LLR> ofec_decode_llr_topk_pruned(const matrix::Matrix<LLR>& llr_mat, const Params& p,
                                        std::vector<TileEarlyStopCounter>* tile_stats = nullptr,
                                        bool normalize_extrinsic = true,
                                        const matrix::Matrix<float>* tx_llr_ref = nullptr);

template <typename LLR>
matrix::Matrix<LLR> ofec_decode_llr_global_pair(const matrix::Matrix<LLR>& llr_mat, const Params& p,
                                                std::vector<TileEarlyStopCounter>* tile_stats = nullptr,
                                                bool normalize_extrinsic = true,
                                                const matrix::Matrix<float>* tx_llr_ref = nullptr);

template <typename LLR>
matrix::Matrix<LLR> ofec_decode_llr_group_minima(const matrix::Matrix<LLR>& llr_mat, const Params& p,
                                                 std::vector<TileEarlyStopCounter>* tile_stats = nullptr,
                                                 bool normalize_extrinsic = true,
                                                 const matrix::Matrix<float>* tx_llr_ref = nullptr);

// 窗口处理：对 work_llr 的 [win_start, win_end] 做一轮（或多轮）tile 扫描（就地修改）
template <typename LLR>
void process_window_plain(matrix::Matrix<LLR>& work_llr,
                          const matrix::Matrix<LLR>& channel_llr,
                          std::size_t win_start, std::size_t win_end, const Params& p,
                          std::size_t tile_height_rows, std::size_t tile_stride_rows, std::size_t TILES_PER_WIN,
                          std::vector<TileEarlyStopCounter>* tile_stats = nullptr,
                          bool normalize_extrinsic = true,
                          const matrix::Matrix<float>* tx_llr_ref = nullptr);

template <typename LLR>
void process_window_ebchPF(matrix::Matrix<LLR>& work_llr,
                           const matrix::Matrix<LLR>& channel_llr,
                           std::size_t win_start, std::size_t win_end, const Params& p,
                           std::size_t tile_height_rows, std::size_t tile_stride_rows, std::size_t TILES_PER_WIN,
                           std::vector<TileEarlyStopCounter>* tile_stats = nullptr,
                           bool normalize_extrinsic = true,
                           const matrix::Matrix<float>* tx_llr_ref = nullptr);

template <typename LLR>
void process_window_topk_pruned(matrix::Matrix<LLR>& work_llr,
                                const matrix::Matrix<LLR>& channel_llr,
                                std::size_t win_start, std::size_t win_end, const Params& p,
                                std::size_t tile_height_rows, std::size_t tile_stride_rows, std::size_t TILES_PER_WIN,
                                std::vector<TileEarlyStopCounter>* tile_stats = nullptr,
                                bool normalize_extrinsic = true,
                                const matrix::Matrix<float>* tx_llr_ref = nullptr);

template <typename LLR>
void process_window_global_pair(matrix::Matrix<LLR>& work_llr,
                                const matrix::Matrix<LLR>& channel_llr,
                                std::size_t win_start, std::size_t win_end, const Params& p,
                                std::size_t tile_height_rows, std::size_t tile_stride_rows, std::size_t TILES_PER_WIN,
                                std::vector<TileEarlyStopCounter>* tile_stats = nullptr,
                                bool normalize_extrinsic = true,
                                const matrix::Matrix<float>* tx_llr_ref = nullptr);

template <typename LLR>
void process_window_group_minima(matrix::Matrix<LLR>& work_llr,
                                 const matrix::Matrix<LLR>& channel_llr,
                                 std::size_t win_start, std::size_t win_end, const Params& p,
                                 std::size_t tile_height_rows, std::size_t tile_stride_rows, std::size_t TILES_PER_WIN,
                                 std::vector<TileEarlyStopCounter>* tile_stats = nullptr,
                                 bool normalize_extrinsic = true,
                                 const matrix::Matrix<float>* tx_llr_ref = nullptr);

template <typename LLR>
inline void process_window(matrix::Matrix<LLR>& work_llr,
                           const matrix::Matrix<LLR>& channel_llr,
                           std::size_t win_start, std::size_t win_end, const Params& p,
                           std::size_t tile_height_rows, std::size_t tile_stride_rows, std::size_t TILES_PER_WIN,
                           std::vector<TileEarlyStopCounter>* tile_stats = nullptr,
                           bool normalize_extrinsic = true,
                           const matrix::Matrix<float>* tx_llr_ref = nullptr)
{
  process_window_plain(work_llr, channel_llr, win_start, win_end, p,
                       tile_height_rows, tile_stride_rows, TILES_PER_WIN,
                       tile_stats, normalize_extrinsic, tx_llr_ref);
}

// Tile 处理：对单个 tile 的 LLR 做一次“Chase 解码/更新”
template <typename LLR>
TileProcessResult<LLR> process_tile_plain(const matrix::Matrix<LLR>& tile_in,
                                         const matrix::Matrix<LLR>& ch_tile,
                                         const Params& p,
                                         std::size_t tile_top_row_global,
                                         bool use_hard_decode,
                                         bool normalize_extrinsic,
                                          const matrix::Matrix<float>* tx_llr_ref = nullptr);

template <typename LLR>
TileProcessResult<LLR> process_tile_ebchPF(const matrix::Matrix<LLR>& tile_in,
                                          const matrix::Matrix<LLR>& ch_tile,
                                          const Params& p,
                                          std::size_t tile_top_row_global,
                                          bool use_hard_decode,
                                          bool normalize_extrinsic,
                                          const matrix::Matrix<float>* tx_llr_ref = nullptr);

template <typename LLR>
TileProcessResult<LLR> process_tile_topk_pruned(const matrix::Matrix<LLR>& tile_in,
                                                const matrix::Matrix<LLR>& ch_tile,
                                                const Params& p,
                                                std::size_t tile_top_row_global,
                                                bool use_hard_decode,
                                                bool normalize_extrinsic,
                                                const matrix::Matrix<float>* tx_llr_ref = nullptr);

template <typename LLR>
TileProcessResult<LLR> process_tile_global_pair(const matrix::Matrix<LLR>& tile_in,
                                                const matrix::Matrix<LLR>& ch_tile,
                                                const Params& p,
                                                std::size_t tile_top_row_global,
                                                bool use_hard_decode,
                                                bool normalize_extrinsic,
                                                const matrix::Matrix<float>* tx_llr_ref = nullptr);

template <typename LLR>
TileProcessResult<LLR> process_tile_group_minima(const matrix::Matrix<LLR>& tile_in,
                                                 const matrix::Matrix<LLR>& ch_tile,
                                                 const Params& p,
                                                 std::size_t tile_top_row_global,
                                                 bool use_hard_decode,
                                                 bool normalize_extrinsic,
                                                 const matrix::Matrix<float>* tx_llr_ref = nullptr);

template <typename LLR>
inline TileProcessResult<LLR> process_tile(const matrix::Matrix<LLR>& tile_in,
                                           const matrix::Matrix<LLR>& ch_tile,
                                           const Params& p,
                                           std::size_t tile_top_row_global,
                                           bool use_hard_decode,
                                           bool normalize_extrinsic,
                                           const matrix::Matrix<float>* tx_llr_ref = nullptr)
{
  return process_tile_plain(tile_in, ch_tile, p, tile_top_row_global,
                            use_hard_decode, normalize_extrinsic, tx_llr_ref);
}

// ===== extern template（减少重复实例化）=====
// 基础类型
extern template matrix::Matrix<float>  ofec_decode_llr_plain<float >(const matrix::Matrix<float>&,  const Params&, std::vector<TileEarlyStopCounter>*, bool, const matrix::Matrix<float>*);

extern template matrix::Matrix<float>  ofec_decode_llr_ebchPF<float >(const matrix::Matrix<float>&,  const Params&, std::vector<TileEarlyStopCounter>*, bool, const matrix::Matrix<float>*);
extern template matrix::Matrix<float>  ofec_decode_llr_topk_pruned<float >(const matrix::Matrix<float>&,  const Params&, std::vector<TileEarlyStopCounter>*, bool, const matrix::Matrix<float>*);
extern template matrix::Matrix<float>  ofec_decode_llr_global_pair<float >(const matrix::Matrix<float>&,  const Params&, std::vector<TileEarlyStopCounter>*, bool, const matrix::Matrix<float>*);
extern template matrix::Matrix<float>  ofec_decode_llr_group_minima<float >(const matrix::Matrix<float>&,  const Params&, std::vector<TileEarlyStopCounter>*, bool, const matrix::Matrix<float>*);

extern template void process_window_plain<float >(matrix::Matrix<float>&,  const matrix::Matrix<float>&,
                                                  std::size_t, std::size_t, const Params&,
                                                  std::size_t, std::size_t, std::size_t,
                                                  std::vector<TileEarlyStopCounter>*, bool,
                                                  const matrix::Matrix<float>*);
extern template void process_window_ebchPF<float >(matrix::Matrix<float>&,  const matrix::Matrix<float>&,
                                                  std::size_t, std::size_t, const Params&,
                                                  std::size_t, std::size_t, std::size_t,
                                                  std::vector<TileEarlyStopCounter>*, bool,
                                                  const matrix::Matrix<float>*);
extern template void process_window_topk_pruned<float >(matrix::Matrix<float>&,  const matrix::Matrix<float>&,
                                                         std::size_t, std::size_t, const Params&,
                                                         std::size_t, std::size_t, std::size_t,
                                                         std::vector<TileEarlyStopCounter>*, bool,
                                                         const matrix::Matrix<float>*);
extern template void process_window_global_pair<float >(matrix::Matrix<float>&,  const matrix::Matrix<float>&,
                                                        std::size_t, std::size_t, const Params&,
                                                        std::size_t, std::size_t, std::size_t,
                                                        std::vector<TileEarlyStopCounter>*, bool,
                                                        const matrix::Matrix<float>*);
extern template void process_window_group_minima<float >(matrix::Matrix<float>&,  const matrix::Matrix<float>&,
                                                         std::size_t, std::size_t, const Params&,
                                                         std::size_t, std::size_t, std::size_t,
                                                         std::vector<TileEarlyStopCounter>*, bool,
                                                         const matrix::Matrix<float>*);
extern template TileProcessResult<float>  process_tile_plain<float >(const matrix::Matrix<float>&,  const matrix::Matrix<float>&,
                                                                     const Params&, std::size_t, bool, bool,
                                                                     const matrix::Matrix<float>*);

extern template TileProcessResult<float>  process_tile_ebchPF<float >(const matrix::Matrix<float>&,  const matrix::Matrix<float>&,
                                                                      const Params&, std::size_t, bool, bool,
                                                                      const matrix::Matrix<float>*);
extern template TileProcessResult<float>  process_tile_topk_pruned<float >(const matrix::Matrix<float>&,  const matrix::Matrix<float>&,
                                                                           const Params&, std::size_t, bool, bool,
                                                                           const matrix::Matrix<float>*);
extern template TileProcessResult<float>  process_tile_global_pair<float >(const matrix::Matrix<float>&,  const matrix::Matrix<float>&,
                                                                           const Params&, std::size_t, bool, bool,
                                                                           const matrix::Matrix<float>*);
extern template TileProcessResult<float>  process_tile_group_minima<float >(const matrix::Matrix<float>&,  const matrix::Matrix<float>&,
                                                                            const Params&, std::size_t, bool, bool,
                                                                            const matrix::Matrix<float>*);


// qfloat 量化类型（按需开启）
extern template matrix::Matrix<qfloat::qfloat<4>> ofec_decode_llr_plain<qfloat::qfloat<4>>(const matrix::Matrix<qfloat::qfloat<4>>&,
                                                                                    const Params&, std::vector<TileEarlyStopCounter>*, bool,
                                                                                    const matrix::Matrix<float>*);
extern template matrix::Matrix<qfloat::qfloat<5>> ofec_decode_llr_plain<qfloat::qfloat<5>>(const matrix::Matrix<qfloat::qfloat<5>>&,
                                                                                    const Params&, std::vector<TileEarlyStopCounter>*, bool,
                                                                                    const matrix::Matrix<float>*);

extern template matrix::Matrix<qfloat::qfloat<4>> ofec_decode_llr_ebchPF<qfloat::qfloat<4>>(const matrix::Matrix<qfloat::qfloat<4>>&,
                                                                                      const Params&, std::vector<TileEarlyStopCounter>*, bool,
                                                                                      const matrix::Matrix<float>*);
extern template matrix::Matrix<qfloat::qfloat<5>> ofec_decode_llr_ebchPF<qfloat::qfloat<5>>(const matrix::Matrix<qfloat::qfloat<5>>&,
                                                                                      const Params&, std::vector<TileEarlyStopCounter>*, bool,
                                                                                      const matrix::Matrix<float>*);
extern template matrix::Matrix<qfloat::qfloat<4>> ofec_decode_llr_topk_pruned<qfloat::qfloat<4>>(const matrix::Matrix<qfloat::qfloat<4>>&,
                                                                                                   const Params&, std::vector<TileEarlyStopCounter>*, bool,
                                                                                                   const matrix::Matrix<float>*);
extern template matrix::Matrix<qfloat::qfloat<5>> ofec_decode_llr_topk_pruned<qfloat::qfloat<5>>(const matrix::Matrix<qfloat::qfloat<5>>&,
                                                                                                   const Params&, std::vector<TileEarlyStopCounter>*, bool,
                                                                                                   const matrix::Matrix<float>*);
extern template matrix::Matrix<qfloat::qfloat<4>> ofec_decode_llr_group_minima<qfloat::qfloat<4>>(const matrix::Matrix<qfloat::qfloat<4>>&,
                                                                                                    const Params&, std::vector<TileEarlyStopCounter>*, bool,
                                                                                                    const matrix::Matrix<float>*);
extern template matrix::Matrix<qfloat::qfloat<5>> ofec_decode_llr_group_minima<qfloat::qfloat<5>>(const matrix::Matrix<qfloat::qfloat<5>>&,
                                                                                                    const Params&, std::vector<TileEarlyStopCounter>*, bool,
                                                                                                    const matrix::Matrix<float>*);


extern template void process_window<qfloat::qfloat<4>>(matrix::Matrix<qfloat::qfloat<4>>&,
                                                        const matrix::Matrix<qfloat::qfloat<4>>&,
                                                        std::size_t, std::size_t, const Params&,
                                                        std::size_t, std::size_t, std::size_t,
                                                        std::vector<TileEarlyStopCounter>*, bool,
                                                        const matrix::Matrix<float>*);
extern template void process_window<qfloat::qfloat<5>>(matrix::Matrix<qfloat::qfloat<5>>&,
                                                        const matrix::Matrix<qfloat::qfloat<5>>&,
                                                        std::size_t, std::size_t, const Params&,
                                                        std::size_t, std::size_t, std::size_t,
                                                        std::vector<TileEarlyStopCounter>*, bool,
                                                        const matrix::Matrix<float>*);
extern template void process_window_plain<qfloat::qfloat<4>>(matrix::Matrix<qfloat::qfloat<4>>&,
                                                              const matrix::Matrix<qfloat::qfloat<4>>&,
                                                              std::size_t, std::size_t, const Params&,
                                                              std::size_t, std::size_t, std::size_t,
                                                              std::vector<TileEarlyStopCounter>*, bool,
                                                              const matrix::Matrix<float>*);
extern template void process_window_plain<qfloat::qfloat<5>>(matrix::Matrix<qfloat::qfloat<5>>&,
                                                              const matrix::Matrix<qfloat::qfloat<5>>&,
                                                              std::size_t, std::size_t, const Params&,
                                                              std::size_t, std::size_t, std::size_t,
                                                              std::vector<TileEarlyStopCounter>*, bool,
                                                              const matrix::Matrix<float>*);
extern template void process_window_plain<qfloat::qfloat<5>>(matrix::Matrix<qfloat::qfloat<5>>&,
                                                              const matrix::Matrix<qfloat::qfloat<5>>&,
                                                              std::size_t, std::size_t, const Params&,
                                                              std::size_t, std::size_t, std::size_t,
                                                              std::vector<TileEarlyStopCounter>*, bool,
                                                              const matrix::Matrix<float>*);
extern template void process_window_plain<qfloat::qfloat<5>>(matrix::Matrix<qfloat::qfloat<5>>&,
                                                              const matrix::Matrix<qfloat::qfloat<5>>&,
                                                              std::size_t, std::size_t, const Params&,
                                                              std::size_t, std::size_t, std::size_t,
                                                              std::vector<TileEarlyStopCounter>*, bool,
                                                              const matrix::Matrix<float>*);
extern template void process_window_ebchPF<qfloat::qfloat<4>>(matrix::Matrix<qfloat::qfloat<4>>&,
                                                              const matrix::Matrix<qfloat::qfloat<4>>&,
                                                              std::size_t, std::size_t, const Params&,
                                                              std::size_t, std::size_t, std::size_t,
                                                              std::vector<TileEarlyStopCounter>*, bool,
                                                              const matrix::Matrix<float>*);

extern template void process_window_ebchPF<qfloat::qfloat<5>>(matrix::Matrix<qfloat::qfloat<5>>&,
                                                              const matrix::Matrix<qfloat::qfloat<5>>&,
                                                              std::size_t, std::size_t, const Params&,
                                                              std::size_t, std::size_t, std::size_t,
                                                              std::vector<TileEarlyStopCounter>*, bool,
                                                              const matrix::Matrix<float>*);
extern template void process_window_topk_pruned<qfloat::qfloat<4>>(matrix::Matrix<qfloat::qfloat<4>>&,
                                                                    const matrix::Matrix<qfloat::qfloat<4>>&,
                                                                    std::size_t, std::size_t, const Params&,
                                                                    std::size_t, std::size_t, std::size_t,
                                                                    std::vector<TileEarlyStopCounter>*, bool,
                                                                    const matrix::Matrix<float>*);
extern template void process_window_topk_pruned<qfloat::qfloat<5>>(matrix::Matrix<qfloat::qfloat<5>>&,
                                                                    const matrix::Matrix<qfloat::qfloat<5>>&,
                                                                    std::size_t, std::size_t, const Params&,
                                                                    std::size_t, std::size_t, std::size_t,
                                                                    std::vector<TileEarlyStopCounter>*, bool,
                                                                    const matrix::Matrix<float>*);
extern template void process_window_group_minima<qfloat::qfloat<4>>(matrix::Matrix<qfloat::qfloat<4>>&,
                                                                     const matrix::Matrix<qfloat::qfloat<4>>&,
                                                                     std::size_t, std::size_t, const Params&,
                                                                     std::size_t, std::size_t, std::size_t,
                                                                     std::vector<TileEarlyStopCounter>*, bool,
                                                                     const matrix::Matrix<float>*);
extern template void process_window_group_minima<qfloat::qfloat<5>>(matrix::Matrix<qfloat::qfloat<5>>&,
                                                                     const matrix::Matrix<qfloat::qfloat<5>>&,
                                                                     std::size_t, std::size_t, const Params&,
                                                                     std::size_t, std::size_t, std::size_t,
                                                                     std::vector<TileEarlyStopCounter>*, bool,
                                                                     const matrix::Matrix<float>*);
extern template TileProcessResult<qfloat::qfloat<4>> process_tile_plain<qfloat::qfloat<4>>(const matrix::Matrix<qfloat::qfloat<4>>&,
                                                                                             const matrix::Matrix<qfloat::qfloat<4>>&,
                                                                                             const Params&, std::size_t, bool, bool,
                                                                                             const matrix::Matrix<float>*);
extern template TileProcessResult<qfloat::qfloat<5>> process_tile_plain<qfloat::qfloat<5>>(const matrix::Matrix<qfloat::qfloat<5>>&,
                                                                                             const matrix::Matrix<qfloat::qfloat<5>>&,
                                                                                             const Params&, std::size_t, bool, bool,
                                                                                             const matrix::Matrix<float>*);
extern template TileProcessResult<qfloat::qfloat<4>> process_tile_ebchPF<qfloat::qfloat<4>>(const matrix::Matrix<qfloat::qfloat<4>>&,
                                                                                              const matrix::Matrix<qfloat::qfloat<4>>&,
                                                                                              const Params&, std::size_t, bool, bool,
                                                                                              const matrix::Matrix<float>*);
extern template TileProcessResult<qfloat::qfloat<5>> process_tile_ebchPF<qfloat::qfloat<5>>(const matrix::Matrix<qfloat::qfloat<5>>&,
                                                                                              const matrix::Matrix<qfloat::qfloat<5>>&,
                                                                                              const Params&, std::size_t, bool, bool,
                                                                                              const matrix::Matrix<float>*);
extern template TileProcessResult<qfloat::qfloat<4>> process_tile_topk_pruned<qfloat::qfloat<4>>(const matrix::Matrix<qfloat::qfloat<4>>&,
                                                                                                   const matrix::Matrix<qfloat::qfloat<4>>&,
                                                                                                   const Params&, std::size_t, bool, bool,
                                                                                                   const matrix::Matrix<float>*);
extern template TileProcessResult<qfloat::qfloat<5>> process_tile_topk_pruned<qfloat::qfloat<5>>(const matrix::Matrix<qfloat::qfloat<5>>&,
                                                                                                   const matrix::Matrix<qfloat::qfloat<5>>&,
                                                                                                   const Params&, std::size_t, bool, bool,
                                                                                                   const matrix::Matrix<float>*);
extern template TileProcessResult<qfloat::qfloat<4>> process_tile_group_minima<qfloat::qfloat<4>>(const matrix::Matrix<qfloat::qfloat<4>>&,
                                                                                                    const matrix::Matrix<qfloat::qfloat<4>>&,
                                                                                                    const Params&, std::size_t, bool, bool,
                                                                                                    const matrix::Matrix<float>*);
extern template TileProcessResult<qfloat::qfloat<5>> process_tile_group_minima<qfloat::qfloat<5>>(const matrix::Matrix<qfloat::qfloat<5>>&,
                                                                                                    const matrix::Matrix<qfloat::qfloat<5>>&,
                                                                                                    const Params&, std::size_t, bool, bool,
                                                                                                    const matrix::Matrix<float>*);
} // namespace newcode
