#include "newcode/two_stream_shared_runner.hpp"

#include <algorithm>
#include <cassert>
#include <cstddef>
#include <cstdint>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <limits>
#include <stdexcept>
#include <utility>
#include <vector>

#include "newcode/channel/awgn.hpp"
#include "newcode/common/interleaver/interleaver.hpp"
#include "newcode/common/matrix/hard_bits_to_llr_matrix.hpp"
#include "newcode/common/matrix/info_extract.hpp"
#include "newcode/common/qfloat/qfloat.hpp"
#include "newcode/llr_known_prefix.hpp"
#include "newcode/ofec/mux/mux_siso_budget.hpp"
#include "newcode/ofec/common/lin_matrix_adapters.hpp"
#include "newcode/ofec/earlystop/tile_early_stop_group_bind.hpp"
#include "newcode/ofec/earlystop/tile_early_stop_stats.hpp"
#include "newcode/ofec_decoder_hard.hpp"
#include "newcode/ofec_llr_matrix.hpp"
#include "newcode/rx/ber/ber.hpp"
#include "newcode/rx/demod/qam_llr.hpp"
#include "newcode/rx/ofec/chase/decoder_core.hpp"
#include "newcode/tx/bitgen/bitgen.hpp"
#include "newcode/tx/mod/qam.hpp"
#include "newcode/tx/ofecencoder/ofec_encoder.hpp"

#include "ofec/detail/ofec_tile_impl.ipp"

namespace newcode::two_stream_shared {
namespace {

using detail::CoreFn;

struct LlrMode {
  LlrFormat format = LlrFormat::Float;
  std::size_t quant_bits = 16;
};

struct FrontendArtifacts {
  matrix::Matrix<float> tx_llr_mat;
  std::vector<uint8_t> tx_info_bits_ref;
  matrix::Matrix<float> channel_llr;
};

struct StreamDecodeArtifacts {
  matrix::Matrix<float> quantized_pre_decoder_llr;
  matrix::Matrix<float> post_decoder_llr;
  StreamQuantizationStats quantization_stats;
};

struct TwoStreamDecodeArtifacts {
  StreamDecodeArtifacts stream_a;
  StreamDecodeArtifacts stream_b;
  SharedCoreAggregateStats shared_core_stats;
};

template <typename LLR>
struct SharedLlrDecodeArtifacts {
  matrix::Matrix<LLR> out_a;
  matrix::Matrix<LLR> out_b;
  SharedCoreAggregateStats shared_core_stats;
};

template <typename LLR>
struct SharedPreparedTile {
  using TilePreparedT = detail::TilePrepared<LLR>;
  using CoreLLR = typename TilePreparedT::CoreLLR;

  struct Slice {
    std::size_t stream_id = 0;
    std::size_t row_count = 0;
    std::vector<std::size_t> merged_rows;
  };

  TilePreparedT merged;
  std::vector<Slice> slices;
};

template <typename LLR>
struct SharedTileResult {
  matrix::Matrix<LLR> tile_out_a;
  matrix::Matrix<LLR> tile_out_b;
  std::size_t rows_early_stop = 0;
  std::size_t rows_total = 0;
  std::size_t rows_hard_finish = 0;
  std::size_t rows_need_siso_before_mux = 0;
  std::size_t rows_unscheduled = 0;
  std::size_t produced_rows = 0;
  std::size_t failed_rows = 0;
  std::size_t produced_rows_a = 0;
  std::size_t produced_rows_b = 0;
  std::size_t failed_rows_a = 0;
  std::size_t failed_rows_b = 0;
};

LlrMode pick_llr_mode(const Params& params) {
  if (params.LLR_BITS == 16) {
    return {LlrFormat::Float, 16};
  }
  if (params.LLR_BITS >= 2 && params.LLR_BITS <= 15) {
    return {LlrFormat::Quantized, params.LLR_BITS};
  }
  throw std::runtime_error("[ERROR] Params::LLR_BITS must be 2..16");
}

FrontendArtifacts build_stream_frontend(const Params& params,
                                        const PipelineConfig& pipeline,
                                        const StreamConfig& stream,
                                        bool quiet) {
  const bool verbose = !quiet;
  if (verbose) {
    std::cout << "[INFO] (" << stream.label
              << ") building frontend with seeds bitgen/channel="
              << params.BITGEN_SEED << "/" << params.CHANNEL_SEED << "\n";
  }

  auto info_bits = bitgen::generate_bits(params);
  auto code_matrix = ofecencoder::ofec_encode(info_bits, params);

  constexpr float kTxRefLlr = 50.0f;
  matrix::Matrix<float> tx_llr_mat =
      hard_bits_to_llr_matrix(code_matrix, kTxRefLlr);
  auto tx_info_bits_ref = matrix::rx_info_from_bit_llr(tx_llr_mat, params);

  auto coded_bits = flatten_row_major(
      code_matrix, [](uint8_t bit) { return static_cast<uint8_t>(bit & 1u); });
  const std::size_t block_dim =
      static_cast<std::size_t>(Params::BITS_PER_SUBBLOCK_DIM);
  auto interleaver = interleaver::Interleaver::build_from_shape(
      code_matrix.rows(), code_matrix.cols(), block_dim, block_dim,
      pipeline.interleaver_name);
  auto coded_bits_itlv = interleaver.interleave_chunks(coded_bits);

  unsigned bits_per_symbol = pipeline.bits_per_symbol;
  if (bits_per_symbol == 0) {
    bits_per_symbol = 2;
  }

  auto tx_syms = mod::qam_modulate(coded_bits_itlv, bits_per_symbol);
  const int n =
      static_cast<int>(params.NUM_SUBBLOCK_COLS * params.BITS_PER_SUBBLOCK_DIM);
  const int k = 239;
  const int take_bits = k - n;
  const float code_rate =
      static_cast<float>(take_bits) / static_cast<float>(n);
  auto rx_syms = channel::add_awgn(
      tx_syms, stream.ebn0_db, bits_per_symbol,
      static_cast<uint32_t>(params.CHANNEL_SEED));
  auto llr = demod::qam_llr_from_ebn0(
      rx_syms, bits_per_symbol, stream.ebn0_db, code_rate);
  auto llr_deint = interleaver.deinterleave_chunks(llr);
  matrix::Matrix<float> llr_mat =
      llr_to_matrix_row_major(llr_deint, code_matrix.rows(), code_matrix.cols());
  apply_known_zero_prefix(llr_mat, params);

  return FrontendArtifacts{
      .tx_llr_mat = std::move(tx_llr_mat),
      .tx_info_bits_ref = std::move(tx_info_bits_ref),
      .channel_llr = std::move(llr_mat),
  };
}

float compute_shared_quant_clip(const FrontendArtifacts& stream_a,
                                const FrontendArtifacts& stream_b,
                                const Params& params) {
  if (params.LLR_BITS == 16) {
    return 0.0f;
  }

  if (params.LLR_CLIP_RATIO > 0.0f) {
    std::vector<float> llr_values;
    llr_values.reserve(stream_a.channel_llr.rows() * stream_a.channel_llr.cols() +
                       stream_b.channel_llr.rows() * stream_b.channel_llr.cols());
    for (std::size_t row = 0; row < stream_a.channel_llr.rows(); ++row) {
      for (std::size_t col = 0; col < stream_a.channel_llr.cols(); ++col) {
        llr_values.push_back(stream_a.channel_llr[row][col]);
      }
    }
    for (std::size_t row = 0; row < stream_b.channel_llr.rows(); ++row) {
      for (std::size_t col = 0; col < stream_b.channel_llr.cols(); ++col) {
        llr_values.push_back(stream_b.channel_llr[row][col]);
      }
    }
    const float clip = qfloat::compute_clip_from_ratio(
        llr_values.begin(), llr_values.end(), params.LLR_CLIP_RATIO);
    if (clip > 0.0f) {
      return clip;
    }
  }

  return params.LLR_CLIP;
}

template <typename T>
matrix::Matrix<T> clone_matrix(const matrix::Matrix<T>& input) {
  matrix::Matrix<T> output(input.rows(), input.cols());
  for (std::size_t row = 0; row < input.rows(); ++row) {
    for (std::size_t col = 0; col < input.cols(); ++col) {
      output[row][col] = input[row][col];
    }
  }
  return output;
}

SharedCoreAggregateStats make_shared_core_aggregate_stats(
    std::size_t num_streams) {
  SharedCoreAggregateStats stats;
  stats.produced_rows_per_stream.assign(num_streams, 0);
  stats.failed_rows_per_stream.assign(num_streams, 0);
  return stats;
}

template <typename LLR>
matrix::Matrix<LLR> make_zero_llr_matrix(std::size_t rows, std::size_t cols) {
  matrix::Matrix<LLR> output(rows, cols);
  for (std::size_t row = 0; row < rows; ++row) {
    for (std::size_t col = 0; col < cols; ++col) {
      output[row][col] = qfloat::llr_from_float<LLR>(0.0f);
    }
  }
  return output;
}

template <typename LLR>
matrix::Matrix<LLR> extract_tile_matrix(const matrix::Matrix<LLR>& input,
                                        std::size_t tile_top_row,
                                        std::size_t tile_height_rows,
                                        std::size_t cols) {
  matrix::Matrix<LLR> tile(tile_height_rows, cols);
  for (std::size_t row = 0; row < tile_height_rows; ++row) {
    const std::size_t global_row = tile_top_row + row;
    for (std::size_t col = 0; col < cols; ++col) {
      tile[row][col] = input[global_row][col];
    }
  }
  return tile;
}

template <typename LLR>
void overwrite_tile_matrix(const matrix::Matrix<LLR>& tile,
                           std::size_t tile_top_row,
                           matrix::Matrix<LLR>* output) {
  for (std::size_t row = 0; row < tile.rows(); ++row) {
    const std::size_t global_row = tile_top_row + row;
    for (std::size_t col = 0; col < tile.cols(); ++col) {
      (*output)[global_row][col] = tile[row][col];
    }
  }
}

template <typename CoreLLR>
CoreFn<CoreLLR> pick_core_fn(const std::string& decoder_name) {
  if (decoder_name == "chase_baseline" ||
      decoder_name == "two_stream_shared_chase_baseline") {
    return &chase::Decoder_Core_plain<CoreLLR>;
  }
  if (decoder_name == "chase_overall_parity_search") {
    return &chase::Decoder_Core_ebchPF<CoreLLR>;
  }
  if (decoder_name == "chase_topk_pruned") {
    return &chase::Decoder_Core_topk_pruned<CoreLLR>;
  }
  if (decoder_name == "chase_global_pair") {
    return &chase::Decoder_Core_global_pair<CoreLLR>;
  }
  if (decoder_name == "chase_group_minima") {
    return &chase::Decoder_Core_group_minima<CoreLLR>;
  }
  throw std::invalid_argument(
      "run_two_stream_shared: unsupported decoder_name '" + decoder_name + "'");
}

Params build_tile_params(const Params& base_params,
                         std::size_t tile_index,
                         std::size_t* chase_invocation_counter) {
  auto pick_float = [](const std::vector<float>& values,
                       std::size_t idx,
                       float fallback) -> float {
    return idx < values.size() ? values[idx] : fallback;
  };
  auto pick_int = [](const std::vector<int>& values,
                     std::size_t idx,
                     int fallback) -> int {
    return idx < values.size() ? values[idx] : fallback;
  };

  Params tile_params = base_params;
  tile_params.beta = pick_float(base_params.beta_list, tile_index, base_params.beta);
  tile_params.EARLY_STOP_ACTION_SIGN_BETA =
      pick_float(base_params.EARLY_STOP_ACTION_SIGN_BETA_LIST,
                 tile_index,
                 base_params.EARLY_STOP_ACTION_SIGN_BETA);
  tile_params.ENABLE_EARLY_STOP =
      pick_int(base_params.EARLY_STOP_ENABLE_LIST,
               tile_index,
               base_params.ENABLE_EARLY_STOP ? 1 : 0) != 0;
  tile_params.EARLY_STOP_CONDITION_MODE =
      pick_int(base_params.EARLY_STOP_CONDITION_MODE_LIST,
               tile_index,
               base_params.EARLY_STOP_CONDITION_MODE);
  tile_params.EARLY_STOP_ACTION_MODE =
      pick_int(base_params.EARLY_STOP_ACTION_MODE_LIST,
               tile_index,
               base_params.EARLY_STOP_ACTION_MODE);
  tile_params.EARLY_STOP_BIND_GROUP_SIZE =
      pick_int(base_params.EARLY_STOP_BIND_GROUP_SIZE_LIST,
               tile_index,
               base_params.EARLY_STOP_BIND_GROUP_SIZE);
  tile_params.HYBRID_ENABLE =
      pick_int(base_params.HYBRID_ENABLE_LIST,
               tile_index,
               base_params.HYBRID_ENABLE ? 1 : 0) != 0;
  tile_params.HYBRID_HARD_LLR_MAG =
      pick_float(base_params.HYBRID_HARD_LLR_MAG_LIST,
                 tile_index,
                 base_params.HYBRID_HARD_LLR_MAG);
  tile_params.ALPHA =
      pick_float(base_params.ALPHA_LIST, tile_index, base_params.ALPHA);
  tile_params.debug_trace.chase_tile_index = static_cast<int>(tile_index);
  tile_params.debug_trace.chase_invocation =
      static_cast<int>(++(*chase_invocation_counter));
  return tile_params;
}

template <typename LLR>
SharedPreparedTile<LLR> merge_prepared_tiles(
    const detail::TilePrepared<LLR>& prep_a,
    const detail::TilePrepared<LLR>& prep_b) {
  using SharedPrepared = SharedPreparedTile<LLR>;
  using CoreLLR = typename SharedPrepared::CoreLLR;

  if (prep_a.lin_matrix.cols() != prep_b.lin_matrix.cols() ||
      prep_a.lch_matrix.cols() != prep_b.lch_matrix.cols()) {
    throw std::invalid_argument(
        "merge_prepared_tiles: stream A/B prepared tile column counts must match");
  }

  const std::size_t cols = prep_a.lin_matrix.cols();
  const std::size_t rows_a = prep_a.lin_matrix.rows();
  const std::size_t rows_b = prep_b.lin_matrix.rows();
  const std::size_t total_rows = rows_a + rows_b;

  SharedPrepared shared;
  shared.merged.lin_matrix = matrix::Matrix<CoreLLR>(total_rows, cols);
  shared.merged.lch_matrix = matrix::Matrix<CoreLLR>(total_rows, cols);
  shared.merged.row_local_lookup.resize(total_rows, 0);
  shared.merged.row_global_lookup.resize(total_rows, 0);
  shared.merged.params_for_core = prep_a.params_for_core;
  shared.merged.trace = prep_a.trace;

  const bool has_expected_bits =
      prep_a.expected_bits && prep_b.expected_bits;
  if (has_expected_bits) {
    shared.merged.expected_bits =
        std::make_shared<std::vector<std::vector<int8_t>>>(
            total_rows, std::vector<int8_t>(cols, -1));
  }

  typename SharedPrepared::Slice slice_a;
  slice_a.stream_id = 0;
  slice_a.row_count = rows_a;
  slice_a.merged_rows.reserve(rows_a);

  typename SharedPrepared::Slice slice_b;
  slice_b.stream_id = 1;
  slice_b.row_count = rows_b;
  slice_b.merged_rows.reserve(rows_b);

  auto append_row = [&](const detail::TilePrepared<LLR>& src,
                        std::size_t src_row,
                        typename SharedPrepared::Slice* slice,
                        std::size_t dst_row) {
    slice->merged_rows.push_back(dst_row);
    shared.merged.row_local_lookup[dst_row] = src.row_local_lookup[src_row];
    shared.merged.row_global_lookup[dst_row] = src.row_global_lookup[src_row];
    for (std::size_t col = 0; col < cols; ++col) {
      shared.merged.lin_matrix[dst_row][col] = src.lin_matrix[src_row][col];
      shared.merged.lch_matrix[dst_row][col] = src.lch_matrix[src_row][col];
    }
    if (has_expected_bits) {
      (*shared.merged.expected_bits)[dst_row] = (*src.expected_bits)[src_row];
    }
  };

  std::size_t dst_row = 0;
  const std::size_t paired_rows = std::min(rows_a, rows_b);
  for (std::size_t row = 0; row < paired_rows; ++row) {
    append_row(prep_a, row, &slice_a, dst_row++);
    append_row(prep_b, row, &slice_b, dst_row++);
  }
  for (std::size_t row = paired_rows; row < rows_a; ++row) {
    append_row(prep_a, row, &slice_a, dst_row++);
  }
  for (std::size_t row = paired_rows; row < rows_b; ++row) {
    append_row(prep_b, row, &slice_b, dst_row++);
  }

  if (dst_row != total_rows) {
    throw std::logic_error(
        "merge_prepared_tiles: merged row count does not match total rows");
  }

  shared.slices.push_back(std::move(slice_a));
  shared.slices.push_back(std::move(slice_b));

  if (shared.merged.expected_bits) {
    shared.merged.params_for_core.debug_trace.chase_expected_bits =
        shared.merged.expected_bits;
  } else {
    shared.merged.params_for_core.debug_trace.chase_expected_bits.reset();
  }
  shared.merged.params_for_core.debug_trace.chase_expected_bits_row = nullptr;
  shared.merged.params_for_core.debug_trace.active_chase_entries.clear();
  shared.merged.params_for_core.debug_trace.chase_decoder_row = -1;
  shared.merged.params_for_core.debug_trace.chase_decoder_col = -1;

  return shared;
}

template <typename LLR>
SharedTileResult<LLR> run_shared_tile(
    const matrix::Matrix<LLR>& tile_in_a,
    const matrix::Matrix<LLR>& tile_in_b,
    const matrix::Matrix<LLR>& ch_tile_a,
    const matrix::Matrix<LLR>& ch_tile_b,
    const Params& tile_params,
    std::size_t tile_top_row,
    int siso_active_for_tile,
    bool capture_last_tile_history,
    bool normalize_extrinsic,
    const matrix::Matrix<float>* tx_llr_ref_a,
    const matrix::Matrix<float>* tx_llr_ref_b,
    matrix::Matrix<float>* last_history_a,
    matrix::Matrix<float>* last_history_b,
    CoreFn<typename detail::TilePrepared<LLR>::CoreLLR> core_fn) {
  const std::size_t rows_to_decode =
      static_cast<std::size_t>(tile_params.CHASE_SBR) *
      Params::BITS_PER_SUBBLOCK_DIM;

  auto prep_a = detail::prepare_tile_inputs(
      tile_in_a, ch_tile_a, tile_params, tile_top_row, tile_params.CHASE_SBR,
      rows_to_decode, tx_llr_ref_a);
  auto prep_b = detail::prepare_tile_inputs(
      tile_in_b, ch_tile_b, tile_params, tile_top_row, tile_params.CHASE_SBR,
      rows_to_decode, tx_llr_ref_b);
  auto shared_prep = merge_prepared_tiles(prep_a, prep_b);

  TileEarlyStopResult early_stop_stats;
  if (tile_params.ENABLE_EARLY_STOP) {
    early_stop_stats =
        detect_tile_early_stop(shared_prep.merged.lin_matrix, tile_params);
    if (tile_params.EARLY_STOP_CONDITION_MODE == 1 &&
        tile_params.EARLY_STOP_BIND_GROUP_SIZE > 1) {
      early_stop_stats = apply_group_bound_early_stop(
          early_stop_stats, tile_params.EARLY_STOP_BIND_GROUP_SIZE);
    }
  } else {
    early_stop_stats.row_passed_flags.assign(shared_prep.merged.lin_matrix.rows(),
                                             false);
    early_stop_stats.row_details.assign(shared_prep.merged.lin_matrix.rows(),
                                        TileEarlyStopRowDetail{});
    early_stop_stats.rows_passed = 0;
    early_stop_stats.rows_total = shared_prep.merged.lin_matrix.rows();
    early_stop_stats.all_rows_passed = false;
  }

  auto dispatch_plan = detail::build_tile_dispatch_plan(
      shared_prep.merged, early_stop_stats, siso_active_for_tile, tile_params);
  const std::size_t rows_need_siso_before_mux =
      detail::count_soft_candidates_before_mux(dispatch_plan);
  detail::run_mux_on_soft_candidates(
      &dispatch_plan, siso_active_for_tile, tile_params, early_stop_stats);
  auto decoder_res = detail::decode_tile_with_plan<LLR>(
      shared_prep.merged, dispatch_plan, normalize_extrinsic, core_fn);

  SharedTileResult<LLR> result;
  result.tile_out_a = tile_in_a;
  result.tile_out_b = tile_in_b;
  result.rows_early_stop = early_stop_stats.rows_passed;
  result.rows_total = early_stop_stats.rows_total;
  result.rows_hard_finish = dispatch_plan.rows_hard_finish;
  result.rows_need_siso_before_mux = rows_need_siso_before_mux;
  result.rows_unscheduled = dispatch_plan.rows_soft_unscheduled;

  for (std::size_t row = 0; row < decoder_res.produced_rows.size(); ++row) {
    if (decoder_res.produced_rows[row]) {
      ++result.produced_rows;
    }
  }
  result.failed_rows = decoder_res.produced_rows.size() - result.produced_rows;

  for (const auto& slice : shared_prep.slices) {
    chase::DecoderCoreResult<typename detail::TilePrepared<LLR>::CoreLLR> sliced;
    sliced.lout =
        matrix::Matrix<float>(slice.row_count, decoder_res.lout.cols());
    sliced.produced_rows.assign(slice.row_count, false);
    for (std::size_t row = 0; row < slice.row_count; ++row) {
      const std::size_t merged_row = slice.merged_rows[row];
      sliced.produced_rows[row] = decoder_res.produced_rows[merged_row];
      if (sliced.produced_rows[row]) {
        if (slice.stream_id == 0) {
          ++result.produced_rows_a;
        } else if (slice.stream_id == 1) {
          ++result.produced_rows_b;
        }
      }
      for (std::size_t col = 0; col < decoder_res.lout.cols(); ++col) {
        sliced.lout[row][col] = decoder_res.lout[merged_row][col];
      }
    }
    const std::size_t slice_failed_rows = slice.row_count -
                                          static_cast<std::size_t>(std::count(
                                              sliced.produced_rows.begin(),
                                              sliced.produced_rows.end(),
                                              true));

    if (slice.stream_id == 0) {
      result.failed_rows_a += slice_failed_rows;
      detail::writeback_tile(
          prep_a, sliced, tile_params, tile_top_row,
          capture_last_tile_history, &result.tile_out_a, last_history_a);
    } else if (slice.stream_id == 1) {
      result.failed_rows_b += slice_failed_rows;
      detail::writeback_tile(
          prep_b, sliced, tile_params, tile_top_row,
          capture_last_tile_history, &result.tile_out_b, last_history_b);
    }
  }

  return result;
}

template <typename LLR>
SharedLlrDecodeArtifacts<LLR> decode_two_stream_shared_llr(
    const matrix::Matrix<LLR>& llr_a,
    const matrix::Matrix<LLR>& llr_b,
    const Params& params,
    const PipelineConfig& pipeline,
    const matrix::Matrix<float>* tx_llr_ref_a,
    const matrix::Matrix<float>* tx_llr_ref_b) {
  using CoreLLR = typename detail::TilePrepared<LLR>::CoreLLR;

  if (llr_a.rows() != llr_b.rows() || llr_a.cols() != llr_b.cols()) {
    throw std::invalid_argument(
        "decode_two_stream_shared_llr: stream A/B matrix shapes must match");
  }

  const std::size_t rows = llr_a.rows();
  const std::size_t cols = llr_a.cols();
  const std::size_t expected_cols =
      Params::NUM_SUBBLOCK_COLS * Params::BITS_PER_SUBBLOCK_DIM;
  if (cols != expected_cols) {
    throw std::invalid_argument(
        "decode_two_stream_shared_llr: input llr_mat cols != N");
  }

  if (rows < params.win_height_rows()) {
    return {
        .out_a = clone_matrix(llr_a),
        .out_b = clone_matrix(llr_b),
        .shared_core_stats = make_shared_core_aggregate_stats(2),
    };
  }

  matrix::Matrix<LLR> channel_a = clone_matrix(llr_a);
  matrix::Matrix<LLR> channel_b = clone_matrix(llr_b);
  matrix::Matrix<LLR> work_a = make_zero_llr_matrix<LLR>(rows, cols);
  matrix::Matrix<LLR> work_b = make_zero_llr_matrix<LLR>(rows, cols);
  matrix::Matrix<float> last_history_a(rows, cols);
  matrix::Matrix<float> last_history_b(rows, cols);
  SharedCoreAggregateStats shared_core_stats =
      make_shared_core_aggregate_stats(2);

  const std::size_t tile_height_rows = params.tile_height_rows();
  const std::size_t tile_stride_rows = params.tile_stride_rows();
  const std::size_t win_height_rows = params.win_height_rows();
  const std::size_t pop_push_rows = params.pop_push_rows();
  const std::size_t tiles_per_window = params.TILES_PER_WIN;
  std::size_t chase_invocation_counter = 0;
  auto core_fn = pick_core_fn<CoreLLR>(pipeline.decoder_name);

  std::vector<bool> hard_tile_mask(tiles_per_window);
  int last_soft_tile_idx = -1;
  auto pick_int = [](const std::vector<int>& values,
                     std::size_t idx,
                     int fallback) -> int {
    return idx < values.size() ? values[idx] : fallback;
  };
  for (std::size_t tile_index = 0; tile_index < tiles_per_window; ++tile_index) {
    const bool is_hard =
        pick_int(params.HARD_TILE_LIST, tile_index,
                 params.HARD_DECODE_DEFAULT ? 1 : 0) != 0;
    hard_tile_mask[tile_index] = is_hard;
    if (!is_hard) {
      last_soft_tile_idx = static_cast<int>(tile_index);
    }
  }

  std::size_t win_start = params.initial_win_start_rows();
  const std::size_t last_ws = rows - win_height_rows;
  while (win_start <= last_ws) {
    const std::size_t win_end = win_start + win_height_rows - 1;
    for (std::size_t tile_index = 0; tile_index < tiles_per_window;
         ++tile_index) {
      const std::size_t tile_bottom_row = win_end - tile_index * tile_stride_rows;
      const std::size_t tile_top_row = tile_bottom_row + 1 - tile_height_rows;

      auto tile_in_a =
          extract_tile_matrix(work_a, tile_top_row, tile_height_rows, cols);
      auto tile_in_b =
          extract_tile_matrix(work_b, tile_top_row, tile_height_rows, cols);
      auto ch_tile_a =
          extract_tile_matrix(channel_a, tile_top_row, tile_height_rows, cols);
      auto ch_tile_b =
          extract_tile_matrix(channel_b, tile_top_row, tile_height_rows, cols);

      const bool use_hard = hard_tile_mask[tile_index];
      if (use_hard) {
        throw std::invalid_argument(
            "decode_two_stream_shared_llr: shared merged-tile path currently "
            "supports soft tiles only; hard-tile compatibility has not been "
            "wired in yet");
      }
      const bool use_history_input =
          use_hard && last_soft_tile_idx >= 0 &&
          static_cast<int>(tile_index) > last_soft_tile_idx;
      if (use_history_input) {
        for (std::size_t row = 0; row < tile_height_rows; ++row) {
          const std::size_t global_row = tile_top_row + row;
          for (std::size_t col = 0; col < cols; ++col) {
            tile_in_a[row][col] = qfloat::llr_from_float<LLR>(
                last_history_a[global_row][col]);
            tile_in_b[row][col] = qfloat::llr_from_float<LLR>(
                last_history_b[global_row][col]);
          }
        }
      }

      Params tile_params =
          build_tile_params(params, tile_index, &chase_invocation_counter);
      const int siso_active_for_tile =
          mux::pick_siso_active_for_tile(params.SISO_ACTIVE_LIST, tile_index);
      const bool capture_last_tile_history =
          last_soft_tile_idx >= 0 &&
          static_cast<int>(tile_index) == last_soft_tile_idx;

      auto tile_result = run_shared_tile<LLR>(
          tile_in_a, tile_in_b, ch_tile_a, ch_tile_b, tile_params, tile_top_row,
          siso_active_for_tile, capture_last_tile_history,
          pipeline.normalize_extrinsic, tx_llr_ref_a, tx_llr_ref_b,
          &last_history_a, &last_history_b, core_fn);

      overwrite_tile_matrix(tile_result.tile_out_a, tile_top_row, &work_a);
      overwrite_tile_matrix(tile_result.tile_out_b, tile_top_row, &work_b);

      shared_core_stats.total_batches += 1;
      shared_core_stats.total_rows += tile_result.rows_total;
      shared_core_stats.produced_rows += tile_result.produced_rows;
      shared_core_stats.failed_rows += tile_result.failed_rows;
      if (shared_core_stats.produced_rows_per_stream.size() < 2) {
        shared_core_stats.produced_rows_per_stream.assign(2, 0);
        shared_core_stats.failed_rows_per_stream.assign(2, 0);
      }
      shared_core_stats.produced_rows_per_stream[0] += tile_result.produced_rows_a;
      shared_core_stats.produced_rows_per_stream[1] += tile_result.produced_rows_b;
      shared_core_stats.failed_rows_per_stream[0] += tile_result.failed_rows_a;
      shared_core_stats.failed_rows_per_stream[1] += tile_result.failed_rows_b;
    }
    win_start += pop_push_rows;
  }

  matrix::Matrix<LLR> out_a(rows, cols);
  matrix::Matrix<LLR> out_b(rows, cols);
  for (std::size_t row = 0; row < rows; ++row) {
    for (std::size_t col = 0; col < cols; ++col) {
      out_a[row][col] = qfloat::llr_from_float<LLR>(
          qfloat::llr_to_float(channel_a[row][col]) + last_history_a[row][col]);
      out_b[row][col] = qfloat::llr_from_float<LLR>(
          qfloat::llr_to_float(channel_b[row][col]) + last_history_b[row][col]);
    }
  }

  return {
      .out_a = std::move(out_a),
      .out_b = std::move(out_b),
      .shared_core_stats = std::move(shared_core_stats),
  };
}

template <int NBITS>
StreamQuantizationStats collect_qfloat_quantization_stats(
    const matrix::Matrix<qfloat::qfloat<NBITS>>& matrix_in,
    float quant_clip) {
  StreamQuantizationStats stats;
  stats.quant_clip = quant_clip;
  stats.total_values = matrix_in.rows() * matrix_in.cols();
  constexpr int kLoCode = -(1 << (NBITS - 1));
  constexpr int kHiCode = (1 << (NBITS - 1)) - 1;
  for (std::size_t row = 0; row < matrix_in.rows(); ++row) {
    for (std::size_t col = 0; col < matrix_in.cols(); ++col) {
      const int code = matrix_in[row][col].code();
      if (code == kLoCode || code == kHiCode) {
        ++stats.saturated_values;
      }
    }
  }
  if (stats.total_values > 0) {
    stats.saturation_ratio =
        static_cast<double>(stats.saturated_values) /
        static_cast<double>(stats.total_values);
  }
  return stats;
}

template <int NBITS>
TwoStreamDecodeArtifacts decode_two_stream_quantized(
    const FrontendArtifacts& stream_a,
    const FrontendArtifacts& stream_b,
    const Config& config,
    float quant_clip) {
  using Q = qfloat::qfloat<NBITS>;
  Q::set_clip(quant_clip);

  auto quantized_a =
      qfloat::quantize_matrix_to_qfloat<NBITS>(stream_a.channel_llr, quant_clip);
  auto quantized_b =
      qfloat::quantize_matrix_to_qfloat<NBITS>(stream_b.channel_llr, quant_clip);
  auto decoded_pair = decode_two_stream_shared_llr(
      quantized_a, quantized_b, config.params, config.pipeline,
      &stream_a.tx_llr_mat, &stream_b.tx_llr_mat);

  return {
      .stream_a =
          StreamDecodeArtifacts{
              .quantized_pre_decoder_llr =
                  qfloat::dequantize_matrix_from_qfloat(quantized_a, quant_clip),
              .post_decoder_llr =
                  qfloat::dequantize_matrix_from_qfloat(decoded_pair.out_a, quant_clip),
              .quantization_stats =
                  collect_qfloat_quantization_stats<NBITS>(quantized_a, quant_clip),
          },
      .stream_b =
          StreamDecodeArtifacts{
              .quantized_pre_decoder_llr =
                  qfloat::dequantize_matrix_from_qfloat(quantized_b, quant_clip),
              .post_decoder_llr =
                  qfloat::dequantize_matrix_from_qfloat(decoded_pair.out_b, quant_clip),
              .quantization_stats =
                  collect_qfloat_quantization_stats<NBITS>(quantized_b, quant_clip),
          },
      .shared_core_stats = std::move(decoded_pair.shared_core_stats),
  };
}

TwoStreamDecodeArtifacts decode_two_stream(const FrontendArtifacts& stream_a,
                                           const FrontendArtifacts& stream_b,
                                           const Config& config,
                                           float quant_clip) {
  const LlrMode llr_mode = pick_llr_mode(config.params);
  switch (llr_mode.format) {
    case LlrFormat::Float: {
      auto decoded_pair = decode_two_stream_shared_llr(
          stream_a.channel_llr, stream_b.channel_llr, config.params,
          config.pipeline, &stream_a.tx_llr_mat, &stream_b.tx_llr_mat);
      return {
          .stream_a =
              StreamDecodeArtifacts{
                  .quantized_pre_decoder_llr = matrix::Matrix<float>{},
                  .post_decoder_llr = std::move(decoded_pair.out_a),
                  .quantization_stats = StreamQuantizationStats{},
              },
          .stream_b =
              StreamDecodeArtifacts{
                  .quantized_pre_decoder_llr = matrix::Matrix<float>{},
                  .post_decoder_llr = std::move(decoded_pair.out_b),
                  .quantization_stats = StreamQuantizationStats{},
              },
          .shared_core_stats = std::move(decoded_pair.shared_core_stats),
      };
    }
    case LlrFormat::Quantized:
      switch (llr_mode.quant_bits) {
        case 2:
          return decode_two_stream_quantized<2>(
              stream_a, stream_b, config, quant_clip);
        case 3:
          return decode_two_stream_quantized<3>(
              stream_a, stream_b, config, quant_clip);
        case 4:
          return decode_two_stream_quantized<4>(
              stream_a, stream_b, config, quant_clip);
        case 5:
          return decode_two_stream_quantized<5>(
              stream_a, stream_b, config, quant_clip);
        case 6:
          return decode_two_stream_quantized<6>(
              stream_a, stream_b, config, quant_clip);
        case 7:
          return decode_two_stream_quantized<7>(
              stream_a, stream_b, config, quant_clip);
        case 8:
          return decode_two_stream_quantized<8>(
              stream_a, stream_b, config, quant_clip);
        case 9:
          return decode_two_stream_quantized<9>(
              stream_a, stream_b, config, quant_clip);
        case 10:
          return decode_two_stream_quantized<10>(
              stream_a, stream_b, config, quant_clip);
        case 11:
          return decode_two_stream_quantized<11>(
              stream_a, stream_b, config, quant_clip);
        case 12:
          return decode_two_stream_quantized<12>(
              stream_a, stream_b, config, quant_clip);
        case 13:
          return decode_two_stream_quantized<13>(
              stream_a, stream_b, config, quant_clip);
        case 14:
          return decode_two_stream_quantized<14>(
              stream_a, stream_b, config, quant_clip);
        case 15:
          return decode_two_stream_quantized<15>(
              stream_a, stream_b, config, quant_clip);
        default:
          break;
      }
      throw std::runtime_error(
          "[ERROR] Unsupported quant_bits for two-stream shared decoder");
  }
  throw std::runtime_error("[ERROR] Unreachable llr mode");
}

PipelineResult build_pipeline_result(const FrontendArtifacts& frontend,
                                     const StreamDecodeArtifacts& decoded,
                                     const Params& params,
                                     const std::string& label,
                                     bool quiet) {
  auto rx_info_bits_pre =
      matrix::rx_info_from_bit_llr(frontend.channel_llr, params);
  std::vector<uint8_t> rx_info_bits_pre_quantized;
  const bool has_pre_quantized =
      decoded.quantized_pre_decoder_llr.rows() > 0 &&
      decoded.quantized_pre_decoder_llr.cols() > 0;
  if (has_pre_quantized) {
    rx_info_bits_pre_quantized =
        matrix::rx_info_from_bit_llr(decoded.quantized_pre_decoder_llr, params);
  }
  auto rx_info_bits_post =
      matrix::rx_info_from_bit_llr(decoded.post_decoder_llr, params);

  PipelineResult result;
  result.ebn0_db = std::numeric_limits<float>::quiet_NaN();
  result.pre_fec_windows = compute_ber_per_window(
      frontend.tx_info_bits_ref, rx_info_bits_pre, params);
  result.pre_fec_tile_windows = compute_ber_per_tile_window(
      frontend.tx_info_bits_ref, rx_info_bits_pre, params);
  result.pre_fec = compute_and_print_ber(
      frontend.tx_info_bits_ref, rx_info_bits_pre,
      (label + " Pre-FEC").c_str(), params, &result.pre_fec_error_positions,
      quiet);

  if (has_pre_quantized) {
    result.pre_fec_quantized_hard_windows = compute_ber_per_window(
        frontend.tx_info_bits_ref, rx_info_bits_pre_quantized, params);
    result.pre_fec_quantized_hard_tile_windows = compute_ber_per_tile_window(
        frontend.tx_info_bits_ref, rx_info_bits_pre_quantized, params);
    result.pre_fec_quantized_hard = compute_and_print_ber(
        frontend.tx_info_bits_ref, rx_info_bits_pre_quantized,
        (label + " Pre-FEC (quantized hard)").c_str(), params,
        &result.pre_fec_quantized_hard_error_positions, quiet);
    result.has_pre_fec_quantized_hard = true;
  }

  result.post_fec_windows = compute_ber_per_window(
      frontend.tx_info_bits_ref, rx_info_bits_post, params);
  result.post_fec_tile_windows = compute_ber_per_tile_window(
      frontend.tx_info_bits_ref, rx_info_bits_post, params);
  result.post_fec = compute_and_print_ber(
      frontend.tx_info_bits_ref, rx_info_bits_post,
      (label + " Post-FEC").c_str(), params,
      &result.post_fec_error_positions, quiet);
  return result;
}

}  // namespace

Result run_two_stream_shared(const Config& config) {
  if (config.stream_a.label.empty() || config.stream_b.label.empty()) {
    throw std::invalid_argument(
        "run_two_stream_shared: stream labels must not be empty");
  }

  Params params_a = config.params;
  Params params_b = config.params;
  params_a.BITGEN_SEED = config.stream_a.bitgen_seed;
  params_a.CHANNEL_SEED = config.stream_a.channel_seed;
  params_b.BITGEN_SEED = config.stream_b.bitgen_seed;
  params_b.CHANNEL_SEED = config.stream_b.channel_seed;

  if (!config.pipeline.quiet) {
    std::cout
        << "[INFO] shared decoder path: shared dispatch plan + shared mux + shared Chase core\n";
  }

  FrontendArtifacts frontend_a = build_stream_frontend(
      params_a, config.pipeline, config.stream_a, config.pipeline.quiet);
  FrontendArtifacts frontend_b = build_stream_frontend(
      params_b, config.pipeline, config.stream_b, config.pipeline.quiet);

  Params decode_params = config.params;
  const float quant_clip =
      compute_shared_quant_clip(frontend_a, frontend_b, decode_params);
  if (pick_llr_mode(decode_params).format == LlrFormat::Quantized) {
    decode_params.LLR_CLIP = quant_clip;
  }

  Config decode_config = config;
  decode_config.params = decode_params;
  auto decoded_pair =
      decode_two_stream(frontend_a, frontend_b, decode_config, quant_clip);

  Result result;
  result.stream_a = build_pipeline_result(frontend_a, decoded_pair.stream_a,
                                          decode_params, config.stream_a.label,
                                          config.pipeline.quiet);
  result.stream_a.ebn0_db = config.stream_a.ebn0_db;
  result.stream_b = build_pipeline_result(frontend_b, decoded_pair.stream_b,
                                          decode_params, config.stream_b.label,
                                          config.pipeline.quiet);
  result.stream_b.ebn0_db = config.stream_b.ebn0_db;
  result.observability.shared_quant_clip = quant_clip;
  result.observability.stream_a_quantization =
      decoded_pair.stream_a.quantization_stats;
  result.observability.stream_b_quantization =
      decoded_pair.stream_b.quantization_stats;
  result.observability.shared_core = std::move(decoded_pair.shared_core_stats);
  return result;
}

}  // namespace newcode::two_stream_shared
