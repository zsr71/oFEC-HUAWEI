#pragma once

#include "ofec_tile_input.ipp"
#include "ofec_tile_extrinsic_normalize.ipp"
#include "newcode/ofec/earlystop/row_early_stop_action.hpp"
#include "newcode/ofec/mux/mux_state_builder.hpp"

#include <cstdint>

namespace newcode {
namespace detail {

template <typename CoreLLR>
chase::DecoderCoreResult<CoreLLR> run_decoder_core(
    const matrix::Matrix<CoreLLR>& lin_matrix,
    const matrix::Matrix<CoreLLR>& lch_matrix,
    bool use_hard_decode,
    const newcode::Params& params_for_core,
    const std::vector<bool>* early_stop_row_flags,
    const std::vector<uint8_t>* mux_state,
    CoreFn<CoreLLR> core_fn) {
  return core_fn(lin_matrix,
                 lch_matrix,
                 use_hard_decode,
                 params_for_core,
                 early_stop_row_flags,
                 mux_state);
}

template <typename Float>
void normalize_selected_rows(matrix::Matrix<Float>& lout,
                             const std::vector<bool>& produced_rows,
                             const std::vector<bool>& normalize_rows,
                             Float beta) {
  // 对输出 lout 做“全局幅度归一化”，但只作用在指定行集合上。
  //
  // 设计目的：
  // - 原先 `normalize_extrinsic_lout(...)` 是对所有 produced rows 做一次全局归一化；
  // - 方案三引入后，tile 内可能同时存在：
  //   - early-stop action 的行
  //   - hard-finish 的行
  //   - soft-decode 的行
  //   其中只有 soft-decode 行的外信息幅度需要按 beta 做归一化，其他行可能不该被“牵连缩放”。
  //
  // 输入语义：
  // - produced_rows[r]=true：这一行本轮产生了输出（lout 有效）
  // - normalize_rows[r]=true：这一行需要参与归一化统计，并会被同一个 scale 缩放
  // - beta：用于识别“fallback ±beta” 的哨兵值（避免把 fallback 常量当成真实幅度统计进去）
  auto is_fallback = [&](Float w) -> bool {
    const Float target = beta;
    const Float diff = std::fabs(std::fabs(w) - target);
    const Float tol = static_cast<Float>(1e-4f) *
                      std::max(static_cast<Float>(1.0f), target);
    return diff <= tol;
  };

  // 第一阶段：计算目标行集合里每个元素的平均绝对幅度 g_alpha。
  // 这里用 “fallback 修正” 的方式，避免 |w|==beta 时把它当作真实外信息幅度。
  double acc = 0.0;
  std::size_t cnt = 0;
  for (std::size_t r = 0; r < lout.rows(); ++r) {
    if (r >= produced_rows.size() || !produced_rows[r]) continue;
    if (r >= normalize_rows.size() || !normalize_rows[r]) continue;
    for (std::size_t j = 0; j < lout.cols(); ++j) {
      const Float w = lout[r][j];
      if (is_fallback(w)) {
        // fallback 位置：把它按 (|w|/beta) 归一到“单位幅度”参与统计，
        // 等价于对 fallback 常量做去幅度化，防止它主导整体 scale。
        acc += std::fabs(w / beta);
        ++cnt;
        continue;
      }
      acc += std::fabs(w);
      ++cnt;
    }
  }

  if (cnt == 0) {
    // 没有任何可归一化的元素：直接返回，保持原值。
    return;
  }

  const Float g_alpha = static_cast<Float>(acc / static_cast<double>(cnt));
  if (g_alpha <= static_cast<Float>(0.0f)) {
    // 防御：避免除零或负数导致的 NaN/inf。
    return;
  }

  // 第二阶段：用同一个 scale 对选中行做统一缩放。
  // 注意：这里是“全局 scale”，不是逐行 scale，这样和旧逻辑保持一致。
  const Float scale = static_cast<Float>(1.0f) / g_alpha;
  for (std::size_t r = 0; r < lout.rows(); ++r) {
    if (r >= produced_rows.size() || !produced_rows[r]) continue;
    if (r >= normalize_rows.size() || !normalize_rows[r]) continue;
    for (std::size_t j = 0; j < lout.cols(); ++j) {
      lout[r][j] *= scale;
    }
  }
}

template <typename LLR>
void postprocess_decoder_result(
    chase::DecoderCoreResult<typename LinMatrixAdapter<LLR>::core_type>* decoder_res,
    bool normalize_extrinsic,
    bool use_hard_decode,
    const newcode::Params& p,
    const std::vector<bool>* normalize_rows = nullptr) {
  // decoder core 输出的统一后处理入口。
  //
  // 处理顺序（与旧逻辑保持一致）：
  // 1. （可选）extrinsic 归一化：
  //    - 仅在 soft 路径（use_hard_decode=false）且 normalize_extrinsic=true 时启用
  //    - normalize_rows==nullptr：对所有 produced rows 做全局归一化（旧行为）
  //    - normalize_rows!=nullptr：只对指定行集合做全局归一化（方案三兼容开关）
  // 2. 对 produced rows 做统一 ALPHA 缩放（p.ALPHA）
  // 3. 对 produced rows 做输出量化（按 LLR 类型走 qfloat / int / float 的 quantize）
  //
  // 注意：
  // - 这里的归一化是“全局 scale”，不是逐行 scale；
  // - normalize_rows 只影响第 1 步（normalize_extrinsic），不会影响 ALPHA/量化步骤。
  if (normalize_extrinsic && !use_hard_decode) {
    if (normalize_rows) {
      normalize_selected_rows(decoder_res->lout,
                              decoder_res->produced_rows,
                              *normalize_rows,
                              p.beta);
    } else {
      normalize_extrinsic_lout(decoder_res->lout, decoder_res->produced_rows, p.beta);
    }
  }

  const std::size_t rcnt = decoder_res->lout.rows();
  const std::size_t ccnt = decoder_res->lout.cols();
  for (std::size_t r = 0; r < rcnt; ++r) {
    if (!decoder_res->produced_rows[r]) continue;
    for (std::size_t j = 0; j < ccnt; ++j) {
      decoder_res->lout[r][j] *= p.ALPHA;
    }
  }

  for (std::size_t r = 0; r < rcnt; ++r) {
    if (!decoder_res->produced_rows[r]) continue;
    for (std::size_t j = 0; j < ccnt; ++j) {
      decoder_res->lout[r][j] =
          ExtrinsicQuantizer<LLR>::quantize(decoder_res->lout[r][j]);
    }
  }
}

inline std::string sanitize_label(const std::string& label) {
  // 输入:
  // - label: 追踪目标的原始标签。
  // 输出:
  // - 适合作为文件名的一段安全字符串。
  // 用途:
  // - Chase 调试导出 CSV 时，避免标签中出现路径或特殊字符。
  std::string safe = label.empty() ? std::string("target") : label;
  for (char& ch : safe) {
    if (!std::isalnum(static_cast<unsigned char>(ch)) && ch != '-') {
      ch = '_';
    }
  }
  return safe;
}

inline std::string quantized_code_csv_cell(float value,
                                           const newcode::Params& params) {
  // 输入:
  // - value: 已经经过公共量化/裁剪后的 extrinsic 浮点值。
  // - params: 当前 tile 参数，提供 LLR 位宽和 clip。
  // 输出:
  // - 若当前是 qfloat 路径，则返回对应的内部 code；否则返回空串。
  // 用途:
  // - 在 target_*.csv 里补充“这一轮实际写回矩阵前，对应的量化码值”。
  const float clip = params.LLR_CLIP;
  switch (params.LLR_BITS) {
    case 2:  return std::to_string(qfloat::qfloat<2>::from_float(value, clip).code());
    case 3:  return std::to_string(qfloat::qfloat<3>::from_float(value, clip).code());
    case 4:  return std::to_string(qfloat::qfloat<4>::from_float(value, clip).code());
    case 5:  return std::to_string(qfloat::qfloat<5>::from_float(value, clip).code());
    case 6:  return std::to_string(qfloat::qfloat<6>::from_float(value, clip).code());
    case 7:  return std::to_string(qfloat::qfloat<7>::from_float(value, clip).code());
    case 8:  return std::to_string(qfloat::qfloat<8>::from_float(value, clip).code());
    case 9:  return std::to_string(qfloat::qfloat<9>::from_float(value, clip).code());
    case 10: return std::to_string(qfloat::qfloat<10>::from_float(value, clip).code());
    case 11: return std::to_string(qfloat::qfloat<11>::from_float(value, clip).code());
    case 12: return std::to_string(qfloat::qfloat<12>::from_float(value, clip).code());
    case 13: return std::to_string(qfloat::qfloat<13>::from_float(value, clip).code());
    case 14: return std::to_string(qfloat::qfloat<14>::from_float(value, clip).code());
    case 15: return std::to_string(qfloat::qfloat<15>::from_float(value, clip).code());
    default: return std::string();
  }
}

inline std::string csv_float_cell(bool has_value, float value) {
  return has_value ? std::to_string(value) : std::string();
}

inline std::string csv_int_cell(bool has_value, int value) {
  return has_value ? std::to_string(value) : std::string();
}

inline void log_target_history(const newcode::Params& params,
                               const matrix::Matrix<float>& lout,
                               const std::vector<bool>& produced_rows) {
  // 输入:
  // - params: 调试配置来源。
  // - lout: 当前 tile core 输出的外信息矩阵。
  // - produced_rows: 每个 decoder row 是否真的产出结果。
  // 输出:
  // - 无返回值；满足条件时向 CSV 文件追加一行追踪记录。
  // 用途:
  // - 将指定坐标/目标比特在每次 Chase 调用中的 extrinsic 变化落盘，便于离线分析。
  const auto& trace = params.debug_trace;
  if (!trace.enable || !trace.dump_chase_csv) return;
  if (trace.active_chase_entries.empty()) return;
  if (trace.chase_tile_index < 0 || trace.chase_invocation < 0) return;
  namespace fs = std::filesystem;
  const fs::path dir = trace.chase_csv_dir.empty()
                           ? fs::path("data/chase_csv")
                           : fs::path(trace.chase_csv_dir);
  std::error_code ec;
  fs::create_directories(dir, ec);
  const std::string csv_header =
      "invocation,tile_index,global_row,global_col,lin_row_index,lin_k,channel_llr_float,channel_llr_code,extrinsic,extrinsic_code,expected_bit";
  for (const auto& entry : trace.active_chase_entries) {
    if (entry.row_index < 0 || entry.k < 0) continue;
    if (entry.row_index >= static_cast<int>(lout.rows())) continue;
    if (static_cast<size_t>(entry.k) >= lout.cols()) continue;
    if (entry.row_index >= static_cast<int>(produced_rows.size())) continue;
    if (!produced_rows[static_cast<size_t>(entry.row_index)]) continue;
    float extrinsic = lout[static_cast<size_t>(entry.row_index)]
                         [static_cast<size_t>(entry.k)];
    std::string label = entry.label;
    if (label.empty()) {
      if (entry.bit_index >= 0) {
        label = "bit" + std::to_string(entry.bit_index);
      } else {
        label = "row" + std::to_string(entry.global_row) + "_col" +
                std::to_string(entry.global_col);
      }
    }
    const std::string safe_label = sanitize_label(label);
    const fs::path file =
        dir / ("target_" + safe_label + ".csv");
    const bool existed = fs::exists(file);
    bool append_mode = existed;
    if (existed) {
      std::ifstream in(file);
      std::string first_line;
      if (!std::getline(in, first_line) || first_line != csv_header) {
        // 旧格式 target.csv 缺少新列时，直接用新表头重建，避免新旧列数混杂。
        append_mode = false;
      }
    }
    std::ofstream out(file, append_mode ? std::ios::app : std::ios::trunc);
    if (!out) continue;
    if (!append_mode) {
      out << csv_header << '\n';
    }
    out << trace.chase_invocation << ',' << trace.chase_tile_index << ','
        << entry.global_row << ',' << entry.global_col << ','
        << entry.row_index << ',' << entry.k << ','
        << csv_float_cell(entry.has_channel_llr_float, entry.channel_llr_float) << ','
        << csv_int_cell(entry.has_channel_llr_code, entry.channel_llr_code) << ','
        << extrinsic << ',' << quantized_code_csv_cell(extrinsic, params)
        << ',' << entry.expected_bit << '\n';
  }
}

template <typename LLR>
chase::DecoderCoreResult<typename LinMatrixAdapter<LLR>::core_type>
decode_tile(const TilePrepared<LLR>& prep,
            bool use_hard_decode,
            bool normalize_extrinsic,
            const newcode::Params& p,
            const std::vector<bool>* early_stop_row_flags,
            const std::vector<uint8_t>* mux_state,
            CoreFn<typename LinMatrixAdapter<LLR>::core_type> core_fn)
{
  // 输入:
  // - prep: 已经排布好的 tile 输入，包括 lin/lch 矩阵和 trace 信息。
  // - use_hard_decode: 是否切换到硬判回退。
  // - normalize_extrinsic: 是否对 soft 路径外信息做归一化。
  // - p: 当前 tile 参数。
  // - early_stop_row_flags: 每个 decoder row 的 early-stop 标志。
  // - mux_state: 当前 tile 的 mux 分配状态。
  // - core_fn: 实际执行 Chase/硬判回退的 core。
  // 输出:
  // - DecoderCoreResult，包含 lout 和 produced_rows。
  // 用途:
  // - 这是 tile 级真正调用 decoder core 的地方，并负责对 core 输出做后处理。
  auto decoder_res = run_decoder_core(prep.lin_matrix,
                                      prep.lch_matrix,
                                      use_hard_decode,
                                      prep.params_for_core,
                                      early_stop_row_flags,
                                      mux_state,
                                      core_fn);

  postprocess_decoder_result<LLR>(&decoder_res,
                                  normalize_extrinsic,
                                  use_hard_decode,
                                  p);

  log_target_history(prep.params_for_core, decoder_res.lout,
                     decoder_res.produced_rows);

  return decoder_res;
}

template <typename LLR>
void materialize_early_stop_rows(
    const TilePrepared<LLR>& prep,
    const TileDispatchPlan<typename TilePrepared<LLR>::CoreLLR>& plan,
    chase::DecoderCoreResult<typename LinMatrixAdapter<LLR>::core_type>* result) {
  using CoreLLR = typename TilePrepared<LLR>::CoreLLR;
  const std::size_t cols = prep.lin_matrix.cols();
  // 早停行不再调用 Chase，而是逐行执行既有 early-stop action，
  // 然后把生成的 y2 直接写进完整结果缓冲。
  for (std::size_t row = 0; row < plan.rows.size(); ++row) {
    if (plan.rows[row].tag != RowDispatchTag::EarlyStopAction) {
      continue;
    }
    std::array<CoreLLR, newcode::Params::BCH_N> lin_vec{};
    std::array<CoreLLR, newcode::Params::BCH_N> lch_vec{};
    std::array<float, newcode::Params::BCH_N> y2{};
    load_row_vectors(prep, row, &lin_vec, &lch_vec);
    const bool produced = newcode::apply_row_early_stop_action(
        lin_vec.data(), lch_vec.data(), y2.data(), prep.params_for_core);
    if (!produced) {
      continue;
    }
    result->produced_rows[row] = true;
    for (std::size_t col = 0; col < cols; ++col) {
      result->lout[row][col] = y2[col];
    }
  }
}

template <typename LLR>
void materialize_hard_finish_rows(
    const TileDispatchPlan<typename TilePrepared<LLR>::CoreLLR>& plan,
    chase::DecoderCoreResult<typename LinMatrixAdapter<LLR>::core_type>* result) {
  const std::size_t cols = result->lout.cols();
  // HardFinish 行的 lout 在 hybrid prepass 阶段已经算完了；
  // decode 阶段只负责把缓存结果灌回完整结果矩阵。
  //
  // 注意：
  // - 这里不会调用 decoder core（不会跑 Chase），只是“回填”；
  // - `hard_finish_valid[row]` 是必要的保护：
  //   plan.hard_finish_lout 是按完整行域分配，但只有少数行真的写入过有效 y2。
  for (std::size_t row = 0; row < plan.rows.size(); ++row) {
    if (plan.rows[row].tag != RowDispatchTag::HardFinish ||
        row >= plan.hard_finish_valid.size() || !plan.hard_finish_valid[row]) {
      continue;
    }
    result->produced_rows[row] = true;
    for (std::size_t col = 0; col < cols; ++col) {
      result->lout[row][col] = plan.hard_finish_lout[row][col];
    }
  }
}

template <typename LLR>
void materialize_classified_clean_rows(
    const TilePrepared<LLR>& prep,
    const TileDispatchPlan<typename TilePrepared<LLR>::CoreLLR>& plan,
    chase::DecoderCoreResult<typename LinMatrixAdapter<LLR>::core_type>* result) {
  using CoreLLR = typename TilePrepared<LLR>::CoreLLR;
  const std::size_t cols = prep.lin_matrix.cols();
  for (std::size_t row = 0; row < plan.rows.size(); ++row) {
    if (plan.rows[row].tag != RowDispatchTag::HardFinish ||
        plan.rows[row].hybrid_class != HybridRowClass::Clean ||
        plan.rows[row].scheduled_for_hard ||
        (row < plan.hard_finish_valid.size() && plan.hard_finish_valid[row])) {
      continue;
    }
    std::array<CoreLLR, newcode::Params::BCH_N> lin_vec{};
    std::array<CoreLLR, newcode::Params::BCH_N> lch_vec{};
    std::array<float, newcode::Params::BCH_N> y2{};
    load_row_vectors(prep, row, &lin_vec, &lch_vec);
    const bool produced = newcode::apply_row_early_stop_action(
        lin_vec.data(), lch_vec.data(), y2.data(), prep.params_for_core);
    if (!produced) {
      continue;
    }
    result->produced_rows[row] = true;
    for (std::size_t col = 0; col < cols; ++col) {
      result->lout[row][col] = y2[col];
    }
  }
}

inline void throw_hard_executor_failure(std::size_t row,
                                        const char* reason) {
  std::ostringstream oss;
  oss << "hybrid hard executor failed at row " << row << ": " << reason;
  throw std::runtime_error(oss.str());
}

template <typename CoreLLR>
void execute_hybrid_hard_class(
    HybridRowClass hard_class,
    const std::array<CoreLLR, newcode::Params::BCH_N>& lin_vec,
    const newcode::Params& p,
    std::array<float, newcode::Params::BCH_N>* y2) {
  auto cw = hard_decision_bits_256(lin_vec);

  switch (hard_class) {
    case HybridRowClass::ParityOnly:
      cw[newcode::Params::BCH_OVERALL_IDX] ^= 1u;
      break;
    case HybridRowClass::OneMain:
    case HybridRowClass::OneMainPlusParity: {
      const auto syndromes = bch::bch_255_239_syndromes_1_4_cw_255(cw.data());
      const int pos = hybrid_fast_gf_log(syndromes[0]);
      if (pos < 0 || pos >= static_cast<int>(newcode::Params::BCH_OVERALL_IDX)) {
        throw std::runtime_error("invalid one-main error position");
      }
      cw[static_cast<std::size_t>(pos)] ^= 1u;
      if (hard_class == HybridRowClass::OneMainPlusParity) {
        cw[newcode::Params::BCH_OVERALL_IDX] ^= 1u;
      }
      break;
    }
    case HybridRowClass::TwoMain: {
      std::array<uint8_t, newcode::Params::BCH_N - 1> decoded{};
      int corrected_errors = 0;
      if (!bch::bch_255_239_decode_hiho_cw_255(cw.data(),
                                               decoded.data(),
                                               &corrected_errors)) {
        throw std::runtime_error("BCH t=2 decode failed");
      }
      if (corrected_errors != 2) {
        throw std::runtime_error("BCH t=2 decode did not correct exactly two errors");
      }
      for (std::size_t i = 0; i < decoded.size(); ++i) {
        cw[i] = decoded[i];
      }
      recompute_overall_parity(&cw);
      break;
    }
    default:
      throw std::runtime_error("unsupported hard class");
  }

  if (!hard_word_valid_256(cw)) {
    throw std::runtime_error("corrected hard word is not a valid BCH+overall codeword");
  }
  materialize_hard_finish_lout(cw, lin_vec, p, y2);
}

template <typename LLR>
void materialize_scheduled_hard_rows(
    const TilePrepared<LLR>& prep,
    const TileDispatchPlan<typename TilePrepared<LLR>::CoreLLR>& plan,
    chase::DecoderCoreResult<typename LinMatrixAdapter<LLR>::core_type>* result) {
  using CoreLLR = typename TilePrepared<LLR>::CoreLLR;
  const std::size_t cols = prep.lin_matrix.cols();
  for (int row_value : plan.hard_scheduled_rows) {
    if (row_value < 0) {
      continue;
    }
    const auto row = static_cast<std::size_t>(row_value);
    if (row >= plan.rows.size() ||
        plan.rows[row].tag != RowDispatchTag::HardFinish ||
        !plan.rows[row].scheduled_for_hard) {
      continue;
    }

    std::array<CoreLLR, newcode::Params::BCH_N> lin_vec{};
    std::array<CoreLLR, newcode::Params::BCH_N> lch_vec{};
    std::array<float, newcode::Params::BCH_N> y2{};
    load_row_vectors(prep, row, &lin_vec, &lch_vec);
    (void)lch_vec;
    try {
      execute_hybrid_hard_class(plan.rows[row].hybrid_class,
                                lin_vec,
                                prep.params_for_core,
                                &y2);
    } catch (const std::exception& ex) {
      throw_hard_executor_failure(row, ex.what());
    }
    result->produced_rows[row] = true;
    for (std::size_t col = 0; col < cols; ++col) {
      result->lout[row][col] = y2[col];
    }
  }
}

template <typename LLR>
std::vector<uint8_t> build_soft_only_mux_state(
    const TileDispatchPlan<typename TilePrepared<LLR>::CoreLLR>& plan) {
  // 不再把 soft 行压缩成 compact batch。
  // 这里保留完整行域，只用 mux_state 告诉 core 哪些行真正要跑 Chase：
  // - scheduled SoftDecode 行：NeedSiso
  // - EarlyStopAction / HardFinish / Unscheduled：Unscheduled，占位但不产出 soft 结果
  std::vector<uint8_t> soft_mux_state(
      plan.rows.size(),
      static_cast<uint8_t>(newcode::mux::StateTag::Unscheduled));
  for (int row : plan.soft_scheduled_rows) {
    if (row < 0) {
      continue;
    }
    const auto row_index = static_cast<std::size_t>(row);
    if (row_index >= plan.rows.size()) {
      continue;
    }
    if (plan.rows[row_index].tag != RowDispatchTag::SoftDecode ||
        !plan.rows[row_index].scheduled_for_soft) {
      continue;
    }
    soft_mux_state[row_index] =
        static_cast<uint8_t>(newcode::mux::StateTag::NeedSiso);
  }
  return soft_mux_state;
}

template <typename LLR>
void merge_soft_decoder_result(
    const TileDispatchPlan<typename TilePrepared<LLR>::CoreLLR>& plan,
    const chase::DecoderCoreResult<typename LinMatrixAdapter<LLR>::core_type>& soft_result,
    chase::DecoderCoreResult<typename LinMatrixAdapter<LLR>::core_type>* full_result) {
  const std::size_t cols = full_result->lout.cols();
  // soft_result 已经是完整行域结果；这里只回填计划中真正 scheduled 的 soft 行，
  // 避免覆盖前面已经 materialize 好的 EarlyStopAction / HardFinish 行。
  for (int row : plan.soft_scheduled_rows) {
    if (row < 0) {
      continue;
    }
    const auto row_index = static_cast<std::size_t>(row);
    if (row_index >= soft_result.produced_rows.size() ||
        !soft_result.produced_rows[row_index]) {
      continue;
    }
    if (row_index >= full_result->produced_rows.size()) {
      continue;
    }
    full_result->produced_rows[row_index] = true;
    for (std::size_t col = 0; col < cols; ++col) {
      full_result->lout[row_index][col] = soft_result.lout[row_index][col];
    }
  }
}

template <typename LLR>
chase::DecoderCoreResult<typename LinMatrixAdapter<LLR>::core_type>
decode_tile_with_plan(
    const TilePrepared<LLR>& prep,
    const TileDispatchPlan<typename TilePrepared<LLR>::CoreLLR>& plan,
    bool normalize_extrinsic,
    CoreFn<typename LinMatrixAdapter<LLR>::core_type> core_fn) {
  using CoreLLR = typename TilePrepared<LLR>::CoreLLR;
  const std::size_t rows = prep.lin_matrix.rows();
  const std::size_t cols = prep.lin_matrix.cols();

  chase::DecoderCoreResult<CoreLLR> merged_result{
      matrix::Matrix<float>(rows, cols),
      std::vector<bool>(rows, false)};

  // 先把不需要 Chase 的两类行落到完整结果里：
  // - EarlyStopAction
  // - classified Clean（复用 early-stop action）
  // - HardFinish（cached 或 hard executor）
  materialize_early_stop_rows(prep, plan, &merged_result);
  materialize_classified_clean_rows(prep, plan, &merged_result);
  materialize_hard_finish_rows<LLR>(plan, &merged_result);
  materialize_scheduled_hard_rows(prep, plan, &merged_result);

  if (!plan.soft_scheduled_rows.empty()) {
    // 保留完整 tile 行域，把是否跑 Chase 交给 soft-only mux_state 控制。
    // 这样 debug trace / expected_bits / 行号 lookup 都继续使用原始 full row 坐标。
    const auto soft_mux_state = build_soft_only_mux_state<LLR>(plan);
    auto soft_result = run_decoder_core(prep.lin_matrix,
                                        prep.lch_matrix,
                                        false,
                                        prep.params_for_core,
                                        nullptr,
                                        &soft_mux_state,
                                        core_fn);
    merge_soft_decoder_result<LLR>(plan, soft_result, &merged_result);
  }

  std::vector<bool> normalize_rows;
  const std::vector<bool>* normalize_rows_ptr = nullptr;
  if (prep.params_for_core.HYBRID_NORMALIZE_SOFT_ONLY) {
    // 兼容性开关：
    // true 只归一化 soft rows；
    // false 保持旧行为，对所有 produced rows 做统一后处理。
    normalize_rows.assign(rows, false);
    for (int row : plan.soft_scheduled_rows) {
      if (row >= 0 && static_cast<std::size_t>(row) < normalize_rows.size()) {
        normalize_rows[static_cast<std::size_t>(row)] = true;
      }
    }
    normalize_rows_ptr = &normalize_rows;
  }

  postprocess_decoder_result<LLR>(&merged_result,
                                  normalize_extrinsic,
                                  false,
                                  prep.params_for_core,
                                  normalize_rows_ptr);
  // 这里记录的是“完整合并后的结果”，这样 trace 看见的是最终语义而不是中间 batch。
  log_target_history(prep.params_for_core,
                     merged_result.lout,
                     merged_result.produced_rows);
  return merged_result;
}

} // 命名空间 detail
} // 命名空间 newcode
