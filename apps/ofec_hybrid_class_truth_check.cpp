// 只读真值核验：复用工程中的 BCH、S0/S1/S3 分类器与 hard-finish 模块，
// 检查“分类名”是否等价于“相对本帧实际发送 BCH 码字的错误数”。
//
// 这不是 BER 仿真，不进入 Level 5/6 调度，也不会改动正式解码流程。
// 它构造一个已知发送 eBCH(256,239) 码字，在其上注入确定性的错误图样，
// 然后记录：输入真实错误数 -> 分类 -> HISO hard-finish 输出真实错误数。

#include <algorithm>
#include <array>
#include <cstdint>
#include <iomanip>
#include <iostream>
#include <map>
#include <random>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

#include "newcode/common/bch/bch_255_239.hpp"
#include "newcode/ofec/hybrid/hybrid_classifier.hpp"
#include "newcode/params.hpp"

namespace {

using Codeword = std::array<uint8_t, newcode::Params::BCH_N>;
using LinVec = std::array<float, newcode::Params::BCH_N>;
using RowClass = newcode::detail::HybridRowClass;

constexpr float kInputMagnitude = 12.0f;
constexpr int kRandomSamplesPerWeight = 25000;
constexpr int kRandomMaxErrorWeight = 32;
constexpr uint32_t kPatternSeed = 20260826u;

struct Counts {
  uint64_t samples = 0;
  uint64_t classify_ok = 0;
  uint64_t hard_finish_ok = 0;
  uint64_t output_matches_tx = 0;
  uint64_t output_differs_from_tx = 0;
  uint64_t two_bch_ok = 0;
  uint64_t two_bch_output_matches_executor = 0;
  // key = HISO hard-finish 输出相对真实 tx 的 Hamming distance。
  // TwoMain 时该输出与直接调用既有 BCH t=2 模块的输出逐 bit 相同。
  std::map<std::size_t, uint64_t> output_error_weight_histogram;
};

struct Example {
  bool present = false;
  int input_weight = 0;
  RowClass cls = RowClass::HardFail;
  std::vector<int> injected_positions;
  std::vector<int> output_error_positions;
  int bch_corrected_errors = -1;
  std::vector<int> bch_output_error_positions;
};

std::string class_name(RowClass cls) {
  switch (cls) {
    case RowClass::Clean: return "Clean";
    case RowClass::ParityOnly: return "ParityOnly";
    case RowClass::OneMain: return "OneMain";
    case RowClass::OneMainPlusParity: return "OneMainPlusParity";
    case RowClass::TwoMain: return "TwoMain";
    case RowClass::HardFail: return "HardFail";
    case RowClass::None: return "None";
    case RowClass::BchHardDecoded: return "BchHardDecoded";
    case RowClass::Suspicious: return "Suspicious";
  }
  return "Unknown";
}

std::string join_positions(const std::vector<int>& positions) {
  std::ostringstream oss;
  oss << '[';
  for (std::size_t i = 0; i < positions.size(); ++i) {
    if (i != 0) oss << ',';
    oss << positions[i];
  }
  oss << ']';
  return oss.str();
}

std::string format_histogram(const std::map<std::size_t, uint64_t>& histogram) {
  if (histogram.empty()) return "-";
  std::ostringstream oss;
  bool first = true;
  for (const auto& [error_weight, count] : histogram) {
    if (!first) oss << ' ';
    first = false;
    // 例如 6:123 表示：输出相对 tx 有 6 错的样本共有 123 个。
    oss << error_weight << ':' << count;
  }
  return oss.str();
}

std::vector<int> differing_positions(const Codeword& a, const Codeword& b) {
  std::vector<int> positions;
  for (std::size_t i = 0; i < a.size(); ++i) {
    if ((a[i] & 1u) != (b[i] & 1u)) positions.push_back(static_cast<int>(i));
  }
  return positions;
}

LinVec make_lin(const Codeword& hard_word) {
  LinVec lin{};
  for (std::size_t i = 0; i < hard_word.size(); ++i) {
    lin[i] = hard_word[i] ? -kInputMagnitude : kInputMagnitude;
  }
  return lin;
}

Codeword recover_hard_finish_word(const LinVec& lin,
                                  const std::array<float, newcode::Params::BCH_N>& y2) {
  Codeword out{};
  for (std::size_t i = 0; i < out.size(); ++i) {
    // materialize_hard_finish_lout() 定义 y2 = lpost - lin；
    // 这里按同一语义恢复 lpost 的硬判，避免引入另一个“自写的 BCH 译码器”。
    out[i] = (lin[i] + y2[i] < 0.0f) ? 1u : 0u;
  }
  return out;
}

Codeword make_tx_codeword() {
  // 固定但非全零的信息，便于证明测试不依赖“零码字”的特殊性。
  std::vector<uint8_t> info(newcode::Params::BCH_K, 0u);
  std::mt19937 rng(20260825u);
  for (auto& bit : info) bit = static_cast<uint8_t>(rng() & 1u);
  const Codeword tx = bch::bch_255_239_encode(info);
  if (!newcode::detail::hard_word_valid_256(tx)) {
    throw std::runtime_error("encoder did not produce a valid eBCH(256,239) word");
  }
  return tx;
}

std::vector<int> sample_positions(int weight, std::mt19937* rng) {
  std::array<int, newcode::Params::BCH_N> all{};
  for (std::size_t i = 0; i < all.size(); ++i) all[i] = static_cast<int>(i);
  std::shuffle(all.begin(), all.end(), *rng);
  std::vector<int> positions(all.begin(), all.begin() + weight);
  std::sort(positions.begin(), positions.end());
  return positions;
}

void print_example(const std::string& title, const Example& ex) {
  if (!ex.present) {
    std::cout << title << ": 在本次扫描范围内未出现。\n";
    return;
  }
  std::cout << title << ":\n"
            << "  输入相对真实发送码字的错误数 = " << ex.input_weight << "\n"
            << "  注入错误位置 = " << join_positions(ex.injected_positions) << "\n"
            << "  分类器结果 = " << class_name(ex.cls) << "\n"
            << "  HISO hard-finish 输出相对真实发送码字的错误数 = "
            << ex.output_error_positions.size() << "\n"
            << "  HISO 输出错误位置 = " << join_positions(ex.output_error_positions) << "\n";
  if (ex.cls == RowClass::TwoMain) {
    std::cout << "  现有 BCH t=2 模块报告 corrected_errors = " << ex.bch_corrected_errors << "\n"
              << "  BCH 输出相对真实发送码字的错误位置 = "
              << join_positions(ex.bch_output_error_positions) << "\n";
  }
}

}  // namespace

int main() {
  try {
    const Codeword tx = make_tx_codeword();
    newcode::Params p;
    p.HYBRID_CLASSIFIER_MODE =
        newcode::HybridClassifierMode::FriendS1S3WithS0Classifier;
    p.HYBRID_USE_FAST_CLASSIFIER = true;
    p.HYBRID_HARD_LLR_MAG = 100.0f;

    std::map<std::pair<int, RowClass>, Counts> table;
    Example first_wrong_onemain;
    Example first_wrong_twomain;
    std::mt19937 rng(kPatternSeed);

    std::cout << "ofec_hybrid_class_truth_check\n"
              << "classifier=FriendS1S3WithS0Classifier\n"
              << "tx=eBCH(256,239), tx_valid="
              << (newcode::detail::hard_word_valid_256(tx) ? "yes" : "no") << "\n"
              << "scan: error weight 0.." << kRandomMaxErrorWeight
              << ", deterministic random patterns per nonzero weight="
              << kRandomSamplesPerWeight << "\n\n";

    for (int weight = 0; weight <= kRandomMaxErrorWeight; ++weight) {
      const int samples = (weight == 0) ? 1 : kRandomSamplesPerWeight;
      for (int sample = 0; sample < samples; ++sample) {
        const std::vector<int> injected = sample_positions(weight, &rng);
        Codeword received = tx;
        for (int pos : injected) received[static_cast<std::size_t>(pos)] ^= 1u;
        const LinVec lin = make_lin(received);

        RowClass cls = RowClass::HardFail;
        const bool classify_ok = newcode::detail::run_selected_hybrid_classifier_classify_only(
            lin, p, &cls);
        Counts& counts = table[{weight, cls}];
        ++counts.samples;
        if (classify_ok) ++counts.classify_ok;

        std::array<float, newcode::Params::BCH_N> y2{};
        RowClass finish_cls = RowClass::HardFail;
        const bool hard_finish_ok = newcode::detail::run_selected_hybrid_classifier_hard_finish(
            lin, &y2, p, &finish_cls);
        if (finish_cls != cls) {
          throw std::runtime_error("classify-only and hard-finish class disagree");
        }
        if (!hard_finish_ok) continue;

        ++counts.hard_finish_ok;
        const Codeword output = recover_hard_finish_word(lin, y2);
        const std::vector<int> output_errors = differing_positions(output, tx);
        ++counts.output_error_weight_histogram[output_errors.size()];
        if (output_errors.empty()) {
          ++counts.output_matches_tx;
        } else {
          ++counts.output_differs_from_tx;
        }

        Example current;
        current.present = true;
        current.input_weight = weight;
        current.cls = cls;
        current.injected_positions = injected;
        current.output_error_positions = output_errors;

        if (cls == RowClass::TwoMain) {
          std::array<uint8_t, newcode::Params::BCH_N - 1> bch_out{};
          int corrected_errors = -1;
          const bool bch_ok = bch::bch_255_239_decode_hiho_cw_255(
              received.data(), bch_out.data(), &corrected_errors);
          if (!bch_ok || corrected_errors != 2) {
            throw std::runtime_error("TwoMain did not reproduce BCH t=2 success/2-correction");
          }
          Codeword direct_output = received;
          for (std::size_t i = 0; i < bch_out.size(); ++i) direct_output[i] = bch_out[i];
          newcode::detail::recompute_overall_parity(&direct_output);
          if (direct_output != output) {
            throw std::runtime_error("TwoMain BCH output differs from HISO hard-finish codeword");
          }
          ++counts.two_bch_ok;
          ++counts.two_bch_output_matches_executor;
          current.bch_corrected_errors = corrected_errors;
          current.bch_output_error_positions = differing_positions(direct_output, tx);
        }

        if (cls == RowClass::OneMain && !output_errors.empty() && !first_wrong_onemain.present) {
          first_wrong_onemain = current;
        }
        if (cls == RowClass::TwoMain && !output_errors.empty() && !first_wrong_twomain.present) {
          first_wrong_twomain = current;
        }
      }
    }

    std::cout << std::left
              << std::setw(10) << "in_err"
              << std::setw(22) << "class"
              << std::setw(12) << "samples"
              << std::setw(14) << "hard_finish"
              << std::setw(14) << "out_eq_tx"
              << std::setw(16) << "out_ne_tx"
              << std::setw(12) << "bch_t2_ok"
              << "output_err_weight:count\n";
    for (const auto& [key, counts] : table) {
      std::cout << std::left
                << std::setw(10) << key.first
                << std::setw(22) << class_name(key.second)
                << std::setw(12) << counts.samples
                << std::setw(14) << counts.hard_finish_ok
                << std::setw(14) << counts.output_matches_tx
                << std::setw(16) << counts.output_differs_from_tx
                << std::setw(12) << counts.two_bch_ok
                << format_histogram(counts.output_error_weight_histogram)
                << '\n';
    }
    std::cout << '\n';
    print_example("OneMain 误纠反例（若存在）", first_wrong_onemain);
    print_example("TwoMain 误纠反例（若存在）", first_wrong_twomain);

    if (!first_wrong_onemain.present) {
      std::cerr << "[FAIL] scan did not find a OneMain miscorrection; enlarge scan before drawing a conclusion\n";
      return 2;
    }
    if (!first_wrong_twomain.present) {
      std::cerr << "[FAIL] scan did not find a TwoMain miscorrection; enlarge scan before drawing a conclusion\n";
      return 2;
    }
    std::cout << "\n[PASS] Found tx-relative miscorrection examples for both OneMain and TwoMain.\n";
  } catch (const std::exception& ex) {
    std::cerr << "[FAIL] hybrid class truth check: " << ex.what() << '\n';
    return 1;
  }
  return 0;
}
