#pragma once

namespace newcode {
namespace detail {

inline const HybridFastGf256Tables& hybrid_fast_gf256_tables() {
  static const HybridFastGf256Tables tables = [] {
    HybridFastGf256Tables t{};
    t.index_of.fill(-1);
    uint16_t a = 1;
    for (int i = 0; i < 255; ++i) {
      t.alpha_to[static_cast<std::size_t>(i)] = static_cast<uint8_t>(a);
      t.index_of[static_cast<std::size_t>(t.alpha_to[static_cast<std::size_t>(i)])] =
          static_cast<int16_t>(i);
      a <<= 1;
      if (a & 0x100u) {
        a ^= 0x11Du;
      }
    }
    return t;
  }();
  return tables;
}

inline uint8_t hybrid_fast_gf_mul(uint8_t x, uint8_t y) {
  if (x == 0u || y == 0u) {
    return 0u;
  }
  const auto& tables = hybrid_fast_gf256_tables();
  const int lx = tables.index_of[static_cast<std::size_t>(x)];
  const int ly = tables.index_of[static_cast<std::size_t>(y)];
  return tables.alpha_to[static_cast<std::size_t>((lx + ly) % 255)];
}

inline uint8_t hybrid_fast_gf_div(uint8_t x, uint8_t y) {
  if (x == 0u) {
    return 0u;
  }
  if (y == 0u) {
    throw std::runtime_error("hybrid_fast_gf_div: division by zero");
  }
  const auto& tables = hybrid_fast_gf256_tables();
  const int lx = tables.index_of[static_cast<std::size_t>(x)];
  const int ly = tables.index_of[static_cast<std::size_t>(y)];
  int diff = lx - ly;
  if (diff < 0) {
    diff += 255;
  }
  return tables.alpha_to[static_cast<std::size_t>(diff)];
}

inline uint8_t hybrid_fast_gf_cube(uint8_t x) {
  return hybrid_fast_gf_mul(hybrid_fast_gf_mul(x, x), x);
}

inline int hybrid_fast_gf_log(uint8_t x) {
  if (x == 0u) {
    return -1;
  }
  return hybrid_fast_gf256_tables().index_of[static_cast<std::size_t>(x)];
}

inline uint8_t hybrid_fast_trace_to_gf2(uint8_t x) {
  uint8_t acc = x;
  uint8_t cur = x;
  for (int i = 1; i < 8; ++i) {
    cur = hybrid_fast_gf_mul(cur, cur);
    acc ^= cur;
  }
  return static_cast<uint8_t>(acc & 0x1u);
}

inline uint8_t overall_parity_syndrome_256(
    const std::array<uint8_t, newcode::Params::BCH_N>& cw) {
  uint8_t parity = 0u;
  for (uint8_t bit : cw) {
    parity ^= static_cast<uint8_t>(bit & 1u);
  }
  return parity;
}

inline void recompute_overall_parity(
    std::array<uint8_t, newcode::Params::BCH_N>* cw) {
  uint8_t parity = 0u;
  for (std::size_t i = 0; i < newcode::Params::BCH_OVERALL_IDX; ++i) {
    parity ^= static_cast<uint8_t>((*cw)[i] & 1u);
  }
  (*cw)[newcode::Params::BCH_OVERALL_IDX] = parity;
}

template <typename CoreLLR>
inline std::array<uint8_t, newcode::Params::BCH_N> hard_decision_bits_256(
    const std::array<CoreLLR, newcode::Params::BCH_N>& lin_vec) {
  std::array<uint8_t, newcode::Params::BCH_N> hard_bits{};
  for (std::size_t i = 0; i < hard_bits.size(); ++i) {
    hard_bits[i] = (qfloat::llr_to_float(lin_vec[i]) < 0.0f) ? 1u : 0u;
  }
  return hard_bits;
}

template <typename CoreLLR>
inline void materialize_hard_finish_lout(
    const std::array<uint8_t, newcode::Params::BCH_N>& cw,
    const std::array<CoreLLR, newcode::Params::BCH_N>& lin_vec,
    const newcode::Params& p,
    std::array<float, newcode::Params::BCH_N>* y2) {
  // hybrid prepass 的 hard-finish 输出单独使用 HYBRID_HARD_LLR_MAG，
  // 这样可以和 legacy hard-decode 的 HARD_LLR_MAG 分开扫参与调节。
  const float hard_mag = std::fabs(p.HYBRID_HARD_LLR_MAG);
  for (std::size_t i = 0; i < cw.size(); ++i) {
    const float sign = cw[i] ? -1.0f : 1.0f;
    const float lpost = sign * hard_mag;
    const float lin = qfloat::llr_to_float(lin_vec[i]);
    (*y2)[i] = lpost - lin;
  }
}

inline bool hard_word_valid_256(
    const std::array<uint8_t, newcode::Params::BCH_N>& cw) {
  return bch::bch_255_239_syndromes_zero_cw_255(cw.data()) &&
         overall_parity_syndrome_256(cw) == 0u;
}

inline newcode::HybridClassifierMode effective_hybrid_classifier_mode(
    const newcode::Params& p) {
  if (p.HYBRID_CLASSIFIER_MODE != newcode::HybridClassifierMode::LegacyHardDecode) {
    return p.HYBRID_CLASSIFIER_MODE;
  }
  if (p.HYBRID_USE_FAST_CLASSIFIER) {
    return newcode::HybridClassifierMode::RepoFastClassifier;
  }
  return newcode::HybridClassifierMode::LegacyHardDecode;
}

template <typename CoreLLR>
bool run_repo_fast_classifier_hard_finish(
    const std::array<CoreLLR, newcode::Params::BCH_N>& lin_vec,
    std::array<float, newcode::Params::BCH_N>* y2,
    const newcode::Params& p,
    HybridRowClass* out_class) {
  // Repo 版分类器的目标是：
  // 1. 直接从 256-bit 硬判码字里识别 clean / parity-only / 单错 / 两错；
  // 2. 一旦能确定最终合法码字，就直接产出 hard-finish 的 y2；
  // 3. 无法确定的情况统一返回 false，让外层继续走 soft path。
  auto cw = hard_decision_bits_256(lin_vec);
  const uint8_t s0 = overall_parity_syndrome_256(cw);
  const auto syndromes = bch::bch_255_239_syndromes_1_4_cw_255(cw.data());
  const uint8_t s1 = syndromes[0];
  const uint8_t s3 = syndromes[2];
  const uint8_t s1_cubed = hybrid_fast_gf_cube(s1);

  auto finish_with = [&](HybridRowClass cls) -> bool {
    // 所有“看起来可 hard-finish”的分支，都必须经过最终合法性校验。
    // 只有 255 位 BCH syndrome 清零且 overall parity 正确时，才允许输出 hard-finish。
    if (!hard_word_valid_256(cw)) {
      *out_class = HybridRowClass::HardFail;
      return false;
    }
    materialize_hard_finish_lout(cw, lin_vec, p, y2);
    *out_class = cls;
    return true;
  };

  if (s0 == 0u && s1 == 0u && s3 == 0u) {
    // 主体 syndrome 与 overall parity 全部为 0，说明整行已经是合法码字。
    return finish_with(HybridRowClass::Clean);
  }
  if (s0 == 1u && s1 == 0u && s3 == 0u) {
    // BCH 主体无错，但 overall parity 位错误，只需翻最后 1 bit。
    cw[newcode::Params::BCH_OVERALL_IDX] ^= 1u;
    return finish_with(HybridRowClass::ParityOnly);
  }
  if (s1 != 0u && s3 == s1_cubed) {
    // 单错模式：利用 S1 的对数定位主体错误位。
    const int pos = hybrid_fast_gf_log(s1);
    if (pos < 0 || pos >= static_cast<int>(newcode::Params::BCH_OVERALL_IDX)) {
      *out_class = HybridRowClass::HardFail;
      return false;
    }
    cw[static_cast<std::size_t>(pos)] ^= 1u;
    if (s0 == 0u) {
      // S0=0 说明除了主体单错外，overall parity 也需要一起翻回去。
      cw[newcode::Params::BCH_OVERALL_IDX] ^= 1u;
      return finish_with(HybridRowClass::OneMainPlusParity);
    }
    return finish_with(HybridRowClass::OneMain);
  }
  if (s0 == 0u && s1 != 0u && s3 != s1_cubed) {
    // Repo 版的两错入口条件是：
    // overall parity 已经匹配“两错应有的奇偶性”，且它又不是单错模式；
    // 这时直接调用现有 BCH t=2 硬解码器做最终确认与纠正。
    std::array<uint8_t, newcode::Params::BCH_N - 1> decoded{};
    int corrected_errors = 0;
    if (!bch::bch_255_239_decode_hiho_cw_255(cw.data(),
                                             decoded.data(),
                                             &corrected_errors) ||
        corrected_errors != 2) {
      *out_class = HybridRowClass::HardFail;
      return false;
    }
    for (std::size_t i = 0; i < decoded.size(); ++i) {
      cw[i] = decoded[i];
    }
    recompute_overall_parity(&cw);
    return finish_with(HybridRowClass::TwoMain);
  }

  *out_class = HybridRowClass::HardFail;
  return false;
}

template <typename CoreLLR>
bool run_friend_fast_classifier_hard_finish(
    const std::array<CoreLLR, newcode::Params::BCH_N>& lin_vec,
    std::array<float, newcode::Params::BCH_N>* y2,
    const newcode::Params& p,
    HybridRowClass* out_class) {
  // 朋友版分类器的主线与 Repo 版不同：
  // 它优先使用 S1 / S3 / Trace(mu) 判断“像 0 错 / 1 错 / 2 错 / 多错”，
  // 再把已识别出的 0/1/2 错模式映射成当前工程里的 hard-finish 语义。
  auto cw = hard_decision_bits_256(lin_vec);
  const uint8_t s0 = overall_parity_syndrome_256(cw);
  const auto syndromes = bch::bch_255_239_syndromes_1_4_cw_255(cw.data());
  const uint8_t s1 = syndromes[0];
  const uint8_t s3 = syndromes[2];

  auto finish_with = [&](HybridRowClass cls) -> bool {
    // 和 Repo 版一样，朋友版即便前级分类命中，也必须通过最终合法性校验，
    // 才能真正把这行升级成 HardFinish。
    if (!hard_word_valid_256(cw)) {
      *out_class = HybridRowClass::HardFail;
      return false;
    }
    materialize_hard_finish_lout(cw, lin_vec, p, y2);
    *out_class = cls;
    return true;
  };

  if (s1 == 0u) {
    if (s3 != 0u) {
      // 按朋友版语义：S1=0 但 S3!=0 时，不可能属于 <=2 bit error。
      *out_class = HybridRowClass::HardFail;
      return false;
    }
    if (s0 == 1u) {
      // BCH 主体无错，仅 overall parity 位错误。
      cw[newcode::Params::BCH_OVERALL_IDX] ^= 1u;
      return finish_with(HybridRowClass::ParityOnly);
    }
    return finish_with(HybridRowClass::Clean);
  }

  const uint8_t s1_sq = hybrid_fast_gf_mul(s1, s1);
  const uint8_t s1_cubed = hybrid_fast_gf_mul(s1_sq, s1);
  if (s3 == s1_cubed) {
    // 单错模式：S3 = S1^3。
    const int pos = hybrid_fast_gf_log(s1);
    if (pos < 0 || pos >= static_cast<int>(newcode::Params::BCH_OVERALL_IDX)) {
      *out_class = HybridRowClass::HardFail;
      return false;
    }
    cw[static_cast<std::size_t>(pos)] ^= 1u;
    if (s0 == 0u) {
      // S0=0 时说明 overall parity 也要一起翻转才能回到合法 256-bit 码字。
      cw[newcode::Params::BCH_OVERALL_IDX] ^= 1u;
      return finish_with(HybridRowClass::OneMainPlusParity);
    }
    return finish_with(HybridRowClass::OneMain);
  }

  // 走到这里表示：
  // 1. S1 != 0
  // 2. 它不是单错模式
  // 于是按朋友版公式计算 mu 与 Trace(mu)。
  const uint8_t numerator = static_cast<uint8_t>(s1_cubed ^ s3);
  const uint8_t mu = hybrid_fast_gf_div(numerator, s1_cubed);
  const uint8_t tr = hybrid_fast_trace_to_gf2(mu);
  if (tr != 0u) {
    // 当前朋友版定义：Trace(mu)!=0 视为多错，不能在 prepass 里硬完成。
    *out_class = HybridRowClass::HardFail;
    return false;
  }

  // 按当前朋友版语义：Trace(mu)==0 视为 2 错，
  // 将该行送入现有 BCH t=2 硬解码器做纠正。
  std::array<uint8_t, newcode::Params::BCH_N - 1> decoded{};
  int corrected_errors = 0;
  if (!bch::bch_255_239_decode_hiho_cw_255(cw.data(),
                                           decoded.data(),
                                           &corrected_errors)) {
    throw std::runtime_error(
        "FriendS1S3Classifier: BCH t=2 decode failed on a row "
        "that already passed the TwoMain prechecks.");
  }
  if (corrected_errors != 2) {
    throw std::runtime_error(
        "FriendS1S3Classifier: BCH decode succeeded but the row "
        "was not corrected as exactly two errors.");
  }
  for (std::size_t i = 0; i < decoded.size(); ++i) {
    cw[i] = decoded[i];
  }
  recompute_overall_parity(&cw);
  return finish_with(HybridRowClass::TwoMain);
}

template <typename CoreLLR>
bool run_friend_fast_classifier_with_s0_hard_finish(
    const std::array<CoreLLR, newcode::Params::BCH_N>& lin_vec,
    std::array<float, newcode::Params::BCH_N>* y2,
    const newcode::Params& p,
    HybridRowClass* out_class) {
  // 第三方案：
  // - 0 错 / 1 错分支保持朋友版思路；
  // - 2 错入口同时要求 Trace(mu)==0 与 S0==0；
  // - 只有“主体像 2 错”且“整体奇偶也像偶数错”时，才进入 BCH t=2 硬解码。
  auto cw = hard_decision_bits_256(lin_vec);
  const uint8_t s0 = overall_parity_syndrome_256(cw);
  const auto syndromes = bch::bch_255_239_syndromes_1_4_cw_255(cw.data());
  const uint8_t s1 = syndromes[0];
  const uint8_t s3 = syndromes[2];

  auto finish_with = [&](HybridRowClass cls) -> bool {
    if (!hard_word_valid_256(cw)) {
      *out_class = HybridRowClass::HardFail;
      return false;
    }
    materialize_hard_finish_lout(cw, lin_vec, p, y2);
    *out_class = cls;
    return true;
  };

  if (s1 == 0u) {
    if (s3 != 0u) {
      *out_class = HybridRowClass::HardFail;
      return false;
    }
    if (s0 == 1u) {
      cw[newcode::Params::BCH_OVERALL_IDX] ^= 1u;
      return finish_with(HybridRowClass::ParityOnly);
    }
    return finish_with(HybridRowClass::Clean);
  }

  const uint8_t s1_sq = hybrid_fast_gf_mul(s1, s1);
  const uint8_t s1_cubed = hybrid_fast_gf_mul(s1_sq, s1);
  if (s3 == s1_cubed) {
    const int pos = hybrid_fast_gf_log(s1);
    if (pos < 0 || pos >= static_cast<int>(newcode::Params::BCH_OVERALL_IDX)) {
      *out_class = HybridRowClass::HardFail;
      return false;
    }
    cw[static_cast<std::size_t>(pos)] ^= 1u;
    if (s0 == 0u) {
      cw[newcode::Params::BCH_OVERALL_IDX] ^= 1u;
      return finish_with(HybridRowClass::OneMainPlusParity);
    }
    return finish_with(HybridRowClass::OneMain);
  }

  const uint8_t numerator = static_cast<uint8_t>(s1_cubed ^ s3);
  const uint8_t mu = hybrid_fast_gf_div(numerator, s1_cubed);
  const uint8_t tr = hybrid_fast_trace_to_gf2(mu);
  if (tr != 0u) {
    *out_class = HybridRowClass::HardFail;
    return false;
  }
  if (s0 != 0u) {
    // 与朋友版相比，这里额外要求 overall parity 也满足偶数错特征。
    *out_class = HybridRowClass::HardFail;
    return false;
  }

  std::array<uint8_t, newcode::Params::BCH_N - 1> decoded{};
  int corrected_errors = 0;
  if (!bch::bch_255_239_decode_hiho_cw_255(cw.data(),
                                           decoded.data(),
                                           &corrected_errors)) {
    throw std::runtime_error(
        "FriendS1S3WithS0Classifier: BCH t=2 decode failed on a row "
        "that already passed the TwoMain prechecks.");
  }
  if (corrected_errors != 2) {
    throw std::runtime_error(
        "FriendS1S3WithS0Classifier: BCH decode succeeded but the row "
        "was not corrected as exactly two errors.");
  }
  for (std::size_t i = 0; i < decoded.size(); ++i) {
    cw[i] = decoded[i];
  }
  recompute_overall_parity(&cw);
  return finish_with(HybridRowClass::TwoMain);
}

template <typename CoreLLR>
bool run_selected_hybrid_classifier_hard_finish(
    const std::array<CoreLLR, newcode::Params::BCH_N>& lin_vec,
    std::array<float, newcode::Params::BCH_N>* y2,
    const newcode::Params& p,
    HybridRowClass* out_class) {
  // 统一分发入口：
  // 外层 prepass 只关心“是否能 hard-finish”以及“属于哪一类 hard path”，
  // 具体使用 Repo 版还是朋友版，由 mode 决定。
  switch (effective_hybrid_classifier_mode(p)) {
    case newcode::HybridClassifierMode::RepoFastClassifier:
      return run_repo_fast_classifier_hard_finish(lin_vec, y2, p, out_class);
    case newcode::HybridClassifierMode::FriendS1S3Classifier:
      return run_friend_fast_classifier_hard_finish(lin_vec, y2, p, out_class);
    case newcode::HybridClassifierMode::FriendS1S3WithS0Classifier:
      return run_friend_fast_classifier_with_s0_hard_finish(
          lin_vec, y2, p, out_class);
    case newcode::HybridClassifierMode::LegacyHardDecode:
    default:
      *out_class = HybridRowClass::HardFail;
      return false;
  }
}

} // namespace detail
} // namespace newcode
