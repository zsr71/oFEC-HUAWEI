#pragma once

#include "newcode/common/bch/bch_255_239.hpp"
#include "newcode/common/qfloat/llr_utils.hpp"
#include "newcode/params.hpp"

#include <array>
#include <cmath>
#include <cstdint>
#include <stdexcept>

namespace newcode {
namespace detail {

enum class HybridRowClass : uint8_t {
  None = 0,
  BchHardDecoded,
  Clean,
  ParityOnly,
  OneMain,
  OneMainPlusParity,
  TwoMain,
  Suspicious,
  HardFail
};

struct HybridFastGf256Tables {
  std::array<uint8_t, 255> alpha_to{};
  std::array<int16_t, 256> index_of{};
};

const HybridFastGf256Tables& hybrid_fast_gf256_tables();
uint8_t hybrid_fast_gf_mul(uint8_t x, uint8_t y);
uint8_t hybrid_fast_gf_div(uint8_t x, uint8_t y);
uint8_t hybrid_fast_gf_cube(uint8_t x);
int hybrid_fast_gf_log(uint8_t x);
uint8_t hybrid_fast_trace_to_gf2(uint8_t x);

uint8_t overall_parity_syndrome_256(
    const std::array<uint8_t, newcode::Params::BCH_N>& cw);
void recompute_overall_parity(
    std::array<uint8_t, newcode::Params::BCH_N>* cw);
bool hard_word_valid_256(
    const std::array<uint8_t, newcode::Params::BCH_N>& cw);

template <typename CoreLLR>
std::array<uint8_t, newcode::Params::BCH_N> hard_decision_bits_256(
    const std::array<CoreLLR, newcode::Params::BCH_N>& lin_vec);

template <typename CoreLLR>
void materialize_hard_finish_lout(
    const std::array<uint8_t, newcode::Params::BCH_N>& cw,
    const std::array<CoreLLR, newcode::Params::BCH_N>& lin_vec,
    const newcode::Params& p,
    std::array<float, newcode::Params::BCH_N>* y2);

template <typename CoreLLR>
void materialize_twomain_parameterized_lout(
    const std::array<uint8_t, newcode::Params::BCH_N>& corrected_cw,
    const std::array<CoreLLR, newcode::Params::BCH_N>& lin_vec,
    const std::array<bool, newcode::Params::BCH_N>& corrected_positions,
    const newcode::Params& p,
    std::array<float, newcode::Params::BCH_N>* y2);

newcode::HybridClassifierMode effective_hybrid_classifier_mode(
    const newcode::Params& p);

template <typename CoreLLR>
bool run_repo_fast_classifier_hard_finish(
    const std::array<CoreLLR, newcode::Params::BCH_N>& lin_vec,
    std::array<float, newcode::Params::BCH_N>* y2,
    const newcode::Params& p,
    HybridRowClass* out_class);

template <typename CoreLLR>
bool run_friend_fast_classifier_hard_finish(
    const std::array<CoreLLR, newcode::Params::BCH_N>& lin_vec,
    std::array<float, newcode::Params::BCH_N>* y2,
    const newcode::Params& p,
    HybridRowClass* out_class);

template <typename CoreLLR>
bool run_friend_fast_classifier_with_s0_hard_finish(
    const std::array<CoreLLR, newcode::Params::BCH_N>& lin_vec,
    std::array<float, newcode::Params::BCH_N>* y2,
    const newcode::Params& p,
    HybridRowClass* out_class);

template <typename CoreLLR>
bool run_selected_hybrid_classifier_hard_finish(
    const std::array<CoreLLR, newcode::Params::BCH_N>& lin_vec,
    std::array<float, newcode::Params::BCH_N>* y2,
    const newcode::Params& p,
    HybridRowClass* out_class);

template <typename CoreLLR>
bool run_selected_hybrid_classifier_classify_only(
    const std::array<CoreLLR, newcode::Params::BCH_N>& lin_vec,
    const newcode::Params& p,
    HybridRowClass* out_class);

} // namespace detail
} // namespace newcode

#include "../../../../src/rx/ofec/hybrid/hybrid_classifier.ipp"
