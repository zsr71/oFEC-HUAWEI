#pragma once

namespace newcode {
namespace detail {

constexpr int BCH_N_TOTAL = 256; // 255 core + 1 overall parity (extended code)
constexpr int BCH_N_CORE  = 255;
constexpr int PAR_IDX     = 255;

} // namespace detail
} // namespace newcode

