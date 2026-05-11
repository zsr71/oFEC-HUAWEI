#include <algorithm>
#include <array>
#include <cstdint>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <random>
#include <sstream>
#include <string>
#include <vector>

#include "newcode/common/bch/bch_255_239.hpp"
#include "newcode/ofec/hybrid/hybrid_classifier.hpp"
#include "newcode/common/qfloat/llr_utils.hpp"
#include "newcode/params.hpp"

namespace {

using Codeword256 = std::array<uint8_t, newcode::Params::BCH_N>;
using LinVec256 = std::array<float, newcode::Params::BCH_N>;
using RepoClass = newcode::detail::HybridRowClass;

constexpr std::size_t kLaneCount = 32u;
constexpr int kInfoBits = static_cast<int>(newcode::Params::BCH_K);
constexpr float kHardDecisionMag = 12.0f;
constexpr int kInfoSeed = 20260511;
constexpr int kErrorSeed = 20260512;

enum class UnifiedClass : uint8_t {
  Clean = 0,
  ParityOnly,
  OneMain,
  OneMainPlusParity,
  TwoMain,
  HardFail
};

enum class FriendClass : uint8_t {
  Invalid = 0,
  ZeroError,
  OneError,
  TwoCandidate
};

struct FriendEvalResult {
  FriendClass cls = FriendClass::Invalid;
  bool flags[3] = {false, false, false};
  bool two_candidate_decoder_success = false;
  int corrected_errors = -1;
  Codeword256 corrected{};
  bool has_corrected = false;
  std::string note;
};

struct RepoEvalResult {
  UnifiedClass cls = UnifiedClass::HardFail;
  bool hard_finish = false;
  Codeword256 corrected{};
  std::string note;
};

struct LaneResult {
  std::size_t lane = 0;
  std::vector<int> flip_positions;
  uint8_t s0 = 0u;
  uint8_t s1 = 0u;
  uint8_t s3 = 0u;
  uint8_t s1_cubed = 0u;
  bool has_tr = false;
  uint8_t tr = 0u;
  FriendClass friend_class = FriendClass::Invalid;
  std::string friend_flags_text;
  int friend_corrected_errors = -1;
  UnifiedClass repo_class = UnifiedClass::HardFail;
  bool same_coarse_class = false;
  bool same_hard_finish_codeword = false;
  bool repo_hard_finish = false;
  std::string friend_note;
  std::string repo_note;
};

uint8_t gf_div(uint8_t x, uint8_t y) {
  if (x == 0u) {
    return 0u;
  }
  if (y == 0u) {
    throw std::runtime_error("gf_div: division by zero");
  }
  const auto& tables = newcode::detail::hybrid_fast_gf256_tables();
  const int lx = tables.index_of[static_cast<std::size_t>(x)];
  const int ly = tables.index_of[static_cast<std::size_t>(y)];
  int diff = lx - ly;
  if (diff < 0) {
    diff += 255;
  }
  return tables.alpha_to[static_cast<std::size_t>(diff)];
}

uint8_t trace_to_gf2(uint8_t x) {
  uint8_t acc = x;
  uint8_t cur = x;
  for (int i = 1; i < 8; ++i) {
    cur = newcode::detail::hybrid_fast_gf_mul(cur, cur);
    acc ^= cur;
  }
  return static_cast<uint8_t>(acc & 0x1u);
}

std::string join_positions(const std::vector<int>& values) {
  std::ostringstream oss;
  oss << '[';
  for (std::size_t i = 0; i < values.size(); ++i) {
    if (i) {
      oss << ' ';
    }
    oss << values[i];
  }
  oss << ']';
  return oss.str();
}

std::string class_name(UnifiedClass cls) {
  switch (cls) {
    case UnifiedClass::Clean: return "Clean";
    case UnifiedClass::ParityOnly: return "ParityOnly";
    case UnifiedClass::OneMain: return "OneMain";
    case UnifiedClass::OneMainPlusParity: return "OneMainPlusParity";
    case UnifiedClass::TwoMain: return "TwoMain";
    case UnifiedClass::HardFail: return "HardFail";
  }
  return "Unknown";
}

std::string friend_class_name(FriendClass cls) {
  switch (cls) {
    case FriendClass::Invalid: return "Invalid";
    case FriendClass::ZeroError: return "ZeroError";
    case FriendClass::OneError: return "OneError";
    case FriendClass::TwoCandidate: return "TwoCandidate";
  }
  return "Unknown";
}

std::string csv_escape(const std::string& s) {
  bool needs_quotes = false;
  for (char ch : s) {
    if (ch == ',' || ch == '"' || ch == '\n') {
      needs_quotes = true;
      break;
    }
  }
  if (!needs_quotes) {
    return s;
  }
  std::string escaped;
  escaped.reserve(s.size() + 8);
  escaped.push_back('"');
  for (char ch : s) {
    if (ch == '"') {
      escaped.push_back('"');
    }
    escaped.push_back(ch);
  }
  escaped.push_back('"');
  return escaped;
}

std::string syndrome_summary(uint8_t s0, uint8_t s1, uint8_t s3, uint8_t s1_cubed) {
  std::ostringstream oss;
  oss << "S0=" << static_cast<int>(s0)
      << " S1=" << static_cast<int>(s1)
      << " S3=" << static_cast<int>(s3)
      << " S1^3=" << static_cast<int>(s1_cubed);
  return oss.str();
}

UnifiedClass map_repo_class(RepoClass cls) {
  switch (cls) {
    case RepoClass::Clean: return UnifiedClass::Clean;
    case RepoClass::ParityOnly: return UnifiedClass::ParityOnly;
    case RepoClass::OneMain: return UnifiedClass::OneMain;
    case RepoClass::OneMainPlusParity: return UnifiedClass::OneMainPlusParity;
    case RepoClass::TwoMain: return UnifiedClass::TwoMain;
    case RepoClass::HardFail: return UnifiedClass::HardFail;
    default: return UnifiedClass::HardFail;
  }
}

Codeword256 decode_codeword_from_y2(const LinVec256& lin_vec,
                                    const std::array<float, newcode::Params::BCH_N>& y2) {
  Codeword256 cw{};
  for (std::size_t i = 0; i < cw.size(); ++i) {
    const float lpost = y2[i] + qfloat::llr_to_float(lin_vec[i]);
    cw[i] = (lpost < 0.0f) ? 1u : 0u;
  }
  return cw;
}

LinVec256 make_lin_vec_from_codeword(const Codeword256& cw) {
  LinVec256 lin{};
  for (std::size_t i = 0; i < cw.size(); ++i) {
    lin[i] = cw[i] ? -kHardDecisionMag : kHardDecisionMag;
  }
  return lin;
}

Codeword256 make_random_codeword(std::mt19937& rng) {
  std::uniform_int_distribution<int> bit_dist(0, 1);
  std::vector<uint8_t> info(static_cast<std::size_t>(kInfoBits), 0u);
  for (auto& bit : info) {
    bit = static_cast<uint8_t>(bit_dist(rng));
  }
  return bch::bch_255_239_encode(info);
}

std::vector<int> inject_random_errors(Codeword256* cw, std::mt19937& rng) {
  std::uniform_int_distribution<int> flip_count_dist(0, 4);
  std::uniform_int_distribution<int> pos_dist(0, static_cast<int>(cw->size() - 1));
  const int flip_count = flip_count_dist(rng);
  std::vector<int> positions;
  positions.reserve(static_cast<std::size_t>(flip_count));
  while (static_cast<int>(positions.size()) < flip_count) {
    const int pos = pos_dist(rng);
    if (std::find(positions.begin(), positions.end(), pos) != positions.end()) {
      continue;
    }
    positions.push_back(pos);
    (*cw)[static_cast<std::size_t>(pos)] ^= 1u;
  }
  std::sort(positions.begin(), positions.end());
  return positions;
}

FriendEvalResult evaluate_friend_classifier(const LinVec256& lin_vec) {
  FriendEvalResult out;
  Codeword256 cw = newcode::detail::hard_decision_bits_256(lin_vec);
  const auto syndromes = bch::bch_255_239_syndromes_1_4_cw_255(cw.data());
  const uint8_t s1 = syndromes[0];
  const uint8_t s3 = syndromes[2];

  if (s1 == 0u) {
    out.flags[0] = (s3 != 0u);
    out.flags[1] = (s3 == 0u);
    out.flags[2] = false;
    out.cls = out.flags[0] ? FriendClass::Invalid : FriendClass::ZeroError;
    out.note = out.flags[0] ? "S1=0,S3!=0" : "zero-error";
    return out;
  }

  const uint8_t s1_sq = newcode::detail::hybrid_fast_gf_mul(s1, s1);
  const uint8_t s1_cubed = newcode::detail::hybrid_fast_gf_mul(s1_sq, s1);
  out.flags[1] = false;
  if (s3 == s1_cubed) {
    out.flags[0] = false;
    out.flags[2] = true;
    out.cls = FriendClass::OneError;
    out.note = "one-error";
    return out;
  }

  const uint8_t numerator = static_cast<uint8_t>(s1_cubed ^ s3);
  const uint8_t mu = gf_div(numerator, s1_cubed);
  const uint8_t tr = trace_to_gf2(mu);
  if (tr != 0u) {
    out.flags[0] = true;
    out.flags[2] = false;
    out.cls = FriendClass::Invalid;
    out.note = "Trace(mu)!=0";
    return out;
  }

  out.flags[0] = false;
  out.flags[2] = false;
  out.cls = FriendClass::TwoCandidate;
  out.note = "two-candidate";

  std::array<uint8_t, newcode::Params::BCH_N - 1> decoded{};
  int corrected_errors = 0;
  if (!bch::bch_255_239_decode_hiho_cw_255(cw.data(),
                                           decoded.data(),
                                           &corrected_errors)) {
    out.note = "two-candidate decoder-fail";
    out.corrected_errors = corrected_errors;
    return out;
  }
  out.corrected_errors = corrected_errors;
  if (corrected_errors == 2) {
    for (std::size_t i = 0; i < decoded.size(); ++i) {
      cw[i] = decoded[i];
    }
    newcode::detail::recompute_overall_parity(&cw);
    out.corrected = cw;
    out.has_corrected = true;
    out.two_candidate_decoder_success = true;
    out.note = "two-candidate decoder-ok";
  } else {
    out.note = "two-candidate corrected_errors!=" + std::to_string(corrected_errors);
  }
  return out;
}

RepoEvalResult evaluate_repo_classifier(const LinVec256& lin_vec,
                                        const newcode::Params& p) {
  RepoEvalResult out;
  std::array<float, newcode::Params::BCH_N> y2{};
  RepoClass repo_class = RepoClass::HardFail;
  const bool hard_ok =
      newcode::detail::run_repo_fast_classifier_hard_finish(lin_vec, &y2, p, &repo_class);
  out.cls = map_repo_class(repo_class);
  out.hard_finish = hard_ok;
  if (hard_ok) {
    out.corrected = decode_codeword_from_y2(lin_vec, y2);
    out.note = "hard-finish";
  } else {
    out.note = "soft-path";
  }
  return out;
}

bool same_codeword(const Codeword256& a, const Codeword256& b) {
  return std::equal(a.begin(), a.end(), b.begin());
}

void print_row(const LaneResult& row) {
  std::cout << std::left
            << std::setw(6) << row.lane
            << std::setw(12) << row.flip_positions.size()
            << std::setw(24) << join_positions(row.flip_positions)
            << std::setw(20) << friend_class_name(row.friend_class)
            << std::setw(20) << class_name(row.repo_class)
            << std::setw(12) << (row.same_coarse_class ? "yes" : "no")
            << std::setw(16)
            << ((row.repo_hard_finish && row.same_hard_finish_codeword)
                    ? "yes"
                    : "no")
            << std::setw(22) << row.friend_note
            << std::setw(22) << row.repo_note
            << '\n';
}

void write_csv_report(const std::vector<LaneResult>& rows) {
  namespace fs = std::filesystem;
  const fs::path out_dir = fs::path("/home/zsr71/projects/newcode/data");
  std::error_code ec;
  fs::create_directories(out_dir, ec);
  const fs::path out_file = out_dir / "ofec_fast_classifier_compare.csv";

  std::ofstream out(out_file, std::ios::trunc);
  if (!out) {
    throw std::runtime_error("failed to open CSV output: " + out_file.string());
  }

  out << "lane,flip_cnt,flip_positions,friend_class,repo_class,same_class,"
         "same_hard_finish_codeword,friend_note,repo_note,friend_flags,"
         "friend_corrected_errors,syndrome_summary,tr\n";
  for (const auto& row : rows) {
    out << row.lane << ','
        << row.flip_positions.size() << ','
        << csv_escape(join_positions(row.flip_positions)) << ','
        << friend_class_name(row.friend_class) << ','
        << class_name(row.repo_class) << ','
        << (row.same_coarse_class ? "yes" : "no") << ','
        << ((row.repo_hard_finish && row.same_hard_finish_codeword)
                ? "yes"
                : "no")
        << ','
        << csv_escape(row.friend_note) << ','
        << csv_escape(row.repo_note) << ','
        << csv_escape(row.friend_flags_text) << ',';
    if (row.friend_corrected_errors >= 0) {
      out << row.friend_corrected_errors;
    }
    out << ','
        << csv_escape(syndrome_summary(row.s0, row.s1, row.s3, row.s1_cubed)) << ',';
    if (row.has_tr) {
      out << static_cast<int>(row.tr);
    }
    out << '\n';
  }
}

}  // namespace

int main() {
  newcode::Params p;
  p.HARD_LLR_MAG = 99.0f;
  p.HYBRID_ENABLE = true;
  p.HYBRID_USE_FAST_CLASSIFIER = true;

  std::mt19937 info_rng(kInfoSeed);
  std::mt19937 error_rng(kErrorSeed);

  std::array<int, 6> repo_counts{};
  std::array<int, 4> friend_class_counts{};
  int class_mismatch_count = 0;
  int hard_finish_codeword_mismatch_count = 0;
  std::vector<LaneResult> rows;
  rows.reserve(kLaneCount);

  std::cout << std::left
            << std::setw(6) << "lane"
            << std::setw(12) << "flip_cnt"
            << std::setw(24) << "flip_positions"
            << std::setw(20) << "friend_class"
            << std::setw(20) << "repo_class"
            << std::setw(12) << "same_map"
            << std::setw(16) << "same_cw"
            << std::setw(22) << "friend_note"
            << std::setw(22) << "repo_note"
            << '\n';

  for (std::size_t lane = 0; lane < kLaneCount; ++lane) {
    Codeword256 injected = make_random_codeword(info_rng);
    const std::vector<int> flips = inject_random_errors(&injected, error_rng);
    const LinVec256 lin_vec = make_lin_vec_from_codeword(injected);
    const Codeword256 hard_bits = newcode::detail::hard_decision_bits_256(lin_vec);
    const uint8_t s0 = newcode::detail::overall_parity_syndrome_256(hard_bits);
    const auto syndromes = bch::bch_255_239_syndromes_1_4_cw_255(hard_bits.data());
    const uint8_t s1 = syndromes[0];
    const uint8_t s3 = syndromes[2];
    const uint8_t s1_cubed = newcode::detail::hybrid_fast_gf_cube(s1);

    const FriendEvalResult friend_result = evaluate_friend_classifier(lin_vec);
    const RepoEvalResult repo_result = evaluate_repo_classifier(lin_vec, p);

    LaneResult row;
    row.lane = lane;
    row.flip_positions = flips;
    row.s0 = s0;
    row.s1 = s1;
    row.s3 = s3;
    row.s1_cubed = s1_cubed;
    if (s1 != 0u && s3 != s1_cubed) {
      const uint8_t numerator = static_cast<uint8_t>(s1_cubed ^ s3);
      const uint8_t mu = gf_div(numerator, s1_cubed);
      row.has_tr = true;
      row.tr = trace_to_gf2(mu);
    }
    row.friend_class = friend_result.cls;
    {
      std::ostringstream flag_oss;
      flag_oss << '['
               << (friend_result.flags[0] ? 1 : 0) << ' '
               << (friend_result.flags[1] ? 1 : 0) << ' '
               << (friend_result.flags[2] ? 1 : 0) << ']';
      row.friend_flags_text = flag_oss.str();
    }
    row.friend_corrected_errors = friend_result.corrected_errors;
    row.repo_class = repo_result.cls;
    auto mapped_friend = UnifiedClass::HardFail;
    switch (friend_result.cls) {
      case FriendClass::Invalid:
        mapped_friend = UnifiedClass::HardFail;
        break;
      case FriendClass::ZeroError:
        mapped_friend = (row.s0 == 0u) ? UnifiedClass::Clean : UnifiedClass::ParityOnly;
        break;
      case FriendClass::OneError:
        mapped_friend = (row.s0 == 0u) ? UnifiedClass::OneMainPlusParity : UnifiedClass::OneMain;
        break;
      case FriendClass::TwoCandidate:
        mapped_friend = UnifiedClass::TwoMain;
        break;
    }
    row.same_coarse_class = (mapped_friend == repo_result.cls);
    row.repo_hard_finish = repo_result.hard_finish;
    row.friend_note = friend_result.note;
    row.repo_note = repo_result.note;
    row.same_hard_finish_codeword =
        friend_result.has_corrected &&
        repo_result.hard_finish &&
        same_codeword(friend_result.corrected, repo_result.corrected);

    ++friend_class_counts[static_cast<std::size_t>(friend_result.cls)];
    ++repo_counts[static_cast<std::size_t>(repo_result.cls)];
    if (!row.same_coarse_class) {
      ++class_mismatch_count;
    }
    if (friend_result.has_corrected && repo_result.hard_finish &&
        !row.same_hard_finish_codeword) {
      ++hard_finish_codeword_mismatch_count;
    }
    rows.push_back(row);
    print_row(row);
  }

  std::cout << "\nSummary\n";
  for (std::size_t i = 0; i < friend_class_counts.size(); ++i) {
    const auto cls = static_cast<FriendClass>(i);
    std::cout << "  friend_" << std::setw(14) << friend_class_name(cls)
              << "=" << std::setw(3) << friend_class_counts[i] << '\n';
  }
  for (std::size_t i = 0; i < repo_counts.size(); ++i) {
    const auto cls = static_cast<UnifiedClass>(i);
    std::cout << "  repo_" << std::setw(16) << class_name(cls)
              << "=" << std::setw(3) << repo_counts[i] << '\n';
  }
  std::cout << "  class_mismatch_count=" << class_mismatch_count << '\n';
  std::cout << "  hard_finish_codeword_mismatch_count="
            << hard_finish_codeword_mismatch_count << '\n';
  write_csv_report(rows);
  std::cout << "  csv_output=/home/zsr71/projects/newcode/data/ofec_fast_classifier_compare.csv\n";
  return 0;
}
