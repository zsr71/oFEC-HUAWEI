#include "newcode/awgn.hpp"
#include "newcode/bch_255_239.hpp"
#include "newcode/chase256.hpp"
#include "newcode/params.hpp"
#include "newcode/linspace.hpp"

#include <algorithm>
#include <array>
#include <complex>
#include <cstdint>
#include <fstream>
#include <future>
#include <iostream>
#include <random>
#include <string>
#include <vector>
#include <chrono>
#include <ctime>
#include <iomanip>
#include <sstream>

namespace {

enum class SoftDecoderKind {
  Plain,
  EbchPF
};

struct DemoConfig {
  std::size_t frames = 50000;
  float ebn0_R_start = 5.5f;
  float ebn0_R_end   = 8.5f;
  std::size_t ebn0_R_points = 7;
  bool run_soft = true;
  SoftDecoderKind soft_decoder = SoftDecoderKind::Plain;
  std::string csv_prefix = "bch_bpsk_demo";
};

uint8_t hard_from_llr(float llr) {
  return (llr < 0.0f) ? 1u : 0u;
}

struct SweepResult {
  float ebn0_R_db = 0.0f;
  float ebn0_db = 0.0f;
  float sigma = 0.0f;
  std::size_t total_bits = 0;
  std::size_t pre_errors = 0;
  std::size_t post_hard_errors = 0;
  std::size_t post_soft_errors = 0;
  std::size_t hard_failures = 0;
};

static void ensure_csv_header(const std::string& path) {
  std::ifstream fin(path);
  if (fin.good() && fin.peek() != std::ifstream::traits_type::eof()) {
    return;
  }
  std::ofstream fout(path, std::ios::out | std::ios::app);
  fout << "ebn0_R_db,ebn0_db,frames,pre_ber,pre_errs,pre_total,"
          "post_hard_ber,post_hard_errs,post_hard_total,"
          "post_soft_ber,post_soft_errs,post_soft_total,"
          "hard_failures\n";
}

static std::string timestamp_stamp()
{
  using clock = std::chrono::system_clock;
  auto t = clock::to_time_t(clock::now());
  std::tm tm{};
#ifdef _WIN32
  localtime_s(&tm, &t);
#else
  localtime_r(&t, &tm);
#endif
  std::ostringstream oss;
  oss << std::put_time(&tm, "%Y%m%d-%H%M%S");
  return oss.str();
}

} // namespace

int main() {
  const DemoConfig cfg{};
  const std::vector<float> ebn0_R_list = newcode::linspace(cfg.ebn0_R_start, cfg.ebn0_R_end, cfg.ebn0_R_points);
  std::vector<float> ebn0_list;
  ebn0_list.reserve(ebn0_R_list.size());
  for (float val : ebn0_R_list) {
    ebn0_list.push_back(val - 0.281422f);
  }
  const std::string csv_path = cfg.csv_prefix + "_" + timestamp_stamp() + ".csv";

  newcode::Params chase_params;
  chase_params.CHASE_L = 2;
  chase_params.CHASE_NTEST = 4;
  chase_params.beta = 1.75f;
  chase_params.ALPHA = 0.8f;

  std::vector<SweepResult> results(ebn0_list.size());
  std::vector<std::future<void>> tasks;
  tasks.reserve(ebn0_list.size());

  for (std::size_t idx = 0; idx < ebn0_list.size(); ++idx) {
    tasks.emplace_back(std::async(std::launch::async, [&, idx]() {
      std::random_device rd_thread;
      std::mt19937 bit_rng_local(rd_thread());
      std::mt19937 awgn_seed_rng_local(rd_thread());
      std::uniform_int_distribution<int> bit_dist_local(0, 1);
      std::uniform_int_distribution<uint32_t> seed_dist_local;

      const float ebn0_R_db = ebn0_R_list[idx];
      const float ebn0_db   = ebn0_list[idx];
      const float sigma = newcode::ebn0_to_sigma(ebn0_db, 1);
      const float inv_sigma_sq = 1.0f / (sigma * sigma);

      SweepResult res;
      res.ebn0_R_db = ebn0_R_db;
      res.ebn0_db = ebn0_db;
      res.sigma = sigma;

      for (std::size_t frame = 0; frame < cfg.frames; ++frame) {
        std::vector<uint8_t> info_bits(239);
        for (auto& b : info_bits) {
          b = static_cast<uint8_t>(bit_dist_local(bit_rng_local));
        }

        const auto codeword = newcode::bch_255_239_encode(info_bits);

        std::vector<std::complex<float>> tx_syms(newcode::Params::BCH_N);
        for (std::size_t i = 0; i < tx_syms.size(); ++i) {
          tx_syms[i] = std::complex<float>(codeword[i] ? -1.0f : +1.0f, 0.0f);
        }

        const uint32_t awgn_seed = seed_dist_local(awgn_seed_rng_local);
        auto rx_syms = newcode::add_awgn(tx_syms, ebn0_db, 1, awgn_seed);

        std::array<float, newcode::Params::BCH_N> channel_llr{};
        std::array<uint8_t, newcode::Params::BCH_N> channel_hard{};
        for (std::size_t i = 0; i < rx_syms.size(); ++i) {
          const float y = rx_syms[i].real();
          const float llr =  y ;
          channel_llr[i] = llr;
          channel_hard[i] = hard_from_llr(llr);
        }

        for (std::size_t i = 0; i < info_bits.size(); ++i) {
          if (channel_hard[i] != (info_bits[i] & 1u)) {
            res.pre_errors++;
          }
        }

        std::array<uint8_t, 255> hard_in{};
        std::array<uint8_t, 255> hard_decoded{};
        std::copy_n(channel_hard.begin(), hard_in.size(), hard_in.begin());
        const bool hard_ok = newcode::bch_255_239_decode_hiho_cw_255(
            hard_in.data(), hard_decoded.data());
        if (!hard_ok) {
          res.hard_failures++;
        }

        for (std::size_t i = 0; i < info_bits.size(); ++i) {
          if (hard_decoded[i] != (info_bits[i] & 1u)) {
            res.post_hard_errors++;
          }
        }

        if (cfg.run_soft) {
          std::array<float, newcode::Params::BCH_N> extrinsic{};
          switch (cfg.soft_decoder) {
            case SoftDecoderKind::Plain:
              newcode::chase_decode_256_plain<float>(
                  channel_llr.data(), channel_llr.data(), extrinsic.data(), chase_params);
              break;
            case SoftDecoderKind::EbchPF:
              newcode::chase_decode_256_ebchPF<float>(
                  channel_llr.data(), channel_llr.data(), extrinsic.data(), chase_params);
              break;
          }
          for (auto& value : extrinsic) {
            value *= chase_params.ALPHA;
          }
          for (std::size_t i = 0; i < info_bits.size(); ++i) {
            const float total_llr = channel_llr[i] + extrinsic[i];
            const uint8_t bit = hard_from_llr(total_llr);
            if (bit != (info_bits[i] & 1u)) {
              res.post_soft_errors++;
            }
          }
        }

        res.total_bits += info_bits.size();
      }

      results[idx] = res;
    }));
  }

  for (auto& task : tasks) task.get();

  for (const auto& res : results) {
    const double pre_ber =
        res.total_bits ? static_cast<double>(res.pre_errors) / res.total_bits : 0.0;
    const double post_hard_ber =
        res.total_bits ? static_cast<double>(res.post_hard_errors) / res.total_bits : 0.0;
    const double post_soft_ber =
        (cfg.run_soft && res.total_bits)
            ? static_cast<double>(res.post_soft_errors) / res.total_bits
            : 0.0;

    std::cout << "Eb/N0: " << res.ebn0_db << " dB\n";
    std::cout << "  Frames: " << cfg.frames << "\n";
    std::cout << "  Sigma: " << res.sigma << "\n";
    std::cout << "  Pre-FEC BER: " << pre_ber << " (" << res.pre_errors << " / "
              << res.total_bits << ")\n";
    std::cout << "  Post-FEC BER (hard): " << post_hard_ber << " ("
              << res.post_hard_errors << " / " << res.total_bits << ")\n";
    if (cfg.run_soft) {
      const char* decoder_label =
          (cfg.soft_decoder == SoftDecoderKind::Plain) ? "plain" : "ebchPF";
      std::cout << "  Post-FEC BER (soft Chase - " << decoder_label << "): "
                << post_soft_ber << " (" << res.post_soft_errors << " / "
                << res.total_bits << ")\n";
    }
    std::cout << "  Hard decoder failures: " << res.hard_failures << "\n\n";
  }

  if (!csv_path.empty()) {
    ensure_csv_header(csv_path);
    std::ofstream fout(csv_path, std::ios::out | std::ios::app);
    for (const auto& res : results) {
      const double pre_ber =
          res.total_bits ? static_cast<double>(res.pre_errors) / res.total_bits : 0.0;
      const double post_hard_ber =
          res.total_bits ? static_cast<double>(res.post_hard_errors) / res.total_bits : 0.0;
      const double post_soft_ber =
          (cfg.run_soft && res.total_bits)
              ? static_cast<double>(res.post_soft_errors) / res.total_bits
              : 0.0;
      fout << res.ebn0_R_db << ','
           << res.ebn0_db << ','
           << cfg.frames << ','
           << pre_ber << ',' << res.pre_errors << ',' << res.total_bits << ','
           << post_hard_ber << ',' << res.post_hard_errors << ',' << res.total_bits << ',';
      if (cfg.run_soft) {
        fout << post_soft_ber << ',' << res.post_soft_errors << ',' << res.total_bits << ',';
      } else {
        fout << "0,0," << res.total_bits << ',';
      }
      fout << res.hard_failures << '\n';
    }
  }

  return 0;
}
