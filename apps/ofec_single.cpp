#include <algorithm>
#include <chrono>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>

#include "newcode/params.hpp"
#include "newcode/pipeline_runner.hpp"  // run_pipeline()

using namespace newcode;
namespace fs = std::filesystem;

// ======== 用户可改区域 ========
// 只需要改这里的常量/列表即可完成一次“单次调试运行”的配置
static constexpr const char* kLabel         = "debug_L6";
static constexpr float       kEbN0_db       = 3.20f;
static constexpr int         kChaseL_override = 6;   // 设为 -1 则沿用 Params 默认

// 方式 A：统一填充值（长度自动取 Params::TILES_PER_WIN）
static constexpr float kAlpha_fill = 0.10f;
static constexpr float kBeta_fill  = 0.80f;

// 方式 B：显式列表（若非空，将覆盖填充值；长度必须等于 TILES_PER_WIN）
static const std::vector<float> kAlpha_explicit = {
  // 例：0.10f, 0.10f, 0.10f, 0.10f, 0.10f
};
static const std::vector<float> kBeta_explicit = {
  // 例：0.80f, 0.80f, 0.80f, 0.80f, 0.80f
};
// =============================

static std::string now_stamp() {
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

static void ensure_dir(const fs::path& p) {
  std::error_code ec;
  fs::create_directories(p, ec);
}

int main() {
  // 1) IO
  const fs::path data_dir = "data";
  ensure_dir(data_dir);
  const std::string log_path = (data_dir / ("run_" + now_stamp() + "_single.log")).string();
  std::ofstream flog(log_path, std::ios::out | std::ios::app);
  auto both = [&](const auto& x) -> void { std::cout << x; if (flog) flog << x; };

  // 2) 组装 Params
  Params p; // 用默认初始化
  if (kChaseL_override >= 0) {
    p.CHASE_L = kChaseL_override;
    p.CHASE_NTEST = 1 << p.CHASE_L;
  }

  const std::size_t T = p.TILES_PER_WIN;

  if (!kAlpha_explicit.empty()) {
    if (kAlpha_explicit.size() != T) {
      both("[ERROR] kAlpha_explicit 长度必须等于 TILES_PER_WIN\n");
      return 2;
    }
    p.ALPHA_LIST = kAlpha_explicit;
  } else {
    p.ALPHA_LIST.assign(T, kAlpha_fill);
  }
  if (!p.ALPHA_LIST.empty()) p.ALPHA = p.ALPHA_LIST.front();

  if (!kBeta_explicit.empty()) {
    if (kBeta_explicit.size() != T) {
      both("[ERROR] kBeta_explicit 长度必须等于 TILES_PER_WIN\n");
      return 2;
    }
    p.beta_list = kBeta_explicit;
  } else {
    p.beta_list.assign(T, kBeta_fill);
  }
  if (!p.beta_list.empty()) p.beta = p.beta_list.front();

  // 3) 运行（串行、单次）
  both("[INFO] run_pipeline(label="); both(kLabel);
  both(", Eb/N0="); both(kEbN0_db);
  both(" dB, CHASE_L="); both(p.CHASE_L); both(")\n");

  PipelineResult r = run_pipeline(p, kLabel, kEbN0_db);

  // 4) 概要
  both("[RESULT] Pre-FEC BER="); both(r.pre_fec.ber);
  both(" (errs="); both(r.pre_fec.errors); both("/"); both(r.pre_fec.total); both(")");
  both(" | Post-FEC BER="); both(r.post_fec.ber);
  both(" (errs="); both(r.post_fec.errors); both("/"); both(r.post_fec.total); both(")\n");

  if (!r.tile_early_stop_pct.empty()) {
    both("[RESULT] EarlyStop hit rates (%): ");
    std::ostringstream oss;
    oss.setf(std::ios::fixed); oss << std::setprecision(1);
    for (size_t i = 0; i < r.tile_early_stop_pct.size(); ++i) {
      oss << r.tile_early_stop_pct[i] << (i + 1 < r.tile_early_stop_pct.size() ? ", " : "");
    }
    both(oss.str()); both("\n");
  }

  both("[INFO] log saved at "); both(log_path); both("\n");
  return 0;
}
