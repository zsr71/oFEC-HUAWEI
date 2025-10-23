#include <algorithm>
#include <chrono>
#include <condition_variable>
#include <ctime>
#include <filesystem>
#include <fstream>
#include <future>
#include <iomanip>
#include <iostream>
#include <limits>
#include <mutex>
#include <sstream>
#include <string>
#include <thread>
#include <vector>

#include "newcode/params.hpp"
#include "newcode/pipeline_runner.hpp"

using namespace newcode;
namespace fs = std::filesystem;

// -------------------- 并发限流：C++17 简单信号量 --------------------
class Semaphore {
  std::mutex m_;
  std::condition_variable cv_;
  size_t count_;

public:
  explicit Semaphore(size_t c) : count_(c) {}
  void acquire() {
    std::unique_lock<std::mutex> lk(m_);
    cv_.wait(lk, [&] { return count_ > 0; });
    --count_;
  }
  void release() {
    std::lock_guard<std::mutex> lk(m_);
    ++count_;
    cv_.notify_one();
  }
};

// -------------------- 双写输出：控制台 + 文件 --------------------
struct DualOut {
  std::ostream& console;
  std::ofstream file;
  DualOut(std::ostream& c, const std::string& filepath)
      : console(c), file(filepath, std::ios::out | std::ios::app) {}
  template <typename T>
  DualOut& operator<<(const T& v) {
    console << v;
    if (file) file << v;
    return *this;
  }
  DualOut& operator<<(std::ostream& (*pf)(std::ostream&)) {
    pf(console);
    if (file) pf(file);
    return *this;
  }
};

// -------------------- 时间戳/工具函数 --------------------
static std::string now_stamp()
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

static void ensure_dir(const fs::path& p)
{
  std::error_code ec;
  fs::create_directories(p, ec);
}

template <typename T>
static std::string join_vec(const std::vector<T>& v, char sep = '|', int prec = 6)
{
  std::ostringstream oss;
  oss.setf(std::ios::fixed);
  oss << std::setprecision(prec);
  for (size_t i = 0; i < v.size(); ++i) {
    oss << v[i];
    if (i + 1 < v.size()) oss << sep;
  }
  return oss.str();
}

static void ensure_csv_header(const std::string& csv_path)
{
  std::ifstream fin(csv_path);
  if (fin.good() && fin.peek() != std::ifstream::traits_type::eof()) return;

  std::ofstream fout(csv_path, std::ios::out | std::ios::app);
  fout << "timestamp,run_id,scenario,alpha,beta,ALPHA_LIST,beta_list,"
          "pre_ber,pre_errs,pre_total,post_ber,post_errs,post_total,"
          "early_stop_mean_pct,early_stop_list\n";
}

static double mean(const std::vector<double>& v)
{
  if (v.empty()) return std::numeric_limits<double>::quiet_NaN();
  long double sum = 0.0L;
  for (double x : v) sum += x;
  return static_cast<double>(sum / v.size());
}

// -------------------- 原始结构 --------------------
struct SweepScenario {
  std::string name;
  std::vector<float> alpha_list;
  std::vector<float> beta_list;
};

static std::vector<SweepScenario> build_scenarios(const Params& base_params,
                                                  const std::vector<float>& alpha_candidates,
                                                  const std::vector<float>& beta_candidates)
{
  std::vector<SweepScenario> scenarios;
  scenarios.push_back({"baseline", base_params.ALPHA_LIST, base_params.beta_list});

  for (float alpha_val : alpha_candidates) {
    for (float beta_val : beta_candidates) {
      SweepScenario sc;
      std::ostringstream oss;
      oss << "alpha" << alpha_val << "_beta" << beta_val;
      sc.name = oss.str();
      sc.alpha_list.assign(base_params.TILES_PER_WIN, alpha_val);
      sc.beta_list.assign(base_params.TILES_PER_WIN, beta_val);
      scenarios.push_back(std::move(sc));
    }
  }
  return scenarios;
}

// 每个任务的返回包（带索引，便于按原顺序汇总）
struct ScenarioOutput {
  std::size_t idx{};
  std::string name;
  std::vector<float> alpha_list;
  std::vector<float> beta_list;
  PipelineResult result;
};

int main()
{
  // ========== 1) IO 准备 ==========
  const fs::path data_dir = "data";
  ensure_dir(data_dir);
  const std::string run_id = now_stamp();
  const std::string log_path = (data_dir / ("run_" + run_id + ".log")).string();
  DualOut out(std::cout, log_path);

  const std::string csv_path = (data_dir / "ofec_sweep_results.csv").string();
  ensure_csv_header(csv_path);

  // ========== 2) 生成 scenario ==========
  Params base_params;

  // 示例候选集，可按需调整
  const std::vector<float> alpha_candidates = {0.3f,0.4f,0.5f};
  const std::vector<float> beta_candidates  = {0.9f,0.8f,0.7f};

  auto scenarios = build_scenarios(base_params, alpha_candidates, beta_candidates);
  const std::size_t NS = scenarios.size();

  out << "[INFO] total scenarios = " << NS << "\n";

  // ========== 3) 并行执行 run_pipeline ==========
  unsigned hw = std::max(1u, std::thread::hardware_concurrency());
  unsigned max_workers = hw;
  out << "[INFO] using up to " << max_workers << " workers\n";
  Semaphore sem(max_workers);

  std::vector<std::future<ScenarioOutput>> futures;
  futures.reserve(NS);

  for (std::size_t idx = 0; idx < NS; ++idx) {
    const auto& scenario = scenarios[idx];

    if (scenario.alpha_list.size() != base_params.TILES_PER_WIN ||
        scenario.beta_list.size()  != base_params.TILES_PER_WIN) {
      out << "[WARN] Scenario '" << scenario.name
          << "' skipped due to list size mismatch (expected "
          << base_params.TILES_PER_WIN << ")\n";
      continue;
    }

    Params params = base_params;
    params.ALPHA_LIST = scenario.alpha_list;
    params.beta_list  = scenario.beta_list;
    if (!params.ALPHA_LIST.empty()) params.ALPHA = params.ALPHA_LIST.front();
    if (!params.beta_list.empty())  params.beta  = params.beta_list.front();

    sem.acquire();
    futures.emplace_back(std::async(std::launch::async,
                                    [idx, scenario, params, &sem]() -> ScenarioOutput {
                                      struct Releaser {
                                        Semaphore& s;
                                        ~Releaser() { s.release(); }
                                      } _releaser{sem};
                                      ScenarioOutput out;
                                      out.idx = idx;
                                      out.name = scenario.name;
                                      out.alpha_list = scenario.alpha_list;
                                      out.beta_list  = scenario.beta_list;
                                      out.result = run_pipeline(params, scenario.name);
                                      return out;
                                    }));
  }

  // ========== 4) 汇总：主线程统一写 CSV/日志，挑选 best ==========
  std::vector<std::string> scenario_summaries;
  scenario_summaries.reserve(NS);

  double best_post_ber = std::numeric_limits<double>::infinity();
  std::size_t best_index = static_cast<std::size_t>(-1);
  PipelineResult best_result{};
  std::vector<double> best_tile_early_stop_pct;

  std::vector<ScenarioOutput> results;
  results.reserve(futures.size());
  for (auto& fut : futures) {
    try {
      results.emplace_back(fut.get());
    } catch (const std::exception& ex) {
      out << "[ERROR] worker threw: " << ex.what() << "\n";
    } catch (...) {
      out << "[ERROR] worker threw unknown exception\n";
    }
  }

  std::sort(results.begin(), results.end(),
            [](const auto& a, const auto& b) { return a.idx < b.idx; });

  std::ofstream csv(csv_path, std::ios::out | std::ios::app);
  csv.setf(std::ios::fixed);
  csv << std::setprecision(8);

  for (const auto& pack : results) {
    const auto& result = pack.result;

    if (!result.tile_early_stop_pct.empty()) {
      out << "[INFO] " << pack.name << " tile early-stop hit rates (%): ";
      out << std::fixed << std::setprecision(1);
      for (size_t i = 0; i < result.tile_early_stop_pct.size(); ++i) {
        out << result.tile_early_stop_pct[i]
            << (i + 1 < result.tile_early_stop_pct.size() ? ", " : "\n");
      }
      out << std::defaultfloat;
    }

    std::string early_stop_summary;
    if (!result.tile_early_stop_pct.empty()) {
      std::ostringstream early_oss;
      early_oss << " | EarlyStop%=[";
      for (size_t i = 0; i < result.tile_early_stop_pct.size(); ++i) {
        early_oss << std::fixed << std::setprecision(1)
                  << result.tile_early_stop_pct[i] << "%";
        if (i + 1 < result.tile_early_stop_pct.size()) early_oss << ", ";
      }
      early_oss << "]";
      early_stop_summary = early_oss.str();
    }

    std::ostringstream summary;
    summary << "[SUMMARY] " << pack.name
            << " Pre-FEC BER=" << result.pre_fec.ber
            << " (errs=" << result.pre_fec.errors << "/" << result.pre_fec.total << ")"
            << " | Post-FEC BER=" << result.post_fec.ber
            << " (errs=" << result.post_fec.errors << "/" << result.post_fec.total << ")"
            << early_stop_summary;
    scenario_summaries.push_back(summary.str());
    out << summary.str() << "\n";

    const double alpha_first = pack.alpha_list.empty()
                                   ? std::numeric_limits<double>::quiet_NaN()
                                   : pack.alpha_list.front();
    const double beta_first  = pack.beta_list.empty()
                                   ? std::numeric_limits<double>::quiet_NaN()
                                   : pack.beta_list.front();
    const double es_mean = mean(result.tile_early_stop_pct);

    csv << now_stamp() << ","
        << run_id << ","
        << pack.name << ","
        << alpha_first << ","
        << beta_first  << ","
        << '"' << join_vec(pack.alpha_list, '|', 6) << "\","
        << '"' << join_vec(pack.beta_list , '|', 6) << "\","
        << result.pre_fec.ber    << ","
        << result.pre_fec.errors << ","
        << result.pre_fec.total  << ","
        << result.post_fec.ber   << ","
        << result.post_fec.errors<< ","
        << result.post_fec.total << ","
        << std::setprecision(3) << es_mean << ","
        << '"' << join_vec(result.tile_early_stop_pct, '|', 1) << "\"\n";
    csv << std::setprecision(8);
    csv.flush();

    if (result.post_fec.total > 0 && result.post_fec.ber < best_post_ber) {
      best_post_ber = result.post_fec.ber;
      best_index = pack.idx;
      best_result = result;
      best_tile_early_stop_pct = result.tile_early_stop_pct;
    }
  }

  if (best_index == static_cast<std::size_t>(-1)) {
    out << "[ERROR] No valid scenarios evaluated.\n";
    return 1;
  }

  out << "\n[SUMMARY] All scenarios:\n";
  for (const auto& line : scenario_summaries) out << "  " << line << '\n';

  const auto& best_scenario = scenarios[best_index];
  out << "\n[RESULT] Best scenario: " << best_scenario.name
      << " with Post-FEC BER=" << best_result.post_fec.ber
      << " (errs=" << best_result.post_fec.errors << "/" << best_result.post_fec.total << ")\n";

  if (!best_tile_early_stop_pct.empty()) {
    out << "[RESULT] Best tile early-stop hit rates (%): ";
    out << std::fixed << std::setprecision(1);
    for (size_t i = 0; i < best_tile_early_stop_pct.size(); ++i) {
      out << best_tile_early_stop_pct[i]
          << (i + 1 < best_tile_early_stop_pct.size() ? ", " : "\n");
    }
    out << std::defaultfloat;
  }

  out << "[RESULT] Best ALPHA_LIST: ";
  for (size_t i = 0; i < best_scenario.alpha_list.size(); ++i) {
    out << best_scenario.alpha_list[i]
        << (i + 1 < best_scenario.alpha_list.size() ? ", " : "\n");
  }

  out << "[RESULT] Best beta_list: ";
  for (size_t i = 0; i < best_scenario.beta_list.size(); ++i) {
    out << best_scenario.beta_list[i]
        << (i + 1 < best_scenario.beta_list.size() ? ", " : "\n");
  }

  return 0;
}

