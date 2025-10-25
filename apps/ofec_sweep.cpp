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
#include <cmath>

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
  fout << "timestamp,run_id,scenario,alpha_start,alpha_step,beta_start,beta_step,chase_L,chase_n_test,ebn0_db,ALPHA_LIST,beta_list,"
          "pre_ber,pre_errs,pre_total,post_ber,post_errs,post_total,"
          "early_stop_mean_pct,early_stop_list\n";
}

static void ensure_ebn0_csv_header(const std::string& csv_path)
{
  std::ifstream fin(csv_path);
  if (fin.good() && fin.peek() != std::ifstream::traits_type::eof()) return;

  std::ofstream fout(csv_path, std::ios::out | std::ios::app);
  fout << "timestamp,run_id,scenario,ebn0_db,pre_ber,pre_errs,pre_total,post_ber,post_errs,post_total,"
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
  float alpha_start = 0.0f;
  float alpha_step  = 0.0f;
  float beta_start  = 0.0f;
  float beta_step   = 0.0f;
  int   chase_L     = 0;
  int   chase_n_test = 0;
  float ebn0_db     = DEFAULT_EBN0_DB;
};

static std::vector<float> generate_sequence(float start, float step, std::size_t length)
{
  std::vector<float> seq(length, start);
  for (std::size_t i = 0; i < length; ++i) {
    seq[i] = start + step * static_cast<float>(i);
  }
  return seq;
}

static float infer_step(const std::vector<float>& values)
{
  return values.size() >= 2 ? values[1] - values[0] : 0.0f;
}

static std::vector<SweepScenario> build_scenarios(const Params& base_params,
                                                  const std::vector<float>& alpha_starts,
                                                  const std::vector<float>& alpha_steps,
                                                  const std::vector<float>& beta_starts,
                                                  const std::vector<float>& beta_steps,
                                                  const std::vector<int>& chase_l_candidates,
                                                  const std::vector<float>& ebn0_candidates)
{
  std::vector<SweepScenario> scenarios;

  SweepScenario baseline;
  baseline.name = "baseline";
  baseline.alpha_list = base_params.ALPHA_LIST;
  baseline.beta_list  = base_params.beta_list;
  if (!baseline.alpha_list.empty()) {
    baseline.alpha_start = baseline.alpha_list.front();
    baseline.alpha_step  = infer_step(baseline.alpha_list);
  }
  if (!baseline.beta_list.empty()) {
    baseline.beta_start = baseline.beta_list.front();
    baseline.beta_step  = infer_step(baseline.beta_list);
  }
  baseline.chase_L = base_params.CHASE_L;
  baseline.chase_n_test = 1 << base_params.CHASE_L;
  baseline.ebn0_db = !ebn0_candidates.empty() ? ebn0_candidates.front() : DEFAULT_EBN0_DB;
  scenarios.push_back(std::move(baseline));

  const std::size_t len = base_params.TILES_PER_WIN;
  for (float a_start : alpha_starts) {
    for (float a_step : alpha_steps) {
      for (float b_start : beta_starts) {
        for (float b_step : beta_steps) {
          for (int chase_L : chase_l_candidates) {
            for (float ebn0_db : ebn0_candidates) {
              SweepScenario sc;
              sc.alpha_start = a_start;
              sc.alpha_step  = a_step;
              sc.beta_start  = b_start;
              sc.beta_step   = b_step;
              sc.chase_L     = chase_L;
              sc.chase_n_test = 1 << chase_L;
              sc.ebn0_db     = ebn0_db;
              sc.alpha_list  = generate_sequence(a_start, a_step, len);
              sc.beta_list   = generate_sequence(b_start, b_step, len);

              const SweepScenario& base = scenarios.front();
              if (sc.alpha_list == base.alpha_list &&
                  sc.beta_list  == base.beta_list &&
                  sc.chase_L    == base.chase_L &&
                  std::fabs(sc.ebn0_db - base.ebn0_db) < 1e-6f) {
                continue; // 与基线重复，跳过
              }

              std::ostringstream oss;
              oss << std::fixed << std::setprecision(3)
                  << "alphaS" << a_start << "_d" << a_step
                  << "_betaS" << b_start << "_d" << b_step
                  << "_chL" << chase_L
                  << "_EbN0_" << ebn0_db;
              sc.name = oss.str();

              scenarios.push_back(std::move(sc));
            }
          }
        }
      }
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
  float alpha_start = 0.0f;
  float alpha_step  = 0.0f;
  float beta_start  = 0.0f;
  float beta_step   = 0.0f;
  int   chase_L     = 0;
  int   chase_n_test = 0;
  float ebn0_db     = DEFAULT_EBN0_DB;
  PipelineResult result;
};

struct EbN0Output {
  float ebn0_db = std::numeric_limits<float>::quiet_NaN();
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

  const std::string csv_path = (data_dir / ("ofec_sweep_results_" + run_id + ".csv")).string();
  ensure_csv_header(csv_path);

  // ========== 2) 生成 scenario ==========
  Params base_params;

  // 示例候选集（起点与步进，可按需调整）
  const std::vector<float> alpha_start_candidates = {0.01f,0.05f,0.1f,0.15f,0.2f,0.25f,0.3f,0.4f};
  const std::vector<float> alpha_step_candidates  = {0.0f,0.1f};
  const std::vector<float> beta_start_candidates  = {0.1f,0.2f,0.3f,0.4f,0.5f,0.6f,0.7f,0.9f,1.2f,1.5f};
  const std::vector<float> beta_step_candidates   = {0.0f,0.1f};
  const std::vector<int>   chase_l_candidates     = {4,6};

  // EbN0 扫描范围配置：起点、终点以及取样点数（均匀分布）
  const float ebn0_start = 3.1f;
  const float ebn0_end   = 3.5f;
  const int   ebn0_points = 8;

  std::vector<float> sweep_ebn0_values;
  if (ebn0_points <= 0) {
    sweep_ebn0_values.push_back(DEFAULT_EBN0_DB);
  } else if (ebn0_points == 1) {
    sweep_ebn0_values.push_back(ebn0_start);
  } else {
    const float eb_step = (ebn0_end - ebn0_start) / static_cast<float>(ebn0_points - 1);
    for (int i = 0; i < ebn0_points; ++i) {
      sweep_ebn0_values.push_back(ebn0_start + eb_step * static_cast<float>(i));
    }
  }

  auto scenarios = build_scenarios(base_params,
                                   alpha_start_candidates, alpha_step_candidates,
                                   beta_start_candidates, beta_step_candidates,
                                   chase_l_candidates, sweep_ebn0_values);
  const std::size_t NS = scenarios.size();

  out << "[INFO] total scenarios = " << NS << "\n";

  // ========== 3) 并行执行 run_pipeline ==========
  unsigned hw = std::max(1u, std::thread::hardware_concurrency());
  unsigned max_workers = hw * 3 / 4;  // 默认值
  if (const char* env = std::getenv("NTHREADS")) {
      max_workers = std::stoi(env);
  }
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
    params.CHASE_L     = scenario.chase_L;
    params.CHASE_NTEST = scenario.chase_n_test;

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
                                      out.alpha_start = scenario.alpha_start;
                                      out.alpha_step  = scenario.alpha_step;
                                      out.beta_start  = scenario.beta_start;
                                      out.beta_step   = scenario.beta_step;
                                      out.chase_L     = scenario.chase_L;
                                      out.chase_n_test = scenario.chase_n_test;
                                      out.ebn0_db     = scenario.ebn0_db;

                                      Params local_params = params;
                                      out.result = run_pipeline(local_params, scenario.name, scenario.ebn0_db);
                                      out.ebn0_db = out.result.ebn0_db;
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
  float best_alpha_start = 0.0f;
  float best_alpha_step  = 0.0f;
  float best_beta_start  = 0.0f;
  float best_beta_step   = 0.0f;
  int   best_chase_L     = base_params.CHASE_L;
  int   best_chase_n_test = 1 << base_params.CHASE_L;
  float best_ebn0_db      = !sweep_ebn0_values.empty() ? sweep_ebn0_values.front() : DEFAULT_EBN0_DB;

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
            << " (errs=" << result.post_fec.errors << "/" << result.post_fec.total << ")";
    summary << std::fixed << std::setprecision(3)
            << " | alpha_start=" << pack.alpha_start
            << " alpha_step=" << pack.alpha_step
            << " | beta_start=" << pack.beta_start
            << " beta_step=" << pack.beta_step
            << " | Eb/N0=" << pack.ebn0_db
            << std::defaultfloat
            << " | CHASE_L=" << pack.chase_L
            << " CHASE_NTEST=" << pack.chase_n_test
            << early_stop_summary;
    scenario_summaries.push_back(summary.str());
    out << summary.str() << "\n";

    const double es_mean = mean(result.tile_early_stop_pct);

    csv << now_stamp() << ","
        << run_id << ","
        << pack.name << ","
        << pack.alpha_start << ","
        << pack.alpha_step  << ","
        << pack.beta_start  << ","
        << pack.beta_step   << ","
        << pack.chase_L     << ","
        << pack.chase_n_test<< ","
        << pack.ebn0_db     << ","
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
      best_alpha_start = pack.alpha_start;
      best_alpha_step  = pack.alpha_step;
      best_beta_start  = pack.beta_start;
      best_beta_step   = pack.beta_step;
      best_chase_L     = pack.chase_L;
      best_chase_n_test = pack.chase_n_test;
      best_ebn0_db      = pack.ebn0_db;
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
  out << std::fixed << std::setprecision(3);
  out << "[RESULT] Best alpha start/step: start=" << best_alpha_start
      << " step=" << best_alpha_step << "\n";
  out << "[RESULT] Best beta start/step: start=" << best_beta_start
      << " step=" << best_beta_step << "\n";
  out << "[RESULT] Best sweep Eb/N0: " << best_ebn0_db << " dB\n";
  out << std::defaultfloat;
  out << "[RESULT] Best CHASE_L/CHASE_NTEST: " << best_chase_L
      << " / " << best_chase_n_test << "\n";

  if (!best_tile_early_stop_pct.empty()) {
    out << "[RESULT] Best tile early-stop hit rates (%): ";
    out << std::fixed << std::setprecision(1);
    for (size_t i = 0; i < best_tile_early_stop_pct.size(); ++i) {
      out << best_tile_early_stop_pct[i]
          << (i + 1 < best_tile_early_stop_pct.size() ? ", " : "\n");
    }
    out << std::defaultfloat;
  }

  // ========== 5) 使用最佳参数进行 EbN0 扫描 ==========
  std::vector<float> ebn0_values = sweep_ebn0_values;
  if (ebn0_values.empty()) {
    ebn0_values.push_back(best_ebn0_db);
  }

  if (!ebn0_values.empty()) {
    const std::string ebn0_log_path = (data_dir / ("run_" + run_id + "_ebn0.log")).string();
    DualOut ebn0_out(std::cout, ebn0_log_path);

    const std::string ebn0_csv_path = (data_dir / ("ofec_ebn0_results_" + run_id + ".csv")).string();
    ensure_ebn0_csv_header(ebn0_csv_path);
    std::ofstream ebn0_csv(ebn0_csv_path, std::ios::out | std::ios::app);
    ebn0_csv.setf(std::ios::fixed);
    ebn0_csv << std::setprecision(8);

    Params best_params = base_params;
    best_params.ALPHA_LIST = best_scenario.alpha_list;
    best_params.beta_list  = best_scenario.beta_list;
    if (!best_params.ALPHA_LIST.empty()) best_params.ALPHA = best_params.ALPHA_LIST.front();
    if (!best_params.beta_list.empty())  best_params.beta  = best_params.beta_list.front();
    best_params.CHASE_L     = best_chase_L;
    best_params.CHASE_NTEST = best_chase_n_test;

    Semaphore ebn0_sem(max_workers);
    std::vector<std::future<EbN0Output>> ebn0_futures;
    ebn0_futures.reserve(ebn0_values.size());

    for (float ebn0_db : ebn0_values) {
      ebn0_sem.acquire();
      ebn0_futures.emplace_back(std::async(std::launch::async,
                                           [best_params, ebn0_db, base_label = best_scenario.name, &ebn0_sem]() mutable -> EbN0Output {
                                             struct Releaser {
                                               Semaphore& s;
                                               ~Releaser() { s.release(); }
                                             } _r{ebn0_sem};
                                             Params params = best_params;
                                             std::ostringstream oss;
                                             oss << base_label << "_EbN0_" << std::fixed << std::setprecision(2) << ebn0_db;
                                             std::string label = oss.str();
                                             PipelineResult result = run_pipeline(params, label, ebn0_db);
                                             return EbN0Output{ebn0_db, std::move(result)};
                                           }));
    }

    std::vector<EbN0Output> ebn0_results;
    ebn0_results.reserve(ebn0_futures.size());
    for (auto& fut : ebn0_futures) {
      try {
        ebn0_results.emplace_back(fut.get());
      } catch (const std::exception& ex) {
        ebn0_out << "[EbN0][ERROR] worker threw: " << ex.what() << "\n";
      } catch (...) {
        ebn0_out << "[EbN0][ERROR] worker threw unknown exception\n";
      }
    }

    std::sort(ebn0_results.begin(), ebn0_results.end(),
              [](const EbN0Output& a, const EbN0Output& b) { return a.ebn0_db < b.ebn0_db; });

    if (!ebn0_results.empty()) {
      ebn0_out << "\n[EbN0] Sweep results for best scenario (" << best_scenario.name << ")\n";
    }

    for (const auto& entry : ebn0_results) {
      const auto& res = entry.result;
      std::vector<double> es_vals(res.tile_early_stop_pct.begin(), res.tile_early_stop_pct.end());
      const double es_mean = mean(es_vals);

      ebn0_out << std::fixed << std::setprecision(2)
               << "[EbN0] " << entry.ebn0_db << " dB | Pre-BER="
               << std::setprecision(8) << res.pre_fec.ber
               << " | Post-BER=" << res.post_fec.ber << "\n";
      if (!es_vals.empty()) {
        ebn0_out << "        EarlyStop%: "
                 << std::fixed << std::setprecision(1) << join_vec(es_vals, ',', 1) << "\n";
      }
      ebn0_out << std::defaultfloat;

      ebn0_csv << now_stamp() << ","
               << run_id << ","
               << best_scenario.name << ","
               << entry.ebn0_db << ","
               << res.pre_fec.ber    << ","
               << res.pre_fec.errors << ","
               << res.pre_fec.total  << ","
               << res.post_fec.ber   << ","
               << res.post_fec.errors<< ","
               << res.post_fec.total << ","
               << std::setprecision(3) << es_mean << ","
               << "\"" << join_vec(es_vals, '|', 1) << "\"\n";
      ebn0_csv << std::setprecision(8);
    }

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
