#include <algorithm>
#include <filesystem>
#include <fstream>
#include <future>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <string>
#include <thread>
#include <vector>

#include "mux_bypass_edges.hpp"
#include "newcode/io/dualwriter.hpp"
#include "newcode/ofec_single_runner.hpp"
#include "newcode/two_stream_shared_runner.hpp"
#include "newcode/utils/now_stamp.hpp"

namespace {

// ======== 用户可调参数区域 ========

// probe 运行标签 / 输出位置
static constexpr const char* kLabelPrefix = "two_stream_ber_window_probe"; // 每个 seed 任务的标签前缀
static constexpr const char* kOutputDir = "data/two_stream_ber_window_probe"; // probe CSV 和日志输出目录
static constexpr bool kQuietPipeline = true;           // true=静默底层 pipeline 日志，false=打印更多过程信息
static constexpr unsigned kMaxParallelEbN0 = 0;        // Eb/N0 维度的最大并行数；0=按硬件线程数自动推断
static constexpr unsigned kMaxParallelSeeds = 0;       // seed 维度的最大并行数；0=按硬件线程数自动推断

// 扫描维度
static const std::vector<float> kEbN0List = {3.05f};   // probe 要运行的 Eb/N0 列表
static constexpr int kSeedCount = 1;                   // 每个 Eb/N0 下运行多少组种子

// 双流各自的随机种子基线
static constexpr int kBitgenSeedBaseA = 20260601;      // Stream A 的信息比特随机种子基线
static constexpr int kChannelSeedBaseA = 6012026;      // Stream A 的信道噪声随机种子基线
static constexpr int kBitgenSeedBaseB = 20260529;      // Stream B 的信息比特随机种子基线
static constexpr int kChannelSeedBaseB = 5292026;      // Stream B 的信道噪声随机种子基线

// 信道 / 前端参数
static constexpr std::size_t kNumInfoBits =
    128u * 132u * 16u * 111u;                           // 每一路独立数据流生成的信息比特总数
static constexpr bool kGenerateRandomBits = true;      // true=发送随机信息比特，false=发送全 0 比特
static constexpr unsigned kBitsPerSymbol = 1;          // 每个调制符号携带的比特数：1=BPSK，偶数=QAM

// Early-stop 参数
static constexpr bool kEnableEarlyStop = true;         // true=启用 early-stop，false=完全关闭
static const std::vector<int> kEarlyStopEnableList = {1, 1, 1, 1, 1, 1}; // 按 tile 覆盖 early-stop 总开关：0=关，非 0=开；空表示全部沿用 kEnableEarlyStop
static constexpr int kEarlyStopConditionMode = 1;      // early-stop 条件模式：1=v1，2=v2
static const std::vector<int> kEarlyStopConditionModeList = {}; // 按 tile 覆盖条件模式；空表示全部沿用 kEarlyStopConditionMode
static constexpr int kEarlyStopActionMode = 1;         // 命中 early-stop 后的动作模式：1=sign beta，2=residual only，3=hard-decode sign LLR，4=sign beta pre-div alpha，5=residual pre-div alpha plus sign beta，6=sign beta without BCH hard-decode
static const std::vector<int> kEarlyStopActionModeList = {}; // 按 tile 覆盖动作模式；空表示全部沿用 kEarlyStopActionMode
static constexpr int kEarlyStopBindGroupSize = 1;      // 条件1专用：多少个 row 绑定为一组；1=逐 row，4=整组都通过才 early-stop
static const std::vector<int> kEarlyStopBindGroupSizeList = {}; // 按 tile 覆盖绑定组大小；空表示全部沿用 kEarlyStopBindGroupSize
static constexpr bool kEarlyStopCondV1RequireBch = true; // 条件1里是否要求 BCH syndrome 为 0
static constexpr bool kEarlyStopCondV1RequireOverall = true; // 条件1里是否要求 overall parity 一致
static constexpr float kEarlyStopV2LlrAbsThreshold = 0.5f; // 条件2里把 bit 视为“不可靠”的 |LLR| 阈值
static constexpr int kEarlyStopV2MaxUnreliableBits = 8; // 条件2里允许的不可靠 bit 数上限
static constexpr bool kEarlyStopCondV2IncludeOverall = true; // 条件2统计不可靠 bit 时是否把 overall bit 算进去

// 量化参数
static constexpr std::size_t kLlrBits = 6;             // LLR 位宽：16=浮点，2~15=qfloat
static constexpr float kQuantClipRatio = 0.5f;         // 动态 clip 比例，0 表示禁用自适应 clip

// Chase / decoder core 参数
static constexpr const char* kInterleaverName = "identity"; // 交织器名称，identity 表示不交织
static constexpr const char* kDecoderName =
    "two_stream_shared_chase_baseline";                // 解码器名称：当前 two-stream shared 主流程入口；底层可映射到 chase_baseline / chase_overall_parity_search / chase_topk_pruned / chase_global_pair / chase_group_minima
static constexpr int kChaseL = 6;                      // Chase L：选择多少个最不可靠位置
static constexpr int kChaseNTestOverride = -1;         // Chase 测试序列数量，-1 表示默认按 2^L 生成
static constexpr int kChaseTopkKeep = 8;               // chase_topk_pruned 中保留参与外信息计算的 Top-K 候选数
static constexpr int kChaseGroupMinimaBits = 3;        // chase_group_minima 按前多少个 test-pattern 位分组
static constexpr bool kNormalizeExtrinsic = false;     // 是否对 Chase 输出的 extrinsic 做归一化
static constexpr bool kNormalizeKnownPrefixTail = false; // 是否对 known-prefix 之后的尾部 LLR 做归一化
static constexpr float kEarlyStopActionResidualDivisor = 1.0f; // 动作2里 residual 的除数
static constexpr float kEarlyStopActionHardLlrMag = 1.0f; // 动作3里硬解成功后输出的固定 |LLR| 幅度

// alpha / beta / shared SISO 预算
static const std::vector<float> kAlphaExplicit = {
    0.428571f, 0.447738f, 0.482782f, 0.528162f, 0.581902f, 0.642857f}; // 每个 tile 的 alpha 显式列表
static const std::vector<float> kBetaExplicit = {
    2.857143f, 6.179301f, 12.253626f, 20.119585f, 29.434408f, 40.0f}; // 每个 tile 的 Chase/fallback beta 显式列表
static const std::vector<float> kEarlyStopActionBetaExplicit = {
    99.857143f, 99.179301f, 99.253626f, 99.119585f, 99.434408f, 99.0f}; // 每个 tile 的 early-stop 动作 beta 显式列表；空表示沿用 beta 配置
static const std::vector<int> kSisoActiveList = {64, 64, 64, 48, 16, 8}; // 每个 tile 在 shared 64-row 域里可参与 SISO 的预算

// MUX 调度参数
static constexpr int kMuxGroupG = 1;                   // MUX 分组数：1=全局池化，>1=按组平均切预算
static constexpr int kMuxSchedulingMode = 0;           // MUX 调度模式：0=legacy 顺序裁剪，1=按 early-stop 细节排序
static constexpr int kMuxPriorityRule = 0;             // MUX 优先级规则：0=harder_first，1=near_threshold_first
static constexpr bool kMuxEnableReconfig = false;      // true=启用 staged reconfig 调度，false=直接预算裁剪
static const std::vector<newcode::mux::MuxEdge> kMuxExtraBypassEdges =
    app_mux::bypass_edges_for_scheme(1);               // reconfig 模式使用的额外旁路边集合；当前 scheme 1

// Hybrid 参数
static constexpr bool kHybridEnable = true;            // true=启用 hybrid prepass，false=全部保留到 soft path
static const std::vector<int> kHybridEnableList = {0, 0, 0, 1, 1, 1}; // 按 tile 覆盖 hybrid 开关：0=关，非 0=开；空表示全部沿用 kHybridEnable
static constexpr float kHybridHardLlrMag = 99.0f;     // hybrid hard-finish 默认输出 |LLR| 幅度
static const std::vector<float> kHybridHardLlrMagList = {}; // 按 tile 覆盖 hybrid hard-finish |LLR| 幅度；空表示全部沿用 kHybridHardLlrMag
static constexpr newcode::HybridClassifierMode kHybridClassifierMode =
    newcode::HybridClassifierMode::FriendS1S3WithS0Classifier;  // hybrid 分类器模式：LegacyHardDecode / RepoFastClassifier / FriendS1S3Classifier / FriendS1S3WithS0Classifier
static constexpr newcode::HybridSisoBackfillMode kHybridSisoBackfillMode =
    newcode::HybridSisoBackfillMode::OneAndTwoErrorPriority;     // hybrid SISO 回填模式：Disabled / TwoErrorOnly / OneAndTwoErrorPriority / ParityOneAndTwoErrorPriority
static constexpr bool kHybridNormalizeSoftOnly = false; // true=只归一化 soft rows，false=对所有 produced rows 保持兼容行为

// ==================================

struct SeedTask {
  int seed_index = 0;
  int bitgen_seed_a = 0;
  int channel_seed_a = 0;
  int bitgen_seed_b = 0;
  int channel_seed_b = 0;
  float ebn0_db = 0.0f;
  std::string label;
  newcode::two_stream_shared::Config config;
};

struct SeedRunResult {
  int seed_index = 0;
  int bitgen_seed_a = 0;
  int channel_seed_a = 0;
  int bitgen_seed_b = 0;
  int channel_seed_b = 0;
  newcode::two_stream_shared::Result result;
};

struct EbN0RunResult {
  float ebn0_db = 0.0f;
  std::vector<SeedRunResult> seed_results;
};

struct RawTileBerRow {
  std::string run_id;
  float ebn0_db = 0.0f;
  int seed_index = 0;
  int chunk_index = 0;
  int bitgen_seed = 0;
  int channel_seed = 0;
  const char* stream_label = "A";
  std::size_t tile_index = 0;
  std::size_t pre_errs = 0;
  std::size_t pre_bits = 0;
  double pre_ber = 0.0;
  std::size_t post_errs = 0;
  std::size_t post_bits = 0;
  double post_ber = 0.0;
};

struct AggregatedTileBerRow {
  std::string run_id;
  float ebn0_db = 0.0f;
  const char* stream_label = "A";
  std::size_t tile_index = 0;
  std::size_t agg_pre_errs = 0;
  std::size_t agg_pre_bits = 0;
  std::size_t agg_post_errs = 0;
  std::size_t agg_post_bits = 0;
};

std::string format_float(float value) {
  std::ostringstream oss;
  oss << std::fixed << std::setprecision(6) << value;
  return oss.str();
}

const char* shared_hybrid_class_name(
    newcode::two_stream_shared::SharedHybridClass row_class) {
  using Class = newcode::two_stream_shared::SharedHybridClass;
  switch (row_class) {
    case Class::None:
      return "none";
    case Class::BchHardDecoded:
      return "bch_hard_decoded";
    case Class::Clean:
      return "clean";
    case Class::ParityOnly:
      return "parity_only";
    case Class::OneMain:
      return "one_main";
    case Class::OneMainPlusParity:
      return "one_main_plus_parity";
    case Class::TwoMain:
      return "two_main";
    case Class::Suspicious:
      return "suspicious";
    case Class::HardFail:
      return "hard_fail";
  }
  return "unknown";
}

const char* shared_row_final_tag_name(
    newcode::two_stream_shared::SharedRowFinalTag tag) {
  using Tag = newcode::two_stream_shared::SharedRowFinalTag;
  switch (tag) {
    case Tag::SoftDecode:
      return "soft_decode";
    case Tag::EarlyStopAction:
      return "early_stop_action";
    case Tag::HardFinish:
      return "hard_finish";
    case Tag::Unscheduled:
      return "unscheduled";
  }
  return "unknown";
}

unsigned resolved_parallel_ebn0() {
  if (kMaxParallelEbN0 > 0) {
    return kMaxParallelEbN0;
  }
  const unsigned hc = std::thread::hardware_concurrency();
  if (hc == 0) {
    return 1u;
  }
  return std::min<unsigned>(hc,
                            static_cast<unsigned>(std::max<std::size_t>(
                                1, kEbN0List.size())));
}

unsigned resolved_parallel_seeds_for_ebn0(unsigned ebn0_parallelism) {
  if (kMaxParallelSeeds > 0) {
    return kMaxParallelSeeds;
  }
  const unsigned hc = std::thread::hardware_concurrency();
  if (hc == 0) {
    return 1u;
  }
  return std::max(1u, hc / std::max(1u, ebn0_parallelism));
}

void write_tile_samples_header(std::ofstream& out) {
  out << "run_id,ebn0_db,seed_index,chunk_index,invocation,tile_index,"
         "stream_rows_a,stream_rows_b,rows_total,"
         "rows_early_stop,rows_not_early_stop,rows_hard_finish,"
         "rows_need_siso_before_mux,rows_soft_scheduled,rows_unscheduled,"
         "produced_rows,failed_rows,produced_rows_a,produced_rows_b,"
         "failed_rows_a,failed_rows_b\n";
}

void write_hybrid_class_counts_header(std::ofstream& out) {
  out << "run_id,ebn0_db,seed_index,chunk_index,invocation,tile_index,"
         "rows_seen_by_hybrid,class_none_count,class_bch_hard_decoded_count,"
         "class_clean_count,class_parity_only_count,class_one_main_count,"
         "class_one_main_plus_parity_count,class_two_main_count,"
         "class_suspicious_count,class_hard_fail_count,"
         "deferred_candidate_count,deferred_priority_0_count,"
         "deferred_priority_1_count,deferred_priority_2_count,"
         "deferred_priority_3_count,"
         "deferred_reclaimed_to_hard_finish_count\n";
}

void write_row_map_header(std::ofstream& out) {
  out << "run_id,ebn0_db,seed_index,chunk_index,invocation,tile_index,"
         "merged_row,stream_id,stream_label,source_local_row,source_global_row,"
         "early_stop_hit,hybrid_class,final_tag,scheduled_for_soft,produced_row\n";
}

void write_seed_summary_header(std::ofstream& out) {
  out << "run_id,ebn0_db,seed_index,chunk_index,"
         "bitgen_seed_a,channel_seed_a,bitgen_seed_b,channel_seed_b,"
         "pre_ber_a,post_ber_a,pre_ber_b,post_ber_b,"
         "tile_sample_rows,hybrid_count_rows,row_map_rows\n";
}

void write_raw_tile_ber_header(std::ofstream& out) {
  out << "run_id,ebn0_db,seed_index,chunk_index,bitgen_seed,channel_seed,"
         "stream_label,tile_index,pre_errs,pre_bits,pre_ber,"
         "post_errs,post_bits,post_ber\n";
}

void write_aggregated_tile_ber_header(std::ofstream& out) {
  out << "run_id,ebn0_db,stream_label,tile_index,"
         "agg_pre_errs,agg_pre_bits,agg_pre_ber,"
         "agg_post_errs,agg_post_bits,agg_post_ber\n";
}

ofec_single::Config make_base_config() {
  ofec_single::Config config{};
  config.label = kLabelPrefix;
  config.ebn0_db = kEbN0List.front();
  config.chaseL_override = kChaseL;
  config.chase_n_test_override = kChaseNTestOverride;
  config.chase_topk_keep = kChaseTopkKeep;
  config.chase_group_minima_bits = kChaseGroupMinimaBits;
  config.normalize_extrinsic = kNormalizeExtrinsic;
  config.bits_per_symbol = kBitsPerSymbol;
  config.bitgen_seed = kBitgenSeedBaseA;
  config.channel_seed = kChannelSeedBaseA;
  config.enable_early_stop = kEnableEarlyStop;
  config.early_stop_enable_list = kEarlyStopEnableList;
  config.early_stop_condition_mode = kEarlyStopConditionMode;
  config.early_stop_condition_mode_list = kEarlyStopConditionModeList;
  config.early_stop_action_mode = kEarlyStopActionMode;
  config.early_stop_action_mode_list = kEarlyStopActionModeList;
  config.early_stop_bind_group_size = kEarlyStopBindGroupSize;
  config.early_stop_bind_group_size_list = kEarlyStopBindGroupSizeList;
  config.early_stop_cond_v1_require_bch = kEarlyStopCondV1RequireBch;
  config.early_stop_cond_v1_require_overall = kEarlyStopCondV1RequireOverall;
  config.early_stop_v2_llr_abs_threshold = kEarlyStopV2LlrAbsThreshold;
  config.early_stop_v2_max_unreliable_bits = kEarlyStopV2MaxUnreliableBits;
  config.early_stop_cond_v2_include_overall = kEarlyStopCondV2IncludeOverall;
  config.alpha_explicit = kAlphaExplicit;
  config.beta_explicit = kBetaExplicit;
  config.early_stop_action_sign_beta_explicit =
      kEarlyStopActionBetaExplicit;
  config.early_stop_action_residual_divisor =
      kEarlyStopActionResidualDivisor;
  config.early_stop_action_hard_llr_mag = kEarlyStopActionHardLlrMag;
  config.siso_active_list = kSisoActiveList;
  config.mux_group_g = kMuxGroupG;
  config.mux_scheduling_mode = kMuxSchedulingMode;
  config.mux_early_stop_priority_rule = kMuxPriorityRule;
  config.mux_enable_reconfig = kMuxEnableReconfig;
  config.mux_extra_bypass_edges = kMuxExtraBypassEdges;
  config.hybrid_enable = kHybridEnable;
  config.hybrid_enable_list = kHybridEnableList;
  config.hybrid_hard_llr_mag = kHybridHardLlrMag;
  config.hybrid_hard_llr_mag_list = kHybridHardLlrMagList;
  config.hybrid_classifier_mode = kHybridClassifierMode;
  config.hybrid_siso_backfill_mode = kHybridSisoBackfillMode;
  config.hybrid_normalize_soft_only = kHybridNormalizeSoftOnly;
  config.interleaver_name = kInterleaverName;
  config.decoder_name = kDecoderName;
  config.generate_random_bits = kGenerateRandomBits;
  config.normalize_known_prefix_tail = kNormalizeKnownPrefixTail;
  config.quant_clip_ratio = kQuantClipRatio;
  config.llr_bits = kLlrBits;
  return config;
}

newcode::two_stream_shared::Config make_seed_config(float ebn0_db,
                                                    int seed_index) {
  ofec_single::Config base_cfg = make_base_config();
  base_cfg.ebn0_db = ebn0_db;
  base_cfg.bitgen_seed = kBitgenSeedBaseA + seed_index;
  base_cfg.channel_seed = kChannelSeedBaseA + seed_index;
  base_cfg.label = std::string(kLabelPrefix) + "_e" + format_float(ebn0_db) +
                   "_s" + std::to_string(seed_index);
  std::ofstream null_file;
  io::DualWriter null_log(null_file);
  auto params_opt = ofec_single::detail::build_params(base_cfg, null_log);
  if (!params_opt.has_value()) {
    throw std::runtime_error("failed to build params for two-stream probe");
  }
  params_opt->NUM_INFO_BITS = kNumInfoBits;

  newcode::PipelineConfig pipeline =
      ofec_single::detail::build_pipeline_config(base_cfg);
  pipeline.quiet = kQuietPipeline;

  return newcode::two_stream_shared::Config{
      .params = *params_opt,
      .pipeline = pipeline,
      .stream_a =
          {
              .label = base_cfg.label + "_A",
              .bitgen_seed = kBitgenSeedBaseA + seed_index,
              .channel_seed = kChannelSeedBaseA + seed_index,
              .ebn0_db = ebn0_db,
          },
      .stream_b =
          {
              .label = base_cfg.label + "_B",
              .bitgen_seed = kBitgenSeedBaseB + seed_index,
              .channel_seed = kChannelSeedBaseB + seed_index,
              .ebn0_db = ebn0_db,
          },
  };
}

std::vector<SeedRunResult> run_seed_tasks(std::vector<SeedTask> tasks,
                                          unsigned max_parallel_seeds) {
  std::vector<SeedRunResult> seed_results;
  seed_results.reserve(tasks.size());
  std::vector<std::future<SeedRunResult>> pending;

  auto launch_task = [](SeedTask task) {
    return std::async(std::launch::async, [task = std::move(task)]() mutable {
      SeedRunResult out;
      out.seed_index = task.seed_index;
      out.bitgen_seed_a = task.bitgen_seed_a;
      out.channel_seed_a = task.channel_seed_a;
      out.bitgen_seed_b = task.bitgen_seed_b;
      out.channel_seed_b = task.channel_seed_b;
      out.result = newcode::two_stream_shared::run_two_stream_shared(task.config);
      return out;
    });
  };

  for (auto& task : tasks) {
    pending.push_back(launch_task(std::move(task)));
    if (pending.size() >= max_parallel_seeds) {
      seed_results.push_back(pending.front().get());
      pending.erase(pending.begin());
    }
  }

  for (auto& future : pending) {
    seed_results.push_back(future.get());
  }

  std::sort(seed_results.begin(), seed_results.end(),
            [](const SeedRunResult& lhs, const SeedRunResult& rhs) {
              return lhs.seed_index < rhs.seed_index;
            });
  return seed_results;
}

EbN0RunResult run_ebn0_probe(float ebn0_db,
                             unsigned seed_parallelism) {
  std::vector<SeedTask> tasks;
  tasks.reserve(kSeedCount);
  for (int seed_index = 0; seed_index < kSeedCount; ++seed_index) {
    auto config = make_seed_config(ebn0_db, seed_index);
    tasks.push_back(SeedTask{
        .seed_index = seed_index,
        .bitgen_seed_a = config.stream_a.bitgen_seed,
        .channel_seed_a = config.stream_a.channel_seed,
        .bitgen_seed_b = config.stream_b.bitgen_seed,
        .channel_seed_b = config.stream_b.channel_seed,
        .ebn0_db = ebn0_db,
        .label = std::string(kLabelPrefix) + "_e" + format_float(ebn0_db) +
                 "_s" + std::to_string(seed_index),
        .config = std::move(config),
    });
  }

  EbN0RunResult out;
  out.ebn0_db = ebn0_db;
  out.seed_results = run_seed_tasks(std::move(tasks), std::max(1u, seed_parallelism));
  return out;
}

}  // namespace

int main() {
  const std::filesystem::path output_dir = kOutputDir;
  std::filesystem::create_directories(output_dir);

  std::ofstream tile_samples_csv(output_dir / "per_invocation_shared_tile_samples.csv");
  std::ofstream hybrid_counts_csv(output_dir / "per_invocation_shared_hybrid_class_counts.csv");
  std::ofstream row_map_csv(output_dir / "per_invocation_shared_row_map.csv");
  std::ofstream seed_summary_csv(output_dir / "per_seed_summary.csv");
  std::ofstream raw_tile_ber_csv(output_dir / "per_seed_per_tile_ber.csv");
  std::ofstream aggregated_tile_ber_csv(output_dir / "aggregated_tile_ber.csv");
  std::ofstream log_file(output_dir / "probe.log");

  write_tile_samples_header(tile_samples_csv);
  write_hybrid_class_counts_header(hybrid_counts_csv);
  write_row_map_header(row_map_csv);
  write_seed_summary_header(seed_summary_csv);
  write_raw_tile_ber_header(raw_tile_ber_csv);
  write_aggregated_tile_ber_header(aggregated_tile_ber_csv);

  io::DualWriter log(log_file);
  log << "[INFO] two-stream shared probe started\n";
  log << "[INFO] seed_count=" << kSeedCount
      << ", max_parallel_ebn0=" << resolved_parallel_ebn0()
      << ", max_parallel_seeds="
      << resolved_parallel_seeds_for_ebn0(resolved_parallel_ebn0()) << '\n';

  const std::string run_id = "two_stream_probe_" + utils::now_stamp();

  const unsigned ebn0_parallel = std::max(1u, resolved_parallel_ebn0());
  const unsigned seed_parallel =
      std::max(1u, resolved_parallel_seeds_for_ebn0(ebn0_parallel));

  std::vector<EbN0RunResult> ebn0_results;
  ebn0_results.reserve(kEbN0List.size());
  std::vector<std::future<EbN0RunResult>> pending_ebn0;

  for (float ebn0_db : kEbN0List) {
    log << "[INFO] schedule Eb/N0=" << ebn0_db
        << " dB, seed_parallelism=" << seed_parallel << '\n';
    pending_ebn0.push_back(std::async(
        std::launch::async,
        [ebn0_db, seed_parallel]() {
          return run_ebn0_probe(ebn0_db, seed_parallel);
        }));
    if (pending_ebn0.size() >= ebn0_parallel) {
      ebn0_results.push_back(pending_ebn0.front().get());
      pending_ebn0.erase(pending_ebn0.begin());
    }
  }
  for (auto& future : pending_ebn0) {
    ebn0_results.push_back(future.get());
  }

  std::sort(ebn0_results.begin(), ebn0_results.end(),
            [](const EbN0RunResult& lhs, const EbN0RunResult& rhs) {
              return lhs.ebn0_db < rhs.ebn0_db;
            });

  std::vector<RawTileBerRow> raw_tile_ber_rows;
  std::vector<AggregatedTileBerRow> aggregated_tile_ber_rows;

  for (const auto& ebn0_result : ebn0_results) {
    std::vector<AggregatedTileBerRow> local_aggregated_tile_rows;
    for (const auto& seed_result : ebn0_result.seed_results) {
      const auto& obs = seed_result.result.observability;
      seed_summary_csv << std::fixed << std::setprecision(6)
                       << run_id << ','
                       << ebn0_result.ebn0_db << ','
                       << seed_result.seed_index << ','
                       << 0 << ','
                       << seed_result.bitgen_seed_a << ','
                       << seed_result.channel_seed_a << ','
                       << seed_result.bitgen_seed_b << ','
                       << seed_result.channel_seed_b << ','
                       << std::setprecision(12)
                       << seed_result.result.stream_a.pre_fec.ber << ','
                       << seed_result.result.stream_a.post_fec.ber << ','
                       << seed_result.result.stream_b.pre_fec.ber << ','
                       << seed_result.result.stream_b.post_fec.ber << ','
                       << obs.tile_samples.size() << ','
                       << obs.hybrid_class_counts.size() << ','
                       << obs.row_map.size() << '\n';

      const auto append_stream_tile_rows =
          [&](const newcode::PipelineResult& stream_result,
              const char* stream_label,
              int bitgen_seed,
              int channel_seed) {
            const std::size_t max_tiles = std::max(
                stream_result.pre_fec_tile_windows.size(),
                stream_result.post_fec_tile_windows.size());
            for (std::size_t tile_index = 0; tile_index < max_tiles; ++tile_index) {
              const auto pre =
                  (tile_index < stream_result.pre_fec_tile_windows.size())
                      ? stream_result.pre_fec_tile_windows[tile_index]
                      : newcode::WindowBerStats{tile_index, 0, 0, 0.0};
              const auto post =
                  (tile_index < stream_result.post_fec_tile_windows.size())
                      ? stream_result.post_fec_tile_windows[tile_index]
                      : newcode::WindowBerStats{tile_index, 0, 0, 0.0};

              raw_tile_ber_rows.push_back(RawTileBerRow{
                  .run_id = run_id,
                  .ebn0_db = ebn0_result.ebn0_db,
                  .seed_index = seed_result.seed_index,
                  .chunk_index = 0,
                  .bitgen_seed = bitgen_seed,
                  .channel_seed = channel_seed,
                  .stream_label = stream_label,
                  .tile_index = tile_index,
                  .pre_errs = pre.errors,
                  .pre_bits = pre.total,
                  .pre_ber = pre.ber,
                  .post_errs = post.errors,
                  .post_bits = post.total,
                  .post_ber = post.ber,
              });
            }
          };

      append_stream_tile_rows(seed_result.result.stream_a,
                              "A",
                              seed_result.bitgen_seed_a,
                              seed_result.channel_seed_a);
      append_stream_tile_rows(seed_result.result.stream_b,
                              "B",
                              seed_result.bitgen_seed_b,
                              seed_result.channel_seed_b);

      for (const auto& sample : obs.tile_samples) {
        tile_samples_csv << run_id << ','
                         << std::fixed << std::setprecision(6)
                         << ebn0_result.ebn0_db << ','
                         << seed_result.seed_index << ','
                         << 0 << ','
                         << sample.invocation << ','
                         << sample.tile_index << ','
                         << sample.stream_rows_a << ','
                         << sample.stream_rows_b << ','
                         << sample.rows_total << ','
                         << sample.rows_early_stop << ','
                         << sample.rows_not_early_stop << ','
                         << sample.rows_hard_finish << ','
                         << sample.rows_need_siso_before_mux << ','
                         << sample.rows_soft_scheduled << ','
                         << sample.rows_unscheduled << ','
                         << sample.produced_rows << ','
                         << sample.failed_rows << ','
                         << sample.produced_rows_a << ','
                         << sample.produced_rows_b << ','
                         << sample.failed_rows_a << ','
                         << sample.failed_rows_b << '\n';
      }

      for (const auto& count : obs.hybrid_class_counts) {
        hybrid_counts_csv << run_id << ','
                          << std::fixed << std::setprecision(6)
                          << ebn0_result.ebn0_db << ','
                          << seed_result.seed_index << ','
                          << 0 << ','
                          << count.invocation << ','
                          << count.tile_index << ','
                          << count.rows_seen_by_hybrid << ','
                          << count.class_none_count << ','
                          << count.class_bch_hard_decoded_count << ','
                          << count.class_clean_count << ','
                          << count.class_parity_only_count << ','
                          << count.class_one_main_count << ','
                          << count.class_one_main_plus_parity_count << ','
                          << count.class_two_main_count << ','
                          << count.class_suspicious_count << ','
                          << count.class_hard_fail_count << ','
                          << count.deferred_candidate_count << ','
                          << count.deferred_priority_0_count << ','
                          << count.deferred_priority_1_count << ','
                          << count.deferred_priority_2_count << ','
                          << count.deferred_priority_3_count << ','
                          << count.deferred_reclaimed_to_hard_finish_count
                          << '\n';
      }

      for (const auto& row : obs.row_map) {
        row_map_csv << run_id << ','
                    << std::fixed << std::setprecision(6)
                    << ebn0_result.ebn0_db << ','
                    << seed_result.seed_index << ','
                    << 0 << ','
                    << row.invocation << ','
                    << row.tile_index << ','
                    << row.merged_row << ','
                    << row.stream_id << ','
                    << (row.stream_id == 0 ? "A" : "B") << ','
                    << row.source_local_row << ','
                    << row.source_global_row << ','
                    << (row.early_stop_hit ? 1 : 0) << ','
                    << shared_hybrid_class_name(row.hybrid_class) << ','
                    << shared_row_final_tag_name(row.final_tag) << ','
                    << (row.scheduled_for_soft ? 1 : 0) << ','
                    << (row.produced_row ? 1 : 0) << '\n';
      }

      log << "[RESULT] Eb/N0=" << ebn0_result.ebn0_db
          << " seed=" << seed_result.seed_index
          << " | A post-FEC BER=" << seed_result.result.stream_a.post_fec.ber
          << " | B post-FEC BER=" << seed_result.result.stream_b.post_fec.ber
          << " | tile_samples=" << obs.tile_samples.size()
          << " | row_map=" << obs.row_map.size() << '\n';
    }

    const auto append_aggregated_stream_tile_rows =
        [&](const char* stream_label, bool stream_a) {
          std::vector<AggregatedTileBerRow> stream_rows;
          for (const auto& seed_result : ebn0_result.seed_results) {
            const auto& stream_result =
                stream_a ? seed_result.result.stream_a : seed_result.result.stream_b;
            const std::size_t max_tiles = std::max(
                stream_result.pre_fec_tile_windows.size(),
                stream_result.post_fec_tile_windows.size());
            if (stream_rows.size() < max_tiles) {
              const std::size_t old_size = stream_rows.size();
              stream_rows.resize(max_tiles);
              for (std::size_t tile_index = old_size; tile_index < max_tiles; ++tile_index) {
                stream_rows[tile_index].run_id = run_id;
                stream_rows[tile_index].ebn0_db = ebn0_result.ebn0_db;
                stream_rows[tile_index].stream_label = stream_label;
                stream_rows[tile_index].tile_index = tile_index;
              }
            }

            for (std::size_t tile_index = 0; tile_index < max_tiles; ++tile_index) {
              const auto pre =
                  (tile_index < stream_result.pre_fec_tile_windows.size())
                      ? stream_result.pre_fec_tile_windows[tile_index]
                      : newcode::WindowBerStats{tile_index, 0, 0, 0.0};
              const auto post =
                  (tile_index < stream_result.post_fec_tile_windows.size())
                      ? stream_result.post_fec_tile_windows[tile_index]
                      : newcode::WindowBerStats{tile_index, 0, 0, 0.0};
              stream_rows[tile_index].agg_pre_errs += pre.errors;
              stream_rows[tile_index].agg_pre_bits += pre.total;
              stream_rows[tile_index].agg_post_errs += post.errors;
              stream_rows[tile_index].agg_post_bits += post.total;
            }
          }
          local_aggregated_tile_rows.insert(local_aggregated_tile_rows.end(),
                                            stream_rows.begin(),
                                            stream_rows.end());
        };

    append_aggregated_stream_tile_rows("A", true);
    append_aggregated_stream_tile_rows("B", false);
    aggregated_tile_ber_rows.insert(aggregated_tile_ber_rows.end(),
                                    local_aggregated_tile_rows.begin(),
                                    local_aggregated_tile_rows.end());
  }

  for (const auto& row : raw_tile_ber_rows) {
    raw_tile_ber_csv << row.run_id << ','
                     << std::fixed << std::setprecision(6)
                     << row.ebn0_db << ','
                     << row.seed_index << ','
                     << row.chunk_index << ','
                     << row.bitgen_seed << ','
                     << row.channel_seed << ','
                     << row.stream_label << ','
                     << row.tile_index << ','
                     << row.pre_errs << ','
                     << row.pre_bits << ','
                     << std::setprecision(12) << row.pre_ber << ','
                     << row.post_errs << ','
                     << row.post_bits << ','
                     << row.post_ber << '\n';
  }

  for (const auto& row : aggregated_tile_ber_rows) {
    const double agg_pre_ber =
        (row.agg_pre_bits == 0)
            ? 0.0
            : static_cast<double>(row.agg_pre_errs) /
                  static_cast<double>(row.agg_pre_bits);
    const double agg_post_ber =
        (row.agg_post_bits == 0)
            ? 0.0
            : static_cast<double>(row.agg_post_errs) /
                  static_cast<double>(row.agg_post_bits);
    aggregated_tile_ber_csv << row.run_id << ','
                            << std::fixed << std::setprecision(6)
                            << row.ebn0_db << ','
                            << row.stream_label << ','
                            << row.tile_index << ','
                            << row.agg_pre_errs << ','
                            << row.agg_pre_bits << ','
                            << std::setprecision(12) << agg_pre_ber << ','
                            << row.agg_post_errs << ','
                            << row.agg_post_bits << ','
                            << agg_post_ber << '\n';
  }

  log << "[INFO] outputs saved under " << output_dir.string() << "\n";
  return 0;
}
