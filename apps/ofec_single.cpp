#include <charconv>
#include <cmath>
#include <cstdlib>
#include <iostream>
#include <string>
#include <vector>

#include "mux_bypass_edges.hpp"
#include "newcode/ofec_single_runner.hpp"

// ======== 用户可改区域 ========
// 只需要改这里的常量/列表即可完成一次“单次调试运行”的配置

// 发射端参数
static constexpr const char* kLabel              = "debug_L6";    // 运行标签：日志名、输出文件名前缀都会带这个名字
static constexpr int         kChaseL_override    = 6;             // Chase L，-1 表示使用 Params 里的默认值
static constexpr int         kBitgenSeed         = 20260319;    // 比特生成随机种子，固定后可复现实验
static constexpr bool        kGenerateRandomBits = true;       // true=发送随机信息比特，false=发送全 0 比特

// 信道参数
static constexpr float       kEbN0_db                      = 3.06f;   // 信道 Eb/N0，单位 dB
static constexpr int         kChannelSeed                  = 3182026; // 信道噪声随机种子，固定后可复现实验
static constexpr unsigned    kBitsPerSymbol                = 1;       // 每个调制符号携带的比特数：1=BPSK，偶数=QAM

//早停参数
static constexpr bool        kEnableEarlyStop              = true;   // true=启用早停，false=完全关闭早停路径
static const std::vector<int> kEarlyStopEnableList         = {1,1,1,1,1,1};      // 按 tile 覆盖早停总开关：0=关，非 0=开；空表示所有 tile 沿用 kEnableEarlyStop
static constexpr int         kEarlyStopConditionMode       = 1;       // 早停条件编号：1=v1，2=v2
static const std::vector<int> kEarlyStopConditionModeList  = {};      // 按 tile 覆盖早停条件模式；空表示所有 tile 沿用 kEarlyStopConditionMode
static constexpr int         kEarlyStopActionMode          = 7;       // 早停动作：1=sign beta，2=residual only，3=硬解输出 ±hard_mag，4=mode1 预除 alpha，5=residual 预除 alpha 后加 sign beta，6=直接 ±beta，7=直接 ±beta 预除 alpha，8=(直接 ±beta-channel) 预除 alpha
static const std::vector<int> kEarlyStopActionModeList     = {};      // 按 tile 覆盖早停动作模式；空表示所有 tile 沿用 kEarlyStopActionMode
static constexpr int         kEarlyStopBindGroupSize       = 1;       // 条件1专用的组绑定大小：1=逐 row；4=每 4 个 row 都通过才整体 early-stop
static const std::vector<int> kEarlyStopBindGroupSizeList  = {};      // 按 tile 覆盖条件1绑定组大小；空表示所有 tile 沿用 kEarlyStopBindGroupSize

static constexpr bool        kEarlyStopCondV1RequireBch    = true;    // 条件1里是否要求 BCH syndrome 为 0
static constexpr bool        kEarlyStopCondV1RequireOverall = true;   // 条件1里是否要求 overall parity 一致

static constexpr float       kEarlyStopV2LlrAbsThreshold   = 26.0f;   // 条件2里把 bit 视为“不可靠”的 |LLR| 阈值
static constexpr int         kEarlyStopV2MaxUnreliableBits = 25;      // 条件2里允许的不可靠 bit 数上限
static constexpr bool        kEarlyStopCondV2IncludeOverall = true;   // 条件2统计不可靠 bit 时是否把 overall bit 算进去

// LLR 量化相关参数
static constexpr std::size_t kLlrBits        = 6;      // LLR 位宽：16=浮点，2~15=qfloat
static constexpr float       kQuantClipRatio = 0.5f;   // 动态裁剪比例，0 表示禁用自适应 clip

// 解码主参数
static constexpr const char* kInterleaverName          = "identity"; // 交织器名称，identity 表示不改变比特顺序
static constexpr const char* kDecoderName              = "chase_baseline"; 
// 解码器名称：chase_baseline=逐 bit 搜索 Cplus/Cminus 的基线 Chase；
//chase_topk_pruned=Top-K 裁剪版（pruned=裁剪，只保留前 K 个 good 候选参与后续外信息计算）
//；chase_global_pair=全局固定一对 best/second 来计算各 bit 可靠度；
//chase_group_minima=分组组内最优版（minima=各组里度量最小/score 最大的代表候选）；
//chase_overall_parity_search=把 overall parity 也纳入搜索的变体
static constexpr int         kChaseNTestOverride       = -1;               // Chase 测试序列数量，-1 表示默认按 2^L 生成
static constexpr int         kChaseTopkKeep            = 24;                // chase_topk_pruned 中保留参与 ML/外信息计算的 Top-K 候选数
static constexpr int         kChaseGroupMinimaBits     = 3;                // chase_group_minima 里按前多少个 test-pattern 位分组；例如取 3 时会分成 2^3=8 组，并在每组里选一个组内最优 good 候选
static constexpr bool        kNormalizeExtrinsic       = false;      // 是否对 Chase 输出的 extrinsic 做归一化
static constexpr bool        kNormalizeKnownPrefixTail = false;      // 是否对 known-prefix 之后的尾部 LLR 做归一化

static constexpr float kEarlyStopActionResidualDivisor = 1.0f;   // 动作2里 residual 的除数
static constexpr float kEarlyStopActionHardLlrMag      = 1.0f;   // 动作3里硬解成功后输出的固定 |LLR| 幅度

// 显式列表（若非空，将直接作为各个 tile 的参数；长度必须等于 TILES_PER_WIN）
static const std::vector<float> kAlpha_explicit = {      // 每个 tile 的 alpha 显式列表
   0.428571,0.447738,0.482782,0.528162,0.581902,0.642857
};
static const std::vector<float> kBeta_explicit = {       // 每个 tile 的 Chase/fallback beta 显式列表
  2.857143,6.179301,12.253626,20.119585,29.434408,40.000000
};
static const std::vector<float> kEarlyStopActionBeta_explicit = { // 每个 tile 的 early-stop 动作 beta 显式列表
  99.857143,99.179301,99.253626,99.119585,99,99
};
static const std::vector<int> kSisoActiveList = {32, 32, 32, 32,32,32}; // 每个 tile 允许参与 SISO 的行数预算
static const std::vector<int> kHiHoActiveList = {32, 32, 32, 32, 32, 32}; // 每个 tile 允许参与 HIHO 硬解码的行数预算
static constexpr int  kMuxGroupG          = 1;                     // MUX 分组粒度，1 表示全局池化
static constexpr int  kMuxSchedulingMode  = 0;                     // MUX 调度模式：0=legacy，1=按 early-stop 细节排序
static constexpr int  kMuxPriorityRule    = 0;                     // 新 MUX 的优先级规则：0=更差优先，1=更接近通过优先
static constexpr bool kMuxEnableReconfig  = false;                 // true 表示启用重配置版 MUX 调度
static constexpr int  kMuxBypassScheme    = 1;                     // 旁路边集合方案编号：1=scheme1，2=scheme2



static constexpr bool kHybridEnable       = true;                 // true=方案三软硬混合前置分流开关
static const std::vector<int> kHybridEnableList = {0, 0, 0, 0, 1, 1}; // 按 tile 覆盖 hybrid 开关：空=沿用 kHybridEnable
static constexpr float kHybridHardLlrMag = 0.0f;                // hybrid hard-finish 默认输出 |LLR| 幅度
static const std::vector<float> kHybridHardLlrMagList = {0,0,0,0,99,99};      // 按 tile 覆盖 hybrid hard-finish |LLR| 幅度；空=沿用 kHybridHardLlrMag
static constexpr newcode::HybridClassifierMode kHybridClassifierMode =
    newcode::HybridClassifierMode::FriendS1S3WithS0Classifier;            // LegacyHardDecode / RepoFastClassifier / FriendS1S3Classifier / FriendS1S3WithS0Classifier
static constexpr newcode::HybridSisoBackfillMode kHybridSisoBackfillMode =
    newcode::HybridSisoBackfillMode::ParityOneAndTwoErrorPriority;                    // Disabled / TwoErrorOnly / OneAndTwoErrorPriority / ParityOneAndTwoErrorPriority
static constexpr bool kHybridNormalizeSoftOnly = false;            // true=只归一化 soft rows，false=保持当前兼容行为
static constexpr bool kLevel56SharedEnable = true;                // true=第五/六级共享 HISO/SISO
static constexpr bool kLevel56TemporalLookaheadEnable = false;    // 旧 t=0/t=1/t=2 三时刻 lookahead；与新 FIFO 互斥
static constexpr int kLevel56TemporalGroupLoadThreshold = 16;     // temporal 分支阈值，默认等价于 K1+K2 < 16-X
static constexpr bool kLevel56BufferedFifoEnable = true;          // 新方案：完整 64-code batch FIFO 与跨时刻 pending
static constexpr std::size_t kLevel56BufferRows = 32;             // R_buf，单位为 block row
static constexpr int kLevel56SharedHisoActive = 8;                 // 第五/六级共享 HISO 容量
static constexpr int kLevel56SharedSisoActive = 8;                 // 第五/六级共享 SISO 容量
static constexpr newcode::Level56PriorityMode kLevel56PriorityMode =
    newcode::Level56PriorityMode::Level5First; // Level5First=同优先级时 Level 5 优先；Level6First=Level 6 优先
static constexpr newcode::Level56ScheduleMode kLevel56ScheduleMode =
    newcode::Level56ScheduleMode::Group4LoadSortedMultiround; // GlobalPriority=全局优先级；Group4LoadSortedMultiround=64 code 固定分为 16 组、按组负载排序并多轮调度
static constexpr newcode::Level56EarlyStopGroupUpdateMode
    kLevel56EarlyStopGroupUpdateMode =
        newcode::Level56EarlyStopGroupUpdateMode::AllGroups;
static constexpr bool kLevel56SingleLevelSelectEnable = false; // true=按 early-stop 命中数动态只解一级
static constexpr bool kLevel56UnselectedEarlyStopActionEnable = true; // true=未选中级仍执行 early-stop action
static constexpr bool kDumpLevel56ScheduleStats = true; // true=导出每次共享调用的轮次与 code 级调度统计
static constexpr const char* kLevel56ScheduleRoundsPath =
    "data/level56_schedule/ofec_single_level56_schedule_rounds.csv";
static constexpr const char* kLevel56ScheduleCodesPath =
    "data/level56_schedule/ofec_single_level56_schedule_codes.csv";
// 5.2 等价性验证：只读逐 code 输入/输出哈希，默认关闭。
static constexpr bool kDumpLevel56EquivalenceObservation = false;
static constexpr const char* kLevel56EquivalenceObservationPath =
    "data/level56_observation/ofec_single_level56_equivalence_observation.csv";

// LLR 导出相关
static constexpr bool        kDumpQuantizedLlr = false;                         // 是否导出量化后的信道 LLR
static constexpr const char* kQuantizedLlrPath = "data/llr/quantized_llr.txt"; // 量化 LLR 输出路径
static constexpr bool        kDumpWorkLlr      = false;                         // 是否保存最终累积得到的 work_llr
static constexpr const char* kWorkLlrPath      = "data/llr/work_llr.txt";      // work_llr 输出路径
static constexpr bool        kDumpTileEarlyStopSamples = true;                 // 是否导出每次进入各个 tile 后的 early-stop 命中码字数样本
static constexpr const char* kTileEarlyStopSamplesPath =
    "data/early_stop_hist/ofec_single_tile_early_stop_samples.csv";             // early-stop 样本 CSV 路径
static constexpr bool        kDumpTileEarlyStopGroupBindDebugSamples = true;   // 是否导出 group-bind 前后的 early-stop bitstring
static constexpr const char* kTileEarlyStopGroupBindDebugSamplesPath =
    "data/early_stop_debug/ofec_single_group_bind_debug.csv";                  // group-bind debug CSV 路径

namespace {
// rx_info_from_bit_llr() 会先跳过全部已知前缀；本配置的前缀为
// 6 tiles × 22 block rows/tile × 16 bit rows/block row = 2112 行。
// 因而 Post-FEC error position 不是从全局矩阵第 0 行开始，而是从
// global row 2112 开始编号。此函数仅供关闭状态下的调试 target 映射使用。
constexpr long BitIndexToRow(long bit_index) { return bit_index / 111 + 2112; }
constexpr long BitIndexToCol(long bit_index) { return bit_index % 111; }        // 把 info bit 编号映射回矩阵列号
struct TraceBitSpec {
  long bit_index;
  const char* label;
};
constexpr TraceBitSpec kTraceBitSpecs[] = { // 需要重点跟踪的目标比特列表
    {677659, "bit677659"},
    {691376, "bit691376"},
    {694548, "bit694548"},
    {706209, "bit706209"},
    {706542, "bit706542"}
};
} 

static const std::vector<newcode::Params::DebugTraceConfig::TraceTarget>
    kDecoderTraceTargets = []() {
      std::vector<newcode::Params::DebugTraceConfig::TraceTarget> targets;
      for (const auto& spec : kTraceBitSpecs) {
        newcode::Params::DebugTraceConfig::TraceTarget entry;
        entry.row = BitIndexToRow(spec.bit_index);
        entry.col = BitIndexToCol(spec.bit_index);
        entry.bit_index = spec.bit_index;
        if (spec.label) entry.label = spec.label;
        targets.push_back(entry);
      }
      return targets;
    }();

// Decoder 调试跟踪配置
// 用于定位单个比特的映射/写回问题：row/col 为全局 work_llr 坐标（0-based）
// enable 与 row/col 同时满足才会输出日志
// log_read_mapping   输出 tile -> global 的读取坐标展开
// log_write_mapping  输出解码回写对应的坐标及 LLR
// log_mismatch       当窗口内多个 tile 写回同一坐标且值不同时报错提示
// log_chase_detail   追踪对应比特在 Chase 内部（chase_baseline / chase_overall_parity_search）的输入/输出
// dump_chase_csv     每次进入 Chase 时导出 256 码字的 LLR/ω/ML/硬判决到 CSV
// chase_csv_dir      CSV 导出目录（可直接用 Excel 打开）

static constexpr bool        kDecoderTraceEnable      = false;            // 是否启用 decoder trace 总开关
static constexpr bool        kDecoderTraceLogRead     = false;            // 是否打印 tile 读取映射
static constexpr bool        kDecoderTraceLogWrite    = false;            // 是否打印 tile 写回映射
static constexpr bool        kDecoderTraceLogMismatch = false;            // 是否打印同坐标写回不一致告警
static constexpr bool        kDecoderTraceLogChase    = false;            // 是否打印 Chase 内部细节
static constexpr bool        kDecoderTraceDumpCsv     = false;            // 是否把 Chase 细节落成 CSV
static constexpr const char* kDecoderTraceCsvDir      = "data/chase_csv"; // Chase CSV 导出目录


// =============================

int main() {
  bool level56_buffered_fifo_enable = kLevel56BufferedFifoEnable;
  std::size_t level56_buffer_rows = kLevel56BufferRows;
  std::size_t level5_buffer_rows =
      std::numeric_limits<std::size_t>::max();
  std::size_t level6_buffer_rows =
      std::numeric_limits<std::size_t>::max();
  std::size_t num_info_bits = 0;
  bool level56_buffered_fifo_drain_at_frame_end = false;
  unsigned level56_hiso_allowed_class_mask = 0x0fu;
  int bitgen_seed = kBitgenSeed;
  int channel_seed = kChannelSeed;
  int level56_shared_hiso_active = kLevel56SharedHisoActive;
  int level56_shared_siso_active = kLevel56SharedSisoActive;
  std::size_t level56_group4_max_entries = 8;
  auto level56_schedule_mode = kLevel56ScheduleMode;
  float ebn0_db = kEbN0_db;
  std::vector<float> hybrid_hard_llr_mag_list = kHybridHardLlrMagList;
  std::string run_label = kLabel;
  std::string level56_schedule_rounds_path = kLevel56ScheduleRoundsPath;
  std::string level56_schedule_codes_path = kLevel56ScheduleCodesPath;
  bool dump_level56_equivalence_observation =
      kDumpLevel56EquivalenceObservation;
  std::string level56_equivalence_observation_path =
      kLevel56EquivalenceObservationPath;
  bool level56_target_trace_enable = false;
  std::size_t level56_target_trace_batch =
      std::numeric_limits<std::size_t>::max();
  int level56_target_trace_level = -1;
  int level56_target_trace_code = -1;
  std::string level56_target_trace_path;
  std::string tile_early_stop_samples_path = kTileEarlyStopSamplesPath;
  std::string tile_early_stop_group_bind_debug_samples_path =
      kTileEarlyStopGroupBindDebugSamplesPath;
  std::string run_suffix;

  if (const char* value = std::getenv("OFEC_EBN0_DB")) {
    const std::string text(value);
    char* end = nullptr;
    const float parsed = std::strtof(text.c_str(), &end);
    if (text.empty() || end != text.c_str() + text.size() ||
        !std::isfinite(parsed)) {
      std::cerr << "[ERROR] OFEC_EBN0_DB must be a finite number\n";
      return 2;
    }
    ebn0_db = parsed;
    run_suffix += "_ebn0" + text;
  }

  if (const char* value = std::getenv("OFEC_NUM_INFO_BITS")) {
    const std::string text(value);
    const auto parse_result = std::from_chars(
        text.data(), text.data() + text.size(), num_info_bits);
    if (text.empty() || parse_result.ec != std::errc{} ||
        parse_result.ptr != text.data() + text.size() || num_info_bits == 0) {
      std::cerr << "[ERROR] OFEC_NUM_INFO_BITS must be a positive integer\n";
      return 2;
    }
    run_suffix += "_nbits" + std::to_string(num_info_bits);
  }

  // 实验覆盖：逗号分隔、长度必须为 6 的 tile 级 HISO hard-finish 幅度。
  // 例如 OFEC_HYBRID_HARD_LLR_MAG_LIST=0,0,0,0,30,30。
  // 未设置时严格保持 kHybridHardLlrMagList 的既有默认配置。
  if (const char* value = std::getenv("OFEC_HYBRID_HARD_LLR_MAG_LIST")) {
    const std::string text(value);
    std::vector<float> parsed;
    std::size_t begin = 0;
    while (begin <= text.size()) {
      const std::size_t end = text.find(',', begin);
      const std::string token = text.substr(
          begin, end == std::string::npos ? std::string::npos : end - begin);
      char* parse_end = nullptr;
      const float magnitude = std::strtof(token.c_str(), &parse_end);
      if (token.empty() || parse_end != token.c_str() + token.size() ||
          !std::isfinite(magnitude)) {
        std::cerr << "[ERROR] OFEC_HYBRID_HARD_LLR_MAG_LIST must be six "
                     "finite comma-separated floats\n";
        return 2;
      }
      parsed.push_back(magnitude);
      if (end == std::string::npos) break;
      begin = end + 1u;
    }
    if (parsed.size() != kHybridHardLlrMagList.size()) {
      std::cerr << "[ERROR] OFEC_HYBRID_HARD_LLR_MAG_LIST must contain exactly "
                << kHybridHardLlrMagList.size() << " values\n";
      return 2;
    }
    hybrid_hard_llr_mag_list = std::move(parsed);
    run_suffix += "_hybridhard" + text;
  }

  const auto parse_seed = [&](const char* env_name,
                              int* destination,
                              const char* suffix_name) -> bool {
    const char* value = std::getenv(env_name);
    if (!value) {
      return true;
    }
    const std::string text(value);
    const auto parse_result = std::from_chars(
        text.data(), text.data() + text.size(), *destination);
    if (text.empty() || parse_result.ec != std::errc{} ||
        parse_result.ptr != text.data() + text.size()) {
      std::cerr << "[ERROR] " << env_name << " must be an integer\n";
      return false;
    }
    run_suffix += std::string("_") + suffix_name +
                  std::to_string(*destination);
    return true;
  };
  if (!parse_seed("OFEC_BITGEN_SEED", &bitgen_seed, "bitseed") ||
      !parse_seed("OFEC_CHANNEL_SEED", &channel_seed, "chseed")) {
    return 2;
  }

  if (const char* value = std::getenv("LEVEL56_BUFFERED_FIFO_ENABLE")) {
    const std::string text(value);
    if (text == "0") {
      level56_buffered_fifo_enable = false;
    } else if (text == "1") {
      level56_buffered_fifo_enable = true;
    } else {
      std::cerr << "[ERROR] LEVEL56_BUFFERED_FIFO_ENABLE must be 0 or 1\n";
      return 2;
    }
    run_suffix += level56_buffered_fifo_enable ? "_fifo1" : "_fifo0";
  }

  if (const char* value = std::getenv("LEVEL56_BUFFER_ROWS")) {
    const std::string text(value);
    const auto parse_result = std::from_chars(
        text.data(), text.data() + text.size(), level56_buffer_rows);
    if (text.empty() || parse_result.ec != std::errc{} ||
        parse_result.ptr != text.data() + text.size()) {
      std::cerr << "[ERROR] LEVEL56_BUFFER_ROWS must be a non-negative integer\n";
      return 2;
    }

    run_suffix += "_rbuf" + std::to_string(level56_buffer_rows);
  }

  const auto parse_split_buffer_rows = [&](const char* name,
                                           std::size_t* target,
                                           const char* suffix) {
    const char* value = std::getenv(name);
    if (!value) return true;
    const std::string text(value);
    const auto parse_result = std::from_chars(
        text.data(), text.data() + text.size(), *target);
    if (text.empty() || parse_result.ec != std::errc{} ||
        parse_result.ptr != text.data() + text.size()) {
      std::cerr << "[ERROR] " << name
                << " must be a non-negative integer\n";
      return false;
    }
    run_suffix += std::string(suffix) + std::to_string(*target);
    return true;
  };
  if (!parse_split_buffer_rows("LEVEL5_BUFFER_ROWS", &level5_buffer_rows,
                               "_rbuf5_") ||
      !parse_split_buffer_rows("LEVEL6_BUFFER_ROWS", &level6_buffer_rows,
                               "_rbuf6_")) {
    return 2;
  }

  if (const char* value = std::getenv("LEVEL56_BUFFERED_FIFO_DRAIN_AT_FRAME_END")) {
    const std::string text(value);
    if (text == "0") {
      level56_buffered_fifo_drain_at_frame_end = false;
    } else if (text == "1") {
      level56_buffered_fifo_drain_at_frame_end = true;
    } else {
      std::cerr << "[ERROR] LEVEL56_BUFFERED_FIFO_DRAIN_AT_FRAME_END must be 0 or 1\n";
      return 2;
    }
    run_suffix += level56_buffered_fifo_drain_at_frame_end ? "_drain1" : "_drain0";
  }

  if (const char* value = std::getenv("LEVEL56_HISO_CLASS_MASK")) {
    const std::string text(value);
    const auto parse_result = std::from_chars(
        text.data(), text.data() + text.size(), level56_hiso_allowed_class_mask,
        10);
    if (text.empty() || parse_result.ec != std::errc{} ||
        parse_result.ptr != text.data() + text.size() ||
        level56_hiso_allowed_class_mask > 0x0fu) {
      std::cerr << "[ERROR] LEVEL56_HISO_CLASS_MASK must be an integer in [0,15] "
                   "(bits: ParityOnly/OneMain/OneMainPlusParity/TwoMain)\\n";
      return 2;
    }
    run_suffix += "_hisomask" + std::to_string(level56_hiso_allowed_class_mask);
  }

  const auto parse_level56_capacity = [&](const char* env_name,
                                          int* destination,
                                          const char* suffix_name) -> bool {
    const char* value = std::getenv(env_name);
    if (!value) {
      return true;
    }
    const std::string text(value);
    const auto parse_result = std::from_chars(
        text.data(), text.data() + text.size(), *destination);
    if (text.empty() || parse_result.ec != std::errc{} ||
        parse_result.ptr != text.data() + text.size() ||
        *destination < 0 || *destination > 64) {
      std::cerr << "[ERROR] " << env_name
                << " must be an integer in [0,64]\n";
      return false;
    }
    run_suffix += std::string("_") + suffix_name +
                  std::to_string(*destination);
    return true;
  };
  if (!parse_level56_capacity("LEVEL56_SHARED_HISO_ACTIVE",
                              &level56_shared_hiso_active, "hiso") ||
      !parse_level56_capacity("LEVEL56_SHARED_SISO_ACTIVE",
                              &level56_shared_siso_active, "siso")) {
    return 2;
  }

  if (const char* value = std::getenv("LEVEL56_GROUP4_MAX_ENTRIES")) {
    const std::string text(value);
    const auto parse_result = std::from_chars(
        text.data(), text.data() + text.size(), level56_group4_max_entries);
    if (text.empty() || parse_result.ec != std::errc{} ||
        parse_result.ptr != text.data() + text.size() ||
        level56_group4_max_entries == 0 || level56_group4_max_entries > 64) {
      std::cerr << "[ERROR] LEVEL56_GROUP4_MAX_ENTRIES must be an integer in [1,64]\n";
      return 2;
    }
    run_suffix += "_g4entries" + std::to_string(level56_group4_max_entries);
  }

  if (const char* value = std::getenv("LEVEL56_SCHEDULE_MODE")) {
    const std::string text(value);
    if (text == "global") {
      level56_schedule_mode = newcode::Level56ScheduleMode::GlobalPriority;
    } else if (text == "group4") {
      level56_schedule_mode =
          newcode::Level56ScheduleMode::Group4LoadSortedMultiround;
    } else {
      std::cerr << "[ERROR] LEVEL56_SCHEDULE_MODE must be global or group4\n";
      return 2;
    }
    run_suffix += "_sched" + text;
  }

  if (const char* value = std::getenv("LEVEL56_EQUIVALENCE_OBSERVATION")) {
    const std::string text(value);
    if (text == "0") {
      dump_level56_equivalence_observation = false;
    } else if (text == "1") {
      dump_level56_equivalence_observation = true;
    } else {
      std::cerr << "[ERROR] LEVEL56_EQUIVALENCE_OBSERVATION must be 0 or 1\n";
      return 2;
    }
    if (dump_level56_equivalence_observation) {
      run_suffix += "_eqobs1";
    }
  }

  if (const char* value = std::getenv("LEVEL56_TARGET_TRACE")) {
    const std::string text(value);
    if (text == "0") {
      level56_target_trace_enable = false;
    } else if (text == "1") {
      level56_target_trace_enable = true;
    } else {
      std::cerr << "[ERROR] LEVEL56_TARGET_TRACE must be 0 or 1, got ["
                << text << "]\n";
      return 2;
    }
  }
  const auto parse_target_trace_number = [&](const char* env_name,
                                             auto* destination) -> bool {
    const char* value = std::getenv(env_name);
    if (!value) return true;
    const std::string text(value);
    const auto parsed = std::from_chars(text.data(), text.data() + text.size(),
                                        *destination);
    if (text.empty() || parsed.ec != std::errc{} ||
        parsed.ptr != text.data() + text.size()) {
      std::cerr << "[ERROR] " << env_name << " must be an integer, got ["
                << text << "]\n";
      return false;
    }
    return true;
  };
  if (!parse_target_trace_number("LEVEL56_TRACE_BATCH",
                                 &level56_target_trace_batch) ||
      !parse_target_trace_number("LEVEL56_TRACE_LEVEL",
                                 &level56_target_trace_level) ||
      !parse_target_trace_number("LEVEL56_TRACE_CODE",
                                 &level56_target_trace_code)) {
    return 2;
  }
  if (level56_target_trace_enable) {
    if (level56_target_trace_batch ==
            std::numeric_limits<std::size_t>::max() ||
        (level56_target_trace_level != 5 && level56_target_trace_level != 6) ||
        level56_target_trace_code < 1 || level56_target_trace_code > 64) {
      std::cerr << "[ERROR] LEVEL56_TARGET_TRACE=1 requires "
                   "LEVEL56_TRACE_BATCH, LEVEL56_TRACE_LEVEL=5|6, and "
                   "LEVEL56_TRACE_CODE=1..64\n";
      return 2;
    }
    run_suffix += "_traceB" + std::to_string(level56_target_trace_batch) +
                  "L" + std::to_string(level56_target_trace_level) +
                  "C" + std::to_string(level56_target_trace_code);
  }

  if (!run_suffix.empty()) {
    run_label += run_suffix;
    level56_schedule_rounds_path =
        "data/level56_schedule/ofec_single" + run_suffix +
        "_level56_schedule_rounds.csv";
    level56_schedule_codes_path =
        "data/level56_schedule/ofec_single" + run_suffix +
        "_level56_schedule_codes.csv";
    level56_equivalence_observation_path =
        "data/level56_observation/ofec_single" + run_suffix +
        "_level56_equivalence_observation.csv";
    if (level56_target_trace_enable) {
      level56_target_trace_path =
          "data/level56_target_trace/ofec_single" + run_suffix +
          "_level56_target_trace.txt";
    }
    tile_early_stop_samples_path =
        "data/early_stop_hist/ofec_single" + run_suffix +
        "_tile_early_stop_samples.csv";
    tile_early_stop_group_bind_debug_samples_path =
        "data/early_stop_debug/ofec_single" + run_suffix +
        "_group_bind_debug.csv";
  }

  const auto& selected_mux_bypass_edges =
      app_mux::bypass_edges_for_scheme(kMuxBypassScheme);
  ofec_single::Config config{
    .label = run_label,
    .ebn0_db = ebn0_db,
    .chaseL_override = kChaseL_override,
    .chase_n_test_override = kChaseNTestOverride,
    .chase_topk_keep = kChaseTopkKeep,
    .chase_group_minima_bits = kChaseGroupMinimaBits,
    .normalize_extrinsic = kNormalizeExtrinsic,
    .bits_per_symbol = kBitsPerSymbol,
    .bitgen_seed = bitgen_seed,
    .num_info_bits = num_info_bits,
    .channel_seed = channel_seed,
    .enable_early_stop = kEnableEarlyStop,
    .early_stop_enable_list = kEarlyStopEnableList,
    .early_stop_condition_mode = kEarlyStopConditionMode,
    .early_stop_condition_mode_list = kEarlyStopConditionModeList,
    .early_stop_action_mode = kEarlyStopActionMode,
    .early_stop_action_mode_list = kEarlyStopActionModeList,
    .early_stop_bind_group_size = kEarlyStopBindGroupSize,
    .early_stop_bind_group_size_list = kEarlyStopBindGroupSizeList,
    .early_stop_cond_v1_require_bch = kEarlyStopCondV1RequireBch,
    .early_stop_cond_v1_require_overall = kEarlyStopCondV1RequireOverall,
    .early_stop_v2_llr_abs_threshold = kEarlyStopV2LlrAbsThreshold,
    .early_stop_v2_max_unreliable_bits = kEarlyStopV2MaxUnreliableBits,
    .early_stop_cond_v2_include_overall = kEarlyStopCondV2IncludeOverall,
    .alpha_explicit = kAlpha_explicit,
    .beta_explicit = kBeta_explicit,
    .early_stop_action_sign_beta_explicit = kEarlyStopActionBeta_explicit,
    .early_stop_action_residual_divisor = kEarlyStopActionResidualDivisor,
    .early_stop_action_hard_llr_mag = kEarlyStopActionHardLlrMag,
    .siso_active_list = kSisoActiveList,
    .hiho_active_list = kHiHoActiveList,
    .mux_group_g = kMuxGroupG,
    .mux_scheduling_mode = kMuxSchedulingMode,
    .mux_early_stop_priority_rule = kMuxPriorityRule,
    .mux_enable_reconfig = kMuxEnableReconfig,
    .mux_extra_bypass_edges = selected_mux_bypass_edges,
    .hybrid_enable = kHybridEnable,
    .hybrid_enable_list = kHybridEnableList,
    .hybrid_hard_llr_mag = kHybridHardLlrMag,
    .hybrid_hard_llr_mag_list = hybrid_hard_llr_mag_list,
    .hybrid_classifier_mode = kHybridClassifierMode,
    .hybrid_siso_backfill_mode = kHybridSisoBackfillMode,
    .hybrid_normalize_soft_only = kHybridNormalizeSoftOnly,
    .level56_shared_enable = kLevel56SharedEnable,
    .level56_temporal_lookahead_enable = kLevel56TemporalLookaheadEnable,
    .level56_temporal_group_load_threshold = kLevel56TemporalGroupLoadThreshold,
    .level56_buffered_fifo_enable = level56_buffered_fifo_enable,
    .level56_buffer_rows = level56_buffer_rows,
    .level5_buffer_rows = level5_buffer_rows,
    .level6_buffer_rows = level6_buffer_rows,
    .level56_buffered_fifo_drain_at_frame_end =
        level56_buffered_fifo_drain_at_frame_end,
    .level56_hiso_allowed_class_mask = level56_hiso_allowed_class_mask,
    .level56_shared_hiso_active = level56_shared_hiso_active,
    .level56_shared_siso_active = level56_shared_siso_active,
    .level56_group4_max_entries = level56_group4_max_entries,
    .level56_priority_mode = kLevel56PriorityMode,
    .level56_schedule_mode = level56_schedule_mode,
    .level56_early_stop_group_update_mode = kLevel56EarlyStopGroupUpdateMode,
    .level56_single_level_select_enable = kLevel56SingleLevelSelectEnable,
    .level56_unselected_early_stop_action_enable =
        kLevel56UnselectedEarlyStopActionEnable,
    .dump_level56_schedule_stats = kDumpLevel56ScheduleStats,
    .level56_schedule_rounds_output_path = level56_schedule_rounds_path,
    .level56_schedule_codes_output_path = level56_schedule_codes_path,
    .dump_level56_equivalence_observation =
        dump_level56_equivalence_observation,
    .level56_equivalence_observation_output_path =
        level56_equivalence_observation_path,
    .level56_target_trace_enable = level56_target_trace_enable,
    .level56_target_trace_batch = level56_target_trace_batch,
    .level56_target_trace_level = level56_target_trace_level,
    .level56_target_trace_code = level56_target_trace_code,
    .level56_target_trace_output_path = level56_target_trace_path,
    .interleaver_name = kInterleaverName,
    .decoder_name = kDecoderName,
    .generate_random_bits = kGenerateRandomBits,
    .normalize_known_prefix_tail = kNormalizeKnownPrefixTail,
    .quant_clip_ratio = kQuantClipRatio,
    .llr_bits = kLlrBits,
    .dump_quantized_llr = kDumpQuantizedLlr,
    .quantized_llr_output_path = kQuantizedLlrPath,
    .dump_work_llr = kDumpWorkLlr,
    .work_llr_output_path = kWorkLlrPath,
    .post_fec_error_positions_output_path =
        std::getenv("POST_FEC_ERROR_POSITIONS_PATH") != nullptr
            ? std::getenv("POST_FEC_ERROR_POSITIONS_PATH")
            : "",
    .dump_tile_early_stop_samples = kDumpTileEarlyStopSamples,
    .tile_early_stop_samples_output_path = tile_early_stop_samples_path,
    .dump_tile_early_stop_group_bind_debug_samples =
        kDumpTileEarlyStopGroupBindDebugSamples,
    .tile_early_stop_group_bind_debug_samples_output_path =
        tile_early_stop_group_bind_debug_samples_path,
    .debug_trace = newcode::Params::DebugTraceConfig{
      .enable = kDecoderTraceEnable,
      .log_read_mapping = kDecoderTraceLogRead,
      .log_write_mapping = kDecoderTraceLogWrite,
      .log_mismatch = kDecoderTraceLogMismatch,
      .log_chase_detail = kDecoderTraceLogChase,
      .dump_chase_csv = kDecoderTraceDumpCsv,
      .row = -1,
      .col = -1,
      .chase_decoder_row = -1,
      .chase_decoder_col = -1,
      .chase_tile_index = -1,
      .chase_invocation = -1,
      .chase_csv_dir = kDecoderTraceCsvDir,
      .chase_expected_bits = {},
      .chase_candidate_s1 = {},
      .chase_candidate_s3 = {},
      .chase_candidate_good = {},
      .chase_candidate_corrected_errors = {},
      .targets = kDecoderTraceTargets,
      .active_chase_entries = {},
    },
  };
  return ofec_single::run_ofec_single(config);
}
