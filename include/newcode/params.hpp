#pragma once
#include <cstddef>
#include <vector>
#include <string>
#include <memory>

#include "newcode/ofec/mux/mux_topology.hpp"

namespace newcode {

struct Params {
  // ===== 编码结构常量（与 oFEC 布局一致） =====
  static constexpr size_t NUM_SUBBLOCK_COLS     = 8;   // 每行包含的子块数量
  static constexpr size_t INFO_SUBROWS_PER_CODE = 16;  // BCH 信息段对应的子块行数
  static constexpr size_t BITS_PER_SUBBLOCK_DIM = 16;  // 子块维度（子块大小 = 16×16 bits）

  // ===== BCH 码参数 =====
  static constexpr size_t BCH_N = 256; // 码长（255 + overall parity）
  static constexpr size_t BCH_K = 239; // 信息长度（对应 Chien/Berlekamp 实现）
  static constexpr size_t BCH_PARITY_BITS = BCH_N - BCH_K - 1; // 本地 parity 数量（256-239-1 = 16）
  static constexpr size_t BCH_OVERALL_IDX = BCH_N - 1;         // overall parity 索引（255）

  // ===== 运行/仿真参数 =====
  size_t NUM_INFO_BITS     =  32 * 110 * 16 * 111; // 信息比特总数
  int    BITGEN_SEED       = 56456;                 // 随机种子
  int    CHANNEL_SEED      = BITGEN_SEED + 656;    // 信道噪声随机种子
  bool   BITGEN_RANDOM_BITS = true;             // true=随机比特，false=全 0
  bool   NORMALIZE_KNOWN_PREFIX_TAIL = false;    // true=对非 known_rows 区域做均值归一化
  size_t NUM_GUARD_SUBROWS = 2;                  // 保护块子行数 G
  bool   DUMP_WORK_LLR = false;                  // 是否保存窗口累积后的 work_llr（解码前）
  std::string WORK_LLR_OUTPUT_PATH;              // work_llr 输出路径（为空则默认命名）

  // ===== 解码组织参数（单位：sub-block rows）=====
  size_t TILES_PER_WIN   = 4;  // 每个 window 含有的 tile 数量（自下而上处理）
  size_t TILE_OVERLAP_BR = 0;  // 相邻 tile 在 sub-block-row 方向的重叠行数
  size_t TILE_HEIGHT_BR  = 22; // 单个 tile 的高度（单位：sub-block-row）
  size_t WINDOW_POP_PUSH = 2;  // window 每次滑动的 sub-block-row 数量（pop/push）

  // ===== LLR 量化参数 =====
  size_t LLR_BITS = 16;   // 量化位宽（2~15 使用 qfloat::qfloat<N>，16=浮点）
  float  LLR_CLIP = 8.0f; // LLR 裁剪幅度（外部计算/配置后写入 qfloat clip）
  float  LLR_CLIP_RATIO = 0.0f; // 若 >0，则按比例动态计算 clip

  // ===== Chase-Pyndiah 控制参数 =====
  int CHASE_L     = 6;  // 选取“最不可靠”位置的数量
  int CHASE_NTEST = 64; // 生成的测试向量数量（<= 2^CHASE_L）
  int CHASE_SBR   = 2;  // 每个 tile 底部解码的子块行数（1 或 2）
  bool ENABLE_EARLY_STOP = true;  // 是否启用 tile/row 级早停判定

  // —— fallback 可靠度系数（Chase(256) 中 L0 的回退幅度）——
  float beta = 0.35f;

  // —— extrinsic scaling 系数 α（eq. (21) 调整外信息幅度）——
  float ALPHA = 1.0f;

  // —— 可按 tile 覆盖的系数表（索引 0..TILES_PER_WIN-1）——
  std::vector<float> ALPHA_LIST = {0.3f, 0.4f, 0.5f, 0.6f};    // 由起点+步进生成的默认列表
  std::vector<float> beta_list  = {0.9f, 1.0f, 1.1f, 1.2f};    // 由起点+步进生成的默认列表
  // 每级 tile 的可用 SISO 数（索引 t 按 bottom->top）
  std::vector<int> SISO_ACTIVE_LIST = {32, 30, 24, 16};
  // MUX 分组数：1=现有 max 全局池化，>1=按组预算裁剪
  int MUX_GROUP_G = 1;
  // 是否启用 scheme B 风格的两阶段重排调度
  bool MUX_ENABLE_RECONFIG = false;
  // 额外允许的跨组旁路线
  std::vector<mux::MuxEdge> MUX_EXTRA_BYPASS_EDGES;

  // —— 每个 tile 是否切换到硬判决译码 —— //
  bool HARD_DECODE_DEFAULT = false;                               // 默认仍使用软判决
  std::vector<int> HARD_TILE_LIST = {0, 0, 0, 0, 1};               // 0=软判决，非 0=硬判决
  float HARD_LLR_MAG = 1.0f;                                      // 硬判决映射的 |LLR| 大小

  struct DebugTraceConfig {
    bool enable = false;
    bool log_read_mapping = false;
    bool log_write_mapping = false;
    bool log_mismatch = false;
    bool log_chase_detail = false;       // 是否在 Chase 内部输出单比特信息
    bool dump_chase_csv = false;         // 是否将 Chase 输入/输出向量写 CSV
    long row = -1;
    long col = -1;
    int chase_decoder_row = -1;          // 在 Chase 输入矩阵中的行（0-based）
    int chase_decoder_col = -1;          // 在 256 码字中的列（0-based）
    int chase_tile_index = -1;           // 当前 tile 索引
    int chase_invocation = -1;           // 本次进入 Chase 的序号
    std::string chase_csv_dir;           // CSV 输出目录（为空则使用默认）
    std::shared_ptr<std::vector<std::vector<int8_t>>> chase_expected_bits;
    const std::vector<int8_t>* chase_expected_bits_row = nullptr;
    struct TraceTarget {
      long row = -1;
      long col = -1;
      long bit_index = -1;
      std::string label;
    };
    struct ChaseTraceEntry {
      int row_index = -1;
      int k = -1;
      long global_row = -1;
      long global_col = -1;
      long bit_index = -1;
      std::string label;
      int expected_bit = -1;
      std::vector<uint8_t> cplus_bits;
      std::vector<uint8_t> cminus_bits;
    };
    std::vector<TraceTarget> targets;
    std::vector<ChaseTraceEntry> active_chase_entries;
  };

  DebugTraceConfig debug_trace{};

  // ===== 便捷派生（统一换算为“比特行 rows”）=====
  constexpr size_t tile_height_rows() const {
    return TILE_HEIGHT_BR * BITS_PER_SUBBLOCK_DIM;
  }
  constexpr size_t tile_stride_rows() const {
    return (TILE_HEIGHT_BR - TILE_OVERLAP_BR) * BITS_PER_SUBBLOCK_DIM;
  }
  constexpr size_t win_height_br() const {
    return TILES_PER_WIN * TILE_HEIGHT_BR - (TILES_PER_WIN - 1) * TILE_OVERLAP_BR;
  }
  constexpr size_t win_height_rows() const {
    return win_height_br() * BITS_PER_SUBBLOCK_DIM;
  }
  constexpr size_t pop_push_rows() const {
    return WINDOW_POP_PUSH * BITS_PER_SUBBLOCK_DIM;
  }
  // 初始 window 起始行（比特行）：跳过 warmup 的 B 行 + 2G 行保护区
  constexpr size_t initial_win_start_rows() const {
    return 0;
  }

  // 基本有效性检查
  constexpr bool valid() const {
    return (TILES_PER_WIN >= 1) &&
           (TILE_OVERLAP_BR < TILE_HEIGHT_BR) &&
           (CHASE_L >= 1) &&
           (CHASE_NTEST >= 1) &&
           (CHASE_SBR == 1 || CHASE_SBR == 2);
  }
};

} // namespace newcode
