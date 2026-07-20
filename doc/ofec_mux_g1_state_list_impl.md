# oFEC MUX（G=1）状态列表实现方案（文档规格版）

## 1. 目标与范围

本说明用于指导后续代码实现，当前只定义规格，不改代码。

目标是基于现有早停判断，构建 `state[i]` 列表来模拟 `SISO` 数量受限时的调度行为，并满足两点：

1. 在 `apps` 入口可配置每一级 tile 的 SISO 数（示例 `[32,30,24,16]`）。
2. `state[i]` 初始化直接复用 `TileEarlyStopResult early_stop_stats`，不重复计算早停。

本阶段范围：

- 仅 `G=1`。
- 覆盖 `apps/ofec_single.cpp` 与 `apps/ofec_sweep.cpp`。
- `state=2` 仅定义为“本轮未分配到 SISO”，不在本规格内确定执行策略。

## 2. 状态定义

对每个 code（或 row）定义状态：

- `state[i] = 1`：该项满足早停条件，不需要 SISO。
- `state[i] = 0`：该项需要 SISO 解码。
- `state[i] = 2`：该项本轮未分配到 SISO。

列表长度等于当前 tile 的 code 数（示例为 32）。

## 3. 与现有结构对齐（直接复用 `TileEarlyStopResult`）

`state` 的初始化直接来自 `TileEarlyStopResult early_stop_stats`：

- `early_stop_stats.row_passed_flags[i] == true` -> `state[i] = 1`
- `early_stop_stats.row_passed_flags[i] == false` -> 初始 `state[i] = 0`

可直接复用的统计字段：

- `early_stop_stats.rows_passed`
- `early_stop_stats.rows_total`
- `early_stop_stats.all_rows_passed`

不再单独维护一份重复的 early-stop flags。

## 4. 每级 tile 的 SISO 列表配置（apps 入口）

建议新增参数名：`SISO_ACTIVE_LIST`，示例：

```cpp
SISO_ACTIVE_LIST = {32, 30, 24, 16};
```

配置约束：

1. 列表长度必须等于 `TILES_PER_WIN`。
2. 索引语义采用自底向上（bottom->top）tile 顺序。
3. `tile t` 使用 `SISO_ACTIVE_LIST[t]` 作为本级可用 SISO 容量。

配置入口范围：

- `apps/ofec_single.cpp`
- `apps/ofec_sweep.cpp`

## 5. G=1 调度规则

`G=1` 表示全可达，无组约束，仅有容量约束。

对每个 tile 的处理流程：

1. 根据 `early_stop_stats.row_passed_flags` 初始化 `state`（1/0）。
2. 统计 `state==0` 的索引集合 `need_indices`。
3. 读取本 tile 容量 `siso_active_for_tile = SISO_ACTIVE_LIST[t]`。
4. 若 `need_indices.size() <= siso_active_for_tile`，所有 `0` 保持不变。
5. 若 `need_indices.size() > siso_active_for_tile`，将超出容量的 `0` 改为 `2`。

裁剪顺序约定：

- 本阶段固定为索引顺序（低索引优先保留为 `0`）。

## 6. 伪代码（输入改为 `early_stop_stats + per-tile siso`）

```cpp
// 输入：
//   early_stop_stats: 来自 tile_early_stop_stats1(...)
//   siso_active_for_tile: 当前 tile 的可用 SISO 数
// 输出：
//   state[i] in {0,1,2}
std::vector<uint8_t> build_state_g1(const TileEarlyStopResult& early_stop_stats,
                                    int siso_active_for_tile) {
  const int n = static_cast<int>(early_stop_stats.row_passed_flags.size());
  std::vector<uint8_t> state(static_cast<size_t>(n), 0);
  std::vector<int> need_indices;
  need_indices.reserve(static_cast<size_t>(n));

  for (int i = 0; i < n; ++i) {
    if (early_stop_stats.row_passed_flags[static_cast<size_t>(i)]) {
      state[static_cast<size_t>(i)] = 1;
    } else {
      state[static_cast<size_t>(i)] = 0;
      need_indices.push_back(i);
    }
  }

  if (static_cast<int>(need_indices.size()) > siso_active_for_tile) {
    for (int k = siso_active_for_tile; k < static_cast<int>(need_indices.size()); ++k) {
      const int idx = need_indices[static_cast<size_t>(k)];
      state[static_cast<size_t>(idx)] = 2;
    }
  }
  return state;
}
```

## 7. `state=2` 的语义约束（本阶段）

`state=2` 在本规格中的定义仅为：

- 该项本轮未分配到 SISO（超出本级容量）。

本阶段不定义其执行策略（例如 `drop` 或 `fallback`），留待后续版本单独决策。

## 8. 工程落点（后续代码实现时参考）

为实现本规格，后续代码改造建议落在以下位置：

1. `apps/ofec_single.cpp`
   - 增加每级 tile 的 SISO 列表配置项。
2. `apps/ofec_sweep.cpp`
   - 增加 sweep 配置中的 SISO 列表透传。
3. `src/ofec_single/ofec_single_params.cpp`
   - 把 `Config` 中的 SISO 列表透传到 `Params`，并做长度校验。
4. `src/rx/ofec/detail/ofec_window_impl.ipp`
   - 按 tile 索引 `t` 读取本级 `siso_active_for_tile`。
5. `src/rx/ofec/detail/ofec_tile_impl.ipp`
   - 基于 `early_stop_stats` 直接构建 `state` 列表。

拟新增接口（规格层）：

1. `Params` 新增：
   - `std::vector<int> SISO_ACTIVE_LIST`
2. `ofec_single::Config` 与 sweep 配置新增：
   - `std::vector<int> siso_active_list`
3. 索引约定：
   - `SISO_ACTIVE_LIST[t]` 的 `t` 为 bottom->top。

## 8.1 新增函数独立文件化（本次新增约束）

在本规格基础上，新增函数不直接堆在现有 `ofec_tile_impl.ipp` / `ofec_window_impl.ipp` 内部，统一拆到独立模块文件：

建议新增文件：

1. `include/newcode/ofec/mux/mux_state_builder.hpp`
2. `src/rx/ofec/mux/mux_state_builder.cpp`
3. `include/newcode/ofec/mux/mux_siso_budget.hpp`
4. `src/rx/ofec/mux/mux_siso_budget.cpp`
5. `include/newcode/ofec/mux/mux_config_validate.hpp`
6. `src/rx/ofec/mux/mux_config_validate.cpp`

文件职责说明：

| 文件 | 作用 | 典型调用方 |
|---|---|---|
| `include/newcode/ofec/mux/mux_state_builder.hpp` | 对外声明 `state` 构建接口与状态枚举常量（0/1/2），供解码流程调用。 | `ofec_tile_impl.ipp` |
| `src/rx/ofec/mux/mux_state_builder.cpp` | 实现“从 `TileEarlyStopResult` 生成初始 `state`（只产生 0/1）”的纯函数逻辑。 | `ofec_tile_impl.ipp` |
| `include/newcode/ofec/mux/mux_siso_budget.hpp` | 对外声明 SISO 容量裁剪与按 tile 读取预算的接口。 | `ofec_window_impl.ipp`、`ofec_tile_impl.ipp` |
| `src/rx/ofec/mux/mux_siso_budget.cpp` | 实现 `G=1` 下预算裁剪（把超预算的 `0` 改为 `2`）和 `SISO_ACTIVE_LIST[t]` 读取逻辑。 | `ofec_window_impl.ipp`、`ofec_tile_impl.ipp` |
| `include/newcode/ofec/mux/mux_config_validate.hpp` | 对外声明 MUX 配置校验接口。 | `ofec_single_params.cpp`、sweep 参数构建逻辑 |
| `src/rx/ofec/mux/mux_config_validate.cpp` | 实现 `SISO_ACTIVE_LIST` 长度与取值合法性校验，统一错误信息。 | `ofec_single_params.cpp`、sweep 参数构建逻辑 |

建议函数归属：

1. `build_state_from_early_stop(const TileEarlyStopResult&)`
   - 放在 `mux_state_builder.*`
   - 功能：仅把 `row_passed_flags` 映射为 `state` 的 `1/0` 初值。
2. `apply_siso_budget_g1(std::vector<uint8_t>& state, int siso_active_for_tile)`
   - 放在 `mux_siso_budget.*`
   - 功能：按索引顺序把超预算的 `0` 改成 `2`。
3. `validate_siso_active_list(const std::vector<int>&, size_t tiles_per_win)`
   - 放在 `mux_config_validate.*`
   - 功能：校验长度与数值合法性（非负、长度匹配）。
4. `pick_siso_active_for_tile(const std::vector<int>&, size_t t)`
   - 放在 `mux_siso_budget.*`
   - 功能：按 bottom->top 语义读取 `SISO_ACTIVE_LIST[t]`。

函数级说明（输入/输出/副作用）：

| 函数 | 输入 | 输出 | 副作用 | 说明 |
|---|---|---|---|---|
| `build_state_from_early_stop(const TileEarlyStopResult& stats)` | `stats.row_passed_flags` | `std::vector<uint8_t> state`（仅 0/1） | 无 | `true->1`，`false->0`，不做预算裁剪。 |
| `apply_siso_budget_g1(std::vector<uint8_t>& state, int siso_active_for_tile)` | 已初始化的 `state`（含 0/1）和本 tile SISO 容量 | 原地更新后的 `state`（可能出现 2） | 修改 `state` | 仅处理 `state==0`，按索引顺序保留前 `siso_active_for_tile` 个，其余改 `2`。 |
| `validate_siso_active_list(const std::vector<int>& list, size_t tiles_per_win)` | `list`、`TILES_PER_WIN` | `bool` 或 `Status`（由实现选型） | 无 | 校验长度等于 `tiles_per_win`，且元素非负；失败返回可打印的错误原因。 |
| `pick_siso_active_for_tile(const std::vector<int>& list, size_t t)` | `list`、tile 索引 `t` | `int siso_active_for_tile` | 无 | 按 bottom->top 语义读取 `list[t]`；越界时按约定报错或回退默认值。 |

建议补充一个轻量统计函数（可选但推荐）：

| 函数 | 输入 | 输出 | 副作用 | 说明 |
|---|---|---|---|---|
| `count_state012(const std::vector<uint8_t>& state)` | `state` | `(n0, n1, n2)` 计数结构 | 无 | 用于日志和验收（确认预算裁剪正确性）。 |

若采用该函数，建议放在：

- `include/newcode/ofec/mux/mux_state_builder.hpp`
- `src/rx/ofec/mux/mux_state_builder.cpp`

现有文件中的角色调整：

1. `src/rx/ofec/detail/ofec_tile_impl.ipp`
   - 只负责调用 `build_state_from_early_stop(...)` 与 `apply_siso_budget_g1(...)`。
2. `src/rx/ofec/detail/ofec_window_impl.ipp`
   - 只负责取 `t`、读取 `siso_active_for_tile`、传参，不内联预算裁剪逻辑。
3. `src/ofec_single/ofec_single_params.cpp`
   - 只负责参数透传，并调用 `validate_siso_active_list(...)`。

## 9. 测试与验收场景

1. 长度校验  
`SISO_ACTIVE_LIST.size() != TILES_PER_WIN` 时应报错，风格与 `ALPHA_LIST/beta_list` 显式长度校验一致。

2. 容量裁剪正确性  
若某 tile 初始 `state==0` 数量为 23，且本级 `SISO=16`，则应有 7 个 `0` 被改为 `2`，`1` 的数量保持不变。

3. 早停复用正确性  
使用 `early_stop_stats.row_passed_flags` 初始化 `state` 的结果，应与独立重复计算 early-stop 的结果一致。

4. 索引语义一致性  
通过日志或注释示例确认 `t=0` 对应 bottom tile，`[32,30,24,16]` 按 bottom->top 解释。

## 10. 默认假设

1. 本文档只定义规格，不直接改代码。
2. 覆盖入口为 `ofec_single + ofec_sweep`，不包含 `ofec_sweep2`。
3. 当前只定义 `G=1`。
4. `state=2` 的执行策略暂不在本文档中定案。
5. 容量裁剪优先级先固定为“索引顺序”。

## 11. 代码组织验收标准（新增）

1. 新增逻辑函数均位于 `include/newcode/ofec/mux/*` 与 `src/rx/ofec/mux/*`。
2. `ofec_tile_impl.ipp` / `ofec_window_impl.ipp` 中不再包含大段状态构建与预算裁剪实现，只保留调用与流程编排。
3. 配置校验逻辑集中在 `mux_config_validate.*`，不在多个入口重复实现。

## 12. 直接可用的建议函数签名（实现时可按此落地）

```cpp
// include/newcode/ofec/mux/mux_state_builder.hpp
#pragma once
#include <cstdint>
#include <vector>
#include "newcode/ofec/earlystop/tile_early_stop_result.hpp"

namespace newcode::mux {

enum class StateTag : uint8_t {
  NeedSiso = 0,      // 0
  EarlyStopped = 1,  // 1
  Unscheduled = 2    // 2
};

std::vector<uint8_t> build_state_from_early_stop(
    const TileEarlyStopResult& stats);

struct StateCount {
  std::size_t n0_need_siso = 0;
  std::size_t n1_early_stopped = 0;
  std::size_t n2_unscheduled = 0;
};

StateCount count_state012(const std::vector<uint8_t>& state);

} // namespace newcode::mux
```

```cpp
// include/newcode/ofec/mux/mux_siso_budget.hpp
#pragma once
#include <cstddef>
#include <cstdint>
#include <vector>

namespace newcode::mux {

int pick_siso_active_for_tile(const std::vector<int>& list, std::size_t t);

void apply_siso_budget_g1(std::vector<uint8_t>& state,
                          int siso_active_for_tile);

} // namespace newcode::mux
```

```cpp
// include/newcode/ofec/mux/mux_config_validate.hpp
#pragma once
#include <cstddef>
#include <string>
#include <vector>

namespace newcode::mux {

struct ValidationResult {
  bool ok = true;
  std::string error;
};

ValidationResult validate_siso_active_list(const std::vector<int>& list,
                                           std::size_t tiles_per_win);

} // namespace newcode::mux
```

## 13. 逐文件改造步骤（按顺序执行）

### 第 1 步：参数结构加字段

1. `include/newcode/params.hpp`
   - 增加：`std::vector<int> SISO_ACTIVE_LIST = {32, 30, 24, 16};`
   - 说明注释：索引 `t` 为 bottom->top。

2. `include/newcode/ofec_single_runner.hpp`
   - `ofec_single::Config` 增加：`std::vector<int> siso_active_list;`

3. `include/newcode/ofec_sweep_runner.hpp`
   - `SweepParameterConfig` 增加：`std::vector<int> siso_active_list;`

### 第 2 步：apps 配置透传

1. `apps/ofec_single.cpp`
   - 新增常量：`kSisoActiveList = {32,30,24,16};`
   - 在 `Config` 初始化中赋值给 `siso_active_list`。

2. `apps/ofec_sweep.cpp`
   - 新增常量：`kSisoActiveList = {32,30,24,16};`
   - 写入 `config.siso_active_list` 与 `config.base_params.SISO_ACTIVE_LIST`。

### 第 3 步：参数校验与拷贝

1. `src/ofec_single/ofec_single_params.cpp`
   - 把 `cfg.siso_active_list` 透传到 `params.SISO_ACTIVE_LIST`（若非空则覆盖默认）。
   - 调用 `validate_siso_active_list(...)`。
   - 校验失败时打印错误并返回 `std::nullopt`（风格对齐 alpha/beta 显式列表）。

2. sweep 参数构建路径（`src/ofec_sweep/*`）
   - 在开始跑场景前执行同样校验。

### 第 4 步：新增 MUX 模块文件

1. 新增：
   - `include/newcode/ofec/mux/mux_state_builder.hpp`
   - `src/rx/ofec/mux/mux_state_builder.cpp`
   - `include/newcode/ofec/mux/mux_siso_budget.hpp`
   - `src/rx/ofec/mux/mux_siso_budget.cpp`
   - `include/newcode/ofec/mux/mux_config_validate.hpp`
   - `src/rx/ofec/mux/mux_config_validate.cpp`

2. `CMakeLists.txt`
   - 把新增 `.cpp` 加入构建目标。

### 第 5 步：接入 tile 主流程

1. `src/rx/ofec/detail/ofec_tile_impl.ipp`
   - 在已有 `TileEarlyStopResult early_stop_stats` 之后：
     - 调用 `build_state_from_early_stop(early_stop_stats)` 得到 `state`。
     - 使用当前 tile 的 `siso_active_for_tile` 调用 `apply_siso_budget_g1(state, ...)`。
   - 把 `state` 传给 row core（新增参数）。

2. `include/newcode/rx/ofec/chase/decoder_core.hpp`
   - `Decoder_Core_plain/ebchPF` 增加参数：`const std::vector<uint8_t>* mux_state`。

3. `src/rx/ofec/ofec_row_decoder_core.cpp`
   - 行分支改为按 `mux_state[row]`：
     - `1` -> 早停路径
     - `0` -> Chase
     - `2` -> 本阶段先“不更新/不 produced”

### 第 6 步：在 window 层按 tile 取预算

1. `src/rx/ofec/detail/ofec_window_impl.ipp`
   - 在每个 tile `t` 开始处理前：
     - `int siso_active_for_tile = pick_siso_active_for_tile(p.SISO_ACTIVE_LIST, t);`
   - 把 `siso_active_for_tile` 传递给 `process_tile_impl(...)`（新增入参）。

## 14. 关键实现细节（避免踩坑）

1. 不要改变已有 early-stop 判定函数  
继续使用 `tile_early_stop_stats1(prep.lin_matrix)`，只复用其结果。

2. 容量裁剪只处理 `state==0`  
`state==1` 绝不改写。

3. `state==2` 在本阶段不引入新算法  
先统一为“本轮不更新该行”，后续再引入 fallback 策略。

4. tile 索引语义统一  
所有读取 `SISO_ACTIVE_LIST[t]` 的位置都用 bottom->top 语义，不允许混用 top->bottom。

5. 长度不匹配必须早失败  
不要 silent fallback，避免调试困难。

## 15. 日志与调试建议（建议但非必须）

建议在 `trace` 打开时每个 tile 打印：

- `tile_index=t`
- `siso_active_for_tile`
- `state` 计数 `(n0,n1,n2)`

这样可直接验证 `[32,30,24,16]` 是否按预期生效。

## 16. 最小回归清单（改完代码后必须通过）

1. 基线一致性  
`SISO_ACTIVE_LIST` 全设为较大值（例如 32）时，结果应接近改造前。

2. 长度校验  
`SISO_ACTIVE_LIST` 长度错误时应立即报错退出。

3. 容量生效  
设置 `[32,30,24,16]` 后，各 tile 的 `n2` 变化应符合直觉（后级更易出现 `2`）。

4. 语义正确  
`n1`（早停数）不应被预算裁剪改变。

5. 可构建  
新增 `mux/*.cpp` 后工程可正常编译链接。
