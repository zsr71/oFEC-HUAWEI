# oFEC MUX 分组 `G` 落地清单（代码改动与新增文件）

## 1. 目标与范围
目标是在现有 `max(G=1)` 逻辑上，支持可配置分组 `G`：
- `G=1` 保持当前行为（全局池化裁剪）。
- `G>1` 时，按组进行 SISO 预算裁剪，组间不可借用。

本文件只定义落地方案，不直接改代码。

## 2. 行为定义（落地后应满足）
输入：
- `state`：长度为 `N_code`，元素含义 `0/1/2`。
- `N_siso`：当前 tile 可用 SISO 数。
- `G`：分组数。

规则：
1. 把 `N_code` 切成 `G` 个 code 组（尽量均分，余数给前几组）。
2. 把 `N_siso` 切成 `G` 个 budget 组（尽量均分，余数给前几组）。
3. 对每个组独立处理：
- 组内 `state==1` 不变。
- 组内 `state==0` 只保留前 `budget[g]` 个。
- 超出的 `0` 改为 `2`。

## 3. 需要修改的现有文件

### 3.1 参数与配置入口
- [params.hpp](/home/zsr71/projects/newcode/include/newcode/params.hpp)  
作用：新增全局参数 `MUX_GROUP_G`（默认 `1`）。
建议新增字段：
```cpp
int MUX_GROUP_G = 1;  // 1=现有max逻辑；>1=分组裁剪
```

- [ofec_single_runner.hpp](/home/zsr71/projects/newcode/include/newcode/ofec_single_runner.hpp)  
作用：`single` 配置入口支持传入 `G`。  
建议新增字段：
```cpp
int mux_group_g = 1;
```

- [ofec_sweep_runner.hpp](/home/zsr71/projects/newcode/include/newcode/ofec_sweep_runner.hpp)  
作用：`sweep` 配置入口支持传入 `G`。  
建议新增字段：
```cpp
int mux_group_g = 1;
// 可选扩展：std::vector<int> mux_group_g_candidates;
```

- [apps/ofec_single.cpp](/home/zsr71/projects/newcode/apps/ofec_single.cpp)  
作用：用户可直接设置 `kMuxGroupG`，并写入 `Config`。

- [apps/ofec_sweep.cpp](/home/zsr71/projects/newcode/apps/ofec_sweep.cpp)  
作用：用户可直接设置 `kMuxGroupG`，并写入 `SweepParameterConfig`。

- [ofec_single_params.cpp](/home/zsr71/projects/newcode/src/ofec_single/ofec_single_params.cpp)  
作用：把 `cfg.mux_group_g` 透传到 `params.MUX_GROUP_G`，并做合法性校验。

- [ofec_sweep_runner.cpp](/home/zsr71/projects/newcode/src/ofec_sweep/ofec_sweep_runner.cpp)  
作用：在 `resolved` 配置阶段透传 `mux_group_g` 到 `base_params`，并统一校验。

### 3.2 解码主流程接线
- [ofec_decode_impl.ipp](/home/zsr71/projects/newcode/src/rx/ofec/detail/ofec_decode_impl.ipp)  
作用：解码入口处统一校验 `G` 配置（防止运行时异常）。

- [ofec_tile_impl.ipp](/home/zsr71/projects/newcode/src/rx/ofec/detail/ofec_tile_impl.ipp)  
作用：把当前
```cpp
apply_siso_budget_g1(mux_state, siso_active_for_tile);
```
替换为分组版本：
```cpp
apply_siso_budget_grouped(mux_state, siso_active_for_tile, p.MUX_GROUP_G);
```

- [CMakeLists.txt](/home/zsr71/projects/newcode/CMakeLists.txt)  
作用：把新增 `.cpp` 源文件加入构建目标。

## 4. 建议新增的代码文件

### 4.1 新增分组预算核心模块
- 新文件：[mux_group_budget.hpp](/home/zsr71/projects/newcode/include/newcode/ofec/mux/mux_group_budget.hpp)  
- 新文件：[mux_group_budget.cpp](/home/zsr71/projects/newcode/src/rx/ofec/mux/mux_group_budget.cpp)  

作用：
- 负责 `G` 分组切片、预算分配、组内裁剪。
- 不和 early-stop 判定耦合，只处理 `state` 与 budget。

建议包含的数据结构：
```cpp
struct GroupRange {
  std::size_t begin = 0;  // 含
  std::size_t end = 0;    // 不含
};
```

建议包含的函数（对外 API）：
```cpp
std::vector<GroupRange> build_groups_even(std::size_t total, int group_g);
std::vector<int> split_budget_even(int total_budget, int group_g);
void apply_siso_budget_grouped(std::vector<uint8_t>& state,
                               int siso_active_for_tile,
                               int group_g);
```

函数职责说明：
- `build_groups_even`：把 `N_code` 按 `G` 均分成区间。
- `split_budget_even`：把 `N_siso` 按 `G` 均分成每组预算。
- `apply_siso_budget_grouped`：按组把超额 `state==0` 改成 `2`。

#### 4.1.1 `mux_group_budget.cpp` 需要包含的函数（详细）
以下是建议的完整函数清单，分为“对外 API”与“文件内 helper”。

1. 对外 API：`build_groups_even(std::size_t total, int group_g)`  
作用：
- 把 `[0, total)` 划分成 `group_g` 个连续区间。
- 尽量均分，余数给前几个组（保证每组差值最多 1）。
输出：
- `std::vector<GroupRange>`，长度为 `group_g`。

2. 对外 API：`split_budget_even(int total_budget, int group_g)`  
作用：
- 把总预算 `total_budget` 分摊到 `group_g` 组。
- 均分 + 余数前置，与 `build_groups_even` 的分配风格一致。
输出：
- `std::vector<int>`，长度为 `group_g`，各组 budget 非负，总和等于 `total_budget`。

3. 对外 API：`apply_siso_budget_grouped(std::vector<uint8_t>& state, int siso_active_for_tile, int group_g)`  
作用：
- 分组预算裁剪的主入口。
- `group_g <= 1` 时直接回退到 `apply_siso_budget_g1`（保持 max 逻辑兼容）。
- `group_g > 1` 时：按组独立裁剪 `state==0`，超额置 `2`。
副作用：
- 原地修改 `state`。

4. 文件内 helper（`static`）：`validate_grouped_budget_args(std::size_t code_count, int siso_active_for_tile, int group_g)`  
作用：
- 对运行时参数做基础防御性检查（例如 `group_g>=1`、`siso_active_for_tile>=0`、`group_g<=code_count`）。
- 非法时抛 `std::invalid_argument`，避免 silent bug。

5. 文件内 helper（`static`）：`collect_need_indices_in_range(const std::vector<uint8_t>& state, std::size_t begin, std::size_t end)`  
作用：
- 收集某个组区间 `[begin,end)` 内所有 `state==0` 的索引。
- 只做收集，不修改数据。
输出：
- `std::vector<std::size_t>`（该组待分配索引列表）。

6. 文件内 helper（`static`）：`trim_need_indices_by_budget(std::vector<uint8_t>& state, const std::vector<std::size_t>& need, int keep_budget)`  
作用：
- 对单组执行“保留前 `keep_budget` 个，剩余改 `2`”。
- 独立成函数可复用，也方便单测覆盖边界（`keep_budget=0`、`keep_budget>=need.size()`）。

7. 文件内 helper（可选）：`apply_group_budget_once(std::vector<uint8_t>& state, GroupRange range, int budget)`  
作用：
- 把“单组处理流程（收集+裁剪）”封装成一个步骤函数。
- 让主函数 `apply_siso_budget_grouped` 更短、更易读。

#### 4.1.2 `mux_group_budget.cpp` 内部调用顺序（建议）
1. `apply_siso_budget_grouped` 先调用 `validate_grouped_budget_args`。  
2. 若 `group_g<=1`，直接调用现有 `apply_siso_budget_g1` 并返回。  
3. 调用 `build_groups_even(state.size(), group_g)` 得到 code 分组。  
4. 调用 `split_budget_even(siso_active_for_tile, group_g)` 得到每组 budget。  
5. 循环每组：
- `collect_need_indices_in_range` 收集组内 `0` 索引；
- `trim_need_indices_by_budget` 把超额索引改成 `2`。

#### 4.1.3 每个函数建议放置位置
- `mux_group_budget.hpp`：仅声明 `GroupRange`、`build_groups_even`、`split_budget_even`、`apply_siso_budget_grouped`。  
- `mux_group_budget.cpp`：实现全部逻辑；helper 函数全部放匿名命名空间（或 `static`）避免泄露符号。  

### 4.2 新增分组配置校验模块（推荐）
- 新文件：[mux_group_config_validate.hpp](/home/zsr71/projects/newcode/include/newcode/ofec/mux/mux_group_config_validate.hpp)  
- 新文件：[mux_group_config_validate.cpp](/home/zsr71/projects/newcode/src/rx/ofec/mux/mux_group_config_validate.cpp)  

作用：
- 统一校验 `G` 的合法性，避免在业务流程中散落判断。

建议包含函数：
```cpp
ValidationResult validate_mux_group_g(int group_g, std::size_t code_count);
ValidationResult validate_mux_group_runtime(int group_g,
                                            int siso_active_for_tile,
                                            std::size_t code_count);
```

检查项建议：
- `group_g >= 1`
- `group_g <= code_count`
- `siso_active_for_tile >= 0`

## 5. 与现有文件的职责边界（避免混乱）
- `mux_state_builder.*`：只负责 `early_stop_stats -> state(0/1)` 与统计，不做预算裁剪。
- `mux_siso_budget.*`：保留当前 `G=1` 实现，作为基线和回归对照。
- `mux_group_budget.*`：只负责 `G>1` 或统一 `G>=1` 的组预算裁剪。
- `mux_config_validate.*`：保留列表长度/非负值校验。
- `mux_group_config_validate.*`：新增 `G` 相关校验。

## 6. 关键函数落点与调用关系
调用链建议：
1. `process_tile_impl` 先构造 `mux_state`。
2. 调用 `apply_siso_budget_grouped(state, siso_active_for_tile, p.MUX_GROUP_G)`。
3. 结果继续传给 row decoder（`state=0/1/2` 语义不变）。

这样 row core 不需要改接口，也不需要知道 `G`。

## 7. 伪代码（落地参考）
```cpp
void apply_siso_budget_grouped(std::vector<uint8_t>& state,
                               int siso_active_for_tile,
                               int group_g) {
  if (group_g <= 1) {
    apply_siso_budget_g1(state, siso_active_for_tile);
    return;
  }
  auto groups = build_groups_even(state.size(), group_g);
  auto budgets = split_budget_even(siso_active_for_tile, group_g);
  for (int g = 0; g < group_g; ++g) {
    std::vector<std::size_t> need;
    for (std::size_t i = groups[g].begin; i < groups[g].end; ++i) {
      if (state[i] == 0u) need.push_back(i);
    }
    for (std::size_t k = static_cast<std::size_t>(budgets[g]); k < need.size(); ++k) {
      state[need[k]] = 2u;
    }
  }
}
```

## 8. 验收点（建议）
1. `G=1` 时结果与当前版本一致。  
2. `N_code=32, N_siso=16, G=2` 时每组最多保留 8 个 `state==0`。  
3. `G` 增大时，总 `state==2` 通常不减（相同输入下）。  
4. 非法配置（如 `G=0`）应在入口直接报错。  

## 9. 推荐实施顺序
1. 新增 `mux_group_budget.*`（先不接主流程，写单元测试/小样例）。  
2. 增加 `Params/Config` 的 `MUX_GROUP_G` 透传与校验。  
3. 在 `ofec_tile_impl.ipp` 接入 grouped budget。  
4. `single/sweep` 增加参数并验证回归（先跑 `G=1` 再跑 `G=2`）。  

## 10. 具体例子：`N_code=32`、`co(=SISO)=24`、`G=4`
这里按“`co` 表示当前 tile 可用 SISO 数”解释，即：
- `state.size() = 32`
- `siso_active_for_tile = 24`
- `group_g = 4`

分组后：
1. code 分组：`32 / 4 = 8`，所以每组 8 个 code。  
组区间为：
- 组0：`[0..7]`
- 组1：`[8..15]`
- 组2：`[16..23]`
- 组3：`[24..31]`

2. budget 分组：`24 / 4 = 6`，所以每组预算 6 个 SISO。  
各组 budget：`[6, 6, 6, 6]`

运行规则（每组独立）：
1. 统计该组内 `state==0` 的数量（记为 `need_g`）。
2. 若 `need_g <= 6`：该组不裁剪。
3. 若 `need_g > 6`：超出的 `need_g-6` 个 `0` 改成 `2`（按索引顺序）。

一个可视化样例：
- 组0 `need_0=8` -> 改 `2` 个（保留 6）
- 组1 `need_1=5` -> 改 `0` 个
- 组2 `need_2=7` -> 改 `1` 个
- 组3 `need_3=4` -> 改 `0` 个

最终：
- 总改动数 `= 2 + 0 + 1 + 0 = 3`
- 总 `state==2` 会在原有基础上再增加 3

对比 `G=1`（max 全局池化）：
- `G=1` 时 24 个 SISO 可以在 32 个 code 上全局共享；
- `G=4` 时被强制切成每组 6，组间不能借用，所以会出现“某组超额但别组有余量”的额外 `state==2`。
