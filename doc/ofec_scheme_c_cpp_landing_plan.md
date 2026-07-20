# 在 C++ 仿真中把 `MUX_ENABLE_RECONFIG=true` 切换为方案 C 的修改方案

## 1. 目标
当前 C++ 中的 MUX 调度分成两条路径：

1. `MUX_ENABLE_RECONFIG=false`
   走原有的按 `G` 分组预算裁剪逻辑。
2. `MUX_ENABLE_RECONFIG=true`
   走当前的方案 B，也就是图匹配式两阶段重排逻辑。

本次修改目标是：

- 保持 `MUX_ENABLE_RECONFIG=false` 不变，继续执行原有 `apply_siso_budget_grouped(...)`
- 把 `MUX_ENABLE_RECONFIG=true` 的执行逻辑，从当前方案 B 改成方案 C

也就是说，修改后的行为应为：

- `false` -> 原始 `G` 分组静态预算
- `true` -> 方案 C 顺序式两阶段调度

## 2. 修改后的行为定义
### 2.1 `MUX_ENABLE_RECONFIG=false`
完全保持当前行为不变：

1. 根据 early-stop 生成 `mux_state`
2. 调用 `apply_siso_budget_grouped(...)`
3. 对超额的 `state==0` 改成 `2`

这条路径不需要改调度算法。

### 2.2 `MUX_ENABLE_RECONFIG=true`
不再执行当前 `scheme_b_reconfig` 的匹配式调度，而改为方案 C：

1. 先调度每组前 `code_per_group - 2` 个普通 code
2. 第一阶段只允许使用本组 SISO
3. 再调度每组最后 `2` 个边界 code
4. 第二阶段边界 code 先尝试本组剩余 SISO，再尝试额外旁路线
5. 未调度上的 code 最终回写为 `state=2`

## 3. 方案 C 的 C++ 语义定义
以当前主要场景 `N_code=32, N_siso=16, G=4` 为例：

- 每组 `8` 个 code
- 每组 `4` 个 SISO
- 每组前 `6` 个 code 是普通 code
- 每组后 `2` 个 code 是边界 code

但为了代码可复用，建议 C++ 实现写成泛化形式：

- `code_per_group = n_code / G`
- `siso_per_group = n_siso / G`
- `normal_count = code_per_group - 2`
- `tail_count = 2`

前提约束：

1. `n_code % G == 0`
2. `n_siso % G == 0`
3. `code_per_group >= 2`

## 4. 需要修改的现有代码位置
### 4.1 tile 调度接线位置
核心改动位置：

- [ofec_tile_impl.ipp](/home/zsr71/projects/newcode/src/rx/ofec/detail/ofec_tile_impl.ipp)

当前这里是：

```cpp
if (p.MUX_ENABLE_RECONFIG) {
  const auto active_codes = collect_active_codes_from_state(mux_state);
  const auto free_siso = build_free_siso_list(siso_active_for_tile);
  const auto schedule = schedule_scheme_b_reconfig_cpp(...);
  apply_schedule_result_to_mux_state(mux_state, schedule.final_code_to_siso);
} else {
  apply_siso_budget_grouped(...);
}
```

要改成：

```cpp
if (p.MUX_ENABLE_RECONFIG) {
  const auto active_codes = collect_active_codes_from_state(mux_state);
  const auto free_siso = build_free_siso_list(siso_active_for_tile);
  const auto schedule = schedule_scheme_c_staged_cpp(...);
  apply_schedule_result_to_mux_state(mux_state, schedule.final_code_to_siso);
} else {
  apply_siso_budget_grouped(...);
}
```

注意：
- `else` 分支保持不变
- 只替换 `true` 分支里的调度函数

### 4.2 参数与校验位置
这些位置仍然保留：

- [params.hpp](/home/zsr71/projects/newcode/include/newcode/params.hpp)
- [ofec_single_runner.hpp](/home/zsr71/projects/newcode/include/newcode/ofec_single_runner.hpp)
- [ofec_sweep_runner.hpp](/home/zsr71/projects/newcode/include/newcode/ofec_sweep_runner.hpp)
- [ofec_single_params.cpp](/home/zsr71/projects/newcode/src/ofec_single/ofec_single_params.cpp)
- [ofec_sweep_runner.cpp](/home/zsr71/projects/newcode/src/ofec_sweep/ofec_sweep_runner.cpp)
- [ofec_decode_impl.ipp](/home/zsr71/projects/newcode/src/rx/ofec/detail/ofec_decode_impl.ipp)

原因是：
- `MUX_ENABLE_RECONFIG` 这个开关仍然需要保留
- `MUX_EXTRA_BYPASS_EDGES` 仍然要被方案 C 使用
- `reconfig` 路径仍然需要 `n_code/n_siso` 可整除 `G`

因此这些参数透传和校验逻辑不需要推翻，只需要确认还能满足方案 C 的前提。

## 5. 建议新增/替换的文件
### 5.1 新增方案 C 头文件
建议新增：

- `include/newcode/ofec/mux/mux_scheme_c_staged.hpp`

作用：
- 声明方案 C 的数据结构和对外函数

建议内容：

```cpp
#pragma once

#include <vector>

#include "newcode/ofec/mux/mux_topology.hpp"

namespace newcode::mux {

struct SchemeCResult {
  std::vector<int> stage1_code_to_siso;
  std::vector<int> final_code_to_siso;
  std::vector<int> waiting_codes;
};

SchemeCResult schedule_scheme_c_staged_cpp(
    int n_code,
    int n_siso,
    int group_g,
    const std::vector<int>& active_codes,
    const std::vector<int>& free_siso,
    const std::vector<MuxEdge>& extra_bypass_edges);

}  // namespace newcode::mux
```

### 5.2 新增方案 C 源文件
建议新增：

- `src/rx/ofec/mux/mux_scheme_c_staged.cpp`

作用：
- 承载方案 C 的核心顺序调度逻辑

## 6. `mux_scheme_c_staged.cpp` 中需要包含的函数
### 6.1 对外主函数
```cpp
SchemeCResult schedule_scheme_c_staged_cpp(
    int n_code,
    int n_siso,
    int group_g,
    const std::vector<int>& active_codes,
    const std::vector<int>& free_siso,
    const std::vector<MuxEdge>& extra_bypass_edges);
```

作用：
- 执行方案 C 的完整两阶段调度
- 输出阶段1结果、最终结果、等待队列

### 6.2 内部辅助函数建议
#### `build_bypass_map(...)`
```cpp
std::vector<std::vector<int>> build_bypass_map(
    int n_code,
    int n_siso,
    const std::vector<MuxEdge>& extra_bypass_edges);
```

作用：
- 把旁路线列表转成 `code -> bypass siso list`
- 供第二阶段直接查询

#### `try_assign_first_free(...)`
```cpp
bool try_assign_first_free(
    int code_idx,
    const std::vector<int>& candidate_sisos,
    const std::vector<uint8_t>& free_siso_mask,
    std::vector<uint8_t>& used_siso,
    std::vector<int>& code_to_siso);
```

作用：
- 按候选 SISO 顺序找到第一个空闲单元
- 如果找到，就完成分配

#### `build_local_siso_list_for_group(...)`
```cpp
std::vector<int> build_local_siso_list_for_group(
    int group_idx,
    int siso_per_group);
```

作用：
- 返回某一组本地 SISO 列表

#### `build_normal_code_list_for_group(...)`
```cpp
std::vector<int> build_normal_code_list_for_group(
    int group_idx,
    int code_per_group);
```

作用：
- 返回某组第一阶段要调度的普通 code 列表

#### `build_tail_code_list_for_group(...)`
```cpp
std::vector<int> build_tail_code_list_for_group(
    int group_idx,
    int code_per_group);
```

作用：
- 返回某组第二阶段要调度的最后两个边界 code

## 7. 方案 C 的推荐实现流程
### 7.1 输入准备
主函数收到：

- `active_codes`
- `free_siso`
- `extra_bypass_edges`

先构造：

- `free_siso_mask`
- `used_siso`
- `code_to_siso`，初始全 `-1`
- `bypass_map`

### 7.2 第一阶段：普通 code 先调度
对每个组依次执行：

1. 取该组前 `code_per_group - 2` 个 code
2. 只尝试该组本地 SISO
3. 按 code 顺序找第一个空闲 SISO
4. 若本组 SISO 全满，则提前结束该组第一阶段

并把这一阶段分配成功的 code 记录到：

- `stage1_code_to_siso`

### 7.3 第二阶段：边界 code 后调度
再对每个组依次执行：

1. 取该组最后 `2` 个 code
2. 若 code 不在 `active_codes` 中，跳过
3. 先尝试本组剩余 SISO
4. 若失败，再尝试 `bypass_map[code]`
5. 若仍失败，则保留未分配状态

### 7.4 结果输出
最后输出：

- `stage1_code_to_siso`
- `final_code_to_siso`
- `waiting_codes`

其中：

- `waiting_codes` = `active_codes` 中仍未分配的 code

## 8. 与现有模块的关系
### 8.1 继续复用的模块
以下模块可以继续复用：

- [mux_state_schedule_apply.cpp](/home/zsr71/projects/newcode/src/rx/ofec/mux/mux_state_schedule_apply.cpp)
  - 继续复用 `collect_active_codes_from_state(...)`
  - 继续复用 `build_free_siso_list(...)`
  - 继续复用 `apply_schedule_result_to_mux_state(...)`

- [mux_group_config_validate.cpp](/home/zsr71/projects/newcode/src/rx/ofec/mux/mux_group_config_validate.cpp)
  - 继续复用 `validate_mux_reconfig_runtime(...)`

### 8.2 不再作为主路径使用的模块
以下模块在 `MUX_ENABLE_RECONFIG=true` 时，不再作为主调度器：

- [mux_scheme_b_reconfig.cpp](/home/zsr71/projects/newcode/src/rx/ofec/mux/mux_scheme_b_reconfig.cpp)
- [mux_matching.cpp](/home/zsr71/projects/newcode/src/rx/ofec/mux/mux_matching.cpp)

注意：
- 不一定要删掉
- 但它们不再是 `true` 分支的默认执行路径

## 9. CMake 需要修改的地方
新增 `mux_scheme_c_staged.cpp` 后，需要把它加入：

- [CMakeLists.txt](/home/zsr71/projects/newcode/CMakeLists.txt)

即加入 `newcode_frontend` 的源文件列表。

## 10. 建议的实施步骤
1. 新增 `mux_scheme_c_staged.hpp/.cpp`
2. 在新文件里实现方案 C 的顺序式两阶段调度
3. 在 `ofec_tile_impl.ipp` 中把 `MUX_ENABLE_RECONFIG=true` 的调用从方案 B 改成方案 C
4. 保持 `MUX_ENABLE_RECONFIG=false` 路径完全不变
5. 编译并验证 `ofec_single` / `ofec_sweep`
6. 用固定 case 做最小回归

## 11. 建议的最小回归用例
### 11.1 保证 `false` 路径不变
用同一组输入验证：

- `MUX_ENABLE_RECONFIG=false`
- 修改前后 `mux_state` 一致

### 11.2 验证方案 C 的阶段性行为
针对 `N_code=32, N_siso=16, G=4`：

1. 第一阶段只允许普通 code 占用本组 SISO
2. 第二阶段边界 code 才允许尝试旁路线
3. 边界 code 若本组尚有空闲，应优先占本组

### 11.3 验证状态回写
最终应满足：

- 原 `state==1` 不变
- 成功调度的 `state==0` 保持 `0`
- 未调度上的 `state==0` 变为 `2`

## 12. 结论
这次修改的本质是：

- `MUX_ENABLE_RECONFIG=false` 保持原有 `G` 分组静态预算
- `MUX_ENABLE_RECONFIG=true` 从“方案 B 图匹配式两阶段重排”切换为“方案 C 顺序式两阶段调度”

因此，这不是单纯替换一个函数，而是把 `reconfig=true` 的语义整体改成新的调度策略。
