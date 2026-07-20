# 在 C++ 项目中实现 `schedule_scheme_b_reconfig` 效果的落地方案

## 1. 目标
当前 C++ 项目的 MUX 仿真逻辑是：
- 先根据 early-stop 生成 `mux_state`
- 再按 `G` 把 `state==0` 做组内预算裁剪
- 超出的 `0` 直接改成 `2`

这对应的是“静态分组预算”模型，不等价于 Python 里的
`schedule_scheme_b_reconfig(...)`。

Python 的目标效果是：
1. 阶段 1：只允许组内边，先得到一个初始匹配。
2. 阶段 2：开放旁路线，在阶段 1 的结果基础上继续增广。
3. 阶段 2 允许重排已有匹配，而不是简单保留阶段 1 结果不动。

因此，如果 C++ 里要实现 `schedule_scheme_b_reconfig` 的效果，核心不是继续改“裁剪规则”，而是要把当前的“计数裁剪”升级成“显式匹配调度”。

## 2. 当前 C++ 逻辑与 Python 逻辑的差异

### 2.1 当前 C++ 做的事情
当前入口在 [ofec_tile_impl.ipp](/home/zsr71/projects/newcode/src/rx/ofec/detail/ofec_tile_impl.ipp)：
1. `tile_early_stop_stats1(...)` 得到 `row_passed_flags`
2. `build_state_from_early_stop(...)` 生成 `mux_state`
3. `apply_siso_budget_grouped(...)` 按 `G` 把超额 `0` 改成 `2`
4. row core 根据 `mux_state` 选择：
- `0` -> 正常 Chase
- `1` -> early-stop
- `2` -> skip

当前组预算实现见 [mux_group_budget.cpp](/home/zsr71/projects/newcode/src/rx/ofec/mux/mux_group_budget.cpp)。

### 2.2 Python `schedule_scheme_b_reconfig` 做的事情
实现见 [group_scheduler.py](/home/zsr71/projects/newcode/python/group_scheduler.py:162)：
1. `build_local_edges(...)` 构造组内边
2. `build_allowed_edges(...)` 构造组内边 + 旁路线
3. `maximum_bipartite_matching(...)` 先做阶段 1 组内最大匹配
4. `augment_from_seed_matching(...)` 再在已有匹配上继续增广

这里的关键差异是：
- Python 是“边级别”的调度；
- C++ 现在只是“每组保留前 K 个 `0`”，没有显式的 code-SISO 匹配关系；
- Python 会因为旁路线和增广路而发生“重排”，C++ 当前不会。

## 3. 在 C++ 里要实现的最小效果
为了在软件仿真中达到 `schedule_scheme_b_reconfig` 的效果，建议目标定义为：

输入：
- `state`
- `N_code`
- `N_siso`
- `G`
- `active_codes`
- `free_siso`
- `extra_bypass_edges`

输出：
- `stage1_match`
- `final_match`
- `waiting_codes`
- 最终回写为 `mux_state`

回写规则建议：
- `state==1` 的行保持 `1`（早停）
- `state==0` 且匹配成功的行保持 `0`（本轮进入 SISO）
- `state==0` 但最终未匹配成功的行改成 `2`

也就是说，`mux_state` 不再由“组内计数裁剪”直接决定，而要由“匹配调度结果”决定。

## 4. 需要修改的现有文件

### 4.1 tile 主流程接线
- [ofec_tile_impl.ipp](/home/zsr71/projects/newcode/src/rx/ofec/detail/ofec_tile_impl.ipp)

当前这里是：
```cpp
std::vector<uint8_t> mux_state =
    newcode::mux::build_state_from_early_stop(early_stop_stats);
newcode::mux::apply_siso_budget_grouped(
    mux_state, siso_active_for_tile, p.MUX_GROUP_G);
```

要改成：
1. 先由 `early_stop_stats` 得到初始 `mux_state`
2. 构造 `active_codes`：收集所有 `state==0` 的 code 索引
3. 构造 `free_siso`：默认是 `0..siso_active_for_tile-1`
4. 调用新的 `schedule_scheme_b_reconfig_cpp(...)`
5. 用 `final_match / waiting_codes` 回写 `mux_state`

### 4.2 参数入口
- [params.hpp](/home/zsr71/projects/newcode/include/newcode/params.hpp)
- [ofec_single_runner.hpp](/home/zsr71/projects/newcode/include/newcode/ofec_single_runner.hpp)
- [ofec_sweep_runner.hpp](/home/zsr71/projects/newcode/include/newcode/ofec_sweep_runner.hpp)
- [ofec_single_params.cpp](/home/zsr71/projects/newcode/src/ofec_single/ofec_single_params.cpp)
- [ofec_sweep_runner.cpp](/home/zsr71/projects/newcode/src/ofec_sweep/ofec_sweep_runner.cpp)

当前这些文件已经支持 `MUX_GROUP_G`。  
如果要支持 Python 方案 B，还应再增加：
- `MUX_ENABLE_RECONFIG`：是否启用两阶段重排调度
- `MUX_EXTRA_BYPASS_EDGES`：旁路线配置

## 5. 建议新增的代码文件

### 5.1 新增：边拓扑构造模块
- `include/newcode/ofec/mux/mux_topology.hpp`
- `src/rx/ofec/mux/mux_topology.cpp`

作用：
- 负责构造组内边和全允许边
- 对应 Python 的 `build_local_edges` / `build_allowed_edges`

建议数据结构：
```cpp
struct MuxEdge {
  int code_idx = -1;
  int siso_idx = -1;
};
```

建议函数：
```cpp
std::vector<MuxEdge> build_local_edges(int n_code, int n_siso, int group_g);
std::vector<MuxEdge> build_allowed_edges(int n_code,
                                         int n_siso,
                                         int group_g,
                                         const std::vector<MuxEdge>& extra_bypass_edges);
std::vector<std::vector<int>> build_adjacency(const std::vector<MuxEdge>& edges,
                                              int n_code);
```

### 5.2 新增：匹配求解模块
- `include/newcode/ofec/mux/mux_matching.hpp`
- `src/rx/ofec/mux/mux_matching.cpp`

作用：
- 负责最大匹配和“从已有匹配继续增广”
- 对应 Python 的 `try_augment` / `maximum_bipartite_matching` / `augment_from_seed_matching`

建议数据结构：
```cpp
struct MatchingResult {
  std::vector<int> code_to_siso;   // size=n_code, -1 表示未匹配
  std::vector<int> waiting_codes;
};
```

建议函数：
```cpp
bool try_augment(int code_idx,
                 const std::vector<std::vector<int>>& adjacency,
                 const std::vector<uint8_t>& free_siso_mask,
                 std::vector<int>& match_siso_to_code,
                 std::vector<uint8_t>& visited_siso);

MatchingResult maximum_bipartite_matching(
    const std::vector<int>& active_codes,
    const std::vector<std::vector<int>>& adjacency,
    const std::vector<int>& free_siso);

MatchingResult augment_from_seed_matching(
    const std::vector<int>& active_codes,
    const std::vector<std::vector<int>>& adjacency,
    const std::vector<int>& free_siso,
    const std::vector<int>& seed_code_to_siso);
```

### 5.3 新增：方案 B 调度模块
- `include/newcode/ofec/mux/mux_scheme_b_reconfig.hpp`
- `src/rx/ofec/mux/mux_scheme_b_reconfig.cpp`

作用：
- 作为 C++ 版 `schedule_scheme_b_reconfig` 的直接承载模块
- 把 topology + matching 串起来

建议数据结构：
```cpp
struct SchemeBResult {
  std::vector<int> stage1_code_to_siso;
  std::vector<int> final_code_to_siso;
  std::vector<int> waiting_codes;
};
```

建议函数：
```cpp
SchemeBResult schedule_scheme_b_reconfig_cpp(
    int n_code,
    int n_siso,
    int group_g,
    const std::vector<int>& active_codes,
    const std::vector<int>& free_siso,
    const std::vector<MuxEdge>& extra_bypass_edges);
```

### 5.4 新增：`mux_state` 回写模块
- `include/newcode/ofec/mux/mux_state_schedule_apply.hpp`
- `src/rx/ofec/mux/mux_state_schedule_apply.cpp`

作用：
- 把 `SchemeBResult` 应用回 `mux_state`
- 这个模块把“调度结果”和“row decoder 所需的 `0/1/2` 语义”隔离开

建议函数：
```cpp
std::vector<int> collect_active_codes_from_state(const std::vector<uint8_t>& state);

std::vector<int> build_free_siso_list(int siso_active_for_tile);

void apply_schedule_result_to_mux_state(std::vector<uint8_t>& state,
                                        const std::vector<int>& final_code_to_siso);
```

语义：
- `state==1` 保持不变
- `final_code_to_siso[code_idx] >= 0` 的 code 保持 `0`
- 其余原本 `state==0` 的 code 改成 `2`

## 6. C++ 主调用链应该怎么改
建议改成下面这样：

1. 在 [ofec_tile_impl.ipp](/home/zsr71/projects/newcode/src/rx/ofec/detail/ofec_tile_impl.ipp) 中得到初始 `mux_state`
2. 如果 `MUX_ENABLE_RECONFIG == false`
   继续走现有 `apply_siso_budget_grouped(...)`
3. 如果 `MUX_ENABLE_RECONFIG == true`
   调用新的 C++ 方案 B：
```cpp
auto active_codes = collect_active_codes_from_state(mux_state);
auto free_siso = build_free_siso_list(siso_active_for_tile);
auto sched = schedule_scheme_b_reconfig_cpp(
    static_cast<int>(mux_state.size()),
    siso_active_for_tile,
    p.MUX_GROUP_G,
    active_codes,
    free_siso,
    p.MUX_EXTRA_BYPASS_EDGES);
apply_schedule_result_to_mux_state(mux_state, sched.final_code_to_siso);
```

## 7. 参数建议
建议在 `Params` 中新增：
```cpp
bool MUX_ENABLE_RECONFIG = false;
std::vector<std::pair<int, int>> MUX_EXTRA_BYPASS_EDGES;
```

说明：
- `MUX_ENABLE_RECONFIG=false` 时，保留当前 grouped budget 行为，方便回归
- `MUX_ENABLE_RECONFIG=true` 时，启用 Python 方案 B 风格的两阶段调度
- `MUX_EXTRA_BYPASS_EDGES` 用来配置 `C_i -> S_j` 的跨组旁路线

## 7.1 兼容性说明
这套改法是“在原有 `G` 分组逻辑旁边增加一条可切换路径”，不是替换原实现。

具体来说：
- 当 `MUX_ENABLE_RECONFIG=false` 时，程序继续走当前已经实现的原始 `G` 分组逻辑：
  - 根据 early-stop 生成 `mux_state`
  - 调用 `apply_siso_budget_grouped(...)`
  - 按组做静态预算裁剪
  - 超额的 `state==0` 改成 `2`
- 当 `MUX_ENABLE_RECONFIG=true` 时，程序才切换到 Python `schedule_scheme_b_reconfig` 风格的两阶段调度逻辑。

因此，这个方案可以随时切回原先的 `G` 分组实现：
- 不需要删除旧代码
- 不需要改回接口
- 只需要把 `MUX_ENABLE_RECONFIG` 设回 `false`

也就是说：
- `G` 本身仍然保留
- 原先 grouped budget 的语义仍然保留
- 新增的只是一个可选的“重排调度模式”

## 8. 为什么不能只改 `apply_siso_budget_grouped(...)`
因为 `schedule_scheme_b_reconfig` 的核心不是“每组保留几个 `0`”，而是：
- 存在显式的 code-SISO 边
- 要做最大匹配
- 阶段 2 允许沿增广路重排已有匹配

当前 `apply_siso_budget_grouped(...)` 只知道“每组有多少个 `0`”，并不知道：
- 哪个 code 可以连哪个 SISO
- 某个 code 未匹配是不是因为组内冲突，还是能通过旁路救回来

所以如果目标是“实现 `schedule_scheme_b_reconfig` 效果”，就必须引入显式匹配模块。

## 9. 建议实施顺序
1. 新增 `mux_topology.*`
2. 新增 `mux_matching.*`
3. 新增 `mux_scheme_b_reconfig.*`
4. 新增 `mux_state_schedule_apply.*`
5. 在 `ofec_tile_impl.ipp` 接入可切换路径
6. `G=1`、`G=2`、`G=4` 用固定 case 做回归

## 10. 最小验收用例
建议至少验证以下几类 case：
1. 无旁路线时，方案 B 退化为纯组内匹配  
2. 有旁路线时，`final_match` 数量 >= `stage1_match` 数量  
3. 某些 code 在阶段 1 未匹配，但在阶段 2 通过重排被匹配成功  
4. 回写后的 `mux_state` 满足：
- 原 `1` 不变
- 匹配到的 `0` 保持 `0`
- 未匹配的 `0` 改成 `2`

这份方案的核心是：把当前“计数裁剪模型”升级为“边约束 + 匹配 + 可重排模型”。只有这样，C++ 仿真结果才真正对应 Python 的 `schedule_scheme_b_reconfig(...)`。
