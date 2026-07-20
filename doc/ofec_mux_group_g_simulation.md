# 在现有 max 逻辑上模拟分组 `G` 影响（设计说明）

## 1. 背景
当前实现本质是 **max 连接自由度**，即 `G=1`：
- `state[i] = 1`：该 code 行已早停。
- `state[i] = 0`：该 code 行需要 SISO。
- `state[i] = 2`：该 code 行本轮未分配到 SISO。

现有预算裁剪是“全局池化”：
- 对整个 `state` 里所有 `0` 统一计数。
- 按 `siso_active_for_tile` 保留前若干个 `0`，其余改 `2`。

这等价于：所有 SISO 都能服务所有 code（全互连 / `G=1`）。

## 2. `G` 的含义
`G` 表示分组数（group 数量）：
- code 被切成 `G` 组；
- SISO 资源也切成 `G` 组；
- 每组 SISO 只能服务本组 code，**不能跨组借用**。

`G` 越大：
- MUX 扇入越小，硬件更省；
- 但调度约束更强，可能出现某组忙、某组闲却不能互借。

## 3. 分组仿真规则（建议）
给定：
- `N_code = state.size()`
- `N_siso = siso_active_for_tile`
- `G >= 1`

先把 code 和 SISO 都按组划分：
- `code_group_size = N_code / G`（若不能整除，余数分配给前几个组）
- `siso_group_budget = N_siso / G`（若不能整除，余数分配给前几个组）

对每个组独立执行裁剪：
1. 组内 `state==1` 保持不变。
2. 统计组内 `state==0` 的位置。
3. 保留前 `siso_group_budget[group]` 个 `0`。
4. 超出的 `0` 改成 `2`。

这就是“组内可选、组间隔离”的 `G>1` 模拟。

## 4. 例子（你提到的场景）
### 例 A：`N_code=32`，`N_siso=16`，`G=2`
- code 分两组：每组 16 个（`[0..15]`, `[16..31]`）
- SISO 分两组：每组 8 个预算

执行后效果：
- 每组最多只能保留 8 个 `state==0`；
- 某组超出的 `0` 变 `2`；
- 另一组就算有空闲，也不能借给它。

### 例 B（如果你希望“每组 4 个 SISO”）
要满足 `G=2` 且每组 4 个预算，则总 SISO 应是 `N_siso=8`。  
即：`32 code / 8 SISO / G=2` -> code 每组 16，SISO 每组 4。

## 5. 与当前 `max(G=1)` 的差异
`G=1`（当前）：
- 全局只裁剪一次。
- `N_siso` 可在全体 code 中自由分配。

`G>1`（分组）：
- 每组单独裁剪。
- 会出现“组内不均衡”导致的额外 `state==2`。

因此，同样的 `N_siso` 下，`G` 越大通常 `state==2` 越多（吞吐更容易受限）。

## 6. 可直接落地的伪代码
```cpp
// 输入：state(0/1/2), N_siso, G
void apply_siso_budget_grouped(std::vector<uint8_t>& state, int N_siso, int G) {
  // 1) 计算每组 code 范围：code_groups[g] = [start, end)
  // 2) 计算每组 siso 预算：siso_budget[g]
  // 3) 对每个组：
  //    - 收集组内所有 state==0 的索引 need
  //    - 保留 need 前 siso_budget[g] 个
  //    - 其余 need[k] -> state=2
}
```

## 7. 建议的仿真输出指标
为了看 `G` 的影响，建议每个 tile 记录：
- `n0_need_siso`、`n1_early_stopped`、`n2_unscheduled`
- 每组的 `need_count`、`budget`、`unscheduled_count`

重点看：
- 总 `n2` 随 `G` 的变化；
- 不同组之间 `unscheduled_count` 是否明显不均衡。

## 8. 建议的参数约束
第一版可先加简单约束，便于分析：
1. `G >= 1`
2. `G <= N_code`
3. `N_code`、`N_siso` 不可为负
4. 若不整除，采用“前几组 +1”的余数分配规则（保证总量守恒）

---

这份文档描述的是“仿真规则层”的改动，不涉及硬件电路实现细节。核心目标是：在软件里先量化 `G` 带来的调度受限影响（主要体现在 `state==2` 的增长和组间不均衡）。
