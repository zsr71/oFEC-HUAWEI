# `G=4` 下重排调度的新方案说明（方案 C）

## 1. 目标
本文档定义一种新的 `G=4, reconfig=true` 调度逻辑，用于替代当前 Python/C++ 中基于“最大匹配 + 增广重排”的方案 B。

为避免命名混淆，本文档把这套新逻辑命名为：

- `方案 C`
- 英文建议命名：`scheme_c_staged`

新的核心思想不是直接在整张允许边图上做一次或两次最大匹配，而是把调度顺序显式拆成两个阶段：

1. 先对每个码字组内的“前 6 个码字”做第一轮调度。
2. 再对每个码字组内的“最后 2 个码字”做第二轮调度。

这样做的目的，是把“组内普通码字优先、组边界特殊码字后调度”的硬件直觉明确写进调度流程中，而不是完全交给图匹配算法自由决定。

## 2. 适用场景
本文档只讨论下面这个固定场景：

- `N_code = 32`
- `N_siso = 16`
- `G = 4`

因此：

- 每组 `8` 个码字
- 每组 `4` 个 SISO
- 一共 `4` 组

码字分组如下：

- 组 0：码字 `0~7`
- 组 1：码字 `8~15`
- 组 2：码字 `16~23`
- 组 3：码字 `24~31`

SISO 分组如下：

- 组 0：SISO `0~3`
- 组 1：SISO `4~7`
- 组 2：SISO `8~11`
- 组 3：SISO `12~15`

## 3. 连接关系定义
### 3.1 组内基础连接
默认情况下，每个码字只连接本组的 4 个 SISO。

例如：

- 码字 `0~7` 默认连接 SISO `0~3`
- 码字 `8~15` 默认连接 SISO `4~7`
- 码字 `16~23` 默认连接 SISO `8~11`
- 码字 `24~31` 默认连接 SISO `12~15`

### 3.2 组间额外旁路线
在基础组内连接之外，额外定义“组边界码字”到其他组 SISO 的旁路线。

这里的“组边界码字”指每组最后两个码字：

- 组 0 的边界码字：`6, 7`
- 组 1 的边界码字：`14, 15`
- 组 2 的边界码字：`22, 23`
- 组 3 的边界码字：`30, 31`

这些码字除了本组 SISO 外，还允许连接预定义的跨组 SISO。具体连哪些边，由额外旁路线配置决定。

### 3.3 当前 Python 中已经定义的旁路线
当前 Python 文件 [group_scheduler.py](/home/zsr71/projects/newcode/python/group_scheduler.py) 中，已经定义了如下旁路线：

```text
(6, 11), (6, 15), (7, 11), (7, 15),
(14, 11), (14, 15), (15, 11), (15, 15),
(22, 3), (22, 7), (23, 3), (23, 7),
(30, 3), (30, 7), (31, 3), (31, 7)
```

按 1-based 写法理解，就是：

- `C7, C8` 额外连接 `S12, S16`
- `C15, C16` 额外连接 `S12, S16`
- `C23, C24` 额外连接 `S4, S8`
- `C31, C32` 额外连接 `S4, S8`

这些边体现的是：

1. 组 0 和组 1 的边界码字，可以跨到组 2 和组 3 的部分 SISO。
2. 组 2 和组 3 的边界码字，可以跨到组 0 和组 1 的部分 SISO。

如果后续要继续扩展方案 C，建议仍然复用这组旁路线定义，不要在不同脚本里各自写一份。

## 4. 新调度逻辑的核心思想
当前你希望的逻辑可以概括为：

1. 每组前 6 个码字先调度。
2. 这一步尽量只消耗本组的 4 个 SISO。
3. 如果前 6 个码字已经把本组 4 个 SISO 占满，那么这一组的第一轮调度结束。
4. 按组顺序把 4 个组都完成第一轮调度。
5. 然后再回过头，依次处理每组最后两个码字。
6. 这最后两个码字优先尝试在本组找剩余 SISO；如果本组没有可用 SISO，再尝试走额外旁路线。

所以，这套方案并不是“所有码字一起做一个全局最大匹配”，而是一个带顺序约束的分阶段分批调度。

## 5. 两个阶段的详细定义
### 5.1 第一阶段：组内普通码字优先调度
对每个码字组，只处理前 6 个码字：

- 组 0：`0~5`
- 组 1：`8~13`
- 组 2：`16~21`
- 组 3：`24~29`

调度规则：

1. 只考虑当前需要解码的码字，即 `state == 0` 的码字。
2. 只允许它们使用本组的 4 个 SISO。
3. 按码字编号顺序调度。
4. 一旦本组 4 个 SISO 全部被占满，这一组第一阶段结束。

这一阶段的目的，是先让组内“普通码字”占用本组资源。

### 5.2 第二阶段：组边界码字后调度
当四个组都完成第一阶段后，再统一处理各组最后两个码字：

- 组 0：`6, 7`
- 组 1：`14, 15`
- 组 2：`22, 23`
- 组 3：`30, 31`

调度顺序建议固定为：

- 先组 0 的 `6, 7`
- 再组 1 的 `14, 15`
- 再组 2 的 `22, 23`
- 再组 3 的 `30, 31`

对每个边界码字，调度规则为：

1. 若该码字不是 `state == 0`，直接跳过。
2. 先尝试本组 SISO 中尚未占用的单元。
3. 若本组没有空闲 SISO，则尝试额外旁路线允许的其他组 SISO。
4. 如果仍找不到可用 SISO，则该码字本轮不调度，最终标成未分配状态。

这一阶段的目标，是把“有跨组能力的特殊码字”放到后面处理，从而避免它们过早占用跨组资源。

## 6. 为什么这样设计
这套顺序化调度和当前最大匹配方案的主要区别是：

1. 当前最大匹配方案是“全局优化”，它只关心最终匹配数。
2. 这里的新方案是“结构优先”，先保证组内普通码字优先使用本组资源。
3. 组边界码字被视为特殊资源申请者，只有在普通码字调度完成后，才使用其更高的连接自由度。

这种做法更接近一种可解释的硬件调度策略：

- 普通码字走普通路由
- 边界码字走补充路由
- 跨组旁路只在后半阶段介入

## 7. 一个具体示例
以组 0 为例：

- 组 0 的普通码字：`0~5`
- 组 0 的边界码字：`6, 7`
- 组 0 的本组 SISO：`0~3`

第一阶段：

1. 先看码字 `0~5` 中哪些需要解码。
2. 按编号顺序把它们分配到 SISO `0~3`。
3. 如果这 4 个 SISO 已被占满，则组 0 第一阶段停止。
4. 码字 `6, 7` 即便当前也需要解码，也先不处理。

第二阶段：

1. 回过头处理码字 `6`。
2. 若组 0 还有空闲 SISO，则优先分配给组 0 的空闲 SISO。
3. 若组 0 没有空闲 SISO，则再检查它的旁路线允许连接的其他 SISO。
4. 如果旁路也没有空闲，则码字 `6` 本轮不调度。
5. 然后再处理码字 `7`。

其他 3 个组同理。

## 8. 需要在实现中固定的细节
如果后续要按本文档改代码，建议把以下细节写死或显式参数化：

1. `G=4` 固定每组 `8` 个码字、`4` 个 SISO。
2. 每组前 `6` 个码字归类为“普通码字”。
3. 每组后 `2` 个码字归类为“边界码字”。
4. 第一阶段只允许普通码字访问本组 SISO。
5. 第二阶段边界码字先尝试本组，再尝试旁路。
6. 第二阶段是否允许边界码字之间再次重排，需要单独决定。

当前按你的描述，更接近“顺序扫描 + 找空闲单元”，而不是再次做全局增广。

## 9. 推荐的数据结构
实现时建议维护以下结构：

- `state[i]`
  当前码字是否需要解码
- `siso_used[j]`
  第 `j` 个 SISO 当前是否已分配
- `code_to_siso[i]`
  第 `i` 个码字最终分配到哪个 SISO，未分配则为 `-1`
- `group_of_code(i)`
  码字所属组
- `local_siso_of_group(g)`
  第 `g` 组本组 SISO 范围
- `bypass_siso_of_code(i)`
  码字 `i` 的额外旁路候选 SISO 列表

## 10. 伪代码
```text
input:
  N_code = 32
  N_siso = 16
  G = 4
  state[0..31]           // 0=need decode, 1=early stop, 2=skip
  bypass_edges           // 额外组间旁路线

output:
  code_to_siso[0..31]    // 未分配则为 -1
  final_state[0..31]     // 成功调度保留 0，未调度改成 2

initialize:
  code_to_siso[i] = -1 for all i
  siso_used[j] = false for all j in [0, 15]

# ---------- Phase 1: 普通码字先调度 ----------
for group g in [0, 1, 2, 3]:
  code_begin = g * 8
  normal_codes = [code_begin + 0, ..., code_begin + 5]
  local_sisos  = [g * 4, ..., g * 4 + 3]

  for code in normal_codes:
    if state[code] != 0:
      continue

    assigned = false
    for siso in local_sisos:
      if siso_used[siso] == false:
        code_to_siso[code] = siso
        siso_used[siso] = true
        assigned = true
        break

    if all local_sisos are used:
      break

# ---------- Phase 2: 边界码字后调度 ----------
for group g in [0, 1, 2, 3]:
  code_begin = g * 8
  tail_codes = [code_begin + 6, code_begin + 7]
  local_sisos = [g * 4, ..., g * 4 + 3]

  for code in tail_codes:
    if state[code] != 0:
      continue

    assigned = false

    # 先尝试本组剩余 SISO
    for siso in local_sisos:
      if siso_used[siso] == false:
        code_to_siso[code] = siso
        siso_used[siso] = true
        assigned = true
        break

    if assigned:
      continue

    # 再尝试额外旁路线
    for siso in bypass_siso_of_code(code, bypass_edges):
      if siso_used[siso] == false:
        code_to_siso[code] = siso
        siso_used[siso] = true
        assigned = true
        break

# ---------- 回写 state ----------
final_state = copy(state)
for code in [0..31]:
  if state[code] != 0:
    continue
  if code_to_siso[code] == -1:
    final_state[code] = 2
  else:
    final_state[code] = 0

return code_to_siso, final_state
```

## 11. 这套方案与当前 `scheme_b_reconfig` 的差异
当前 Python/C++ 中的 `scheme_b_reconfig` 是：

1. 阶段 1：组内最大匹配
2. 阶段 2：开放旁路后继续增广

而本文档的新方案是：

1. 阶段 1：每组前 6 个普通码字按顺序占本组 SISO
2. 阶段 2：每组最后 2 个边界码字按顺序尝试本组剩余 + 旁路线

所以本文档实际上定义的是一种新的“顺序式两阶段调度”，而不是原始的“图匹配式两阶段重排”。

后续如果按本文档改代码，建议新建单独函数，不要直接复用当前 `schedule_scheme_b_reconfig` 的实现名字，以免语义混淆。

## 12. Python 中应如何修改
如果要先在 Python 里落地方案 C，建议改两处文件。

### 12.1 修改 `group_scheduler.py`
在 [group_scheduler.py](/home/zsr71/projects/newcode/python/group_scheduler.py) 中新增：

1. `schedule_scheme_c_staged(...)`
   作用：实现本文档定义的顺序式两阶段调度。

2. `_build_bypass_map(...)`
   作用：把 `EXTRA_BYPASS_EDGES` 转成 `code -> bypass_siso_list` 的映射，方便第二阶段直接查询。

3. `run_scheme_c(...)`
   作用：仿照 `run_scheme_b(...)` 提供一个对外入口，便于画图脚本或统计脚本直接调用。

建议接口形式：

```text
schedule_scheme_c_staged(
    N_code,
    N_siso,
    G,
    active_codes,
    free_siso
) -> (final_match, waiting_codes)
```

其中：

- `final_match`：`code -> siso`
- `waiting_codes`：未调度上的 code 列表

### 12.2 修改 `compare_mux_siso_utilization.py`
在 [compare_mux_siso_utilization.py](/home/zsr71/projects/newcode/python/compare_mux_siso_utilization.py) 中新增：

1. `simulate_scheme_c(...)`
   作用：调用 `schedule_scheme_c_staged(...)`，统计方案 C 的已调度数与未调度数。

2. 在对比表中增加方案 C 的列：
   - `方案C最终利用率`
   - `方案C平均已调度`
   - `方案C平均未调度`

这样脚本就可以同时比较：

1. 静态方案 `G=2, reconfig=false`
2. 静态方案 `G=4, reconfig=false`
3. 方案 B `G=4, reconfig=true`
4. 方案 C `G=4, reconfig=true`

## 13. 方案 C 的建议命名
为保证后续文档、脚本和 C++ 实现一致，建议统一命名如下：

- 文档名：`方案 C`
- Python 函数：`schedule_scheme_c_staged`
- Python 运行入口：`run_scheme_c`
- C++ 未来函数名：`schedule_scheme_c_staged_cpp`

这样可以和当前已有的方案 B 明确区分：

- `scheme_b_reconfig`：图匹配式两阶段重排
- `scheme_c_staged`：顺序式两阶段调度
