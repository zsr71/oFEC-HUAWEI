# oFEC 早停场景下的分组 MUX 实现说明

## 1. 背景与核心问题

在面向低功耗、低面积的 oFEC/译码器架构中，引入早停后通常会减少同一时刻启用的 SISO 数量。  
这会打破 `code <-> SISO` 的静态一一映射，原先“近似直连”的连接方式变为“动态分配”问题。

动态分配在硬件上的直接代价是选择网络：

- 每个 SISO 需要在多个 code 之间做选择；
- 需要在 SISO 输入侧（或 code 输出侧）增加 MUX/交叉开关；
- 连接自由度越高，MUX 扇入越大、布线越复杂、面积和功耗越高。

因此，早停带来的算力节省并不会自动等价为面积节省，连接网络可能成为新的主开销。

## 2. 无早停与早停下的连接差异

### 2.1 无早停（或固定并行度）场景

- 典型做法：每个 SISO 绑定固定 code 位置；
- 连接关系稳定，接近直连；
- 选择逻辑极少，MUX 成本可忽略。

### 2.2 引入早停后的动态场景

以 `32 code`、`有效启用 16 个 SISO` 为例：

- 每个 SISO 不再固定服务单个 code；
- SISO 需要从多个候选 code 中选择当前译码对象；
- 连接关系从静态映射变为动态调度；
- MUX 网络成为必需模块。

## 3. 全相连方案与其代价

若采用全相连（Fully Connected）拓扑：

- 每个 SISO 可访问全部 code；
- 调度自由度最高，负载均衡能力最强；
- 但每个 SISO 需要大扇入选择器（例：`32:1`）。

在一阶近似下，MUX 网络规模与扇入成正相关。  
当 code 数增大时，面积、布线拥塞和动态功耗会快速上升，时序也更难收敛。

## 4. 方案一：组相连（Grouped Connectivity）

### 4.1 基本思想

将“全局任意选择”约束为“组内选择”：

1. 把 `C_total` 个 code 划分为 `G` 组；
2. 每组仅配置本组可用的 SISO；
3. 每个 SISO 只能在组内 code 中选择；
4. 组与组之间不互连（不跨组选 code）。

### 4.2 直观效果

- 每个 MUX 的输入路数显著下降；
- 连线范围收敛到组内，布线更局部；
- 负载和时序压力降低；
- 代价是跨组资源复用能力下降。

### 4.3 示例

`32 code` 划分为 `4` 组，每组 `8 code`：

- 全相连：每个 SISO 需要 `32:1` 选择器；
- 组相连：每个 SISO 仅需 `8:1` 选择器。

即 MUX 扇入从 32 降到 8，选择网络规模明显缩小，面积和互连复杂度同步下降。

## 5. 分组数 G 的权衡

`G` 是该方案的关键设计参数。

### 5.1 小 G（组少、组大）

- 更接近全相连；
- 调度灵活、性能风险小；
- 但 MUX 降幅有限，面积收益有限。

### 5.2 大 G（组多、组小）

- MUX 扇入明显下降，面积收益更好；
- 但调度约束更强，可能出现组间负载不均：
  - 某组忙、某组闲；
  - 闲置 SISO 无法跨组支援；
  - 利用率下降，未被服务任务需要按策略处理（丢弃或走早停路径），吞吐稳定性受影响。

因此，不能只追求最大分组数，需要在“面积收益”和“性能退化”之间寻找折中点。

## 6. 推荐仿真评估方法（面向程序实现）

建议围绕分组数 `G` 做参数扫描，至少记录以下指标：

1. **面积代理指标**  
   - 单 SISO MUX 扇入：`C_total / G`  
   - 全芯片 MUX 复杂度代理：`S_active * (C_total / G)` 或门级估算值
2. **利用率指标**  
   - 每组 SISO busy 比例
   - 未获分配任务占比（按组统计）
3. **性能指标**  
   - 吞吐（码块/周期或 bit/s）
   - 平均时延与尾时延（P95/P99）
   - BER/FER（确认分组不会引入不可接受性能损失）
4. **早停相关指标**  
   - 组内早停率分布
   - 组间负载方差（反映“早停不均衡”导致的资源失配）

建议的扫描维度：

- `G = 1, 2, 4, 8, ...`
- 在每个 `G` 下联合扫描 SNR、早停率、code 负载分布。

## 7. 工程结论（用于方案决策）

该方案的目标不是简单替换连接拓扑，而是通过分组数优化得到一个可行点：

- 相比全相连，显著压缩 MUX 面积与布线复杂度；
- 在吞吐、时延与误码性能约束下，性能退化可控；
- 最终为早停后的低功耗、低面积 oFEC 架构提供可实现的连接方案。

## 8. 在现有软件中实现 MUX 仿真需要修改的代码

下面按“参数定义 -> 调度执行 -> 指标回传”的实现路径列出改动点。

### 8.1 参数与配置入口（先把 MUX 模式配置接进来）

建议在 `Params` 中新增 MUX 仿真参数，例如：

- `MUX_ENABLE`：是否启用 MUX 仿真
- `MUX_MODE`：`full` / `grouped`
- `MUX_GROUP_COUNT (G)`：分组数
- `MUX_ACTIVE_SISO`：可并行 SISO 数
- `MUX_SCHED_POLICY`：组内调度策略（round-robin / oldest-first）
- `MUX_UNSERVED_POLICY`：未分配任务处理策略（`drop` / `earlystop_fallback`）

需要修改：

- `include/newcode/params.hpp`
- 若有参数加载/命令行透传，也要同步修改对应入口（`apps/` 或 sweep runner 的参数解析）

### 8.2 主要调度逻辑（MUX 约束真正生效的地方）

当前窗口处理是按 tile 顺序调用 `process_tile_impl()`，要在这里加入“资源受限 + 连接约束”：

- 位置：`src/rx/ofec/detail/ofec_window_impl.ipp`

建议改动：

1. 在处理每个 tile 前，先构建本次待处理 code 列表（或 row 任务列表）。
2. 当 `MUX_MODE=full` 时，允许任意 SISO 选择任意 code（作为基线）。
3. 当 `MUX_MODE=grouped` 时，按 `G` 把 code 划组，只允许组内分配。
4. 加入 `MUX_ACTIVE_SISO` 限制：每个调度步只发放不超过该数量的任务。
5. 对未分配任务按 `MUX_UNSERVED_POLICY` 处理：直接放弃译码，或走早停回退路径。
6. 统计空闲率、未分配比例、回退触发次数。

说明：软件里虽然不是 RTL，但这里可以抽象出“每一步可服务任务数”和“可连接集合”，并结合未分配任务处理策略，等效仿真 MUX 约束对吞吐与 BER/FER 的影响。

### 8.3 tile/row 解码执行路径（把调度结果传递到解码器）

当前 tile 内部直接把 `row_passed_flags` 交给 row core，row core逐行分支执行。
要支持 MUX 仿真，需要把“本步被分配到 SISO 的行/码字”显式传下去。

需要修改：

- `src/rx/ofec/detail/ofec_tile_impl.ipp`
- `src/rx/ofec/detail/ofec_tile_decode.ipp`
- `src/rx/ofec/ofec_row_decoder_core.cpp`
- 相关函数声明：`include/newcode/rx/ofec/chase/decoder_core.hpp`

建议接口新增：

- `scheduled_row_mask` 或 `scheduled_code_indices`

执行语义：

- 被调度到的对象：按现有流程跑 Chase 或 early-stop 分支；
- 未被调度到的对象：按策略二选一
  - `drop`：直接放弃本次译码更新；
  - `earlystop_fallback`：走早停对应流程（例如 `row_early_stop_process_2` 路径）。

### 8.4 结果统计结构（支持 MUX 评估）

当前统计主要是 early-stop 百分比，需新增 MUX 指标字段。

建议修改：

- `include/newcode/ofec_decoder.hpp`（扩展 `TileEarlyStopCounter` 或新增 `MuxStats`）
- `include/newcode/decoder_api.hpp`（`DecodeResult` 增加 MUX 统计）
- `include/newcode/pipeline_runner.hpp`（`PipelineResult` 增加 MUX 指标）
- `src/common/pipeline/pipeline_runner.cpp`（汇总与输出）

建议新增指标：

- `mux_avg_fanin`（理论扇入）
- `mux_sched_utilization`（SISO 利用率）
- `mux_group_load_variance`（组间负载方差）
- `mux_unscheduled_count`（未获分配次数）
- `mux_drop_count`（直接放弃译码次数）
- `mux_earlystop_fallback_count`（走早停回退次数）

### 8.5 Sweep 与结果导出（用于分组数扫描）

为了做 `G` 扫描并导出对比，需要把新指标接入 sweep CSV。

建议修改：

- `src/ofec_sweep/ofec_sweep_runner.cpp`
- `src/ofec_sweep/ofec_sweep_io.cpp`

建议最少导出：

- `group_count`
- `mux_mode`
- `avg_fanin`
- `utilization`
- `unscheduled_count`
- `drop_count`
- `earlystop_fallback_count`
- BER/FER 与吞吐指标

### 8.6 推荐实施顺序

1. 先实现 `full` 模式（逻辑上等价当前行为）并加统计，确保回归一致。  
2. 再实现 `grouped` 模式，只加连接约束，不改解码算法。  
3. 最后接入 sweep 做 `G` 参数扫描，验证面积代理收益与性能退化拐点。
