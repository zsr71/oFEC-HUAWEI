# oFEC 增加 MUX 逻辑的代码改动清单与功能需求

本文给出一份面向实现的改造清单，目标是在现有软件中引入 MUX 约束仿真（`full` / `grouped`），并支持“未分配任务”策略（`drop` / `earlystop_fallback`）。

## 1. 改造目标

1. 在不改 oFEC 基本译码算法（Chase/early-stop 判定）的前提下，引入“连接受限”的调度仿真。
2. 支持两种连接拓扑：
   - `full`：任意 code 可分配给任意可用 SISO。
   - `grouped`：code 仅能分配给组内 SISO。
3. 支持未分配任务处理策略：
   - `drop`：本轮不更新该任务。
   - `earlystop_fallback`：走早停回退路径。
4. 在 pipeline/sweep 输出可分析指标（利用率、未分配比例、回退次数、BER/FER/吞吐）。

## 2. 需要修改的代码位置

## 2.1 参数层（配置入口）

文件：
- `include/newcode/params.hpp`
- 参数透传相关入口（按项目实际可能在 `apps/*` 与 `src/ofec_sweep/*`）

需要新增：
- `MUX_ENABLE`（bool）
- `MUX_MODE`（enum/string：`full`/`grouped`）
- `MUX_GROUP_COUNT`（int）
- `MUX_ACTIVE_SISO`（int）
- `MUX_UNSERVED_POLICY`（enum/string：`drop`/`earlystop_fallback`）
- 可选：`MUX_SCHED_POLICY`（`round_robin`/`oldest_first`）

## 2.2 调度核心（MUX 约束生效点）

文件：
- `src/rx/ofec/detail/ofec_window_impl.ipp`

需要改动：
1. 在 tile 处理循环中增加“任务集 -> 分配集”的调度步骤。
2. 根据 `MUX_MODE` 构建可连接集合：
   - `full`：全可达
   - `grouped`：仅组内可达
3. 根据 `MUX_ACTIVE_SISO` 截断本轮可服务任务数量。
4. 输出本轮“被分配任务掩码”与“未分配任务掩码”。
5. 未分配任务按 `MUX_UNSERVED_POLICY` 处理，不进入等待队列模型。

## 2.3 tile/row 执行链路（把调度结果传下去）

文件：
- `src/rx/ofec/detail/ofec_tile_impl.ipp`
- `src/rx/ofec/detail/ofec_tile_decode.ipp`
- `src/rx/ofec/ofec_row_decoder_core.cpp`
- `include/newcode/rx/ofec/chase/decoder_core.hpp`

需要改动：
1. 在 core 接口增加调度输入参数（例如 `scheduled_row_mask` 或 `scheduled_code_indices`）。
2. 行处理逻辑改为三分支：
   - 已分配：按现有逻辑（early-stop 行走 `row_early_stop_process_2`，其余跑 Chase）
   - 未分配 + `drop`：保持原值，不更新
   - 未分配 + `earlystop_fallback`：执行早停回退路径（建议复用 `row_early_stop_process_2`）
3. 保持现有 `early_stop_row_flags` 判定和统计路径不破坏。

## 2.4 统计结构与结果回传

文件：
- `include/newcode/ofec_decoder.hpp`
- `include/newcode/decoder_api.hpp`
- `include/newcode/pipeline_runner.hpp`
- `src/common/pipeline/pipeline_runner.cpp`

需要新增统计字段（建议）：
- `mux_sched_utilization`
- `mux_unscheduled_count`
- `mux_drop_count`
- `mux_earlystop_fallback_count`
- `mux_group_load_variance`
- `mux_avg_fanin`（理论值）

需要改动：
1. 在 decode 过程中累计 MUX 指标。
2. 在 pipeline 结果中暴露新指标。
3. 保持原 `tile_early_stop_pct`、`tile_row_early_stop_pct` 兼容。

## 2.5 扫描实验与导出

文件：
- `src/ofec_sweep/ofec_sweep_runner.cpp`
- `src/ofec_sweep/ofec_sweep_io.cpp`

需要改动：
1. sweep 参数增加 `MUX_GROUP_COUNT`、`MUX_MODE`、`MUX_UNSERVED_POLICY`。
2. CSV 增加 MUX 指标列，至少包括：
   - `mux_mode`
   - `group_count`
   - `utilization`
   - `unscheduled_count`
   - `drop_count`
   - `earlystop_fallback_count`
   - BER/FER 与吞吐字段

## 3. 需要新增的功能模块

建议新增以下轻量模块，避免逻辑散落：

1. `MuxConfig`：集中保存 MUX 配置参数。
2. `MuxScheduler`：输入任务集合，输出分配结果（已分配/未分配）。
3. `MuxStats`：统一计数与汇总接口。
4. `apply_unscheduled_policy(...)`：对未分配任务执行 `drop` 或 `earlystop_fallback`。

建议放置：
- 头文件：`include/newcode/ofec/mux/*`
- 实现：`src/rx/ofec/mux/*`

## 4. 推荐实现顺序

1. 先接参数与统计结构（不改变行为，默认 `MUX_ENABLE=false`）。
2. 实现 `MUX_MODE=full`，验证与当前结果一致（回归基线）。
3. 接入 `grouped` 分配约束。
4. 接入 `drop`/`earlystop_fallback` 未分配策略。
5. 接入 sweep 导出并做 `G` 扫描。

## 5. 验证点（最小必做）

1. 功能回归：
   - `MUX_ENABLE=false` 时结果与现网一致。
   - `MUX_MODE=full` 时与无 MUX 约束近似一致。
2. 约束生效：
   - `grouped` 模式下不同组不可交叉分配。
3. 策略正确：
   - `drop` 与 `earlystop_fallback` 的计数与行为一致。
4. 指标可用：
   - pipeline 与 sweep 输出含新增字段，且数值随 `G` 有合理变化趋势。

## 6. 交付标准

1. 编译通过，核心入口可运行。
2. 新增参数有默认值，老脚本不需要改也能跑。
3. 文档与 CSV 字段一致。
4. 至少一组 `full vs grouped` 对比结果可复现。

