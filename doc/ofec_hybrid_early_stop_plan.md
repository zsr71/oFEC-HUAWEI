# oFEC 软硬结合早停策略改造说明

本文说明如何把当前的 `detect_tile_early_stop_v2` 改造成一种“基于 LLR 可靠性的早停判据”：

- 在每次进入行级软解码前，统计当前 256 维 LLR 中“不可靠”比特的数量
- 当不可靠比特数量降低到阈值以下时，认为继续做软迭代的收益已经较小
- 第一阶段只修改“早停判据”，不修改早停命中后的后处理路径

这份文档只说明推荐改法，不直接修改算法细节。

## 1. 当前代码现状

当前主流程里有两类“早停相关”的机制：

1. `detect_tile_early_stop_v1 / v2`
   - 位置：
     - `src/rx/ofec/earlystop/tile_early_stop_stats.ipp`
     - `src/rx/ofec/earlystop/tile_early_stop_stats2.ipp`
   - 当前 `v1` 和 `v2` 实现内容相同
   - 它们都在做“硬判后 BCH + overall parity 检查”

2. `row_early_stop_process_1 / 2`
   - 位置：
     - `src/rx/ofec/earlystop/row_early_stop_process_1.ipp`
     - `src/rx/ofec/earlystop/row_early_stop_process_2.ipp`
   - 它们并不是硬判代数译码，而是“早停命中后如何构造外信息”

另外，当前主流程已经支持在 app 层选择 `v1 / v2`：

- `apps/ofec_single.cpp`
- `apps/ofec_sweep.cpp`
- `src/rx/ofec/detail/ofec_tile_impl.ipp`

但注意：

- 当前 `detect_tile_early_stop_v2` 还不是“不可靠比特阈值判据”
- 当前“早停后处理”默认走的是 `row_early_stop_process_1(...)`
- 本文当前版本只准备改 `v2` 的判据，不改后处理方式

## 2. 目标行为

希望新增的 `v2` 行为是：

1. 对每一行的 `LinVec[0..255]` 统计不可靠比特数
2. 不可靠比特定义为：
   - `abs(LLR) < unreliable_llr_threshold`
3. 若某行满足：
   - `unreliable_count <= max_unreliable_bits`
4. 则认为：
   - 继续做软迭代的收益已经较小
   - 本行可以提前命中早停判据
   - 第一阶段仍沿用当前早停后的外信息路径

第一阶段只把 `v2` 改造成“LLR 可靠性驱动的早停判据”，不改当前早停后的输出构造方式。

也就是说，新版 `v2` 的语义是：

- 这行虽然不一定已经通过 BCH 校验
- 但它的 LLR 已经足够可靠，可以直接按“已早停”路径处理

这和当前 `v1` 的“已经满足 BCH 校验，所以直接早停”不同。

## 3. 关键设计决定

这里最重要的一点是：

**第一阶段可以继续复用当前 `row_passed_flags=true` 的语义。**

虽然从严格语义上讲：

- `v1` 的 `row_passed=true` 表示“这行已经满足 BCH + overall parity”
- 新版 `v2` 的 `row_passed=true` 更像“这行的 LLR 已经足够可靠，可以提前停”

但如果当前阶段只改判据、不改后处理，那么继续复用现有结构是最小改法：

- `row_passed_flags[r] = true`
- `build_state_from_early_stop(...)` 继续把它映射成 `StateTag::EarlyStopped`
- row decoder 继续调用 `row_early_stop_process_1(...)`

这样可以把改动范围限制在：

- `detect_tile_early_stop_v2(...)`
- 新增阈值参数
- 入口参数配置

而不用同时改：

- `TileEarlyStopResult`
- `mux_state` 状态定义
- row decoder 的分支逻辑

## 4. `detect_tile_early_stop_v2` 应该怎么改

### 4.1 新增可配置参数

建议在 `Params` 里增加至少两个参数：

```cpp
float EARLY_STOP_V2_LLR_ABS_THRESHOLD = 0.5f;
int   EARLY_STOP_V2_MAX_UNRELIABLE_BITS = 8;
```

同时在：

- `include/newcode/ofec_single_runner.hpp`
- `include/newcode/ofec_sweep_runner.hpp`
- `apps/ofec_single.cpp`
- `apps/ofec_sweep.cpp`

增加入口配置和赋值，和现在 `kEarlyStopDetectMode` 的接法一致。

### 4.2 `v2` 的核心判据

把 `detect_tile_early_stop_v2(...)` 改成如下逻辑：

对每一行：

1. 遍历 `lin_matrix[r][0..255]`
2. 统计

```cpp
unreliable_count = number of j such that abs(llr(j)) < threshold
```

3. 若 `unreliable_count <= max_unreliable_bits`
   - 标记这一行为“通过早停判据”
   - 即 `row_passed_flags[r] = true`
4. 否则
   - 保持 `row_passed_flags[r] = false`

推荐第一版不要把 BCH syndrome 校验混进 `v2`，保持语义清晰：

- `v1`：码字合法性驱动的 early stop
- `v2`：LLR 可靠性驱动的 early stop 判据

### 4.3 `all_rows_passed` 的定义

对 `v2` 来说，第一阶段仍然保持和当前接口一致：

```cpp
all_rows_passed = all rows satisfy unreliable_count <= max_unreliable_bits
```

也就是：

- 只要每一行都满足不可靠比特阈值条件
- 就认为这整个 tile 已经“判据命中”

## 5. 第一阶段不修改 row decoder 后处理

第一阶段不改以下逻辑：

- `src/rx/ofec/ofec_row_decoder_core.cpp`

也就是说，`mux_tag == 1` 仍然继续走：

```cpp
row_early_stop_process_1(...)
```

不会改成 `perform_hard_decode(...)`。

这样做的目的，是先单独验证：

- 仅替换早停判据后
- 对 early-stop 命中率、调度比例、误码性能会产生什么影响

等这个阶段验证稳定后，再决定是否进入第二阶段：

- 把“早停后处理”从 `row_early_stop_process_1(...)`
- 改成真正的 HIHO 收尾

## 6. 第二阶段（后续可选）：再接 HIHO 收尾

如果后续要继续实现你最初描述的“软硬结合”完整策略，可以再单独做第二阶段改造：

- 扩展 `TileEarlyStopResult`
- 扩展 `mux_state`
- 在 row decoder 里增加专门的 `HardDecodeFinish` 分支
- 复用 `perform_hard_decode(...)`

但这些内容不属于当前阶段的最小改造范围。

## 7. 需要修改的文件

### 7.1 参数与入口

- `include/newcode/params.hpp`
- `include/newcode/ofec_single_runner.hpp`
- `include/newcode/ofec_sweep_runner.hpp`
- `src/ofec_single/ofec_single_params.cpp`
- `src/ofec_sweep/ofec_sweep_runner.cpp`
- `apps/ofec_single.cpp`
- `apps/ofec_sweep.cpp`

### 7.2 `v2` 判据实现

- `src/rx/ofec/earlystop/tile_early_stop_stats2.ipp`

### 7.3 tile 主流程

- `src/rx/ofec/detail/ofec_tile_impl.ipp`

## 8. 推荐的最小落地顺序

1. 先把 `Params` 和 app 层参数接好
2. 把 `detect_tile_early_stop_v2(...)` 改成“不可靠比特计数判据”
3. 保持 `TileEarlyStopResult` 和 `build_state_from_early_stop(...)` 不变
4. 保持 row decoder 后处理路径不变
5. 最后补统计和日志

## 9. 最小验收用例

建议至少验证以下场景：

1. `detect_mode=1`
   - 行为必须与当前 `v1` 完全一致

2. `detect_mode=2`，阈值很小
   - 基本不触发早停
   - 行为接近当前“全走软解码”

3. `detect_mode=2`，阈值较大 / `max_unreliable_bits` 较宽松
   - 部分行会更早地被标记成 `EarlyStopped`
   - row decoder 仍然走 `row_early_stop_process_1(...)`

4. `mux_state` 统计
   - 第一阶段仍然只会出现 `0/1/2`
   - 不引入新的状态码

## 10. 一句话总结

当前阶段的目标是：

- 只把 `detect_tile_early_stop_v2(...)` 改成“不可靠比特计数驱动的早停判据”
- 保持 `TileEarlyStopResult`、`mux_state`、`row_early_stop_process_1(...)` 这条后处理路径不变

也就是说，这一版先只验证“判据改了会怎样”，而不同时引入 HIHO 收尾这第二个变量。
