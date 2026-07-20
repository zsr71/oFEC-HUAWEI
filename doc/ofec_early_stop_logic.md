# oFEC 当前早停逻辑说明

本文只说明当前代码中的 oFEC early-stop 行为，不讨论理想算法或硬件目标行为。

## 1. 早停判定入口（按 tile、按行）

在每个 tile 进入核心解码前，会先计算 `TileEarlyStopResult`：

- 位置：`src/rx/ofec/detail/ofec_tile_impl.ipp`
- 调用：`tile_early_stop_stats1(prep.lin_matrix)`

`tile_early_stop_stats1()` 的判定规则：

1. 对 `lin_matrix` 每一行做硬判决（`LLR < 0 -> 1`，否则 `0`）。
2. 检查该行 `BCH(255,239)` syndrome 是否为 0。
3. 检查 overall parity 是否一致。
4. 得到：
- `row_passed_flags[r]`：该行是否通过
- `rows_passed / rows_total`
- `all_rows_passed`

实现位置：`src/rx/ofec/earlystop/tile_early_stop_stats.ipp`。

## 2. 真正的“早停动作”发生在行分支

当前实现并不会因为 `all_rows_passed=true` 直接跳过整个 tile 解码；  
而是把 `row_passed_flags` 传给行解码器，在每行处做分支：

- 位置：`src/rx/ofec/ofec_row_decoder_core.cpp`

分支逻辑：

- 若 `early_stop_row_flags[row] == true`：
  - 不调用 Chase
  - 调用 `row_early_stop_process_2(LinVec, LchVec, Y2, row_params)`
- 否则：
  - 正常调用 `chase_fn(...)`

也就是说，当前 oFEC 的 early-stop 是**行级分流**，不是**tile 级硬停止**。

## 3. 早停行的外信息构造方式（当前使用 process_2）

当前代码使用的是 `row_early_stop_process_2`（`process_1` 在调用处被注释）：

- 位置：`src/rx/ofec/earlystop/row_early_stop_process_2.ipp`
- 公式：`y2[j] = (lin[j] - lch[j]) / ALPHA`

之后在 tile decode 流程里，对所有 `produced_rows` 会继续执行：

1. 可选归一化（`normalize_extrinsic_lout`，受 `normalize_extrinsic` 开关控制）
2. 统一乘 `ALPHA`
3. 量化回目标 LLR 精度

位置：`src/rx/ofec/detail/ofec_tile_decode.ipp`。

## 4. 统计口径（对外看到的 early-stop 百分比）

每个 tile 处理后更新两个计数：

- `triggered / total`：tile 级命中（`all_rows_passed`）
- `row_triggered / row_total`：行级命中（`row_passed_flags` 累计）

计数更新位置：`src/rx/ofec/detail/ofec_window_impl.ipp`。  
最终转换为百分比并返回 pipeline 结果：

- `tile_early_stop_pct`
- `tile_row_early_stop_pct`

位置：`src/common/pipeline/pipeline_runner.cpp`。

## 5. 一句话总结

当前 oFEC early-stop 的本质是：

- 先用 BCH+overall parity 做**行级可早停判定**
- 再在行解码处做**早停行/正常行二选一分支**
- tile 级 `all_rows_passed` 目前主要用于统计，不直接短路 tile 解码流程
