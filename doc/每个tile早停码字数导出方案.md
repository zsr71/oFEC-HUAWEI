# 每个 tile 早停码字数导出方案

## 目标
记录“每次进入某个 tile 做早停检查后，有多少个 row/码字命中 early-stop”，并把这组原始样本导出成 CSV，方便后续用 Matlab 画分布直方图。

这里要记录的是逐次样本，不是最后的平均值。

也就是说，最终希望拿到类似下面这样的数据：

- tile 0: `3, 5, 0, 7, 2, ...`
- tile 1: `16, 12, 9, 14, ...`
- tile 2: `1, 0, 4, 3, ...`

这样后面就可以直接按 `tile_index` 分组画直方图。

## 为什么单独做这件事
现在已有的 early-stop 输出更偏向汇总统计，例如：

- 每个 tile 的平均 early-stop 百分比
- 每个 tile 的行级 early-stop 百分比

这些适合看整体趋势，但不适合看分布。

如果后面想回答下面这类问题，就需要逐次样本：

- 某个 tile 是“经常中等数量 early-stop”，还是“要么几乎全停、要么几乎不停”？
- 同样的平均值下，分布形状是否不同？
- 新旧 MUX、不同 early-stop 规则下，tile 级 early-stop 数量的分布是否发生了明显变化？

## 推荐记录的最小样本字段
第一版建议每条样本至少记录这些字段：

- `invocation`
- `tile_index`
- `rows_total`
- `rows_passed`

含义如下：

- `invocation`
  - 第几次进入 tile early-stop 检查
  - 用于恢复时间顺序

- `tile_index`
  - 当前是第几个 tile
  - 后续 Matlab 可按 tile 分组

- `rows_total`
  - 这次参与 early-stop 检查的总 row 数
  - 通常会等于当前 tile 的 `rows_to_decode`

- `rows_passed`
  - 这次检查后命中 early-stop 的 row 数
  - 这是画直方图最核心的一列

## 如果以后要扩展，可以再加的字段
第一版不一定要做，但后面如果需要更细分析，可以再加：

- `window_start_row`
- `decoder_name`
- `Eb/N0`
- `bitgen_seed`
- `channel_seed`
- `bind_group_size`
- `mux_scheduling_mode`
- `mux_priority_rule`
- `rows_need_siso_before_mux`
- `rows_unscheduled`

这些字段的价值在于：

- 可以把早停分布和具体实验配置关联起来
- 可以把早停数量分布和后续 MUX 裁剪联系起来分析

## 推荐输出文件形式
建议单独导出一个 CSV，不要混在现有 summary CSV 里。

推荐文件名形式：

- `data/early_stop_hist/ofec_single_tile_early_stop_samples.csv`
- `data/early_stop_hist/ofec_sweep_<run_id>_tile_early_stop_samples.csv`

推荐表头：

```csv
invocation,tile_index,rows_total,rows_passed
```

这样 Matlab 里会很直接，例如：

```matlab
T = readtable("ofec_single_tile_early_stop_samples.csv");
histogram(T.rows_passed(T.tile_index == 0));
```

## 推荐的实现位置
### 1. 在 tile 级拿到单次样本
最自然的位置是：

- `src/rx/ofec/detail/ofec_tile_impl.ipp`

原因是这里已经有：

- `early_stop_stats.rows_passed`
- `early_stop_stats.rows_total`

也就是记录这次样本所需的核心数据已经在这里拿到了。

### 2. 在 window / pipeline 层累计样本
推荐在下面这条链上往上传：

- `ofec_tile_impl`
- `ofec_window_impl`
- `pipeline_runner`

也就是说，每处理完一个 tile，就把这次样本 append 到一个样本列表里。

建议新增一个结构，例如：

- `TileEarlyStopSample`
  - `invocation`
  - `tile_index`
  - `rows_total`
  - `rows_passed`

### 3. 在运行结束时统一导出 CSV
推荐在：

- `ofec_single`
- `ofec_sweep`

结束时统一落盘。

这样有几个好处：

- tile 层只负责采样，不负责 I/O
- 文件输出集中在 runner/summary 层，结构更清楚
- 后面如果同时想支持 `ofec_single` 和 `ofec_sweep`，也更容易复用

## 新增的开关参数
为了避免普通运行时总是多生成文件，加两个控制项：

- `DUMP_TILE_EARLY_STOP_SAMPLES = false`
- `TILE_EARLY_STOP_SAMPLES_OUTPUT_PATH`

建议 `ofec_single` 和 `ofec_sweep` 顶层都暴露。

这样：

- 不需要画图时，保持关闭
- 要分析分布时，再打开导出

## 第一版建议只做什么
第一版建议尽量收敛，先只做：

- 记录每次 tile 检查后的 `rows_passed`
- 记录对应的 `rows_total`
- 导出逐次 CSV

暂时不要一上来就把下面这些都塞进去：

- syndrome 明细
- raw/bound 双版本 rows_passed
- unscheduled 联动统计
- 多种额外标签

因为你这次的主要目标只是：

- 看“每次有多少码字 early-stop”
- 后面画直方图

这个最小版本已经足够支持。

## 后续可以自然扩展的方向
如果第一版跑通，后面可以继续扩展成：

- `raw_rows_passed`
- `bound_rows_passed`
- `rows_need_siso_before_mux`
- `rows_unscheduled`

这样后面还可以分析：

- 组绑定前后 early-stop 数量分布差异
- early-stop 数量和 unscheduled 数量之间的关系
- 不同 MUX 模式下的分布变化

## 推荐实现顺序
建议按下面顺序做：

1. 新增 `TileEarlyStopSample` 结构
2. 在 tile 级采集每次样本
3. 在 pipeline 结果里挂上样本列表
4. 新增 CSV 导出函数
5. 先在 `ofec_single` 接通
6. 再接到 `ofec_sweep`

这样做的好处是：

- 先用 `ofec_single` 验证输出格式
- 再推广到 `ofec_sweep`
- 风险最小，调试最方便
