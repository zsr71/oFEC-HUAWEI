# ofec_sweep CSV 画早停率说明

本文说明 `apps/ofec_sweep.cpp` 生成的 `ofec_sweep_results_*.csv`，主要给 MATLAB 画图使用。

## 文件位置

运行：

```bash
./build/apps/ofec_sweep
```

输出文件通常位于：

```text
data/ofec_sweep_results_<run_id>.csv
```

每一行对应一个 sweep scenario，也就是某个 Eb/N0、seed、decoder 参数、alpha/beta、early-stop 配置组合下的一次 pipeline 结果。

## 最常用画图列

画 Eb/N0 下每个 tile 的早停率曲线时，优先使用这些列：

```text
ebn0_db
early_stop_bind_group_size
early_stop_condition_mode
early_stop_action_mode
enable_early_stop
early_stop_row_list
early_stop_row_mean_pct
scenario
run_id
```

其中最关键的是：

```text
ebn0_db
early_stop_row_list
early_stop_bind_group_size
```

`early_stop_row_list` 是每个 tile 的 row-level early-stop percentage，单位是百分比。例如：

```text
"12.3|45.6|78.9|90.1"
```

表示：

```text
tile0 = 12.3%
tile1 = 45.6%
tile2 = 78.9%
tile3 = 90.1%
```

如果 `early_stop_bind_group_size = 4`，并且 `early_stop_condition_mode = 1`，那么这里的 `early_stop_row_list` 已经是 `apply_group_bound_early_stop` 生效之后的早停率。也就是说，它不是 raw early-stop rate，而是“四个 row 绑在一起”之后的 bound early-stop rate。

## early_stop_row_list 和 early_stop_list 的区别

画每个 tile 的早停率曲线时，通常应该使用：

```text
early_stop_row_list
```

不要优先使用：

```text
early_stop_list
```

二者含义不同：

```text
early_stop_row_list
```

每个 tile 里，row 级别 early-stop 的比例。这个是你要画“各 tile 早停率随 Eb/N0 变化”时最直接的列。

```text
early_stop_list
```

每个 tile 里，“整个 tile 的所有 row 都 early-stop”的比例。这个条件更严格，数值通常会明显不同，不等价于 row-level 早停率。

对应的均值列：

```text
early_stop_row_mean_pct
```

是 `early_stop_row_list` 对所有 tile 求平均后的百分比。

```text
early_stop_mean_pct
```

是 `early_stop_list` 对所有 tile 求平均后的百分比。

## early-stop 配置列

这些列用于判断当前行到底是哪一种 early-stop 实验：

```text
enable_early_stop
early_stop_condition_mode
early_stop_action_mode
early_stop_bind_group_size
early_stop_cond_v1_require_bch
early_stop_cond_v1_require_overall
early_stop_v2_llr_abs_threshold
early_stop_v2_max_unreliable_bits
early_stop_cond_v2_include_overall
```

常见配置含义：

```text
enable_early_stop = 1
```

表示 early-stop 总开关开启。

```text
early_stop_condition_mode = 1
```

表示使用条件 1，也就是 BCH syndrome / overall parity 相关的早停检测。

```text
early_stop_bind_group_size = 4
```

表示条件 1 下启用 4-row group binding：同一组 4 个 row 都通过 raw early-stop 检测时，这 4 个 row 才会在 bound 结果中 early-stop。

注意：`early_stop_bind_group_size` 主要对 `early_stop_condition_mode = 1` 有意义。当前代码中 group-bound 逻辑的触发条件是：

```text
ENABLE_EARLY_STOP = true
EARLY_STOP_CONDITION_MODE = 1
EARLY_STOP_BIND_GROUP_SIZE > 1
```

## BER 相关列

如果想顺便看这一行的 BER，可以使用：

```text
pre_ber
pre_errs
pre_total
post_ber
post_errs
post_total
```

其中：

```text
pre_ber
```

译码前 BER。

```text
post_ber
```

译码后 BER。

不过 `ofec_sweep` 通常更适合快速看 early-stop 率和参数趋势。如果要跑低 BER 曲线，建议使用 `ofec_sweep3`。

## alpha/beta 和 tile 参数列

这些列也是 `|` 分隔的 per-tile 列表：

```text
ALPHA_LIST
beta_list
early_stop_beta_list
siso_active_list
```

例如：

```text
"0.342857|0.387439|0.435806|0.485714"
```

表示 tile0 到 tile3 的 alpha。

```text
siso_active_list
```

表示每个 tile 的 SISO 行预算。

这些列主要用于确认不同曲线是不是同一套配置，不一定每次画图都需要。

## unscheduled 相关列

如果要分析 MUX / SISO 预算导致的未调度情况，可以看：

```text
unscheduled_count_list
unscheduled_list
unscheduled_need_list
unscheduled_mean_pct
unscheduled_need_mean_pct
```

含义：

```text
unscheduled_count_list
```

每个 tile 中未被调度的 row 数量列表。

```text
unscheduled_list
```

每个 tile 的未调度比例，单位是百分比。

```text
unscheduled_need_list
```

在“需要 SISO 的 row”里面，每个 tile 未被调度的比例，单位是百分比。

## MATLAB 读取和解析示例

下面示例读取 `ofec_sweep_results_*.csv`，筛选 `early_stop_bind_group_size = 4` 的行，然后画每个 tile 的 `early_stop_row_list` 曲线。

```matlab
csvPath = "data/ofec_sweep_results_xxxxx.csv";
T = readtable(csvPath, "TextType", "string");

% 按需要筛选实验条件
idx = T.enable_early_stop == 1 & ...
      T.early_stop_condition_mode == 1 & ...
      T.early_stop_bind_group_size == 4;

S = T(idx, :);

% 按 Eb/N0 排序，保证曲线横轴有序
S = sortrows(S, "ebn0_db");

% 将 "12.3|45.6|78.9|90.1" 转成数值矩阵
n = height(S);
firstList = split(S.early_stop_row_list(1), "|");
numTiles = numel(firstList);
earlyStopPct = zeros(n, numTiles);

for i = 1:n
    parts = split(S.early_stop_row_list(i), "|");
    earlyStopPct(i, :) = str2double(parts).';
end

figure;
plot(S.ebn0_db, earlyStopPct, "-o", "LineWidth", 1.5);
grid on;
xlabel("Eb/N0 (dB)");
ylabel("Row early-stop rate after group binding (%)");
legend(compose("tile%d", 0:numTiles-1), "Location", "best");
title("Per-tile group-bound early-stop rate vs Eb/N0");
```

如果同一个 Eb/N0 下有多个 seed 或多个 scenario，需要先分组平均。一个简单做法是按 `ebn0_db` 聚合：

```matlab
ebVals = unique(S.ebn0_db);
meanPct = zeros(numel(ebVals), numTiles);

for k = 1:numel(ebVals)
    rows = S.ebn0_db == ebVals(k);

    tmp = zeros(sum(rows), numTiles);
    sub = S(rows, :);
    for i = 1:height(sub)
        parts = split(sub.early_stop_row_list(i), "|");
        tmp(i, :) = str2double(parts).';
    end

    meanPct(k, :) = mean(tmp, 1, "omitnan");
end

figure;
plot(ebVals, meanPct, "-o", "LineWidth", 1.5);
grid on;
xlabel("Eb/N0 (dB)");
ylabel("Mean row early-stop rate after group binding (%)");
legend(compose("tile%d", 0:numTiles-1), "Location", "best");
```

## 建议画图过滤条件

为了避免把不同实验混在一张曲线上，建议至少用这些列做过滤：

```text
decoder_name
chase_L
chase_n_test
chase_topk_keep
chase_group_minima_bits
mux_group_g
mux_scheduling_mode
enable_early_stop
early_stop_condition_mode
early_stop_action_mode
early_stop_bind_group_size
ALPHA_LIST
beta_list
early_stop_beta_list
siso_active_list
```

如果只是比较 `bind_group_size = 1` 和 `bind_group_size = 4`，建议保持其他列完全一致，只改变：

```text
early_stop_bind_group_size
```

这样两条曲线才是干净对照。

## 列表字段格式总结

CSV 中这些字段是字符串列表，元素之间用 `|` 分隔：

```text
ALPHA_LIST
beta_list
early_stop_beta_list
early_stop_list
early_stop_row_list
unscheduled_count_list
unscheduled_list
unscheduled_need_list
siso_active_list
```

MATLAB 中可以统一用：

```matlab
parts = split(T.some_list_column(i), "|");
values = str2double(parts).';
```

转成数值向量。

