# Matlab 对齐任务说明

目标是把 `ofec_ber_window_probe` 导出的两类 CSV 在 Matlab 里做对齐分析，重点看：

- 最后一个 tile 的 `need_decode = rows_total - rows_passed` 在时间上的峰值变化
- 这个峰值变化和 BER 是否相关

## 现有 CSV

当前主要看这两个文件：

- `build/data/ber_window_probe/per_seed_tile_early_stop_samples.csv`
- `build/data/ber_window_probe/per_seed_per_tile.csv`

### 1. `per_seed_tile_early_stop_samples.csv`

这是时域样本表，每行对应：

- `ebn0_db`
- `seed_index`
- `bitgen_seed`
- `channel_seed`
- `invocation`
- `tile_index`
- `rows_total`
- `rows_passed`

这里最关键的派生量是：

- `need_decode = rows_total - rows_passed`

它表示当前这次 tile early-stop 检查后，仍然需要解码的 code 数。

### 2. `per_seed_per_tile.csv`

这是 tile 级 BER 表，每行对应：

- `ebn0_db`
- `seed_index`
- `bitgen_seed`
- `channel_seed`
- `tile_idx`
- `pre_errs`
- `pre_bits`
- `pre_ber`
- `post_errs`
- `post_bits`
- `post_ber`

这份表是 tile 汇总结果，不是 invocation 级时序数据。

## 重点说明

不要直接把 `per_seed_tile_early_stop_samples.csv` 和 `per_seed_per_tile.csv` 按 `invocation` 一一硬对齐。

原因是：

- `per_seed_tile_early_stop_samples.csv` 是 invocation 级时序样本
- `per_seed_per_tile.csv` 是 tile 汇总结果

它们粒度不同，不能直接逐点配对。

## 推荐的 Matlab 处理流程

### 第一步：画时域曲线

对 `per_seed_tile_early_stop_samples.csv`：

1. 先按 `(ebn0_db, seed_index, tile_index)` 分组
2. 再按 `invocation` 排序
3. 计算：
   - `need_decode = rows_total - rows_passed`
4. 画出 `need_decode` 随 `invocation` 的曲线

这一步的目标是找出：

- `need_decode` 是否会出现峰值
- 峰值持续多久
- 峰值前后是否有明显变化

### 第二步：看 tile 级 BER

对 `per_seed_per_tile.csv`：

1. 按 `(ebn0_db, seed_index, tile_idx)` 分组
2. 直接看 `post_ber`
3. 画出 tile 级 BER 分布或曲线

这一步的目标是看：

- 某个 tile 的整体 BER 水平
- tile 之间 BER 是否有明显差异

### 第三步：做特征压缩后再对齐

如果要分析 `need_decode` 峰值和 BER 的关系，建议不要用原始 invocation 逐点对齐，而是先把时域曲线压缩成 tile 级特征，再和 `post_ber` 做 join。

建议的特征包括：

- `mean_need_decode`
- `peak_need_decode`
- `first_N_invocations_mean_need_decode`
- `area_under_curve`

然后构造一个新的分析表，每行一个：

- `(ebn0_db, seed_index, tile_idx)`

字段可以包括：

- `mean_need_decode`
- `peak_need_decode`
- `post_ber`
- `pre_ber`

这样就可以直接做散点图或者相关性分析。

## 建议的 Matlab 输出图

建议至少画三类图：

1. `need_decode` 的时域曲线
2. `post_ber` 的 tile 级分布图
3. `need_decode` 特征和 `post_ber` 的散点图

如果 `need_decode` 的峰值和 `post_ber` 明显相关，就说明这个时域异常段确实影响译码性能。

## 如果要进一步增强分析

如果后面想做严格的“随时间变化的 tile BER 曲线”，当前这两份 CSV 还不够，因为 `per_seed_per_tile.csv` 没有 invocation 级时序信息。

那时需要再让 probe 输出一份：

- `tile + invocation` 粒度的 BER 样本

但这一步不是第一版必须做的。

## 给 Codex 的执行提示

如果你要把这份说明转成任务，建议直接要求：

1. 读取 `per_seed_tile_early_stop_samples.csv`
2. 计算 `need_decode = rows_total - rows_passed`
3. 读取 `per_seed_per_tile.csv`
4. 按 `(ebn0_db, seed_index, tile_idx)` 做汇总
5. 生成：
   - 时域曲线
   - tile BER 曲线
   - 特征相关性散点图

不要一开始就要求它改主程序逻辑。第一版只做 Matlab 分析就够了。
