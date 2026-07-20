# `ofec_two_stream_ber_window_probe` tile BER CSV 说明

本文说明 `apps/ofec_two_stream_ber_window_probe.cpp` 目前新增导出的两份 tile 长度 BER CSV：

1. `per_seed_per_tile_ber.csv`
2. `aggregated_tile_ber.csv`

这两份表的目标不是替代 shared 内部调试三表：

- `per_invocation_shared_tile_samples.csv`
- `per_invocation_shared_hybrid_class_counts.csv`
- `per_invocation_shared_row_map.csv`

而是补一层“按 tile 位置观察 BER 分布”的视角，用来回答：

- Stream A 和 Stream B 在不同 `tile_index` 上的 BER 是否不同
- 某些 tile 是否天然更难，或者更容易出 shared 不对称
- 某个 tile 的 BER 异常，是否能和 shared 内部的 `early-stop / hybrid / MUX / unscheduled` 异常对应起来

---

## 1. 两份表分别看什么

### 1.1 `per_seed_per_tile_ber.csv`

这是逐 seed、逐 tile 的原始 BER 表。

一行代表：

- 某次具体运行中的一个 `seed`
- 某一路 `stream_label`
- 某个 `tile_index`
- 在这个 tile 对应的比特段上，统计出来的 `pre-FEC BER` 和 `post-FEC BER`

这个表适合回答的问题是：

- 在同一个 seed 下，A/B 两路在哪些 tile 上开始分叉
- 某个 seed 的 BER 异常，是只集中在几个 tile，还是沿整条 frame 都有问题
- 某个 `tile_index` 的问题，是偶发的，还是多个 seed 都反复出现

### 1.2 `aggregated_tile_ber.csv`

这是按 `Eb/N0 + stream_label + tile_index` 聚合后的 BER 表。

一行代表：

- 同一个 `Eb/N0`
- 同一路 `stream_label`
- 同一个 `tile_index`
- 把该点下所有 seed 的该 tile 误差和比特数累加后得到的 BER

这个表适合回答的问题是：

- 从总体统计上看，A/B 两路在哪些 tile 上 BER 差异最大
- 某个 tile 的异常是不是 seed 偶然性，还是总体结构性问题
- shared 调度导致的问题，是不是集中出现在特定 tile

---

## 2. “tile BER” 到底是怎么统计的

这里的 “tile BER” 不是：

- shared 内部某次 invocation 的 `produced/failed` 比例
- 也不是某个 tile 里 64 个 row 的调度结果占比

这里的 “tile BER” 指的是：

- 把一整条发送信息比特序列与接收信息比特序列
- 按“一个 tile 覆盖的比特长度”做非重叠切段
- 每一段单独统计误码数和 BER

当前底层逻辑沿用仓库里已有的：

- `compute_ber_per_tile_window(...)`

它的段长定义是：

- `tile_bits = p.tile_height_rows() * 111`

也就是说，一个 `tile_index` 对应的是：

- 第 `tile_index` 个 tile 长度的信息比特段

它是“固定长度分段统计”，不是重叠式滑动窗口统计。

---

## 3. A/B 两路是怎么放进同一张表里的

这次的设计是：

- A/B 两路共用同一张 CSV
- 用 `stream_label` 区分当前这行属于哪一路

约定如下：

- `stream_label = A`
  - 表示这一行对应 Stream A
- `stream_label = B`
  - 表示这一行对应 Stream B

这样做的好处是：

- 后处理时可以直接按 `tile_index` 对齐比较 A/B
- 不需要额外 merge 两张不同文件
- 很适合和 shared debug 三表一起做 join 或对照

---

## 4. 第一份表：`per_seed_per_tile_ber.csv`

## 4.1 一行代表什么

每一行代表：

- 一个 `run_id`
- 一个 `Eb/N0`
- 一个 `seed`
- 一路 `stream_label`
- 一个 `tile_index`
- 在这个 tile 长度比特段上的 `pre/post BER`

因此，如果：

- 有 1 个 `Eb/N0`
- 有 1 个 `seed`
- 有 384 个 tile 段
- 有 2 路流 A/B

那么这张表会有：

- `1 * 1 * 384 * 2 = 768` 行数据

加上一行表头，总行数就是 769。

## 4.2 列定义

```text
run_id,ebn0_db,seed_index,chunk_index,bitgen_seed,channel_seed,
stream_label,tile_index,pre_errs,pre_bits,pre_ber,
post_errs,post_bits,post_ber
```

## 4.3 各列含义

### 标识列

`run_id`

- 当前这次 probe 运行的编号
- 用来区分不同批次导出的数据

`ebn0_db`

- 当前这行数据所属的 `Eb/N0`

`seed_index`

- 当前这行数据所属的 seed 编号

`chunk_index`

- 预留给和其他 sweep / chunk 化流程对齐的字段
- 当前 `apps/ofec_two_stream_ber_window_probe.cpp` 中固定写成 `0`

### 种子列

`bitgen_seed`

- 当前这一路生成发送信息比特时使用的随机种子
- 对于 `stream_label=A`，它对应 Stream A 的 bitgen seed
- 对于 `stream_label=B`，它对应 Stream B 的 bitgen seed

`channel_seed`

- 当前这一路生成信道噪声时使用的随机种子
- 对于 `stream_label=A`，它对应 Stream A 的 channel seed
- 对于 `stream_label=B`，它对应 Stream B 的 channel seed

### 路别 / 位置列

`stream_label`

- 当前这行属于哪一路
- 取值约定：
  - `A`
  - `B`

`tile_index`

- 当前统计的是第几个 tile 长度比特段
- 它对应的是 BER 分段编号，不是 shared 内部某次 invocation 的 `tile_index`
- 这一点很重要：
  - shared debug 三表里的 `tile_index` 指的是 decoder 流程里当前 window 内的 tile 位置
  - 本文这两张 BER 表里的 `tile_index` 指的是沿整条 frame 线性切分后的第几个 tile 比特段

### pre-FEC 统计列

`pre_errs`

- 当前 tile 比特段内，pre-FEC 的误码数
- 这里的 pre-FEC 指的是：
  - 用信道 LLR 直接硬判得到的信息比特
  - 和发送端参考信息比特比较后的误差

`pre_bits`

- 当前 tile 比特段内，参与统计的总比特数

`pre_ber`

- 当前 tile 比特段的 pre-FEC BER
- 计算方式为：
  - `pre_errs / pre_bits`

### post-FEC 统计列

`post_errs`

- 当前 tile 比特段内，post-FEC 的误码数
- 这里的 post-FEC 指的是：
  - decoder 最终输出 LLR 经过硬判后的信息比特
  - 与发送端参考信息比特比较后的误差

`post_bits`

- 当前 tile 比特段内，参与统计的总比特数

`post_ber`

- 当前 tile 比特段的 post-FEC BER
- 计算方式为：
  - `post_errs / post_bits`

## 4.4 适合怎么用

这张表最适合做：

- 同一 seed 下 A/B 两路的 `tile_index -> post_ber` 曲线
- 某个 seed 的热点 tile 排序
- 某个 `tile_index` 在不同 seed 上的离散程度分析

最直接的看法是：

- 先固定一个 `seed_index`
- 再比较：
  - `stream_label=A`
  - `stream_label=B`
- 看两路在哪些 `tile_index` 上 `post_ber` 差得最明显

---

## 5. 第二份表：`aggregated_tile_ber.csv`

## 5.1 一行代表什么

每一行代表：

- 一个 `run_id`
- 一个 `Eb/N0`
- 一路 `stream_label`
- 一个 `tile_index`
- 把所有 seed 在这个 tile 上的误码和比特数累加后得到的聚合 BER

它不是简单平均每个 seed 的 BER，而是：

- 先累加 `errs`
- 再累加 `bits`
- 最后用总误码数除以总比特数

这种做法更稳定，也和仓库里其他 sweep / BER 聚合逻辑保持一致。

## 5.2 列定义

```text
run_id,ebn0_db,stream_label,tile_index,
agg_pre_errs,agg_pre_bits,agg_pre_ber,
agg_post_errs,agg_post_bits,agg_post_ber
```

## 5.3 各列含义

### 标识列

`run_id`

- 当前这次 probe 运行的编号

`ebn0_db`

- 当前聚合结果所属的 `Eb/N0`

`stream_label`

- 当前聚合结果属于哪一路
- 取值为：
  - `A`
  - `B`

`tile_index`

- 当前聚合的是第几个 tile 长度比特段

### 聚合后的 pre-FEC 列

`agg_pre_errs`

- 所有 seed 在该 `tile_index` 上的 pre-FEC 总误码数之和

`agg_pre_bits`

- 所有 seed 在该 `tile_index` 上的 pre-FEC 总统计比特数之和

`agg_pre_ber`

- 聚合后的 pre-FEC BER
- 计算方式为：
  - `agg_pre_errs / agg_pre_bits`

### 聚合后的 post-FEC 列

`agg_post_errs`

- 所有 seed 在该 `tile_index` 上的 post-FEC 总误码数之和

`agg_post_bits`

- 所有 seed 在该 `tile_index` 上的 post-FEC 总统计比特数之和

`agg_post_ber`

- 聚合后的 post-FEC BER
- 计算方式为：
  - `agg_post_errs / agg_post_bits`

## 5.4 适合怎么用

这张表最适合做：

- A/B 两路的总体 `tile_index -> agg_post_ber` 曲线
- 排查哪些 tile 在总体上对 A/B 不对称最敏感
- 把 BER 热点 tile 和 shared debug 三表里的：
  - `rows_unscheduled`
  - `class_hard_fail_count`
  - `class_two_main_count`
  - `deferred_candidate_count`
  做对应

最常见的用法是：

- 先看 `aggregated_tile_ber.csv`
- 找到 `A/B` 的 `agg_post_ber` 差距最大的几个 `tile_index`
- 再回到同一批 `run_id` 的 shared debug 三表里，看这些 tile 是否同时伴随：
  - 更高的 `rows_unscheduled`
  - 更高的 `hard_fail`
  - 更高的 `two_main`

---

## 6. 它和 shared 内部调试三表是什么关系

这两份 tile BER 表回答的是：

- “哪个 tile 的 BER 高”
- “A/B 两路哪个 tile 的 BER 差异大”

shared 内部调试三表回答的是：

- “这个 tile 在 decoder 内部到底经历了什么”

可以这样分工理解：

- `per_seed_per_tile_ber.csv`
  - 负责告诉你：
    - 某个 seed 下，问题出现在哪些 tile

- `aggregated_tile_ber.csv`
  - 负责告诉你：
    - 总体上，哪些 tile 最容易有问题

- `per_invocation_shared_tile_samples.csv`
  - 负责告诉你：
    - 这个 invocation / tile 上，early-stop / hard-finish / MUX 的数量关系是什么

- `per_invocation_shared_hybrid_class_counts.csv`
  - 负责告诉你：
    - 这个 invocation / tile 上，hybrid 分类结构是什么

- `per_invocation_shared_row_map.csv`
  - 负责告诉你：
    - 这个 invocation / tile 上，每个 row 最后去了哪

所以推荐的排查顺序通常是：

1. 先看 `aggregated_tile_ber.csv`
2. 再看 `per_seed_per_tile_ber.csv`
3. 最后去 shared debug 三表里解释原因

---

## 7. 一个例子

假设某次运行里，`aggregated_tile_ber.csv` 出现下面两行：

```csv
run_id,ebn0_db,stream_label,tile_index,agg_pre_errs,agg_pre_bits,agg_pre_ber,agg_post_errs,agg_post_bits,agg_post_ber
two_stream_probe_20260604-191215,3.050000,A,120,900,39072,0.023034397164,3,39072,0.000076781221
two_stream_probe_20260604-191215,3.050000,B,120,915,39072,0.023418304034,0,39072,0.000000000000
```

这表示：

- 在 `tile_index = 120` 这个 tile 比特段上
- A/B 两路的 pre-FEC 难度差不多
- 但 post-FEC 结果已经不同：
  - A 路这里还有 3 个 post-FEC 错误
  - B 路这里是 0 个 post-FEC 错误

这时下一步就应该：

- 去 `per_seed_per_tile_ber.csv` 看这是单个 seed 偶发，还是多个 seed 都存在
- 再去 shared debug 三表里找同批 `run_id` 下相关 tile 的内部调度信息

---

## 8. 一句话总结

这两份新表的职责可以概括成：

- `per_seed_per_tile_ber.csv`
  - 看单个 seed、单个 tile 的 BER 明细

- `aggregated_tile_ber.csv`
  - 看总体上每个 tile 的 BER 分布与 A/B 差异

它们和 shared 内部调试三表配合起来，才能把：

- “哪个 tile BER 高”
- “A/B 为什么在这个 tile 上不一样”

这两类问题真正串起来。
