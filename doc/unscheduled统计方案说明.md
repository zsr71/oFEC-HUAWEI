# unscheduled 统计方案说明

本文说明在当前 oFEC + MUX + early-stop 机制下，如何定义和统计 `unscheduled`，以及推荐把统计逻辑放在什么位置。

## 1. 为什么需要这个统计

当前一个 soft tile 内的 row，在进入行级执行前会被分成三类：

- `EarlyStopped`
- `NeedSiso`
- `Unscheduled`

如果只看最终 BER，很难直接判断：

- 有多少 row 原本需要正常解码，但因为预算不够被裁掉了
- 不同 MUX 规则下，被裁掉的 row 比例有没有变化
- early-stop 是否缓解了 SISO 资源压力

所以增加 `unscheduled` 统计，可以更清楚地回答：

- 本轮到底有多少 row 被 MUX 丢掉了
- 被丢掉的比例是多少
- 新旧 MUX 对资源分配的影响是什么

## 2. 推荐的定义

这里推荐把 `unscheduled` 定义为：

- **本轮原本需要正常 SISO/Chase 解码，但最终没有被保留下来的 row**

在当前机制下，这正对应：

- `mux_state == Unscheduled`

所以最直接的定义是：

```text
rows_unscheduled = count(mux_state == Unscheduled)
```

## 3. 为什么 `EarlyStopped` 不算 unscheduled

`EarlyStopped` 的 row 虽然不占正常 SISO 预算，
但它本轮已经被处理了：

- 已经做了 early-stop 判定
- 已经被决定走 early-stop action

所以它不能算成 `unscheduled`。

`unscheduled` 只应该包含：

- 本轮完全没被处理、保持旧值的 row

## 4. 推荐统计的量

建议每个 tile 至少统计下面这些量：

- `rows_total`
- `rows_need_siso_before_mux`
- `rows_unscheduled`
- `unscheduled_pct`

其中：

```text
unscheduled_pct = rows_unscheduled / rows_total
```

如果还想更直接分析 MUX 对正常解码路径的影响，再加一个：

- `unscheduled_among_need_pct = rows_unscheduled / rows_need_siso_before_mux`

这个值更能直接回答：

- 在原本需要 SISO 的那批 row 里，有多少被裁掉了

## 5. 当前代码里最适合统计的位置

最推荐的统计位置是：

- **tile 级**
- **在最终 `mux_state` 确定之后**

原因是此时三类状态已经都明确了：

- `EarlyStopped`
- `NeedSiso`
- `Unscheduled`

不需要再去反推 row 到底经历了什么。

## 6. 为什么推荐基于 `mux_state` 统计

因为当前流程中：

- `EarlyStopped` 已经在 early-stop 阶段确定
- `NeedSiso` 是否被保留，由 MUX 决定
- `Unscheduled` 是 MUX 最终裁掉的结果

所以基于最终 `mux_state` 做统计，语义最清楚。

可以直接理解为：

- `NeedSiso`：本轮保留下来的正常 SISO row
- `EarlyStopped`：本轮走 early-stop action 的 row
- `Unscheduled`：本轮完全没处理的 row

## 7. 不推荐只看 `produced_rows`

虽然也可以在 row decoder 后，通过 `produced_rows=false` 去找“没产出输出的 row”，
但这不是最好的主统计口径。

原因是：

- `produced_rows=false` 既可能来自 `Unscheduled`
- 也可能来自某些特殊动作失败
- 不如直接在 tile 层根据 `mux_state` 统计

所以：

- `produced_rows` 更适合作为辅助校验
- `mux_state` 更适合作为主统计口径

## 8. 推荐的 tile 级定义

假设当前 tile 最终 `mux_state` 已经确定，则可以这样定义：

- `rows_total = state.size()`
- `rows_need_siso_before_mux = count(初始状态为 NeedSiso)`
- `rows_unscheduled = count(state == Unscheduled)`
- `unscheduled_pct = rows_unscheduled / rows_total`
- `unscheduled_among_need_pct = rows_unscheduled / rows_need_siso_before_mux`

这里有一个细节：

- `rows_need_siso_before_mux` 需要在 MUX 裁剪之前统计
- 因为裁剪之后，那些被裁掉的 row 已经不再是 `NeedSiso`

所以推荐做法是：

1. 先统计一遍 MUX 前的 `NeedSiso` 数量
2. 再做 MUX 裁剪
3. 最后统计 `Unscheduled`

## 9. 推荐的数据流

推荐沿着下面这条链做：

### 第一步：tile 级产生统计

在 `TileProcessResult` 里增加：

- `rows_need_siso_before_mux`
- `rows_unscheduled`

### 第二步：window 级按 tile 累加

像现在累计 early-stop 统计一样，
把每个 tile 的 `unscheduled` 统计累加到 window 级数组。

### 第三步：single / sweep 输出

在最终输出里增加：

- 每个 tile 的 `unscheduled_pct`
- 必要时加 `unscheduled_among_need_pct`

## 10. 推荐先输出哪些

如果想控制输出复杂度，建议优先输出这两个：

- `rows_unscheduled`
- `unscheduled_pct`

如果还想更好地解释资源压力，再加：

- `rows_need_siso_before_mux`
- `unscheduled_among_need_pct`

这样就能区分：

- 是因为总 row 多
- 还是因为真正需要 SISO 的 row 太多

## 11. 和新旧 MUX 对比时怎么用

如果后面你想比较：

- 旧 MUX
- 新 MUX

那么 `unscheduled` 统计非常适合作为解释 BER 差异的辅助量。

例如可以一起看：

- `post_ber`
- `tile_early_stop_pct`
- `tile_unscheduled_pct`
- `tile_unscheduled_among_need_pct`

这样就能分清：

- BER 变化是因为 early-stop 判得多了
- 还是因为 MUX 更会挑 row 了
- 还是因为被裁掉的 row 变少了

## 12. 一句话总结

最推荐的定义是：

- `unscheduled = 本轮原本需要正常 SISO，但最终被 MUX 裁掉、完全没被处理的 row`

最推荐的统计位置是：

- tile 级、最终 `mux_state` 确定之后

最推荐先记录的量是：

- `rows_unscheduled`
- `unscheduled_pct`

如果还想更细：

- `rows_need_siso_before_mux`
- `unscheduled_among_need_pct`

这样最方便后面做新旧 MUX 的性能对比。
