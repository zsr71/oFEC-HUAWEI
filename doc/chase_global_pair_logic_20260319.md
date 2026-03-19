# `chase_global_pair` 的候选选择逻辑

这份说明对应：

- `decoder_name = "chase_global_pair"`

重点说明它和 `chase_baseline` 的差别：它不再对每个 bit 单独找一对最优竞争码字，而是先固定一对“全局 best / second”，再复用这对码字去算各 bit 的可靠度。

## 1. 顶层参数入口

### `ofec_single`

在单次运行里，相关顶层参数是：

- [apps/ofec_single.cpp](/home/zsr71/projects/newcode/apps/ofec_single.cpp)
  - `kDecoderName`
  - `kChaseL_override`
  - `kChaseNTestOverride`

这里没有单独属于 `global_pair` 的专用参数。

它复用的是通用 Chase 参数：

- `CHASE_L`
- `CHASE_NTEST`
- `beta`

也就是说，在顶层要切到这个方法，主要改的是：

- `kDecoderName = "chase_global_pair"`

### `ofec_sweep`

在 sweep 里，这个 decoder 同样没有专属扫描参数。

相关入口是：

- [apps/ofec_sweep.cpp](/home/zsr71/projects/newcode/apps/ofec_sweep.cpp)
  - `kDecoderName`
  - `kDecoderNameCandidates`
  - `kChaseLCandidates`
  - `kChaseNTestCandidates`

所以 sweep 下如果要比较它，通常是：

- 把 `decoder_name` 扫描到 `chase_global_pair`
- 再配合通用的 `CHASE_L / CHASE_NTEST / beta` 扫描

## 2. 先生成全部候选

真正实现位于：

- [src/rx/ofec/plain/detail/chase256_global_pair_impl.ipp](/home/zsr71/projects/newcode/src/rx/ofec/plain/detail/chase256_global_pair_impl.ipp)

它前半段和 baseline 基本一样：

1. 从 `Lin256` 构造 `y`
2. 得到信道硬判 `hard_ch`
3. 找最不可靠位
4. 生成全部 `NTEST` 个 test pattern
5. 每个 pattern 都跑 BCH 硬解
6. 构造候选码字 `CW_all`
7. 计算每个候选的 `score = -dist`

所以它不是提前减少候选数，而是：

- 先完整生成全部候选
- 再在这些候选上选择全局 best / second

## 3. 全局 best / second 只从 `good` 候选里选

当前实现已经统一口径了：

- `global best`
- `global second`

都只从 `good == true` 的候选里选。

这里的 `good` 含义仍然是：

- BCH 硬解成功

实现里先做的是：

- 从 `comps` 中筛出全部 `good` 候选，得到 `ranked_good`
- 再对 `ranked_good` 排序

排序规则是：

- 先按 `score` 降序
- 若分数相同，则按 `idx` 升序

因此：

- `global_best_idx = ranked_good.front().idx`
- `global_second_idx = ranked_good[1].idx`

前提是：

- 至少有 2 个 `good` 候选

## 4. 如果 `good` 候选不足两个，会发生什么

### 情况 A：只有 1 个 `good`

这时：

- `global_best_idx` 有值
- `global_second_idx = -1`

结果是：

- 全局 best / second 这对竞争码字不完整
- 所有 bit 都无法用“全局对公式”计算
- 每个 bit 最终都会回退到 `beta * sign(ML[j])`

### 情况 B：一个 `good` 都没有

这时：

- `ML` 也选不出来

会退回到：

- `channel hard decision + overall parity`

然后每个 bit 继续走：

- `beta * sign(ML[j])`

所以这个 decoder 在极端情况下并不会崩，而是会自然退回到和 baseline 一样的 fallback 口径。

## 5. `ML` 是怎么选的

`ML` 的选择和“全局 best / second”不是同一件事。

当前实现里：

- `global best / second`：从 `ranked_good` 里取前 2 个
- `ML`：重新从全部 `good` 候选里选分数最高的那个

在当前排序规则下，这两个结果实际上会一致到：

- `ML` 等于 `global best`

但代码层面它们还是分开写的，语义也分得很清楚：

- `ML` 负责给 fallback 提供符号基准
- `global best / second` 负责给“全局对”提供竞争码字

## 6. 每个 bit 的 `omega[j]` 怎么算

这是 `global_pair` 和 baseline 最大的差别。

### baseline

baseline 是：

- 对每个 bit `j`
- 都单独找一对最优的 `Cplus / Cminus`

所以不同 bit 往往会用到不同的竞争码字对。

### global_pair

global_pair 是：

- 先固定全局 `best` 和 `second`
- 对每个 bit `j`，只看这两条码字在 `j` 位上是否不同

然后分两种情况：

#### 情况 1：这两条全局码字在 `j` 位不同

这时：

- bit=0 的那条当 `Cplus`
- bit=1` 的那条当 `Cminus`

然后直接复用 baseline 同样的求和公式来算 `omega[j]`。

也就是说，这个 bit 会真正使用：

- 同一对固定的全局竞争码字

#### 情况 2：这两条全局码字在 `j` 位相同

这时：

- 不再继续给该 bit 另找别的候选对
- 直接回退到 `beta * sign(ML[j])`

所以 `global_pair` 的核心不是“所有 bit 都一定用同一对码字算出 `omega`”，而是：

- 所有 bit 都先试图使用同一对全局码字
- 如果这对码字在某个 bit 上没有形成 0/1 分歧，该 bit 就直接 fallback

## 7. 哪些 bit 更可能走 fallback

如果全局 best / second 两条码字彼此很像，那么：

- 它们在大多数 bit 上都相同

此时：

- 很多 bit 都不会走“全局对公式”
- 而会直接落到 `beta * sign(ML[j])`

所以这个 decoder 的行为和两条全局码字之间的“分歧位数量”关系很大。

分歧位越少：

- 真正用公式算 `omega` 的 bit 越少
- fallback 的 bit 越多

## 8. trace 里看到的 `Cplus / Cminus` 是什么意思

在 `global_pair` 下，trace CSV 里：

- 若某个 bit 使用了全局 best / second
  - 会把这对码字写进 `cplus_bits / cminus_bits`
- 若某个 bit 走的是 fallback
  - 不会伪造一对新的竞争码字

所以 trace 口径和算法本身是一致的：

- 只有真正用了全局对的 bit，才会看到对应的 `Cplus / Cminus`

## 9. 一个具体例子

假设这次 `good` 候选按分数排序后前两名是：

- `A`
- `B`

那么：

- `global best = A`
- `global second = B`

对每个 bit：

- 如果 `A[j] != B[j]`
  - 用 `A/B` 这一对算 `omega[j]`
- 如果 `A[j] == B[j]`
  - 不再找别的候选
  - 直接用 `beta * sign(ML[j])`

所以 `global_pair` 实际上是把 baseline 的“逐 bit 搜索”改成了：

- “单次全局选对 + 按 bit 判断能不能用这对”

## 10. 总结

`chase_global_pair` 的候选逻辑可以概括成一句话：

- 先从全部 `good` 候选里选出全局 best / second；
- 然后所有 bit 都优先尝试用这一对固定竞争码字计算可靠度；
- 若这对码字在某个 bit 上没有形成分歧，该 bit 直接回退到 `beta * sign(ML[j])`。

所以它和 baseline 的根本区别不是：

- 候选怎么生成

而是：

- 外信息阶段是否还为每个 bit 单独重新搜索一对竞争码字
