# `chase_group_minima` 的候选选择逻辑

这份说明对应：

- `decoder_name = "chase_group_minima"`

重点说明它在解码流程里是怎么选候选、怎么分组、以及哪些候选会真正参与 `ML` 和外信息计算。

## 1. 顶层参数入口

### `ofec_single`

在单次运行里，相关顶层参数是：

- [apps/ofec_single.cpp](/home/zsr71/projects/newcode/apps/ofec_single.cpp)
  - `kDecoderName`
  - `kChaseGroupMinimaBits`
  - `kChaseL_override`
  - `kChaseNTestOverride`

其中：

- `kDecoderName = "chase_group_minima"` 时才会进入这条 decoder
- `kChaseGroupMinimaBits` 决定“按前多少个 test-pattern 位分组”

这些参数会进入：

- [include/newcode/ofec_single_runner.hpp](/home/zsr71/projects/newcode/include/newcode/ofec_single_runner.hpp)
  - `ofec_single::Config::decoder_name`
  - `ofec_single::Config::chase_group_minima_bits`

然后在构造 `Params` 时写进：

- [src/ofec_single/ofec_single_params.cpp](/home/zsr71/projects/newcode/src/ofec_single/ofec_single_params.cpp)
  - `params.CHASE_GROUP_MINIMA_BITS = cfg.chase_group_minima_bits;`

### `ofec_sweep`

在 sweep 里，这个 decoder 的相关入口是：

- [apps/ofec_sweep.cpp](/home/zsr71/projects/newcode/apps/ofec_sweep.cpp)
  - `kDecoderName`
  - `kDecoderNameCandidates`
  - `kChaseGroupMinimaBits`
  - `kChaseGroupMinimaBitsCandidates`

规则是：

- 只有当当前 decoder 是 `chase_group_minima` 时
  - `kChaseGroupMinimaBitsCandidates` 才会真正展开
- 否则只保留一个固定值，不参与场景扩展

优先级仍然是：

1. `kChaseGroupMinimaBitsCandidates`
2. `kChaseGroupMinimaBits`

也就是说：

- 候选列表非空时，优先扫描候选列表
- 候选列表为空时，才回退到单值

## 2. 先生成全部候选

真正实现位于：

- [src/rx/ofec/plain/detail/chase256_group_minima_impl.ipp](/home/zsr71/projects/newcode/src/rx/ofec/plain/detail/chase256_group_minima_impl.ipp)

它前半段和 `chase_baseline` 基本一致：

1. 从 `Lin256` 得到 `y` 和信道硬判 `hard_ch`
2. 找最不可靠位
3. 生成全部 `NTEST` 个 test pattern
4. 每个 pattern 都跑一次 BCH 硬解
5. 构造候选码字 `CW_all`
6. 计算每个候选的 `score = -dist`

所以一开始它并不是直接“先分组再生成”，而是：

- 先把全部候选都生成出来
- 再在这些候选上做分组筛选

## 3. 只在 `good` 候选里参加分组竞争

和现在其它几个 decoder 一样，`chase_group_minima` 不会在全部候选里直接挑组内代表，而是：

- 只有 `good == true` 的候选才允许参加组内竞争

这里的 `good` 含义是：

- 该候选在 BCH 硬解码时成功

因此：

- `good = false` 的候选虽然有 `score`
- 但不会被选为任何一组的代表

## 4. 分组是按 `test pattern` 分，不是按解后的码字分

这是这个 decoder 最关键的地方。

组键不是看：

- BCH 解后的最终码字 `CW`

而是看：

- 原始 `test pattern` 在前若干位上的取值

具体实现里，分组位数是：

```cpp
const int group_bits = std::min(std::max(0, p.CHASE_GROUP_MINIMA_BITS), L_eff);
```

也就是说：

- 先取你设置的 `CHASE_GROUP_MINIMA_BITS`
- 再截到 `[0, L_eff]` 之间

所以实际用于分组的位数是：

```text
实际分组位数 = min(max(0, CHASE_GROUP_MINIMA_BITS), L_eff)
```

然后组数就是：

```text
group_count = 2 ^ group_bits
```

例如：

- `group_bits = 3` 时，有 `8` 组
- `group_bits = 2` 时，有 `4` 组
- `group_bits = 0` 时，只有 `1` 组

这里给一个具体例子。

假设：

- `L_eff = 6`
- `CHASE_GROUP_MINIMA_BITS = 3`

那么：

- 实际 `group_bits = 3`
- 一共 `2^3 = 8` 组

现在有 4 个 test pattern：

```text
候选 A 的 pattern = 0 1 1 0 0 1
候选 B 的 pattern = 0 1 1 1 1 0
候选 C 的 pattern = 1 0 0 0 1 1
候选 D 的 pattern = 0 0 1 1 0 0
```

分组时只看前 `group_bits=3` 位，所以它们的组键分别是：

```text
A: 前3位 = 0 1 1
B: 前3位 = 0 1 1
C: 前3位 = 1 0 0
D: 前3位 = 0 0 1
```

也就是说：

- A 和 B 会被分到同一组
- C 在另一组
- D 在第三组

注意这里看的是：

- “原始 test pattern 的前 3 位是否相同”

不是看：

- “BCH 解完以后得到的最终码字是否相同”

所以即使 A 和 B 最后解出来的 `CW` 很不一样，它们也还是同一组；  
反过来，即使 A 和 C 最后恰好解出了很像的码字，只要它们的 test pattern 前 3 位不同，它们就不在同一组。

再换一种更接近代码的写法看，组键就是把前 `group_bits` 位当成二进制数来拼：

```text
0 1 1 -> group_key = 6
1 0 0 -> group_key = 1
0 0 1 -> group_key = 4
```

这里的位权顺序按当前实现是：

- 第 0 位 pattern 决定最低位
- 第 1 位 pattern 决定次低位
- 第 2 位 pattern 决定第三位

所以 `0 1 1` 在当前代码里对应的 `group_key` 其实是：

```text
0*(2^0) + 1*(2^1) + 1*(2^2) = 6
```

这就是“按 test pattern 分组”的真正含义：

- 先看翻的是哪几个最不可靠位中的前若干位
- 用这些翻转位模式给候选分桶
- 每个桶里再选一个组内最优 good 候选

## 5. 每组怎么选“组内最优”

对每个 `good` 候选，代码会：

1. 读取该候选对应的 test pattern
2. 取这个 pattern 的前 `group_bits` 位
3. 把这几位拼成一个 `group_key`
4. 放进对应组里竞争

组内竞争规则是：

- 先比较 `score`
- `score` 更高的保留
- 若 `score` 相同，则 `idx` 更小的保留

因此，每一组最多保留一个代表候选。

这就是 `group_minima` 里 `minima` 的含义：

- 每组里保留一个“组内最优”代表
- 这里“最优”按当前实现等价于“`score` 最大”

## 6. 某一组没有 `good` 候选怎么办

允许为空。

如果某一组里一个 `good` 候选都没有：

- 该组不会补人
- 最后 `kept_comps` 里就没有这组的代表

所以最后真正保留下来的候选数不是固定的 `2^group_bits`，而是：

- 最多 `2^group_bits`
- 最少 `0`

也就是说：

```text
实际保留数 <= 2^group_bits
```

## 7. 组内最优挑完之后，后面怎么用

所有非空组的代表会组成：

- `kept_comps`

然后这批代表会按 `score` 再排一次序，但这一步只是整理顺序，不会再放回被裁掉的候选。

后续两个关键步骤都只在 `kept_comps` 上做：

- `ML` 选择
- 每个 bit 的 `Cplus / Cminus` 搜索

所以 `group_minima` 的 pruning 影响的是整条后续软输出路径，不只是最后一小步。

## 8. `ML` 是怎么选的

`ML` 直接在 `kept_comps` 里选分数最高的那个。

因为 `kept_comps` 本身就只包含组内最优的 `good` 候选，所以这里的 `ML` 等价于：

- 在“每组最多保留一个代表”的前提下，再选全局得分最好的代表

如果 `kept_comps` 为空：

- 说明所有组都没有 `good` 候选
- 这时会回退到 `channel hard decision + overall parity`

## 9. 每个 bit 的 `Cplus / Cminus` 怎么找

对每个 bit `j`：

- 只在 `kept_comps` 里找
- 看哪个代表在 `j` 位上是 `0`
  - 从里面选 `score` 最高的，作为 `Cplus`
- 看哪个代表在 `j` 位上是 `1`
  - 从里面选 `score` 最高的，作为 `Cminus`

如果两边都找到了，就按 baseline 同样的公式算 `omega[j]`。

如果某一边缺失：

- `omega[j]` 会保持 `NaN`
- 然后退回到 `beta * sign(ML[j])`

所以这里的 fallback 触发条件是：

- 不是“没有代表候选”
- 而是“虽然有代表候选，但在这个 bit 上没有形成 `0/1` 两边竞争”

## 10. 一般化支持别的 `L`

当前实现不是只写死给 `L=6` 用的。

因为它的分组位数是：

- `min(CHASE_GROUP_MINIMA_BITS, L_eff)`

所以：

- `L_eff < CHASE_GROUP_MINIMA_BITS` 时，会自动缩小分组位数
- `L_eff >= CHASE_GROUP_MINIMA_BITS` 时，就按你给定的分组位数来

这就是为什么它可以支持别的 `L`，而不是只能支持“64 候选分 8 组”这一种情况。

## 11. 一个具体例子

假设：

- `CHASE_L = 6`
- `CHASE_NTEST = 64`
- `CHASE_GROUP_MINIMA_BITS = 3`

那么：

1. 会先生成 64 个 test pattern
2. 用前 3 个 pattern 位做组键
3. 一共对应 8 组
4. 每组只在本组的 `good` 候选里保留一个分数最高的代表
5. 最后最多得到 8 个代表
6. `ML` 和逐 bit 的 `Cplus/Cminus` 都只在这最多 8 个代表上做

如果某两组没有 `good` 候选，那么最后就只剩：

- 6 个代表

## 12. 总结

`chase_group_minima` 的候选逻辑可以概括成一句话：

- 先生成全部候选；
- 再按 test pattern 的前若干位分组；
- 每组只在 `good` 候选里选一个组内最优；
- 然后只用这些组代表去做 `ML` 和逐 bit 外信息搜索。

所以它更准确的含义是：

- “按 test pattern 分组后的组内最优剪枝版 Chase”

而不是：

- “简单的全局 Top-K”
