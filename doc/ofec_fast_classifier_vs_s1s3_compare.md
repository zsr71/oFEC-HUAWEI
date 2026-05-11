# `run_fast_classifier_hard_finish(...)` 与朋友版 `S1/S3` 分类器对比

本文对比两套逻辑：

- 当前仓库实现：`run_fast_classifier_hard_finish(...)`
  - 代码位置：[src/rx/ofec/detail/ofec_tile_dispatch.ipp](/home/zsr71/projects/newcode/src/rx/ofec/detail/ofec_tile_dispatch.ipp:200)
- 朋友提供的代码：仅基于 `S1/S3` 和 `Trace(mu)` 的快速分类器

目标问题：

- 两者是否存在逻辑差异？
- 如果有，差异属于“功能范围不同”还是“判断矛盾”？

## 结论

存在逻辑差异，但不是简单的“谁对谁错”。

更准确地说：

- 朋友版逻辑更像是 `BCH(255,239)` 主体上的一级 syndrome 分类器。
- 当前仓库逻辑是面向扩展 `256` 位码字的“可直接 HardFinish 的前置分流器”。

因此两者在目标上就不完全相同：

- 朋友版回答的是：`0 错 / 1 错 / 2 错候选 / 非法`
- 当前仓库回答的是：`这行现在能不能直接硬完成；如果能，属于哪一类`

如果只看 `255` 位 BCH 主体，二者在 `S1/S3` 的基础判据上没有根本冲突；但一旦把扩展 parity 位和“最终是否允许 HardFinish”算进去，当前仓库比朋友版多了几层判断。

## 当前仓库逻辑

当前实现使用三个量：

- `S0`：256 位整体 parity syndrome
- `S1`
- `S3`

其中 `S0` 的定义在 [src/rx/ofec/detail/ofec_tile_dispatch.ipp](/home/zsr71/projects/newcode/src/rx/ofec/detail/ofec_tile_dispatch.ipp:142)：

- 对 256 位码字逐位异或
- `S0 = 0` 表示整体 parity 正确
- `S0 = 1` 表示整体 parity 不正确

当前实现的判定分支如下：

1. `S0=0, S1=0, S3=0`
- 判为 `Clean`
- 直接 `HardFinish`

2. `S0=1, S1=0, S3=0`
- 判为 `ParityOnly`
- 翻 overall parity 位后 `HardFinish`

3. `S1!=0, S3=S1^3`
- 先按单错模型通过 `log(S1)` 定位 BCH 主体错误位
- 翻该主体位
- 若原始 `S0=0`，再翻 overall parity 位，分类为 `OneMainPlusParity`
- 若原始 `S0=1`，不翻 overall parity 位，分类为 `OneMain`

4. `S0=0, S1!=0, S3!=S1^3`
- 进入两错分支
- 但不是直接认定 `TwoMain`
- 而是调用现有 BCH `t=2` hard decoder 做确认
- 只有 `decode` 成功且 `corrected_errors == 2`，才 `HardFinish`

5. 其余情况
- `HardFail`
- 交回 soft path

## 朋友版逻辑

朋友版代码只用 `S1` 和 `S3`，逻辑可以整理成：

1. `S1=0`
- `S3=0`：判为 0 错
- `S3!=0`：判为非法

2. `S1!=0`
- 若 `S3=S1^3`：判为 1 错
- 否则计算
  - `mu = (S1^3 + S3) / S1^3`
  - 若 `Trace(mu) = 0`：判为 2 错候选
  - 若 `Trace(mu) != 0`：判为非法

这段逻辑没有显式使用 `S0`，也没有对 two-error candidate 做最终确认。

## 差异一：当前仓库显式处理了扩展 parity，朋友版没有

这是两者最本质的差异。

朋友版只看 `S1/S3`，因此无法区分下面两类情况：

- 只有 overall parity 位错
- BCH 主体 1 位错，同时 overall parity 位也错

而当前仓库通过 `S0` 把它们拆开了：

- `S0=1, S1=0, S3=0` -> `ParityOnly`
- `S1!=0, S3=S1^3, S0=0` -> `OneMainPlusParity`
- `S1!=0, S3=S1^3, S0=1` -> `OneMain`

因此：

- 对扩展 256 位码字来说，朋友版的分类信息不够完整。
- 对 255 位 BCH 主体来说，朋友版可以看作是“忽略整体 parity 维度”的简化版本。

## 差异二：two-error 分支的确认强度不同

朋友版在 `S1!=0` 且 `S3!=S1^3` 时，使用

- `mu = (S1^3 + S3) / S1^3`
- `Trace(mu) == 0`

来判断“是否是双错候选”。

这一步本质上只是筛选 `TwoCandidate`，不是最终确认“双错一定成立”。

当前仓库没有把这个候选条件直接作为 `TwoMain` 的充分条件，而是更保守：

- 先要求 `S0=0`
- 再调用现有 `bch_255_239_decode_hiho_cw_255(...)`
- 只有返回成功且 `corrected_errors == 2` 才算 `TwoMain`

因此两者在两错分支上的关系是：

- 朋友版：做“候选分类”
- 当前仓库：做“最终可执行 hard-finish 的确认”

从工程语义上讲，当前仓库更严格，误判风险更低。

## 差异三：输出语义不同

朋友版的输出是 `flags[]`，更像 syndrome 分类结果。

当前仓库的输出是：

- `bool hard_ok`
- `HybridRowClass out_class`
- 若 `hard_ok=true`，还要顺带产出可直接写回的 `y2`

也就是说，当前仓库不是单纯在回答“它像几错”，而是在回答：

- 这行现在能否不经过 Chase，直接作为 `HardFinish` 收尾？

因此当前仓库在每个分支末尾还会调用 [hard_word_valid_256(...)](/home/zsr71/projects/newcode/src/rx/ofec/detail/ofec_tile_dispatch.ipp:185) 做最终合法性验证。朋友版没有这一步。

这意味着即便某条路径在 syndrome 层面“看起来像对”，只要最终码字不合法，当前仓库仍会把它打回 `HardFail`。

## 差异四：`Trace(mu)==0` 在当前仓库里没有直接落地

朋友版在 two-error candidate 上采用的是解析判据：

- `Trace(mu)==0`

当前仓库并没有实现这条闭式判据，而是直接复用现有 BCH `t=2` decoder。

所以从实现层面说：

- 朋友版更接近算法分类器
- 当前仓库更接近“拿现成 hard decoder 做最终确认”的落地版

这不一定表示当前仓库逻辑更“数学正确”，而是表示它更偏工程保守实现。

## 是否存在直接矛盾

### 1. `0 错 / 1 错` 主体判据

这里没有本质矛盾。

两边都认可：

- `S1=0, S3=0` 对应 0 错主体模型
- `S3=S1^3` 对应 1 错主体模型

差异只在于当前仓库还额外引入了 `S0`，把扩展 parity 的情况分开处理。

### 2. `2 错` 判据

这里有实现策略差异，但不是数学上的直接冲突。

朋友版说的是：

- `Trace(mu)==0` -> 2 错候选

当前仓库说的是：

- 不直接相信候选
- 继续让 BCH hard decoder 证实
- 成功且确认为 2 错，才 `HardFinish`

因此当前仓库相当于：

- 接受 “two-error 需要进一步确认” 这个大前提
- 只是没有采用朋友版的 `Trace(mu)` 作为一级筛选条件

### 3. `OneMainPlusParity` 分支

这里要特别说明一下，因为代码初读时容易误会。

当前实现中：

- `S0=1, S1=0, S3=0` -> `ParityOnly`
- `S1!=0, S3=S1^3, S0=1` -> `OneMain`
- `S1!=0, S3=S1^3, S0=0` -> `OneMainPlusParity`

这和扩展 parity 码的总错误个数奇偶性是一致的：

- 仅 BCH 主体 1 位错：总错误数为奇数，因此 `S0=1`
- BCH 主体 1 位错 + overall parity 位也错：总错误数为偶数，因此 `S0=0`

所以这个分支虽然和“直觉上的 parity 错就 `S0=1`”容易混淆，但从扩展码的总 parity syndrome 角度看是自洽的。

## 哪些情况下可以认为两者“等价”

如果满足下面两个前提，可以把朋友版看作当前仓库的简化前级分类器：

1. 暂时忽略扩展 parity 位，只讨论 255 位 BCH 主体
2. 把 `Trace(mu)==0` 仅理解为 `TwoCandidate`，而不是最终 `TwoMain`

在这两个前提下：

- 朋友版提供的是一级候选分类
- 当前仓库可以看作在这个基础上又加了：
  - `S0` 维度
  - 最终合法性校验
  - BCH `t=2` hard decoder 确认

## 哪些情况下不能认为两者“等价”

下面这些场景下，二者不能直接视为同一逻辑：

1. 你希望识别 `ParityOnly`
- 朋友版做不到

2. 你希望区分 `OneMain` 和 `OneMainPlusParity`
- 朋友版做不到

3. 你希望 `TwoMain` 是“已经确认能硬纠完成”
- 朋友版不够，最多只能给 `TwoCandidate`

4. 你希望分类结果可以直接驱动 `HardFinish`
- 朋友版还缺最终合法性验证和输出 materialization

## 建议

如果你的目标是“比较算法判据是否一致”，建议这样理解：

- 朋友版适合当作 `255` 位 BCH 主体上的一级快速分类器参考
- 当前仓库适合当作 `256` 位扩展码字上的工程落地版前置 hard-finish 分流器

如果后面要把朋友版思路合入当前仓库，比较合理的方式不是直接替换现有逻辑，而是：

1. 保留当前 `S0` 处理
2. 将 `Trace(mu)==0` 引入为 `TwoCandidate` 的前级快速筛选
3. 后面仍保留 BCH `t=2` decoder 或等价根搜索做最终确认
4. 最后继续保留 `hard_word_valid_256(...)` 作为 `HardFinish` 出口校验

这样可以同时保留：

- 朋友版的低成本双错候选筛选
- 当前仓库的扩展 parity 语义
- 当前仓库的工程安全性

## 一句话结论

两者有逻辑差异，但主要是“分层不同、目标不同”：

- 朋友版是 `S1/S3` 的候选分类器
- 当前仓库是包含 `S0`、最终合法性验证和 `t=2` hard decoder 确认的 `HardFinish` 决策器

所以不能直接说两者等价；更不能把朋友版里的 `Trace(mu)==0` 直接当成当前仓库里的 `TwoMain`。
