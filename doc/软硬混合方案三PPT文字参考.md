# 软硬混合方案三 PPT 文字参考

本文用于后续制作 PPT，重点介绍方案三的主体思想、整体流程和硬件时序。本文不展开具体代码接口、文件路径或实现细节。

## 1. 方案三的核心目标

方案三的核心目标是把 oFEC 解码过程中的行级任务拆成更清晰的几类，让每一行在进入软解码资源竞争之前，就尽可能完成前置判断。

传统流程里，所有没有早停的行通常都被视为需要进入 SISO / Chase 软解码的候选。这样做虽然结构简单，但会把一些本来可以用硬判决快速完成的行也送入软解码资源池，造成 Chase 资源消耗和调度压力。

方案三的思路是：

> 在 MUX / Chase 之前增加一层 row 级分流，把已经 early-stop 的行、可以硬完成的行、仍需软解码的行明确拆开。

最终每一行只会落入下面几类之一：

- `EarlyStopAction`：已经满足早停条件，直接走早停动作。
- `HardFinish`：通过前置硬判决或 syndrome 分类确认可以硬完成，不再进入 Chase。
- `SoftDecode`：无法在前端确认完成，需要继续参与 MUX 和 Chase。
- `Unscheduled`：由于 SISO 预算限制，本轮没有获得软解码资源。

这种拆分的意义是：MUX 不再关心一行为什么需要处理，只负责在真正需要软解码的行里分配 SISO 资源。

## 2. 为什么不是简单增加一个状态

早停、硬完成、软解码、未调度这四类行在系统中的语义不同。

早停行表示当前判据已经认为可以提前停止，但它的输出通常由早停动作生成。硬完成行表示前端已经完成了 BCH 硬纠或等价确认，它应该输出硬纠后的结果。未调度行则表示本轮没有获得软解码资源，不应该产生有效输出。

如果只是在原来的 MUX 状态里继续增加状态，很容易让一个状态同时承担两种语义：

- 是否参与 SISO 资源竞争。
- 最终执行哪一种输出动作。

方案三的本质是把这两件事拆开：

- 资源竞争只针对 `SoftDecode` 行。
- 执行动作由每一行的最终分类决定。

这样系统结构更清楚，也更适合后续扩展硬件资源约束、可靠度保护和更复杂的分类策略。

## 3. 高层处理流程

方案三可以从更高层概括成三步：

```text
第一拍：获得 32 个码字的校验式
第二拍：根据 syndrome 对 32 个码字分类
第三拍：按照分类结果分别生成硬解码外信息，或送入 Chase 软解码
```

这比直接从软件流程描述更适合 PPT 表达，因为它把方案三的核心动作抽象成了硬件流水中的三个阶段。

第一拍关注的是“观察当前 32 行码字的状态”。系统并行计算每个码字的 `S0/S1/S3`，得到后续分类所需的 syndrome 信息。

第二拍关注的是“判断每个码字属于哪一类”。根据 syndrome，可以把码字分成 clean、parity-only、single-error、two-error candidate、hard-fail candidate 等类别。

第三拍关注的是“按类别分流”。能够在前端确认的码字直接生成硬解码外信息；无法确认或需要软辅助的码字进入 Chase 解码队列，由 MUX / SISO 调度选择真正参与 Chase 的行。

因此方案三的关键变化是：

> Chase 不再面对全部未早停码字，而是只处理 syndrome 分类后仍然需要软解码的码字。

## 4. 前置硬完成的分类思想

方案三中的前置硬完成可以有不同实现方式。最直接的方式是对当前行做 BCH 硬译码，成功则认为这行可以硬完成，失败则送入软解码。

更适合硬件的方式是先使用 syndrome 做快速分类。对于 eBCH(256,239)，可以基于三个关键 syndrome：

- `S0`：整体 parity syndrome。
- `S1`：BCH 主体的一阶 syndrome。
- `S3`：BCH 主体的三阶 syndrome。

通过 `S0/S1/S3` 可以把行快速分成：

- `Clean`：没有错误，可以直接完成。
- `ParityOnly`：只有整体 parity 位错误，可以直接修正。
- `OneMain`：BCH 主体存在 1 位错误，可以直接定位并修正。
- `OneMainPlusParity`：BCH 主体 1 位错误，同时 parity 位也需要修正。
- `TwoCandidate`：疑似 BCH 主体 2 位错误，需要进一步确认。
- `HardFailCandidate`：前端无法解释，直接送软解码。

前四类可以在快速分类阶段基本确定为硬完成。真正需要注意的是 `TwoCandidate`。

## 5. TwoCandidate 的特殊性

`TwoCandidate` 不是最终结果，而是一个中间候选状态。

当 syndrome 满足疑似两错条件时，只能说明这行可能是 BCH 主体中的 2 位错误，但还不能直接认为它可以硬完成。它必须进一步经过 BCH hard decoder，或者等价的 two-error locator 与 Chien search，确认是否真的存在两个合法错误位置。

确认成功后，这行才可以变成 `TwoMain / HardFinish`。

确认失败后，这行必须回到 `HardFail / SoftDecode`，进入后续 MUX 和 Chase。

因此在方案三里：

```text
TwoCandidate 不是 HardFinish
TwoCandidate 也不是 Chase
TwoCandidate 是需要延迟确认的中间状态
```

这个点对硬件时序非常重要。

## 6. 三拍硬件时序

如果把方案三映射到硬件流水，可以按三拍理解。

第 1 拍：获得 32 个校验式。

```text
输入 32 个码字
并行计算每个码字的 S0 / S1 / S3
```

这一拍的输出不是最终调度结果，而是 32 组 syndrome。它们描述每个码字当前是否满足 clean、单错、两错候选或失败类判据。

第 2 拍：进行 syndrome 分类。

```text
Clean / ParityOnly / OneMain / OneMainPlusParity
  -> 归入可硬解码外信息生成类

明显 HardFailCandidate
  -> 归入 Chase 软解码候选类

TwoCandidate
  -> 送入 BCH hard decoder 或 two-error confirm 模块继续确认
```

这一拍的输出是每个码字的分类结果。对于 clean、parity-only、single-error 类，分类已经足够决定后续硬处理。对于明显 hard-fail 类，可以直接进入 Chase 候选。对于 `TwoCandidate`，还需要等待确认结果。

第 3 拍：按照分类结果分别处理。

```text
Clean / ParityOnly / OneMain / OneMainPlusParity
  -> 生成硬解码外信息

TwoCandidate confirm 成功
  -> 生成硬解码外信息

TwoCandidate confirm 失败
  -> 进入 Chase 解码

HardFailCandidate
  -> 进入 Chase 解码
```

第三拍结束后，可以形成完整的两类输出：

```text
hard-output rows:
  直接生成硬解码外信息

chase rows:
  进入 MUX / SISO 调度，再执行 Chase 软解码
```

从更抽象的角度看，方案三的第三拍就是一次分流：

```text
可硬完成码字 -> 硬解码外信息
不可硬完成码字 -> Chase 软解码
```

## 7. 为什么 Chase 调度建议放在第三拍之后

从硬件调度角度看，第二拍已经可以知道一部分码字一定可以硬处理，也可以知道一部分码字一定需要 Chase。但由于 `TwoCandidate` 还没有确认完成，第二拍还不能得到完整的 Chase 候选集合。

如果第二拍就先对已确定的 Chase 行做 MUX 调度，会遇到一个问题：

```text
第二拍可能已经把 SISO 预算用完
第三拍 TwoCandidate confirm 失败后，又产生新的 Chase 候选
```

这时系统必须支持抢占、重排、额外缓冲，或者接受某些 two-fail 行被推迟处理。这会显著增加 MUX 和调度控制复杂度。

因此第一版硬件实现更建议采用保守时序：

> 等 TwoCandidate 确认完成后，再统一生成完整 Chase candidate mask，并统一执行 MUX / SISO 调度。

这样虽然 Chase 调度晚一拍，但控制逻辑更清晰，也不会错误占用或遗漏 Chase 资源。对于 PPT 表达，可以概括为：

> 第三拍再统一决定“哪些给硬解码外信息，哪些进入 Chase”，可以避免跨拍调度回流。

## 8. 行级输出合并

方案三最终仍然需要输出一个完整 tile 的结果。不同类型的行输出来源不同：

- 早停行：由 early-stop action 生成输出。
- 硬完成行：由硬纠确认后的码字生成输出。
- 软解码行：由 Chase / soft decoder 生成输出。
- 未调度行：本轮不产生有效输出。

这些结果最后重新合并到完整行域中，再统一写回。

这种设计的好处是：前端可以做复杂分流，但后端看到的仍然是一份完整 tile 输出，不需要把整个系统拆成多个不一致的数据域。

## 9. 方案三的价值

方案三的主要价值有三点。

第一，减少软解码资源压力。能够 early-stop 或 hard-finish 的行不会再参与 Chase 资源竞争，SISO 预算集中给真正困难的行。

第二，结构语义清楚。早停、硬完成、软解码、未调度被明确拆开，避免一个 MUX 状态同时表达调度和执行动作。

第三，便于硬件扩展。后续可以自然加入独立 hard prepass 预算、TwoCandidate 队列、LLR 可靠度保护、Suspicious 分类和更复杂的调度策略。

## 10. PPT 可用总结

方案三可以总结为：

> 方案三把 32 个码字的处理前移到 syndrome 分类阶段：第一拍获得 32 个校验式，第二拍完成 syndrome 分类，第三拍按分类结果决定生成硬解码外信息还是进入 Chase 软解码。这样可以把明显可硬完成的码字提前移出 Chase 资源池，让 SISO 资源集中服务真正困难的码字。

一句话版本：

> 方案三的本质是：先分类，再分流；可硬完成的生成外信息，不可硬完成的进入 Chase。
