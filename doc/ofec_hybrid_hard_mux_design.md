# oFEC hybrid hard MUX 后续实现设计

## 1. 背景

当前 hybrid prepass 已经接入 `TileDispatchPlan`，主流程大致是：

```text
early-stop
  -> run_hybrid_prepass
  -> run_mux_on_soft_candidates
  -> decode_tile_with_plan
```

在 `FriendS1S3WithS0Classifier` 模式下，当前分类器不是只分类。它会在分类阶段直接完成部分 hard-finish 操作：

- `Clean`：当前代码里校验合法后直接生成 `y2`；新设计里改为复用 early-stop 处理路径
- `ParityOnly`：翻 overall parity 后生成 `y2`
- `OneMain`：根据 `S1` 定位并翻 1 个主体 bit 后生成 `y2`
- `OneMainPlusParity`：翻 1 个主体 bit 和 overall parity 后生成 `y2`
- `TwoMain`：调用 `bch_255_239_decode_hiho_cw_255(...)`，要求 `corrected_errors == 2`，再生成 `y2`

如果要模拟“每个 tile 只有 `n1` 个 hard decoder / hard-finish lane”，这些动作不能继续放在分类器里。新的设计要把职责拆开：

```text
分类器：只判断 row 属于哪一类
现有 deferred/backfill 优先级：继续决定哪些 hard candidate 需要从 soft path 收回
hard MUX：只限制 hard path 资源数量，不在 MUX 内做优先级
decode_tile_with_plan：对 hard MUX 选中的 row 执行真正 hard decode / hard finish
```

目标流程：

```text
early-stop
  -> classify-only hybrid classifier
  -> existing deferred/backfill priority
  -> hard MUX
  -> soft MUX
  -> decode_tile_with_plan
       -> EarlyStopAction materialize
       -> Clean/free materialize
       -> Hard executor materialize
       -> scheduled SoftDecode run Chase
```

这里有两个关键点：

- hard executor 更适合放在 `decode_tile_with_plan` 里面，因为它本质上是执行阶段的一部分，和现有 `materialize_early_stop_rows(...)`、`materialize_hard_finish_rows(...)` 语义一致。
- hard MUX 只做资源门控，不做 class priority。class priority 复用当前 deferred/backfill 逻辑，放在 hard MUX 前面。

## 2. 总体资源模型

一个 tile 内有两类受限资源：

```text
n1 = hard decoder / hard-finish lane 数量
n2 = soft decoder / SISO 数量
```

经过 early-stop 和 classify-only 后，每个 row 先进入以下几类：

```text
EarlyStopAction
HardCandidate
SoftDecode
```

含义：

- `EarlyStopAction`：早停命中，不进 deferred/backfill priority，不进 hard MUX，不进 soft MUX。
- `HardCandidate`：分类器认为这行可以走 hard path，但还没有真正翻 bit、BCH decode 或生成 `y2`。
- `SoftDecode`：分类器认为这行不能安全走 hard path，需要走 soft path。

之后先复用现有 deferred/backfill 优先级逻辑：

```text
HardCandidate 先暂时留在 SoftDecode 集合里。
如果当前 SoftDecode 数量超过 n2，
就按现有 deferred priority 从 HardCandidate 里收回一部分，
形成需要进入 hard path 的 ordered reclaim list。
```

再做 hard MUX：

```text
ordered reclaim list 进入 hard MUX。
hard MUX 最多放行 n1 个。
hard MUX 内不做优先级排序，只按 existing priority 已经排好的顺序截断。
```

最后做 soft MUX：

```text
所有 SoftDecode 进入现有 soft MUX。
抢到 SISO 的 row 后续跑 Chase。
没抢到的 row 变成 Unscheduled。
```

外部 app 配置层应和现有 `kSisoActiveList` 对齐，新增一张并列表：

```cpp
static const std::vector<int> kHiHoActiveList = {32, 32, 32, 16, 8, 4};
```

参数层对应字段建议命名为：

```cpp
std::vector<int> HIHO_ACTIVE_LIST;
```

两张表的语义分别是：

```text
kSisoActiveList / SISO_ACTIVE_LIST:
  每个 tile 可用的 soft decoder / SISO 数量，即 n2。

kHiHoActiveList / HIHO_ACTIVE_LIST:
  每个 tile 可用的 hard decoder / hard-finish lane 数量，即 n1。
```

对于 two-stream shared，`kHiHoActiveList` 的口径也应和 `kSisoActiveList` 一样使用 shared row 域。例如 shared tile 一次合并 A/B 两路后有 64 个 decoder row，则 `kHiHoActiveList[t]` 表示这个 64-row shared tile 内可用的 hard lane 数量。

## 3. 三个核心改动

### 3.1 分类器只分类

分类器职责：

```text
输入：lin_vec
输出：HybridRowClass / 分类标签
```

分类器不再做：

- 不翻 bit
- 不调用 `bch_255_239_decode_hiho_cw_255(...)`
- 不调用 `recompute_overall_parity(...)`
- 不调用 `materialize_hard_finish_lout(...)`
- 不生成 `y2`

分类器只输出：

```text
Clean
ParityOnly
OneMain
OneMainPlusParity
TwoMain
HardFail
```

建议新增 classify-only helper：

```cpp
template <typename CoreLLR>
HybridRowClass classify_friend_s1s3_with_s0(
    const std::array<CoreLLR, newcode::Params::BCH_N>& lin_vec);
```

正式实现继续复用现有 `HybridRowClass`，不新增分类 enum。需要在注释中明确：在 classify-only 阶段，`HybridRowClass` 只是分类标签，不表示已经 hard-finish 成功。

### 3.2 复用现有 deferred/backfill 优先级

这一层位于 classify-only 和 hard MUX 之间，不是新增一套优先级。它应该复用当前 `run_hybrid_prepass(...)` 里已经存在的 deferred/backfill 逻辑。

当前代码里的核心语义是：

```text
所有 deferred hard candidates 先留在 SoftDecode 集合里。
如果 soft 候选数量超过 n2，就按 priority 从 deferred candidates 中收回一部分。
被收回的候选后续走 hard path。
没被收回的候选继续留在 SoftDecode，进入 soft MUX。
```

输入：

```text
hard_candidate_rows
soft_candidate_rows
n2
```

计算：

```text
soft_rows_before_reclaim = soft_candidate_rows.size() + hard_candidate_rows.size()
siso_budget = n2
need_hard_reclaim = max(soft_rows_before_reclaim - siso_budget, 0)
hard_reclaim_count = min(need_hard_reclaim, hard_candidate_rows.size())
```

然后复用现有 priority bucket，从 `HardCandidate` 中选 `hard_reclaim_count` 个，形成 ordered reclaim list。

当前 `ParityOneAndTwoErrorPriority` 的顺序保持不变：

```text
ParityOnly -> OneMain -> OneMainPlusParity -> TwoMain
```

这个顺序就是当前代码里 `hybrid_siso_backfill_reclaim_priority(...)` 的语义。后续实现时应尽量直接复用这段 priority bucket 逻辑，避免再写一套新的排序。

`Clean` 的处理方式固定为：

- `Clean` 不进入 hard MUX，不消耗 `HIHO_ACTIVE_LIST[t]`。
- `Clean` 按类似 early-stop 的免费处理路径直接 materialize。
- `Clean` 不参与 deferred/backfill priority。

deferred/backfill priority 之后：

```text
被收回的 hard candidate -> ordered reclaim list，等待 hard MUX
没被收回的 hard candidate -> SoftDecode，进入 soft MUX
```

### 3.3 hard MUX 只做资源门控

hard MUX 的输入是上一层 deferred/backfill priority 产生的 ordered reclaim list。

hard MUX 的职责只有一个：

```text
最多放行 n1 个 hard candidate
```

hard MUX 不做：

- 不计算 class priority
- 不翻 bit
- 不调用 BCH t=2 decode
- 不生成 `y2`

hard MUX 输出：

```text
hard_scheduled_rows
hard_not_scheduled_rows
```

含义：

- `hard_scheduled_rows`：拿到 hard 资源，后续在 `decode_tile_with_plan` 里执行 hard executor。
- `hard_not_scheduled_rows`：没拿到 hard 资源，转成 `SoftDecode`，后续进入 soft MUX。

第一版 hard MUX 可以按输入顺序截断：

```text
hiho_budget = HIHO_ACTIVE_LIST[t]
hard_scheduled_rows = first hiho_budget rows from ordered reclaim list
hard_not_scheduled_rows = rest
```

等价地，用前文符号表示：

```text
hard_scheduled_rows = first n1 rows from ordered reclaim list
hard_not_scheduled_rows = rest
```

如果之后需要改 class-based hard order，也应该改上一层 deferred/backfill priority，而不是把优先级写进 hard MUX 本身。

## 4. hard executor 放在 decode_tile_with_plan

hard executor 不建议放在 dispatch 阶段。原因是 dispatch 阶段应该只决定 row 的执行计划，真正执行 hard decode / hard finish 更像 decode 阶段的工作。

目标结构：

```text
decode_tile_with_plan
  -> materialize_early_stop_rows
  -> materialize_clean_rows
  -> materialize_hard_scheduled_rows
       -> hard executor
  -> run soft decode for scheduled SoftDecode rows
  -> merge result
  -> postprocess
```

当前 `materialize_hard_finish_rows(...)` 是直接回填 prepass 已经算好的 `hard_finish_lout`。新设计里可以改成：

```text
materialize_hard_scheduled_rows(...)
```

`materialize_clean_rows(...)` 处理 `HybridRowClass::Clean`，语义直接对齐 early-stop：不走 hard MUX，不消耗 `HIHO_ACTIVE_LIST[t]`，复用 `apply_row_early_stop_action(...)` 或同等 early-stop materialize helper 生成输出并写入 `merged_result`。

`materialize_hard_scheduled_rows(...)` 只遍历 hard MUX 选中的 row，对每行调用 hard executor，生成 `y2` 后写入 `merged_result`。

hard executor 失败时，正式实现选择直接 `throw`：

- classify-only 与 hard executor 应该严格一致。
- 如果 hard MUX 选中的 row 在 executor 阶段失败，说明分类逻辑、位置映射、BCH decode 或 parity 口径存在不一致。
- 这种情况不回退 soft path，直接抛异常，方便第一时间定位。

## 5. clean/free materialize 与 hard executor 具体逻辑

hard executor 输入：

```text
lin_vec
hybrid_class
params
```

输出：

```text
y2
bool produced
```

### 5.1 Clean/free materialize

```text
复用 EarlyStopAction 的 materialize 逻辑：
apply_row_early_stop_action(lin_vec, lch_vec, y2, params)
```

`Clean` 不属于 hard executor 资源范围，不经过 hard MUX，不消耗 `HIHO_ACTIVE_LIST[t]`。它按 early-stop 的免费路径直接 materialize。

### 5.2 hard executor: ParityOnly

```text
cw = hard_decision_bits_256(lin_vec)
cw[255] ^= 1
检查 hard_word_valid_256(cw)
合法则 materialize_hard_finish_lout(cw, lin_vec, p, y2)
```

只翻 overall parity。

### 5.3 hard executor: OneMain

```text
cw = hard_decision_bits_256(lin_vec)
s1 = syndrome S1
pos = hybrid_fast_gf_log(s1)
如果 pos 不合法，hard executor 失败
cw[pos] ^= 1
检查 hard_word_valid_256(cw)
合法则 materialize_hard_finish_lout(cw, lin_vec, p, y2)
```

只翻 1 个主体 bit。

### 5.4 hard executor: OneMainPlusParity

```text
cw = hard_decision_bits_256(lin_vec)
s1 = syndrome S1
pos = hybrid_fast_gf_log(s1)
如果 pos 不合法，hard executor 失败
cw[pos] ^= 1
cw[255] ^= 1
检查 hard_word_valid_256(cw)
合法则 materialize_hard_finish_lout(cw, lin_vec, p, y2)
```

翻 1 个主体 bit，再翻 overall parity。

### 5.5 hard executor: TwoMain

```text
cw = hard_decision_bits_256(lin_vec)
调用 bch_255_239_decode_hiho_cw_255(cw.data(), decoded.data(), &corrected_errors)
要求 decode 成功
要求 corrected_errors == 2
decoded[0..254] 拷回 cw[0..254]
recompute_overall_parity(&cw)
检查 hard_word_valid_256(cw)
合法则 materialize_hard_finish_lout(cw, lin_vec, p, y2)
```

注意：这段 BCH t=2 decode 必须从分类器里挪出来，只能在 hard MUX 选中后由 hard executor 执行。这样 `TwoMain` 才真正受 `n1` 的 hard 资源限制。

## 6. FriendS1S3WithS0 classify-only 逻辑

这一节只描述分类，不做翻转和硬解码。

输入：

```text
lin_vec
```

先做硬判：

```text
cw = hard_decision_bits_256(lin_vec)
```

计算：

```text
s0 = overall_parity_syndrome_256(cw)
syndromes = bch_255_239_syndromes_1_4_cw_255(cw.data())
s1 = syndromes[0]
s3 = syndromes[2]
```

### 6.1 S1 为 0

```text
if s1 == 0:
```

如果：

```text
s3 != 0
```

分类为：

```text
HardFail
```

如果：

```text
s3 == 0 && s0 == 0
```

分类为：

```text
Clean
```

如果：

```text
s3 == 0 && s0 == 1
```

分类为：

```text
ParityOnly
```

### 6.2 单主体错误

当：

```text
s1 != 0
```

计算：

```text
s1_cubed = s1^3
```

如果：

```text
s3 == s1_cubed
```

则进入单错分类。

计算：

```text
pos = hybrid_fast_gf_log(s1)
```

如果 `pos` 不合法：

```text
HardFail
```

如果 `pos` 合法，并且：

```text
s0 == 1
```

分类为：

```text
OneMain
```

含义：主体 255 位中有 1 个错误，整体错误数为奇数。

如果 `pos` 合法，并且：

```text
s0 == 0
```

分类为：

```text
OneMainPlusParity
```

含义：主体 255 位中有 1 个错误，同时 overall parity 位也错了，整体错误数为偶数。

### 6.3 两主体错误

如果：

```text
s1 != 0
s3 != s1^3
```

则不是单错模式，进入两错候选判定。

计算：

```text
numerator = s1^3 ^ s3
mu = numerator / s1^3
tr = Trace(mu)
```

如果：

```text
tr != 0
```

分类为：

```text
HardFail
```

如果：

```text
tr == 0 && s0 != 0
```

分类为：

```text
HardFail
```

这是 `FriendS1S3WithS0Classifier` 相比 `FriendS1S3Classifier` 的额外约束：两主体错误是偶数错，所以 overall parity syndrome 必须满足 `s0==0`。

如果：

```text
tr == 0 && s0 == 0
```

分类为：

```text
TwoMain
```

classify-only 阶段到这里结束。不能调用 `bch_255_239_decode_hiho_cw_255(...)`，不能翻 bit，不能生成 `y2`。

## 7. 新旧流程对比

当前流程：

```text
early-stop
  -> run_hybrid_prepass
       classifier 内部直接翻转 / BCH decode / 生成 y2
       deferred candidates 根据 SISO 缺口回收成 HardFinish
  -> soft MUX
  -> decode_tile_with_plan
       直接回填 cached hard_finish_lout
```

目标流程：

```text
early-stop
  -> classify-only
       只产出 HybridRowClass
  -> existing deferred/backfill priority
       复用当前 priority bucket，决定哪些 hard candidates 从 soft path 收回
  -> hard MUX
       只按 n1 做资源门控
  -> soft MUX
       只处理最终 SoftDecode
  -> decode_tile_with_plan
       hard executor 对 hard_scheduled_rows 生成 y2
       scheduled SoftDecode 跑 Chase
```

## 8. 建议落点

建议优先改以下文件：

- `src/rx/ofec/hybrid/hybrid_classifier.ipp`

  - 新增 classify-only helper。
  - `FriendS1S3WithS0Classifier` 的 TwoMain 分支不再调用 BCH t=2 decode。
- `include/newcode/ofec/hybrid/hybrid_classifier.hpp`

  - 声明 classify-only helper 和 hard executor helper。
- `src/rx/ofec/detail/ofec_tile_dispatch.ipp`

  - 将当前 `run_hybrid_prepass(...)` 拆成 classify-only、existing deferred/backfill priority、hard MUX 三段。
  - 继续保留 `rebuild_soft_candidate_rows(...)` 和现有 soft MUX。
- `src/rx/ofec/detail/ofec_tile_decode.ipp`

  - 将当前 `materialize_hard_finish_rows(...)` 改成执行型 materialize。
  - 在 `decode_tile_with_plan(...)` 内对 hard MUX 选中的 row 调 hard executor。
- `include/newcode/params.hpp`

  - 增加 hard 资源参数：

```cpp
std::vector<int> HIHO_ACTIVE_LIST = {32, 32, 32, 32, 32, 32};
```

- `HIHO_ACTIVE_LIST[t]` 控制第 `t` 个 tile 的 hard MUX / hard executor 资源数量。
- app 入口，例如 `apps/ofec_ber_window_probe.cpp`、`apps/ofec_two_stream_shared_sweep.cpp`

  - 在 `kSisoActiveList` 附近新增外部配置：

```cpp
static const std::vector<int> kHiHoActiveList = {32, 32, 32, 16, 8, 4};
```

- build config / build params 时把它写入：

```cpp
config.hiho_active_list = kHiHoActiveList;
params.HIHO_ACTIVE_LIST = kHiHoActiveList;
```

- 如果中间配置结构暂时没有 `hiho_active_list` 字段，也需要和 `siso_active_list` 同步补齐。

## 9. 统计建议

为了验证新 hard MUX 行为，建议新增统计：

```text
rows_hard_candidate
rows_hard_reclaim_by_priority
rows_hard_scheduled
rows_hard_not_scheduled
rows_hard_executor_success
rows_hard_executor_fail
rows_soft_candidate_after_backfill_priority
rows_soft_candidate_before_soft_mux
hiho_active_for_tile
hiho_active_list
```

现有统计：

```text
rows_hard_finish
rows_need_siso_before_mux
rows_unscheduled
class_two_main_count
deferred_reclaimed_to_hard_finish_count
```

可以继续保留，但含义需要同步更新：

- `rows_hard_finish`：hard executor 成功并最终进入 `HardFinish` 的行数。
- `rows_need_siso_before_mux`：deferred/backfill priority、hard MUX 之后，剩余要进 soft MUX 的行数。
- `deferred_reclaimed_to_hard_finish_count`：如果保留 deferred 命名，它应表示 hard MUX 选中并执行成功的 hard candidates 数量。

## 10. 第一版实现原则

第一版建议尽量少改外层流程：

```text
build_tile_dispatch_plan
  -> classify-only
  -> existing deferred/backfill priority
  -> hard MUX
  -> rebuild_soft_candidate_rows

run_mux_on_soft_candidates
  -> 保持现有 soft MUX 逻辑

decode_tile_with_plan
  -> EarlyStopAction materialize
  -> Clean/free materialize
  -> hard executor materialize
  -> scheduled SoftDecode run Chase
```

这样 single-stream 和 two-stream shared 都可以继续复用 `TileDispatchPlan`，避免在两个入口里写两套调度逻辑。
