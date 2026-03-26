# 早停动作与公共 Alpha 的两种改造方案比较

## 1. 背景

当前解码流程里，早停命中的行虽然不再走正常 Chase，但它生成的 `y2` 仍然会进入 tile 级公共后处理。

这一步公共后处理会对所有 `produced=true` 的行统一乘当前 tile 的 `ALPHA`。

因此，当前早停动作 1 的实际写回值不是：

`final = (lin - lch) + sign(lin) * beta`

而是：

`final = ALPHA * ((lin - lch) + sign(lin) * beta)`

如果某个 tile 的 `ALPHA < 1`，那么早停动作即使本来想“增强当前方向”，最终也可能被统一缩小，甚至写回值比原来的 `lin-lch` 还弱。

这和“早停命中后应当让该位更可靠”这个目标存在冲突。

---

## 2. 方案一：早停动作输出跳过公共 Alpha

### 2.1 核心思路

保持早停动作函数本身不变：

`y2 = (lin - lch) + sign(lin) * beta`

但在 tile 公共后处理阶段，对“来自 early-stop 的行”不再乘 `ALPHA`。

这样最终写回值就是：

`final = y2`

而正常 Chase 行仍然保持：

`final = ALPHA * y2`

### 2.2 语义

这个方案的语义最直接：

- 正常 Chase 输出仍然受公共 `ALPHA` 控制
- 早停动作输出被视为“已经是最终希望写回的外信息”
- 不再对它做额外缩放

### 2.3 优点

- 语义清楚，最符合“早停命中后应该增强可靠度”的直觉
- 日志、trace、文档都容易解释
- 后续如果还要增加别的早停动作，也可以明确决定它们是否跳过 `ALPHA`

### 2.4 风险和代价

- 需要让 tile 级后处理知道“哪些 produced 行来自 early-stop”
- 也就是说，现有 `produced_rows` 这个布尔标记不够用了，还要再传一层“行来源”
- 改动面会比表面看起来稍大，因为它会穿过 decoder core 到 tile 后处理

### 2.5 可能需要改动的代码位置

最核心的是这几处：

- [ofec_row_decoder_core.cpp](/home/zsr71/projects/newcode/src/rx/ofec/ofec_row_decoder_core.cpp)
  - 当前这里只记录某行有没有 `produced`
  - 需要再额外记录这行是不是通过 early-stop action 产生的

- [DecoderCoreResult 定义处](/home/zsr71/projects/newcode/src/rx/ofec/ofec_row_decoder_core.cpp)
  - 需要增加一个与 `produced_rows` 同长度的标记，例如：
  - `from_early_stop_rows`
  - 或者更通用一点的 `row_output_kind`

- [ofec_tile_decode.ipp](/home/zsr71/projects/newcode/src/rx/ofec/detail/ofec_tile_decode.ipp)
  - 当前这里对所有 `produced_rows[r]` 都统一乘 `ALPHA`
  - 需要改成：
  - 正常 Chase 行乘 `ALPHA`
  - early-stop 行跳过这一步

- [ofec_tile_writeback.ipp](/home/zsr71/projects/newcode/src/rx/ofec/detail/ofec_tile_writeback.ipp)
  - 通常不用改公式
  - 但最好重新核对 trace/log 输出，确认记录的是“跳过 Alpha 后的最终写回值”

### 2.6 复杂度判断

这个方案不是改一行就完事。

真正的工作量主要不在“乘不乘 Alpha”本身，而在：

- 给 row decoder 输出结构补充“来源标记”
- 把这个标记一路传到 tile 公共后处理

所以它是结构上更正确，但改动会稍大一些的方案。

---

## 3. 方案二：新增一个早停动作，在动作内部先除以 Alpha

### 3.1 核心思路

保持 tile 公共后处理完全不变，仍然对所有 `produced` 行统一乘 `ALPHA`。

但是新增一个早停动作，例如“动作 4”，内部先算：

`y2_raw = (lin - lch) + sign(lin) * beta`

然后输出：

`y2 = y2_raw / ALPHA`

后面公共后处理再乘回去：

`final = ALPHA * y2 = ALPHA * (y2_raw / ALPHA) = y2_raw`

于是最终等效于“早停动作输出不乘 Alpha”。

### 3.2 语义

这个方案的本质是：

- 不改公共后处理
- 在动作内部做一次预补偿

所以动作函数输出的并不是“直观上的最终外信息”，而是“为了通过公共后处理后得到目标值而预先补偿过的值”。

### 3.3 优点

- 改动范围更小
- 不需要给 `DecoderCoreResult` 增加额外的行来源标记
- 不需要改 tile 级公共后处理的分支结构

### 3.4 风险和代价

- 语义比较绕
- 看代码时，动作输出值不是最终真实想写回的值
- 如果后面又插入别的公共处理环节，这个“先除再乘”的等效关系更容易被破坏
- 还要处理 `ALPHA=0` 或极小值的边界

### 3.5 可能需要改动的代码位置

这个方案主要只动早停动作分发这条链：

- [row_early_stop_action.ipp](/home/zsr71/projects/newcode/src/rx/ofec/earlystop/row_early_stop_action.ipp)
  - 增加一个新的 `action_mode`

- [row_early_stop_process_1.ipp](/home/zsr71/projects/newcode/src/rx/ofec/earlystop/row_early_stop_process_1.ipp)
  - 不建议直接改原动作 1
  - 更合理的是新增一个新动作文件，例如：
  - `row_early_stop_process_4.ipp`

- [params.hpp](/home/zsr71/projects/newcode/include/newcode/params.hpp)
  - 放开新的 `EARLY_STOP_ACTION_MODE` 合法值

- [ofec_single_runner.hpp](/home/zsr71/projects/newcode/include/newcode/ofec_single_runner.hpp)
- [ofec_single_params.cpp](/home/zsr71/projects/newcode/src/ofec_single/ofec_single_params.cpp)
- [apps/ofec_single.cpp](/home/zsr71/projects/newcode/apps/ofec_single.cpp)
  - 顶层参数和注释需要同步补上新动作模式

- [ofec_sweep_runner.hpp](/home/zsr71/projects/newcode/include/newcode/ofec_sweep_runner.hpp)
- [ofec_sweep_scenarios.cpp](/home/zsr71/projects/newcode/src/ofec_sweep/ofec_sweep_scenarios.cpp)
- [apps/ofec_sweep.cpp](/home/zsr71/projects/newcode/apps/ofec_sweep.cpp)
  - 如果 sweep 里也要扫描新动作，需要同步接线

### 3.6 复杂度判断

这个方案实现速度会更快，改动点更集中。

从纯编码工作量看，它通常比方案一小。

---

## 4. 两种方案的直接对比

### 4.1 结构正确性

- 方案一更像真正的结构修正
- 方案二更像为了适配现有公共后处理而做的补偿

### 4.2 改动范围

- 方案一改动更广
- 方案二改动更集中

### 4.3 代码可解释性

- 方案一更好解释
- 方案二更容易让人困惑：“为什么动作里要先除 Alpha”

### 4.4 实验推进速度

- 如果目的是尽快验证“早停动作不乘 Alpha 会不会更好”，方案二更快
- 如果目的是后面长期保留这套机制，方案一更稳

---

## 5. 建议

如果目标是：

### 5.1 先快速做实验验证

优先考虑方案二。

原因是：

- 改动少
- 上手快
- 很适合先验证“去掉早停动作上的 Alpha 缩放”到底有没有收益

### 5.2 后面打算长期保留

优先考虑方案一。

原因是：

- 语义最清楚
- 后续维护成本更低
- 也更适合继续扩展更多早停动作

---

## 6. 我当前的判断

如果只看你现在这个阶段，我更建议：

1. 先用方案二快速做实验
2. 如果结果确实明显更合理，再决定是否升级到方案一

因为你现在最关心的是：

- 早停动作是不是被 `ALPHA` 缩弱了
- 去掉这个缩弱后 BER 会不会改善

这个问题先用小改动验证，性价比最高。
