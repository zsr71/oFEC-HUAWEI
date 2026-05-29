# two-stream shared 中 A/B 交换后结果不对称问题说明

## 1. 问题背景

当前 two-stream shared 主流程入口是：

- [apps/ofec_two_stream_shared_decoder.cpp](/home/zsr71/projects/newcode_two_stream_shared_chase/apps/ofec_two_stream_shared_decoder.cpp)

底层共享主链路在：

- [src/two_stream_shared/two_stream_shared_runner.cpp](/home/zsr71/projects/newcode_two_stream_shared_chase/src/two_stream_shared/two_stream_shared_runner.cpp)

当前主流程的关键语义是：

1. 两路原生独立数据流各自生成发端码字
2. 各自经过独立信道并完成量化
3. 各自先做 `prepare_tile_inputs(...)`
4. 从 `detect_tile_early_stop(...)` 开始，把两边 tile 合并成一个 shared 64-row 域
5. 在这个 merged 域上统一做：
   - early-stop
   - `build_tile_dispatch_plan(...)`
   - `run_mux_on_soft_candidates(...)`
   - `decode_tile_with_plan(...)`
6. 最后再把结果切回 A/B 两路，分别写回并统计 BER

这个结构的目标，不是做两路各自独立的 tile-local 调度，而是让两路在合并后的 64 code 域里统一竞争 shared 资源。

---

## 2. 当前观察到的问题

在当前配置下，交换 Stream A / Stream B 的输入 seed 之后，前端结果会严格对调，但 post-FEC 结果不会只做简单对调，而是会发生明显变化。

当前相关配置包括：

- `kDecoderName = "two_stream_shared_chase_baseline"`
- `kSisoActiveList = {64, 64, 64, 32, 16, 8}`
- `kMuxGroupG = 1`
- `kMuxSchedulingMode = 0`
- `kMuxPriorityRule = 0`
- `kMuxEnableReconfig = false`
- `kHybridEnable = true`

位置见：

- [apps/ofec_two_stream_shared_decoder.cpp](/home/zsr71/projects/newcode_two_stream_shared_chase/apps/ofec_two_stream_shared_decoder.cpp:40)
- [apps/ofec_two_stream_shared_decoder.cpp](/home/zsr71/projects/newcode_two_stream_shared_chase/apps/ofec_two_stream_shared_decoder.cpp:68)
- [apps/ofec_two_stream_shared_decoder.cpp](/home/zsr71/projects/newcode_two_stream_shared_chase/apps/ofec_two_stream_shared_decoder.cpp:71)

---

## 3. 复现实验现象

### 3.1 第一次运行

- Stream A seeds:
  - bitgen=`20260319`
  - channel=`3182026`
- Stream B seeds:
  - bitgen=`20260320`
  - channel=`3182027`

结果摘要：

- A `Pre-FEC BER = 0.022304`
- B `Pre-FEC BER = 0.0222865`
- A `Post-FEC BER = 0.00010389338`，错误数 `1437`
- B `Post-FEC BER = 0.00013411428`，错误数 `1855`

shared core 摘要：

- total rows = `1622400`
- produced = `1622129`
- failed = `271`
- stream A: produced = `811084`, failed = `116`
- stream B: produced = `811045`, failed = `155`

### 3.2 第二次运行

交换两路 seed：

- Stream A seeds:
  - bitgen=`20260320`
  - channel=`3182027`
- Stream B seeds:
  - bitgen=`20260319`
  - channel=`3182026`

结果摘要：

- A `Pre-FEC BER = 0.0222865`
- B `Pre-FEC BER = 0.022304`
- A `Post-FEC BER = 0.00018226528`，错误数 `2521`
- B `Post-FEC BER = 7.8588797e-05`，错误数 `1087`

shared core 摘要：

- total rows = `1622400`
- produced = `1622086`
- failed = `314`
- stream A: produced = `811004`, failed = `196`
- stream B: produced = `811082`, failed = `118`

---

## 4. 为什么这个现象可疑

## 4.1 前端结果是严格对调的

两次运行里，下面这些前端结果会随着 seed 交换而严格互换：

- `Pre-FEC BER`
- `Pre-FEC (quantized hard) BER`
- shared quant clip
- A/B 各自的量化饱和统计

这说明：

- 发端比特生成是稳定的
- 独立信道是稳定的
- 量化口径是稳定的
- A/B seed 的交换是正确生效的

也就是说，问题不在前端。

## 4.2 post-FEC 结果没有只做标签互换

如果 shared 主流程对 A/B 是完全对称的，那么交换两路输入之后，理论上应接近：

- 第一次 A 的 post-FEC 结果，变成第二次 B 的结果
- 第一次 B 的 post-FEC 结果，变成第二次 A 的结果

但实际并不是这样。

同一组 seed 在不同槽位下的结果是：

- `20260319/3182026`
  - 放在 A 时：`1437` errors
  - 放在 B 时：`1087` errors

- `20260320/3182027`
  - 放在 B 时：`1855` errors
  - 放在 A 时：`2521` errors

这说明：

> 同一条物理输入流，仅仅因为它被放在 A 槽还是 B 槽，最终 post-FEC 结果就明显不同。

这不是理想对称 shared 结构应有的现象。

## 4.3 这不是早期那种“明显饿死一边”的问题

从 shared core 统计看，当前实现已经不像更早版本那样出现：

- A 全部 produced
- B 几乎全部 failed

这种极端偏流。

现在的 produced/failed 已经接近平衡，这说明：

- 交织 merged 顺序 `A0,B0,A1,B1,...` 的修正是有效的
- 但系统仍然没有达到“交换 A/B 后结果只做标签互换”的程度

因此当前问题更准确的描述不是：

- “某一路完全饿死”

而是：

- “系统仍然存在顺序敏感或槽位敏感”

---

## 5. 当前最可能的问题性质

当前最可能的问题不是：

- 前端错了
- seed 没生效
- 写回完全错位

而是：

> 合并成 64-row shared tile 之后，后续调度/预算裁剪路径仍然依赖 merged row 的先后顺序。

这会导致：

- 同一组物理码字，如果落在不同 merged row 位置
- 即使其“质量”本身没有变化
- 也可能在 MUX / scheduling 阶段遭遇不同命运

最终反映到：

- shared core 的 produced/failed 数不同
- post-FEC BER 不再只随物理输入而定，而是也受 A/B 槽位影响

---

## 6. 当前最可疑的顺序敏感链条

## 6.1 merged row 顺序本身是显式构造出来的

在 [src/two_stream_shared/two_stream_shared_runner.cpp](/home/zsr71/projects/newcode_two_stream_shared_chase/src/two_stream_shared/two_stream_shared_runner.cpp:333) 的 `merge_prepared_tiles(...)` 里，当前 merged 顺序是：

- `A0, B0, A1, B1, A2, B2, ...`

对应关键位置：

- [src/two_stream_shared/two_stream_shared_runner.cpp](/home/zsr71/projects/newcode_two_stream_shared_chase/src/two_stream_shared/two_stream_shared_runner.cpp:392)
- [src/two_stream_shared/two_stream_shared_runner.cpp](/home/zsr71/projects/newcode_two_stream_shared_chase/src/two_stream_shared/two_stream_shared_runner.cpp:395)

这一步本身不是 bug。  
它的作用是避免早期 `AAAA...BBBB...` 排布造成的明显偏流。

但它只是“换了一种顺序”，并没有提供真正与顺序无关的全局排序依据。

## 6.2 dispatch plan 的 soft candidate 顺序就是 row 顺序

在 [src/rx/ofec/detail/ofec_tile_dispatch.ipp](/home/zsr71/projects/newcode_two_stream_shared_chase/src/rx/ofec/detail/ofec_tile_dispatch.ipp:298) 的 `build_tile_dispatch_plan(...)` 中：

- 会顺着 row `0..N-1` 扫
- 未命中 early-stop 的行先标成 `SoftDecode`
- 并按这个顺序 `push_back` 到 `soft_candidate_rows`

关键位置：

- [src/rx/ofec/detail/ofec_tile_dispatch.ipp](/home/zsr71/projects/newcode_two_stream_shared_chase/src/rx/ofec/detail/ofec_tile_dispatch.ipp:323)
- [src/rx/ofec/detail/ofec_tile_dispatch.ipp](/home/zsr71/projects/newcode_two_stream_shared_chase/src/rx/ofec/detail/ofec_tile_dispatch.ipp:337)

因此：

- `soft_candidate_rows` 的天然顺序
- 就是 merged row 顺序

## 6.3 MUX 输入状态也是按 row 顺序构造的

在 [src/rx/ofec/detail/ofec_tile_dispatch.ipp](/home/zsr71/projects/newcode_two_stream_shared_chase/src/rx/ofec/detail/ofec_tile_dispatch.ipp:365) 的 `run_mux_on_soft_candidates(...)` 中：

- 会按 row 顺序生成 `mux_candidate_state`
- `SoftDecode` 行标为 `NeedSiso`
- 其他行标为 `EarlyStopped`

关键位置：

- [src/rx/ofec/detail/ofec_tile_dispatch.ipp](/home/zsr71/projects/newcode_two_stream_shared_chase/src/rx/ofec/detail/ofec_tile_dispatch.ipp:379)
- [src/rx/ofec/detail/ofec_tile_dispatch.ipp](/home/zsr71/projects/newcode_two_stream_shared_chase/src/rx/ofec/detail/ofec_tile_dispatch.ipp:382)

这意味着：

- merged row 的先后顺序
- 会直接传进 MUX 预算裁剪器

## 6.4 当前 `MUX_SCHEDULING_MODE=0` 时，本质是顺序裁剪

当前 app 配置里：

- `kMuxSchedulingMode = 0`
- `kMuxGroupG = 1`
- `kMuxEnableReconfig = false`

位置：

- [apps/ofec_two_stream_shared_decoder.cpp](/home/zsr71/projects/newcode_two_stream_shared_chase/apps/ofec_two_stream_shared_decoder.cpp:71)
- [apps/ofec_two_stream_shared_decoder.cpp](/home/zsr71/projects/newcode_two_stream_shared_chase/apps/ofec_two_stream_shared_decoder.cpp:72)
- [apps/ofec_two_stream_shared_decoder.cpp](/home/zsr71/projects/newcode_two_stream_shared_chase/apps/ofec_two_stream_shared_decoder.cpp:74)

在这种配置下，`run_mux_on_soft_candidates(...)` 会走普通 grouped budget 路径。

而 [src/rx/ofec/mux/mux_siso_budget.cpp](/home/zsr71/projects/newcode_two_stream_shared_chase/src/rx/ofec/mux/mux_siso_budget.cpp:76) 的 `apply_siso_budget_g1(...)` 逻辑是：

1. 按下标顺序收集所有 `NeedSiso`
2. 保留前 `budget` 个
3. 后面的全部标成 `Unscheduled`

关键位置：

- [src/rx/ofec/mux/mux_siso_budget.cpp](/home/zsr71/projects/newcode_two_stream_shared_chase/src/rx/ofec/mux/mux_siso_budget.cpp:82)
- [src/rx/ofec/mux/mux_siso_budget.cpp](/home/zsr71/projects/newcode_two_stream_shared_chase/src/rx/ofec/mux/mux_siso_budget.cpp:95)

这本质上就是：

> 保前裁后

如果输入顺序改变，谁被保留、谁被裁掉，就可能变化。

## 6.5 priority 路径也仍然保留 tie 下的顺序偏置

即使后面切到 `MUX_SCHEDULING_MODE=1`，问题也不一定完全消失。

原因是 [src/rx/ofec/mux/mux_siso_budget.cpp](/home/zsr71/projects/newcode_two_stream_shared_chase/src/rx/ofec/mux/mux_siso_budget.cpp:35) 与
[src/rx/ofec/mux/mux_group_budget.cpp](/home/zsr71/projects/newcode_two_stream_shared_chase/src/rx/ofec/mux/mux_group_budget.cpp:84)
里的 `trim_need_indices_by_priority(...)` 使用的是：

- `std::stable_sort(...)`

比较器只看：

- `early_stop_priority_score(...)`

关键位置：

- [src/rx/ofec/mux/mux_siso_budget.cpp](/home/zsr71/projects/newcode_two_stream_shared_chase/src/rx/ofec/mux/mux_siso_budget.cpp:48)
- [src/rx/ofec/mux/mux_group_budget.cpp](/home/zsr71/projects/newcode_two_stream_shared_chase/src/rx/ofec/mux/mux_group_budget.cpp:94)

这意味着：

- 如果很多 row 的 priority score 相同
- `stable_sort` 会保留它们原本的相对顺序
- 原本的相对顺序仍然是 merged row 顺序

因此：

- priority MUX 只能减弱部分顺序问题
- 但不能自动消灭 tie 情况下的槽位偏置

## 6.6 slice back 暂时不像第一嫌疑点

在 [src/two_stream_shared/two_stream_shared_runner.cpp](/home/zsr71/projects/newcode_two_stream_shared_chase/src/two_stream_shared/two_stream_shared_runner.cpp:499) 之后，当前逻辑是：

- 根据 `slice.merged_rows` 把 `decoder_res` 切回 A/B
- 再分别调用 `writeback_tile(...)`

从当前实现看，这一层更像是：

- 按既定 merged row 映射做搬运

而不像前面的 MUX 那样天然带有“保前裁后”的顺序偏置。

因此当前优先怀疑：

- merged 后的调度与预算裁剪

而不是：

- slice back 本身先写错

---

## 7. 当前问题的更准确表述

当前问题不应简单表述为：

- “A/B 不公平”

因为与早期版本相比，现在已经没有��现“一边全 produced，一边几乎全 failed”的极端失衡。

更准确的表述应该是：

> 当前 two-stream shared 实现仍然存在 A/B 槽位敏感性。  
> 即同一组物理输入流，在交换到另一侧后，post-FEC 结果不会只做简单标签互换。

或者更技术一点地说：

> 当前 merged-64-row shared path 尚未达到 A/B 交换不变性。

---

## 8. 这个问题对当前结果解读的影响

这意味着当前版本的 BER 结果可以用于：

- 看当前实现是否基本跑通
- 看量级是否合理
- 看是否存在极端退化

但还不能直接当作“理想对称 shared 结构”的 BER 基线。

因为当前扫出来的 BER 曲线实际上代表的是：

- 当前这版实现
- 加上当前 merged 顺序与 MUX 顺序语义

而不是：

- 一个对 A/B 完全对称的 shared 资源竞争系统

---

## 9. 后续排查建议

当前最值得优先验证的，不是继续跑更大样本，而是做小规模、可定位的顺序敏感诊断。

建议优先检查：

1. merged row index
2. 每行是否进入 `SoftDecode`
3. 每行是否被 MUX 保留或裁掉
4. 每行最终是否 `produced`
5. 同一物理 row 在 A 槽和 B 槽下，是否只是因为 merged index 改变，就进入了不同调度结果

最小可行的诊断点应集中在：

- `merge_prepared_tiles(...)`
- `build_tile_dispatch_plan(...)`
- `run_mux_on_soft_candidates(...)`

而不是先从 BER 汇总层面继续猜。

---

## 10. 总结

当前 two-stream shared 主流程已经解决了最早的明显偏流问题，但仍未达到理想的 A/B 对称性。

当前最可能的问题根因是：

- merged 64-row 域中的 row 顺序
- 被后续 dispatch / MUX 预算裁剪路径继续当成隐式 tie-break
- 从而让“交换 A/B 输入”不再等价于“只交换标签”

因此当前问题的核心不是：

- 前端不稳定
- seed 不对
- BER 随机波动过大

而是：

> shared 调度路径中仍然存在顺序敏感语义，这种语义已经足以影响最终 post-FEC BER。
