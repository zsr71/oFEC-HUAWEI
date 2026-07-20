# `ofec_two_stream_shared_sweep` CSV 列说明

本文说明
[apps/ofec_two_stream_shared_sweep.cpp](/home/zsr71/projects/newcode_two_stream_shared_chase/apps/ofec_two_stream_shared_sweep.cpp)
运行后生成的结果 CSV 中，各列分别表示什么。

当前这个 sweep app 的定位，是对 two-stream shared 主流程做低 BER 扫描，因此 CSV 的重点是：

- 每个 `Eb/N0` 点上，A/B 两路各自的 BER 统计
- 当前这次扫描对应的关键参数快照

第一版不会在 CSV 里聚合 shared 内部可观测性，例如 shared core 的 produced/failed、量化饱和数等。这些仍建议回到单点入口
[apps/ofec_two_stream_shared_decoder.cpp](/home/zsr71/projects/newcode_two_stream_shared_chase/apps/ofec_two_stream_shared_decoder.cpp)
里看。

---

## 1. 一行代表什么

这个 CSV 里每一行代表：

- 一次 sweep 运行中的一个 `Eb/N0` 点
- 并且这一行已经是该点若干个 chunk 聚合后的最终结果

也就是说，它不是单个 chunk 的结果，而是：

- 多个 chunk 累加后的 A 路统计
- 多个 chunk 累加后的 B 路统计
- 再附带这次点级结果对应的参数快照

---

## 2. 基础标识列

### `timestamp`

- 写入这一行 CSV 时的时间戳
- 使用 `utils::now_stamp()` 生成
- 主要用于区分文件里各行的实际落盘时间

### `run_id`

- 这次 sweep 运行的统一编号
- 同一次运行生成的所有行，`run_id` 都相同
- 可以用它来区分不同批次的 sweep 结果

### `ebn0_db`

- 当前这一行对应的 `Eb/N0`，单位是 dB
- 这是 sweep 的主扫描轴

---

## 3. Stream A 结果列

### `pre_ber_a`

- A 路聚合后的 pre-FEC BER
- 分子是 `pre_errs_a`
- 分母是 `pre_total_a`

### `pre_errs_a`

- A 路 pre-FEC 累计错误比特数

### `pre_total_a`

- A 路 pre-FEC 累计比较的总比特数

### `pre_quant_ber_a`

- A 路“量化后立刻做硬判”的 BER
- 对应你平时看到的 `Pre-FEC (quantized hard) BER`
- 如果当前流程没有生成这项统计，则该列为空

### `pre_quant_errs_a`

- A 路量化硬判口径下的累计错误比特数
- 如果没有该统计，则为空

### `pre_quant_total_a`

- A 路量化硬判口径下的累计比较总比特数
- 如果没有该统计，则为空

### `post_ber_a`

- A 路聚合后的 post-FEC BER
- 分子是 `post_errs_a`
- 分母是 `post_total_a`

### `post_errs_a`

- A 路 post-FEC 累计错误比特数

### `post_total_a`

- A 路 post-FEC 累计比较总比特数

---

## 4. Stream B 结果列

### `pre_ber_b`

- B 路聚合后的 pre-FEC BER

### `pre_errs_b`

- B 路 pre-FEC 累计错误比特数

### `pre_total_b`

- B 路 pre-FEC 累计比较总比特数

### `pre_quant_ber_b`

- B 路量化硬判口径下的 BER
- 如果当前流程没有生成该统计，则该列为空

### `pre_quant_errs_b`

- B 路量化硬判口径下的累计错误比特数
- 如果没有该统计，则为空

### `pre_quant_total_b`

- B 路量化硬判口径下的累计比较总比特数
- 如果没有该统计，则为空

### `post_ber_b`

- B 路聚合后的 post-FEC BER

### `post_errs_b`

- B 路 post-FEC 累计错误比特数

### `post_total_b`

- B 路 post-FEC 累计比较总比特数

---

## 5. 聚合规模与停止条件列

### `chunk_num_info_bits`

- 每个 chunk、每一路发端生成的信息比特数
- 对应代码常量 `kChunkNumInfoBits`
- 值越大，单个 chunk 更重，但统计收敛更快

### `chunks_completed`

- 当前这个 `Eb/N0` 点最终一共跑完了多少个 chunk
- 注意这是点级聚合后的 chunk 数，不是全局总 chunk 数

### `target_post_errors`

- post-FEC 错误数停止门限
- 对应 `kTargetPostErrors`
- 当前 sweep 是 A/B 两路都满足各自停止条件之后，这个点才停止

### `max_post_fec_total_bits`

- post-FEC 总比较比特数上限
- 对应 `kMaxPostFecTotalBits`
- 用来避免极低 BER 点无限跑下去

### `confidence_level`

- 零错上界停止时使用的置信水平
- 对应 `kConfidenceLevel`
- 例如 `0.95` 表示 95% 置信水平

---

## 6. A/B 路停止原因与零错上界列

这几列是为了说明：

- 这一点为什么停
- 如果 post-FEC 一直没有错误，那么当前记录到底是“真 0”，还是“零错上界意义下的提前停止”

### `stop_reason_a`

- A 路停止原因
- 可能值包括：
  - `target_post_errors`
  - `max_post_fec_total_bits`
  - `zero_error_upper_bound`
  - `no_data`

其中：

- `target_post_errors` 表示 A 路累计错误数已达到目标
- `max_post_fec_total_bits` 表示 A 路累计总比特已达到上限
- `zero_error_upper_bound` 表示 A 路虽然仍是零错，但按当前置信水平算出来的 BER 上界已经足够低，因此允许提前停止

### `post_ber_is_upper_bound_a`

- `1` 表示这一行的 A 路 post-FEC 结果应理解为“上界意义下的零错结果”
- `0` 表示不是这种情况

更具体地说：

- 如果 A 路 `post_errs_a == 0`
- 且停止原因是 `zero_error_upper_bound`
- 那么这里会写 `1`

### `post_ber_upper_bound_a`

- A 路在零错条件下，根据 `confidence_level` 计算出来的 BER 上界
- 只有在：
  - `post_errs_a == 0`
  - 并且 `post_total_a > 0`
  时才有意义
- 否则可能为空

### `stop_reason_b`

- B 路停止原因
- 含义与 `stop_reason_a` 相同

### `post_ber_is_upper_bound_b`

- B 路是否是“零错上界意义下的结果”
- 含义与 A 路对应列相同

### `post_ber_upper_bound_b`

- B 路零错条件下的 BER 上界
- 含义与 A 路对应列相同

---

## 7. 运行耗时列

### `elapsed_seconds`

- 当前这个 `Eb/N0` 点从开始到结束的总耗时，单位是秒
- 是点级耗时，不是整个 sweep 总耗时

---

## 8. 基础 seed 快照列

这几列保存的是“基础 seed”，不是每个 chunk 派生出来的最终 seed。

### `bitgen_seed_a_base`

- A 路基础信息比特 seed
- 对应 `kBitgenSeedA`

### `channel_seed_a_base`

- A 路基础信道 seed
- 对应 `kChannelSeedA`

### `bitgen_seed_b_base`

- B 路基础信息比特 seed
- 对应 `kBitgenSeedB`

### `channel_seed_b_base`

- B 路基础信道 seed
- 对应 `kChannelSeedB`

实际每个 chunk 会在这些基础 seed 的基础上，再结合 `chunk_index` 派生出 chunk 级 seed。

---

## 9. alpha / beta / SISO 预算列

### `alpha_list`

- 当前运行使用的每个 tile 的 alpha 显式列表
- 用 `|` 分隔保存
- 例如：`0.428571|0.447738|...`

### `beta_list`

- 当前运行使用的每个 tile 的 beta 显式列表
- 同样用 `|` 分隔保存

### `siso_active_list`

- 当前运行使用的每个 tile 的 shared SISO 预算列表
- 用 `|` 分隔保存
- 在 two-stream shared 里，它表示 merged 64-code 域上每个 tile 允许参与 soft path 的预算

---

## 10. Chase core 参数列

### `chase_L`

- Chase L 参数
- 表示选取多少个最不可靠位参与 Chase 候选展开

### `chase_n_test`

- 当前实际使用的 Chase 测试 pattern 数
- 如果代码里 `kChaseNTestOverride < 0`，这里会写实际采用的 `2^L`
- 如果设置了 override，则这里写 override 后的值

### `chase_topk_keep`

- `topk/pruned` 类路径保留的候选数
- 这个 baseline sweep 里即便不一定直接用到，也会一并记录

### `chase_group_minima_bits`

- `group_minima` 类路径的分组 bit 数
- 同样作为参数快照保留下来

---

## 11. MUX 调度参数列

### `mux_group_g`

- MUX 分组数
- `1` 表示全局共享预算池

### `mux_scheduling_mode`

- MUX 调度模式
- 当前常见值：
  - `0`：legacy 顺序裁剪
  - `1`：按 early-stop 细节排序再裁剪

### `mux_early_stop_priority_rule`

- 在新调度模式下使用的优先级规则
- 当前常见值：
  - `0`
  - `1`

具体语义建议结合当前 MUX 设计文档一起看。

### `mux_enable_reconfig`

- 是否启用 reconfig 调度
- `1` 表示启用
- `0` 表示关闭

---

## 12. Hybrid 参数列

### `hybrid_enable`

- hybrid 总开关
- `1` 表示启用
- `0` 表示关闭

### `hybrid_enable_list`

- 每个 tile 的 hybrid 开关覆盖列表
- 用 `|` 分隔

### `hybrid_classifier_mode`

- hybrid 分类器模式名称
- 例如：
  - `legacy_hard_decode`
  - `repo_fast_classifier`
  - `friend_s1s3_classifier`
  - `friend_s1s3_with_s0_classifier`

### `hybrid_siso_backfill_mode`

- hybrid 的 SISO 回填模式名称
- 例如：
  - `disabled`
  - `two_error_only`
  - `one_and_two_error_priority`

### `hybrid_normalize_soft_only`

- 是否只对 soft rows 做归一化
- `1` 表示是
- `0` 表示否

---

## 13. Early-stop 参数列

### `early_stop_enable`

- early-stop 总开关
- `1` 表示启用
- `0` 表示关闭

### `early_stop_enable_list`

- 每个 tile 的 early-stop 开关覆盖列表
- 用 `|` 分隔

### `early_stop_condition_mode`

- early-stop 条件模式
- 当前常见值：
  - `1`：v1
  - `2`：v2

### `early_stop_action_mode`

- early-stop 动作模式
- 当前代码里会直接以数字保存

### `early_stop_bind_group_size`

- early-stop 组绑定大小
- `1` 表示逐 row 判断
- `4` 表示按四个 row 绑定判断

### `early_stop_cond_v1_require_bch`

- 条件 v1 是否要求 BCH syndrome 通过
- `1` 表示要求
- `0` 表示不要求

### `early_stop_cond_v1_require_overall`

- 条件 v1 是否要求 overall parity 通过
- `1` 表示要求
- `0` 表示不要求

### `early_stop_v2_llr_abs_threshold`

- 条件 v2 里用于定义“不可靠位”的 `|LLR|` 阈值

### `early_stop_v2_max_unreliable_bits`

- 条件 v2 允许的不可靠 bit 数上限

### `early_stop_cond_v2_include_overall`

- 条件 v2 是否把 overall parity bit 纳入统计
- `1` 表示纳入
- `0` 表示不纳入

---

## 14. 实际使用时怎么读这个 CSV

如果你的目标是看 BER 曲线，通常最先关心的是：

- `ebn0_db`
- `post_ber_a`
- `post_ber_b`
- `post_errs_a`
- `post_errs_b`
- `post_total_a`
- `post_total_b`

如果你的目标是判断某个点是否“统计得够不够稳”，重点看：

- `chunks_completed`
- `stop_reason_a`
- `stop_reason_b`
- `post_ber_is_upper_bound_a`
- `post_ber_is_upper_bound_b`
- `post_ber_upper_bound_a`
- `post_ber_upper_bound_b`

如果你的目标是确认不同 CSV 文件是不是配置一致，重点看参数快照列：

- `alpha_list`
- `beta_list`
- `siso_active_list`
- `chase_*`
- `mux_*`
- `hybrid_*`
- `early_stop_*`

---

## 15. 当前边界

当前这个 CSV 是“BER 曲线导向”的结果文件，不是 shared 内部调试 dump。

因此它不会回答下面这些问题：

- 某个 tile 里 shared core 实际 produced 了多少行
- failed 了多少行
- A/B 是否在 shared 预算竞争里偏流
- 某个 chunk 的内部调度细节是什么

如果要看这些，建议回到：

- [apps/ofec_two_stream_shared_decoder.cpp](/home/zsr71/projects/newcode_two_stream_shared_chase/apps/ofec_two_stream_shared_decoder.cpp)

在单点模式下配合控制台输出或后续专门的观测项去看。
