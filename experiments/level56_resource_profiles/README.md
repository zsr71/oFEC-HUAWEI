# 双窗口组内资源配置对照

保留每个物理时刻的 8 个 group-entry；仅扩展单个 entry 的软/硬解码机会。

| 配置 | SISO 总容量 | HISO 总容量 | 单 entry |
|---|---:|---:|---|
| 1S1H（默认） | 8 | 8 | 最多 1 SISO + 1 HISO |
| 2S1H | 16 | 8 | 最多 2 SISO + 1 HISO |
| 1S2H | 8 | 16 | 最多 1 SISO + 2 HISO |

两个应用都读取 `LEVEL56_SHARED_SISO_ACTIVE`、`LEVEL56_SHARED_HISO_ACTIVE` 和
`LEVEL56_BUFFER_ROWS`。例如在项目根目录运行：

```bash
LEVEL56_SHARED_SISO_ACTIVE=16 LEVEL56_SHARED_HISO_ACTIVE=8 \
LEVEL56_BUFFER_ROWS=64 ./build_prior_fix_release/ofec_single
```

将程序改为 `ofec_sweep3` 即可对扫描使用相同资源配置；它仍按自身信噪比、
种子及累计错误停止规则运行，不等于本次固定单帧实验。

## 组内分配

先填满 SISO：仅软候选优先，可软可硬候选按原分类顺序补足。
随后从剩余可软可硬候选填 HISO。同一 code 不重复分配。
组负载排序、同负载 L5 优先、多轮共享预算、分类冻结、窗口运动均保持原规则。
解码单元编号为 `entry_slot * width + lane`，SISO/HISO 各自编号。

`2S1H` 不是“先按 1S1H 分配，再补一个 SISO”；它先分配两个 SISO。
因此 A 仅软、B/C/D 可软可硬的一组中：1S1H 服务 A(S)+B(H)，
2S1H 服务 A/B(S)+C(H)，1S2H 服务 A(S)+B/C(H)。
未用组内机会不能跨组借用，未用 entry 仍可按原规则再次服务有 pending 的组。

## 可重复实验

```bash
cmake --build build_prior_fix_release --target \
  ofec_single ofec_sweep3 ofec_level56_shared_regression_check -j2
./build_prior_fix_release/ofec_level56_shared_regression_check
python3 experiments/level56_resource_profiles/run_matrix.py \
  --output data/resource_profiles_long_20261014
python3 experiments/level56_resource_profiles/analyze.py \
  data/resource_profiles_long_20261014
```

完整矩阵：3.06/3.07 dB × Rbuf=0/32/64/128 × 三种配置，共 24 组。
每帧 60,014,592 信息位（256 window）；固定信息位种子 20260319、噪声种子
3182026。最多并行两组。帧尾规则、量化、BER 裁剪沿用修复后的当前基线。
每组保留独立日志、有效环境配置、程序哈希和代码差异，不覆盖现有输出。
`--short` 使用 10-window 短帧、3.06 dB、Rbuf=0/8/16/32，共 12 组；
`--resume` 只跳过同一程序哈希下已经正常结束的组。

## 统计口径

时间 CSV 新增 `hiso_calls/siso_calls/pending5/pending6`：前两项是本 t
实际分配并执行的 code 数；后两项是本 t 服务结束后的当前阶段残留 code 数。
`pending_code_times` 是残留数按 t 求和，不是互不重复的 code 数。
组调度 CSV 的 `consumed_entries` 为本次调度新消耗的 entry，不能把累计
`total_group_entries` 在微轮之间重复求和。分析脚本还验证时间/微步窗口合法性、
entry 与单元上限、阶段守恒及真实输入/输出一致。

资源利用率以全部物理 t 为分母，包括无普通解码的 t，不按调度器调用次数计。
HISO/SISO 调用次数分别比较；两者每次操作成本不同，不把调用次数简单相加
当成等成本复杂度。单种子有限帧的 BER 用于配对对照，不作为统计显著性证明。
