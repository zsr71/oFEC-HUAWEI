# TwoMain 单次运行清单

每次实验由统一运行器自动生成 `manifest.json` 和 `result.json`。它们位于 `runs/twomain/<experiment_id>/<run_id>/`，不进入普通 Git 历史。

## 一、运行身份

- 实验编号和方案编号；
- 配置文件路径及其 SHA-256；
- 运行编号；
- 开始时间、结束时间和耗时；
- 主机名、操作系统和 Python 版本；
- Git 分支、提交 SHA、标签；
- 工作区和暂存区是否存在差异；
- 使用的可执行文件路径和 SHA-256；
- 编译器和 CMake 版本。

## 二、完整参数

- 信噪比；
- 信息 bit 数；
- bit 生成和信道随机种子；
- LLR 位宽、量化裁剪比例；
- EarlyStop 条件、动作和逐级开关；
- 五六级共享、FIFO、`R_buf` 和帧尾排空；
- HISO、SISO 和 group-entry 数量；
- 调度模式和优先级；
- HISO 类别掩码；
- TwoMain 准入和输出策略；
- TwoMain 输出幅度或门限；
- BER 的固定首尾裁剪窗口数；
- 调度统计、错误位置和详细跟踪开关。

## 三、结果摘要

- 程序返回值；
- Pre-FEC 和 Post-FEC 错误数、比较 bit 数和 BER；
- Post-FEC 错误位置文件 SHA-256；
- Level 5 和 Level 6 的 EarlyStop、HISO、SISO、Unscheduled 数量；
- 各 HISO 类别的 HISO、SISO、Unscheduled 数量；
- FIFO 服务时刻数、完成 batch 数、全 EarlyStop batch 数；
- forced eviction 数量；
- 最大 FIFO 深度和最小窗口位置；
- 与配置内冻结签名逐项比较的结果。

## 四、原始输出

- `stdout.log`：程序完整标准输出和标准错误；
- `post_fec_error_positions.txt`：完整 Post-FEC 错误位置；
- `data/level56_schedule/`：调度轮次、逐 code 动作和 FIFO 时刻表；
- `data/early_stop_hist/`：EarlyStop 样本；
- `data/early_stop_debug/`：分组绑定观测；
- 按配置显式开启的其他 trace。

大型输出只保存在运行目录或单独归档。Git 中只保存配置、小型汇总和最终结论。
