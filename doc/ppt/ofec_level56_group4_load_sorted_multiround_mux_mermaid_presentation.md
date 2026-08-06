---
marp: true
paginate: true
size: 16:9
theme: default
style: |
  section {
    font-family: "Noto Sans CJK SC", "Microsoft YaHei", "PingFang SC", sans-serif;
    color: #111111;
    background: #ffffff;
    padding: 54px 72px;
    font-size: 25px;
  }
  section::after {
    color: #444444;
    font-size: 16px;
    content: 'oFEC · Level 5/6 共享译码流程（Mermaid）  |  ' attr(data-marpit-pagination);
  }
  h1 { color: #111111; font-size: 52px; margin: 0 0 18px; }
  h2 { color: #111111; font-size: 38px; margin: 0 0 18px; border-bottom: 3px solid #222222; padding-bottom: 10px; }
  h3 { color: #111111; font-size: 28px; margin: 8px 0; }
  strong { color: #c00000; font-weight: 700; }
  code { color: #111111; background: #eeeeee; padding: 2px 6px; border-radius: 4px; }
  table { font-size: 20px; margin: 12px auto; }
  th { background: #e3e3e3; color: #111111; }
  th, td { padding: 8px 12px; border-color: #bdbdbd; }
  blockquote { border-left: 6px solid #c00000; background: #f3f3f3; padding: 12px 18px; margin: 14px 0; }
  .lead { background: #eeeeee; color: #111111; }
  .lead h1, .lead h2 { color: #111111; }
  .lead::after { color: #444444; }
  .columns { display: grid; grid-template-columns: 1fr 1fr; gap: 36px; }
  .card { background: #f3f3f3; border-radius: 0; box-shadow: none; padding: 16px 22px; }
  .diagram { display: block; width: 100%; max-height: 420px; margin: 24px auto 18px; }
  .diagram-wide { display: block; width: 100%; max-height: 340px; margin: 24px auto 18px; }
---

<!-- _class: lead -->

# oFEC 第五 / 六级共享译码流程

### Mermaid 流程图版本

`64 code` · `16 组 × 4 code` · `8 SISO + 8 HISO`

---

## 解码目标：有限共享资源优先服务更高负载组

<div class="columns">
<div class="card">

### 输入与资源

- Level 5 和 Level 6 合计 **64 个 code**
- 固定分为 **16 个组**，每组 4 个 code
- 共享 **8 路 SISO** 和 **8 路 HISO**

</div>
<div class="card">

### 核心策略

- 先识别无需共享译码的 code
- 统计每个组的待译码负载
- 优先调度负载更高的组
- 空余预算可再次分给仍有待译码 code 的组

</div>
</div>

> 关键目标：在最多 **8 次组级译码机会** 内，尽可能覆盖需要译码的 code。

---

## 总体解码流程

![w:1120](mermaid/01_overview.svg)

流程先完成全体 code 的判断与分类，再统一决定共享资源如何分配。

---

## 步骤一：EarlyStop 与 Hybrid 分类

![w:1120](mermaid/02_preprocess.svg)

- **EarlyStop**：已满足停止条件的 code 直接输出，不再竞争共享译码资源。
- **Hybrid 分类**：为其余 code 判断译码优先级，以及适合 HISO、SISO 或两者均可的路径。
- 分类本身不产生最终译码结果；它为后续资源分配提供依据。

---

## 步骤二：固定四行分组并统计负载

| 固定分组 | Group 1 | Group 2 | … | Group 16 |
| --- | --- | --- | --- | --- |
| code 范围 | 1–4 | 5–8 | … | 61–64 |
| 组内资源 | 4-to-1 MUX | 4-to-1 MUX | … | 4-to-1 MUX |

每个组的负载 = 其中**未命中 EarlyStop 的 code 数量**。

| 组负载 | 含义 |
| ---: | --- |
| 0 | 该组无需共享译码 |
| 1–4 | 该组仍有 1–4 个 code 等待共享译码 |

> 固定分组保证每个 4-to-1 MUX 只在本组内选择 code。

---

## 步骤三：根据非空组数选择调度方式

![w:1120](mermaid/03_schedule.svg)

**负载高的组优先**，使有限资源优先覆盖更多待译码 code。

---

## 步骤四：低负载时的多轮调度

当 `0 < K < 8`，第一次调度后仍有空余资源：

![w:1080](mermaid/04_multiround.svg)

重复直到 **所有待译码 code 都已获得机会**，或 **8 次组级机会用完**。

---

## 步骤五：组内 HISO / SISO 选择

![w:1100](mermaid/05_mux_output.svg)

> 一次组进入最多产生一个 HISO 译码和一个 SISO 译码；两路始终选择不同的 code，并优先服务更难译码的可选 code。

---

## 输出与未覆盖 code 的处理

<div class="columns">
<div class="card">

### 已获得资源

- EarlyStop 命中 → EarlyStop 输出
- HISO 资源 → HISO 译码输出
- SISO 资源 → SISO 译码输出

</div>
<div class="card">

### 未获得资源

- 非 EarlyStop code 保持原值
- 不产生新的译码输出
- 不占用额外共享资源

</div>
</div>

> 同一个 code 在整个流程中只会走一条输出路径。

---

<!-- _class: lead -->

# 解码流程总结

1. **EarlyStop 先行**：已满足条件的 code 直接输出，不占共享资源。  
2. **固定四行分组**：64 个 code 形成 16 个稳定的 4-to-1 MUX 候选组。  
3. **负载驱动调度**：优先选择待译码 code 更多的组。  
4. **多轮补充覆盖**：非空组少于 8 个时，剩余资源继续服务仍有负载的组。  
5. **HISO/SISO 分工**：每次组进入最多服务两个不同 code。  

### 在严格的 8 组共享资源预算内，获得更高的待译码 code 覆盖率。
