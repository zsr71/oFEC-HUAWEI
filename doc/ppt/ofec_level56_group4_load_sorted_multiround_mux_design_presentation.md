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
    content: 'oFEC · Level 5/6 共享译码流程  |  ' attr(data-marpit-pagination);
  }
  h1 { color: #111111; font-size: 52px; margin: 0 0 18px; }
  h2 { color: #111111; font-size: 38px; margin: 0 0 18px; border-bottom: 3px solid #222222; padding-bottom: 10px; }
  h3 { color: #111111; font-size: 28px; margin: 8px 0; }
  strong { color: #c00000; font-weight: 700; }
  code { color: #111111; background: #eeeeee; padding: 2px 6px; border-radius: 4px; }
  pre { font-size: 19px; line-height: 1.3; background: #f0f0f0; color: #111111; padding: 16px 20px; border-radius: 0; }
  table { font-size: 20px; margin: 12px auto; }
  th { background: #e3e3e3; color: #111111; }
  th, td { padding: 8px 12px; border-color: #bdbdbd; }
  blockquote { border-left: 6px solid #c00000; background: #f3f3f3; padding: 12px 18px; margin: 14px 0; }
  .lead { background: #eeeeee; color: #111111; }
  .lead h1, .lead h2 { color: #111111; }
  .lead::after { color: #444444; }
  .columns { display: grid; grid-template-columns: 1fr 1fr; gap: 36px; }
  .card { background: #f3f3f3; border-radius: 0; box-shadow: none; padding: 16px 22px; }
  .flow { display: flex; align-items: stretch; justify-content: center; gap: 12px; margin: 40px 0 28px; }
  .flow-card { flex: 1; min-height: 120px; background: #f2f2f2; border-top: 4px solid #111111; padding: 18px 14px; text-align: center; display: flex; flex-direction: column; justify-content: center; }
  .flow-card.red { border-top-color: #c00000; }
  .flow-card p { margin: 0; font-size: 22px; line-height: 1.35; }
  .flow-arrow { align-self: center; color: #555555; font-size: 34px; font-weight: 300; }
  .branch-flow { display: grid; grid-template-columns: 1fr 70px 1fr; gap: 14px; align-items: center; margin: 26px 0; }
  .branch-source, .branch-result { background: #f2f2f2; border-top: 4px solid #111111; padding: 18px; text-align: center; }
  .branch-result.red { border-top-color: #c00000; }
  .branch-arrow { text-align: center; font-size: 30px; color: #555555; }
  .branch-label { color: #c00000; font-size: 19px; font-weight: 700; margin: 4px 0; }
  .group-row { display: flex; gap: 16px; margin: 34px 0 22px; }
  .group { flex: 1; background: #f2f2f2; border-top: 4px solid #111111; padding: 14px; text-align: center; font-size: 20px; }
  .group b { display: block; font-size: 24px; margin-bottom: 6px; }
  .group-dots { align-self: center; font-size: 30px; color: #555555; }
  .loop { border: 2px solid #bdbdbd; background: #fafafa; padding: 18px 24px; margin: 24px 0; text-align: center; }
  .output-grid { display: grid; grid-template-columns: 1fr 1fr; gap: 18px; margin: 28px 0; }
  .output-row { display: grid; grid-template-columns: 1fr 52px 1fr; align-items: center; }
  .output-card { background: #f2f2f2; border-top: 4px solid #111111; padding: 16px; text-align: center; }
  .output-card.red { border-top-color: #c00000; }
  .center { text-align: center; }
---

<!-- _class: lead -->

# oFEC 第五 / 六级共享译码流程

### 四行分组、负载排序与多轮 MUX

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

<div class="flow">
<div class="flow-card"><p>Level 5 / 6<br><strong>64 个 code</strong></p></div>
<div class="flow-arrow">→</div>
<div class="flow-card red"><p><strong>EarlyStop</strong><br>预判</p></div>
<div class="flow-arrow">→</div>
<div class="flow-card"><p>Hybrid<br><strong>分类</strong></p></div>
<div class="flow-arrow">→</div>
<div class="flow-card"><p>固定分组<br><strong>统计负载</strong></p></div>
<div class="flow-arrow">→</div>
<div class="flow-card red"><p>排序与<br><strong>资源调度</strong></p></div>
<div class="flow-arrow">→</div>
<div class="flow-card"><p>HISO / SISO<br><strong>译码输出</strong></p></div>
</div>

流程先完成全体 code 的判断与分类，再统一决定共享资源如何分配。

---

## 步骤一：EarlyStop 与 Hybrid 分类

<div class="branch-flow">
<div class="branch-source"><strong>每个输入 code</strong><br>执行 EarlyStop 判断</div>
<div class="branch-arrow">→</div>
<div class="branch-result red"><strong>命中</strong><br>直接 EarlyStop 输出</div>
<div></div>
<div class="branch-arrow">↓<div class="branch-label">未命中</div></div>
<div class="branch-result"><strong>Hybrid 分类</strong><br>确定难度与可用路径</div>
</div>

- **EarlyStop**：已满足停止条件的 code 不再竞争共享译码资源。
- **Hybrid 分类**：为其余 code 判断译码优先级，以及适合 HISO、SISO 或两者均可的路径。
- 分类本身不产生最终译码结果；它为后续资源分配提供依据。

---

## 步骤二：固定四行分组并统计负载

<div class="group-row">
<div class="group"><b>Group 1</b>code 1–4<br>4-to-1 MUX</div>
<div class="group"><b>Group 2</b>code 5–8<br>4-to-1 MUX</div>
<div class="group-dots">···</div>
<div class="group"><b>Group 16</b>code 61–64<br>4-to-1 MUX</div>
</div>

每个组的负载 = 其中**未命中 EarlyStop 的 code 数量**。

| 组负载 | 含义 |
| ---: | --- |
| 0 | 该组无需共享译码 |
| 1–4 | 该组仍有 1–4 个 code 等待共享译码 |

> 固定分组保证每个 4-to-1 MUX 只在本组内选择 code。

---

## 步骤三：根据非空组数选择调度方式

令 **K** 为负载大于 0 的组数：

| 条件 | 调度方式 |
| --- | --- |
| `K = 0` | 所有 code 已 EarlyStop，直接结束 |
| `K > 8` | 按组负载从高到低，选择前 8 组 |
| `K = 8` | 8 个非空组全部进入译码 |
| `0 < K < 8` | 先让全部非空组进入，再用剩余资源进行多轮调度 |

<div class="flow">
<div class="flow-card"><p>统计<br><strong>组负载</strong></p></div>
<div class="flow-arrow">→</div>
<div class="flow-card"><p>计算<br><strong>K</strong></p></div>
<div class="flow-arrow">→</div>
<div class="flow-card red"><p>一次选择<br>或多轮选择</p></div>
<div class="flow-arrow">→</div>
<div class="flow-card"><p>最多 <strong>8</strong> 次<br>组级译码机会</p></div>
</div>

**负载高的组优先**，使有限资源优先覆盖更多待译码 code。

---

## 步骤四：低负载时的多轮调度

当 `0 < K < 8`，第一次调度后仍有空余资源：

<div class="flow">
<div class="flow-card"><p>第一轮<br><strong>全部非空组</strong>进入</p></div>
<div class="flow-arrow">→</div>
<div class="flow-card"><p>每组最多服务<br><strong>两个不同 code</strong></p></div>
<div class="flow-arrow">→</div>
<div class="flow-card"><p>重新计算<br><strong>剩余负载</strong></p></div>
</div>

<div class="loop"><strong>检查：</strong>是否仍有待译码 code，且组级译码机会尚未用完？</div>

<div class="columns">
<div class="card center"><strong>否</strong><br>停止调度</div>
<div class="card center"><strong>是</strong><br>按剩余负载重排<br>高负载组再次优先进入</div>
</div>

重复直到 **所有待译码 code 都已获得机会**，或 **8 次组级机会用完**。

---

## 步骤五：组内 HISO / SISO 选择

每次一个组进入共享译码时：

<div class="flow">
<div class="flow-card"><p>组内<br><strong>待译码 code</strong></p></div>
<div class="flow-arrow">→</div>
<div class="flow-card red"><p>SISO 优先服务<br><strong>SISO 专用 code</strong></p></div>
<div class="flow-arrow">→</div>
<div class="flow-card"><p>其余 code<br><strong>按难度排序</strong></p></div>
</div>

<div class="columns">
<div class="card center"><strong>SISO MUX</strong><br>选择 1 个 code</div>
<div class="card center"><strong>HISO MUX</strong><br>从剩余 code 中选择 1 个不同 code</div>
</div>

> 一次组进入最多产生一个 HISO 译码和一个 SISO 译码；两路始终选择不同的 code，并优先服务更难译码的可选 code。

---

## 步骤六：输出与未覆盖 code 的处理

<div class="output-grid">
<div class="output-row"><div class="output-card red"><strong>EarlyStop</strong><br>命中 code</div><div class="branch-arrow">→</div><div class="output-card">EarlyStop<br><strong>输出</strong></div></div>
<div class="output-row"><div class="output-card red">获得 <strong>HISO</strong><br>资源的 code</div><div class="branch-arrow">→</div><div class="output-card">HISO<br><strong>译码输出</strong></div></div>
<div class="output-row"><div class="output-card red">获得 <strong>SISO</strong><br>资源的 code</div><div class="branch-arrow">→</div><div class="output-card">SISO<br><strong>译码输出</strong></div></div>
</div>

<div class="loop"><strong>未获得资源的非 EarlyStop code</strong>：保持原值，不产生新的译码输出。</div>

- 同一个 code 在整个流程中只会走一条输出路径。
- 未获得 HISO/SISO 机会的 code 不产生新的译码结果。
- 因此，组负载排序与多轮调度直接决定有限共享资源的最终覆盖范围。

---

<!-- _class: lead -->

# 解码流程总结

1. **EarlyStop 先行**：已满足条件的 code 直接输出，不占共享资源。  
2. **固定四行分组**：64 个 code 形成 16 个稳定的 4-to-1 MUX 候选组。  
3. **负载驱动调度**：优先选择待译码 code 更多的组。  
4. **多轮补充覆盖**：非空组少于 8 个时，剩余资源继续服务仍有负载的组。  
5. **HISO/SISO 分工**：每次组进入最多服务两个不同 code。  

### 在严格的 8 组共享资源预算内，获得更高的待译码 code 覆盖率。
