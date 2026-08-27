#!/usr/bin/env python3
"""Summarize Level56 buffered FIFO batch service lifetimes.

The source of truth for configuration is the run log matched by the CSV
run_id.  Service counts are counts of ordinary_batch appearances.  Full
EarlyStop batches therefore have zero ordinary services by definition.
Only retired batches are used for lifetime distribution statistics; batches
still resident at the end of a non-drained run are reported as censored.
"""

from __future__ import annotations

import csv
import glob
import math
import os
import re
from collections import Counter, defaultdict
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, Iterable, List, Optional, Tuple


ROOT = Path(__file__).resolve().parent
DATA = ROOT / "data"
TIMES = DATA / "level56_schedule"
OUT = TIMES / "fifo_service_analysis.md"


@dataclass
class Meta:
    run_id: str
    log_path: str
    label: str
    ebn0: str = "?"
    fifo: str = "?"
    rbuf: str = "?"
    drain: str = "OFF"
    hiso: str = "?"
    siso: str = "?"
    schedule: str = "?"
    entries: str = "default"
    nbits: str = "default"
    seed: str = "?"
    pre_errors: str = "?"
    post_errors: str = "?"
    ber_bits: str = "?"
    pre_ber: str = "?"
    post_ber: str = "?"


def first(pattern: str, text: str, default: str = "?") -> str:
    m = re.search(pattern, text)
    return m.group(1) if m else default


def build_meta() -> Dict[str, Meta]:
    result: Dict[str, Meta] = {}
    for p in sorted(DATA.glob("run_*_single.log")):
        text = p.read_text(errors="ignore")
        m = re.search(r"run_(\d{8}-\d{6})_single\.log", p.name)
        if not m:
            continue
        run_id = m.group(1).replace("-", "_")
        label = first(r"run_pipeline\(label=([^,]+)", text)
        hs = re.search(r"HISO/SISO = (\d+)/(\d+)", text)
        nb = re.search(r"nbits(\d+)", label)
        ber = re.search(
            r"Pre-FEC BER=([^ ]+) \(errs=(\d+)/(\d+)\).*?"
            r"Post-FEC BER=([^ ]+) \(errs=(\d+)/(\d+)\)", text)
        result[run_id] = Meta(
            run_id=run_id,
            log_path=str(p),
            label=label,
            ebn0=first(r"Eb/N0=([0-9.]+)", text),
            fifo=first(r"buffered_fifo = (ON|OFF)", text),
            rbuf=first(r"buffer_rows = (\d+)", text),
            drain=first(r"drain_at_frame_end = (ON|OFF)", text, "OFF"),
            hiso=hs.group(1) if hs else "?",
            siso=hs.group(2) if hs else "?",
            schedule=("Group4" if "schedule_mode = 1" in text else
                      "GlobalPriority" if "schedule_mode = 0" in text else "?"),
            entries=first(r"group4_max_entries = (\d+)", text, "default"),
            nbits=nb.group(1) if nb else "default",
            seed=first(r"RNG seeds \(bitgen/channel\) = ([^\\n]+)", text),
            pre_ber=ber.group(1) if ber else "?",
            pre_errors=ber.group(2) if ber else "?",
            ber_bits=ber.group(3) if ber else "?",
            post_ber=ber.group(4) if ber else "?",
            post_errors=ber.group(5) if ber else "?",
        )
    return result


def parse_retirements(value: str) -> Iterable[Tuple[str, str]]:
    for token in re.split(r"[|;]", value or ""):
        if ":" in token:
            batch, reason = token.split(":", 1)
            if batch:
                yield batch, reason


def percentile(values: List[float], q: float) -> float:
    if not values:
        return float("nan")
    a = sorted(values)
    if len(a) == 1:
        return a[0]
    pos = (len(a) - 1) * q
    lo = int(math.floor(pos))
    hi = int(math.ceil(pos))
    if lo == hi:
        return a[lo]
    return a[lo] + (a[hi] - a[lo]) * (pos - lo)


def fmt(x: float, digits: int = 2) -> str:
    if isinstance(x, float) and math.isnan(x):
        return "—"
    if abs(x - round(x)) < 1e-10:
        return str(int(round(x)))
    return f"{x:.{digits}f}"


def config_key(m: Meta) -> str:
    return (f"Eb/N0={m.ebn0} dB; FIFO={m.fifo}; R_buf={m.rbuf}; "
            f"HISO/SISO={m.hiso}/{m.siso}; {m.schedule}; "
            f"entries={m.entries}; drain={m.drain}")


def selected_files() -> List[Path]:
    # One canonical file per run/configuration.  eqobs1 is retained where it
    # is the only detailed run; exact duplicate reruns are intentionally not
    # counted as independent samples.
    names = [
        # Complete/near-complete primary comparisons.
        "debug_L6_ebn03.065_nbits60014592_fifo1_rbuf16000_drain1_hiso0_siso8_schedgroup4_level56_buffered_times.csv",
        "debug_L6_ebn03.065_nbits60014592_fifo1_rbuf16000_drain1_hiso8_siso8_schedgroup4_level56_buffered_times.csv",
        "debug_L6_ebn03.09_fifo1_rbuf1600_drain1_hiso0_siso8_schedgroup4_eqobs1_level56_buffered_times.csv",
        "debug_L6_ebn03.1_rbuf64_level56_buffered_times.csv",
        # 3.09 dB R_buf sweep, 0/8 Group4.
        "debug_L6_ebn03.09_fifo1_rbuf138_hiso0_siso8_schedgroup4_eqobs1_level56_buffered_times.csv",
        "debug_L6_ebn03.09_fifo1_rbuf512_hiso0_siso8_schedgroup4_eqobs1_level56_buffered_times.csv",
        "debug_L6_ebn03.09_fifo1_rbuf1200_hiso0_siso8_schedgroup4_eqobs1_level56_buffered_times.csv",
        "debug_L6_ebn03.09_fifo1_rbuf1600_hiso0_siso8_schedgroup4_eqobs1_level56_buffered_times.csv",
        # Low-SNR pressure/drain runs.
        "debug_L6_ebn03.06_fifo1_rbuf1600_drain1_hiso0_siso8_schedgroup4_eqobs1_level56_buffered_times.csv",
        "debug_L6_ebn03.06_fifo1_rbuf3200_drain1_hiso0_siso8_schedgroup4_eqobs1_level56_buffered_times.csv",
        "debug_L6_ebn03.065_nbits60014592_fifo1_rbuf6400_drain1_hiso0_siso8_schedgroup4_level56_buffered_times.csv",
        # 8/8 Group4 Eb/N0/R_buf sweep files.
        "debug_L6_fifo1_rbuf0_level56_buffered_times.csv",
        "debug_L6_fifo1_rbuf8_level56_buffered_times.csv",
        "debug_L6_fifo1_rbuf16_level56_buffered_times.csv",
        "debug_L6_fifo1_rbuf32_level56_buffered_times.csv",
        "debug_L6_ebn03.07_fifo1_rbuf0_level56_buffered_times.csv",
        "debug_L6_ebn03.07_fifo1_rbuf8_level56_buffered_times.csv",
        "debug_L6_ebn03.07_fifo1_rbuf16_level56_buffered_times.csv",
        "debug_L6_ebn03.07_fifo1_rbuf32_level56_buffered_times.csv",
        "debug_L6_ebn03.07_rbuf64_level56_buffered_times.csv",
        "debug_L6_ebn03.07_rbuf80_level56_buffered_times.csv",
        "debug_L6_ebn03.07_rbuf96_level56_buffered_times.csv",
        "debug_L6_ebn03.08_fifo1_rbuf0_level56_buffered_times.csv",
        "debug_L6_ebn03.08_fifo1_rbuf8_level56_buffered_times.csv",
        "debug_L6_ebn03.08_fifo1_rbuf16_level56_buffered_times.csv",
        "debug_L6_ebn03.08_fifo1_rbuf32_level56_buffered_times.csv",
        "debug_L6_ebn03.08_rbuf64_level56_buffered_times.csv",
        "debug_L6_ebn03.09_rbuf64_level56_buffered_times.csv",
        "debug_L6_ebn03.09_rbuf80_level56_buffered_times.csv",
        "debug_L6_ebn03.09_rbuf128_level56_buffered_times.csv",
        "debug_L6_ebn03.09_rbuf136_level56_buffered_times.csv",
        "debug_L6_ebn03.09_rbuf138_level56_buffered_times.csv",
        "debug_L6_ebn03.09_rbuf140_level56_buffered_times.csv",
        "debug_L6_ebn03.09_rbuf144_level56_buffered_times.csv",
        "debug_L6_ebn03.1_rbuf32_level56_buffered_times.csv",
        # FIFO GlobalPriority control.
        "debug_L6_ebn03.09_fifo1_rbuf138_hiso32_siso32_schedglobal_eqobs1_level56_buffered_times.csv",
    ]
    return [TIMES / n for n in names if (TIMES / n).exists()]


def analyze(path: Path, metas: Dict[str, Meta]):
    with path.open(newline="") as f:
        rows = list(csv.DictReader(f))
    if not rows:
        return None
    run_id = rows[0]["run_id"]
    meta = metas.get(run_id)
    if meta is None:
        return None

    arrivals = set()
    retire_reason: Dict[str, str] = {}
    retired = set()
    for r in rows:
        b = r.get("arrived_batch", "")
        if b and b != "B18446744073709551615":
            arrivals.add(b)
        for batch, reason in parse_retirements(r.get("retirements", "")):
            # A batch must retire once; repeated records would indicate a bad
            # export, so preserve first occurrence and expose duplicates.
            if batch not in retire_reason:
                retire_reason[batch] = reason
            retired.add(batch)

    ordinary: Counter[str] = Counter()
    events: Dict[str, List[dict]] = defaultdict(list)
    for r in rows:
        b = r.get("ordinary_batch", "")
        if b:
            ordinary[b] += 1
            events[b].append(r)
        hb = r.get("head_batch_before", "")
        if hb and hb != "":
            # Include full-EarlyStop head visits for queue/S_t relations.
            if hb not in events:
                events[hb] = []
            if not b or b != hb:
                events[hb].append(r)
        # A single time row may retire the original ordinary head and then
        # one or more full-EarlyStop heads.  Associate the row with every
        # retirement batch so S_t/FIFO/pending relations are attributed to
        # the correct batch rather than only to head_batch_before.
        for rb, _reason in parse_retirements(r.get("retirements", "")):
            if not events[rb] or events[rb][-1] is not r:
                events[rb].append(r)

    # Service distribution is over retired batches only.  A retired batch not
    # present in ordinary_batch is a full-EarlyStop fast-path batch (count=0).
    by_reason: Dict[str, List[dict]] = defaultdict(list)
    for b, reason in retire_reason.items():
        ev = events.get(b, [])
        # A FullEarlyStop retirement is a fast-path completion and does not
        # consume ordinary HISO/SISO service.  In a time row that retires a
        # normal head and then one or more fast-path heads, the exported
        # ordinary_batch field can coexist with a FullEarlyStop retirement;
        # never attribute that field to the fast-path batch.
        service = 0 if reason == "FullEarlyStop" else ordinary.get(b, 0)
        pending_before = [int(x["pending_before"] or 0) for x in ev]
        fifo_depth = [int(x["fifo_depth_before"] or 0) for x in ev]
        svals = [int(x["S_t"] or 0) for x in ev]
        last_ord = [x for x in ev if x.get("ordinary_batch") == b]
        last_pending = int(last_ord[-1]["pending_before"] or 0) if last_ord else 0
        by_reason[reason].append({
            "batch": b,
            "service": service,
            "max_pending": max(pending_before) if pending_before else 0,
            "last_pending_before": last_pending,
            "max_fifo": max(fifo_depth) if fifo_depth else 0,
            "min_s": min(svals) if svals else 0,
        })

    # Include all ordinary batches that somehow lack retirement as censored;
    # they are not used in lifetime quantiles but are counted in metadata.
    censored = arrivals - retired
    return meta, rows, arrivals, retired, censored, by_reason


def dist(values: List[int]) -> str:
    if not values:
        return "—"
    return "/".join(fmt(percentile([float(v) for v in values], q), 2)
                    for q in (0, .5, .9, .95, .99, 1.0))


def render_table(reason_data: Dict[str, List[dict]]) -> str:
    lines = []
    lines.append("| 退休原因 | batch数 | max服务 | mean服务 | 服务1次 | 服务2次 | 服务3次 | 服务4次以上 | mean max_pending | P95 max_pending |")
    lines.append("|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|")
    for reason in ("Normal", "FullEarlyStop", "ForcedEvicted"):
        vals = reason_data.get(reason, [])
        services = [x["service"] for x in vals]
        counts = Counter(services)
        if vals:
            mean = sum(services) / len(services)
            maxv = max(services)
            mp = [x["max_pending"] for x in vals]
            rel = (sum(mp)/len(mp), percentile([float(x) for x in mp], .95))
        else:
            mean = maxv = float("nan")
            rel = (float("nan"),)*2
        lines.append("| {} | {} | {} | {} | {} | {} | {} | {} | {} | {} |".format(
            reason, len(vals), fmt(maxv), fmt(mean),
            counts.get(1, 0), counts.get(2, 0), counts.get(3, 0),
            sum(v for k, v in counts.items() if k >= 4),
            fmt(rel[0]), fmt(rel[1])))
    return "\n".join(lines)


def render_ber_trend_analyses(analyses) -> List[str]:
    """Render requested BER trend analyses from the selected run records."""
    records = []
    for p, a in analyses:
        meta, rows, arrivals, retired, censored, by_reason = a
        if meta.pre_ber == "?" or meta.post_ber == "?":
            continue
        try:
            records.append({
                "meta": meta,
                "file": p.name,
                "pre_errors": int(meta.pre_errors),
                "post_errors": int(meta.post_errors),
                "pre_ber": float(meta.pre_ber),
                "post_ber": float(meta.post_ber),
                "bits": int(meta.ber_bits),
            })
        except ValueError:
            continue

    def same_group(r, ebn0, hiso, siso, schedule="Group4"):
        m = r["meta"]
        return (m.ebn0 == ebn0 and m.hiso == str(hiso) and
                m.siso == str(siso) and m.schedule == schedule and
                m.fifo == "ON")

    def trend_table(items):
        lines = [
            "| Eb/N0 (dB) | R_buf | drain | Pre-FEC BER | Post-FEC错误数 | Post-FEC BER | 比较总比特数 |",
            "|---:|---:|---|---:|---:|---:|---:|",
        ]
        for r in sorted(items, key=lambda x: (float(x["meta"].ebn0), int(x["meta"].rbuf))):
            m = r["meta"]
            lines.append(
                f"| {m.ebn0} | {m.rbuf} | {m.drain} | {m.pre_ber} | "
                f"{r['post_errors']} | {m.post_ber} | {r['bits']} |"
            )
        return lines

    lines = [
        "",
        "## BER 随 R_buf / HISO-SISO 的变化分析",
        "",
        "> 下面的趋势只使用上方 BER 汇总表中实际存在的运行。不同 Eb/N0、帧长、drain 状态或随机输入不混合计算；`R_buf` 变化趋势优先在固定 Eb/N0、资源和调度模式下观察。",
        "",
        "### 1. 8/8 配置下，Post-FEC BER 随 R_buf 的变化",
        "",
        "这里的 8/8 指 `HISO/SISO=8/8`，调度模式为 Group4。",
    ]
    for eb in sorted({r["meta"].ebn0 for r in records}, key=float):
        items = [r for r in records if same_group(r, eb, 8, 8)]
        if not items:
            continue
        lines += ["", f"#### Eb/N0={eb} dB", ""]
        lines += trend_table(items)
        # Report monotonicity only when there are at least two R_buf values.
        ordered = sorted(items, key=lambda x: int(x["meta"].rbuf))
        if len(ordered) >= 2:
            vals = [x["post_ber"] for x in ordered]
            if all(vals[i] >= vals[i+1] for i in range(len(vals)-1)):
                conclusion = "Post-FEC BER 随 R_buf 增大整体不升（单调下降或持平）。"
            elif all(vals[i] <= vals[i+1] for i in range(len(vals)-1)):
                conclusion = "Post-FEC BER 随 R_buf 增大整体不降（单调上升或持平）。"
            else:
                conclusion = "Post-FEC BER 随 R_buf 不是严格单调，存在运行噪声或其他状态因素。"
            lines += ["", f"结论：{conclusion}"]
        if any(x["meta"].drain == "OFF" for x in items):
            lines += ["> 注意：该组包含 drain=OFF 运行；BER 可能受到尾部未完成/ForcedEvicted 影响，不能只归因于 R_buf。"]

    lines += [
        "",
        "### 2. 0/8 配置下，Post-FEC BER 随 R_buf 的变化",
        "",
        "这里的 0/8 指 `HISO/SISO=0/8`，即纯 SISO 的受限 FIFO 配置。",
    ]
    for eb in sorted({r["meta"].ebn0 for r in records}, key=float):
        items = [r for r in records if same_group(r, eb, 0, 8)]
        if not items:
            continue
        lines += ["", f"#### Eb/N0={eb} dB", ""]
        lines += trend_table(items)
        ordered = sorted(items, key=lambda x: int(x["meta"].rbuf))
        if len(ordered) >= 2:
            vals = [x["post_ber"] for x in ordered]
            if all(vals[i] >= vals[i+1] for i in range(len(vals)-1)):
                conclusion = "Post-FEC BER 随 R_buf 增大整体不升（单调下降或持平）。"
            elif all(vals[i] <= vals[i+1] for i in range(len(vals)-1)):
                conclusion = "Post-FEC BER 随 R_buf 增大整体不降（单调上升或持平）。"
            else:
                conclusion = "Post-FEC BER 随 R_buf 不是严格单调，存在明显非单调变化。"
            lines += ["", f"结论：{conclusion}"]
        if any(x["meta"].drain == "OFF" for x in items):
            lines += ["> 注意：该组包含 drain=OFF 运行；请结合 censored/ForcedEvicted 数量阅读。"]

    lines += [
        "",
        "### 3. 同 R_buf 下，8/8 与 0/8 的 Post-FEC BER 对比",
        "",
        "只列出当前数据中 Eb/N0、帧长/运行样本和 R_buf 都能对齐的配对；没有对应数据的 R_buf 不做推断。",
        "",
        "| Eb/N0 (dB) | R_buf | 8/8 Pre-FEC BER | 0/8 Pre-FEC BER | 8/8 Post-FEC BER | 0/8 Post-FEC BER | 8/8错误数 | 0/8错误数 | 说明 |",
        "|---:|---:|---:|---:|---:|---:|---:|---:|---|",
    ]
    # Match the clean long-frame pair and the 3.09 short-frame pair.  More
    # pairs can be added automatically if future selected files provide them.
    pairs = []
    for a in records:
        ma = a["meta"]
        if ma.hiso != "8" or ma.siso != "8" or ma.schedule != "Group4":
            continue
        for b in records:
            mb = b["meta"]
            if (mb.hiso == "0" and mb.siso == "8" and mb.schedule == "Group4" and
                    ma.ebn0 == mb.ebn0 and ma.rbuf == mb.rbuf and
                    ma.nbits == mb.nbits):
                # Prefer same drain state; if unavailable, still report the
                # pair and expose the drain difference in the note.
                pairs.append((a, b))
    seen = set()
    for a, b in sorted(pairs, key=lambda x: (float(x[0]["meta"].ebn0), int(x[0]["meta"].rbuf))):
        key = (a["meta"].ebn0, a["meta"].rbuf, a["meta"].nbits)
        if key in seen:
            continue
        seen.add(key)
        ma, mb = a["meta"], b["meta"]
        note = "同 Eb/N0、R_buf、帧长"
        if ma.drain != mb.drain:
            note += f"，但 drain={ma.drain}/{mb.drain}"
        ratio = (a["post_ber"] / b["post_ber"]
                 if b["post_ber"] > 0 else float("inf") if a["post_ber"] > 0 else 1.0)
        if math.isinf(ratio):
            note += "；0/8 为 0，8/8 非 0"
        else:
            note += f"；BER比值(8/8 ÷ 0/8)={ratio:.3g}"
        lines.append(
            f"| {ma.ebn0} | {ma.rbuf} | {ma.pre_ber} | {mb.pre_ber} | "
            f"{ma.post_ber} | {mb.post_ber} | {a['post_errors']} | "
            f"{b['post_errors']} | {note} |"
        )
    if not pairs:
        lines.append("| — | — | — | — | — | — | — | — | 当前选定数据没有同 R_buf 的 8/8 与 0/8 配对 |")

    lines += [
        "",
        "#### 解读",
        "",
        "- `R_buf` 增大主要改变 pending 的保存深度和 ForcedEvicted 风险；它不直接改变单个 code 的 HISO/SISO 解码算法。",
        "- 因此在同一资源配置下，BER 随 `R_buf` 下降通常表示 buffer 减少了边界淘汰或尾部未完成影响；若不单调，应结合 drain、ForcedEvicted 和 censored 数量判断。",
        "- 8/8 与 0/8 的差异不能只解释为服务次数差异：HISO 会改变 code 的输出和后续 SRAM/history 输入，因此即使两组都没有 ForcedEvicted，BER 也可能不同。",
        "- 当前同 R_buf 配对数量有限，以上是已有样本的直接比较，不代表对所有随机帧和 Eb/N0 的统计规律。",
        "",
        "#### 为什么长帧和短帧的结论方向相反",
        "",
        "- 长帧 `Eb/N0=3.065, R_buf=16000, drain=ON` 配对中，两组都是 16897/16897 个 arrival batch 已退休，`ForcedEvicted=0`，因此比较的是完整生命周期。0/8 的 Post-FEC BER 为 `1.87688e-7`，8/8 为 `2.71294e-6`。",
        "- 短帧 `Eb/N0=3.09, R_buf=138, drain=OFF` 配对中，8/8 仍有 24 个 censored batch，0/8 有 67 个 censored batch，并且 0/8 有 855 个 `ForcedEvicted`，8/8 为 0。0/8 的 `1.19228e-5` 与 8/8 的 `2.79713e-7` 不能视为只由 HISO/SISO 导致的差异。",
        "- 短帧的两个运行虽然 Pre-FEC BER 相同，但 buffer 压力和退休路径不同：0/8 没有 HISO 服务，SISO 压力更集中，更容易出现 FIFO 边界淘汰；8/8 的 HISO 服务会改变历史状态写回和后续窗口输入。两种机制叠加后，方向完全可能与长帧不同。",
        "- 另外，长帧和短帧的 `R_buf`、drain 状态、观测长度不同；短帧 Post-FEC 错误数只有 8（8/8）或 341（0/8），低错误数的相对统计不确定性也更大。因此这两个配对应作为两个不同实验解读，不能拼成一条“8/8 一定优于/劣于 0/8”的规律。",
        "- 要判断 HISO/SISO 的稳定影响，建议固定 `Eb/N0`、帧长、R_buf、drain、seeds 和 equivalence-observation 设置，至少各重复多次，并优先使用 `drain=ON`、`ForcedEvicted=0`、`censored=0` 的运行。",
    ]
    return lines


def main() -> None:
    metas = build_meta()
    analyses = []
    for p in selected_files():
        a = analyze(p, metas)
        if a:
            analyses.append((p, a))
    lines = [
        "# Level 5/6 Buffered FIFO：不同解码配置下 batch 服务次数汇总",
        "",
        "> 生成时间：2026-08-26。服务次数定义为 `ordinary_batch` 在时间 CSV 中出现的次数；`FullEarlyStop` 不占普通 HISO/SISO 预算，因此普通服务次数为 0。每张表只对已经退休的 batch 统计生命周期分布；未 drain 运行中仍留在 FIFO 的 batch 计入 `censored`，不计入分位数。百分位采用线性插值。",
        "",
        "## 字段说明",
        "",
        "- `max_pending`：该 batch 被观察到的最大 `pending_before`；比退休后的 `final pending` 更有信息，因为退休时 pending 会被清零。",
        "- `P95 max_pending`：先对每个 batch 取其生命周期内的最大 pending，再对这些 batch 级最大值取 95 分位；它不是某一时刻所有 pending 值的 95 分位。",
        "- `服务4次以上` 是服务次数 `>=4` 的 batch 数。",
        "",
        "## 各配置 Pre-FEC / Post-FEC 错误汇总",
        "",
        "| 配置 | run_id | Pre-FEC错误数 | Pre-FEC BER | Post-FEC错误数 | Post-FEC BER | BER比较总比特数 |",
        "|---|---|---:|---:|---:|---:|---:|",
    ]
    for p, a in analyses:
        meta, rows, arrivals, retired, censored, by_reason = a
        lines.append(
            f"| {config_key(meta)} | `{meta.run_id}` | {meta.pre_errors} | "
            f"{meta.pre_ber} | {meta.post_errors} | {meta.post_ber} | "
            f"{meta.ber_bits} |"
        )
    lines += render_ber_trend_analyses(analyses)
    lines += [
        "",
        "> Pre-FEC/Post-FEC 错误数和比较总比特数直接取对应 run log 的 `[RESULT]` 行；`?` 表示该日志没有可解析的 BER 结果。",
        "",
        "## 数据源清单",
        "",
        "| 配置 | CSV | run_id | arrival | retired | censored | rows |",
        "|---|---|---|---:|---:|---:|---:|",
    ]
    for p, a in analyses:
        meta, rows, arrivals, retired, censored, by_reason = a
        lines.append(f"| {config_key(meta)} | `{p.name}` | `{meta.run_id}` | {len(arrivals)} | {len(retired)} | {len(censored)} | {len(rows)} |")
    for idx, (p, a) in enumerate(analyses, 1):
        meta, rows, arrivals, retired, censored, by_reason = a
        lines += ["", f"## {idx}. {config_key(meta)}", "", f"数据源：`{p}`；run_id：`{meta.run_id}`。", "", render_table(by_reason)]
        if censored:
            lines += ["", f"> 注意：该运行有 {len(censored)} 个 arrival batch 在观测结束时尚未退休（右截尾）；上表的服务次数分布只包含 {len(retired)} 个已退休 batch。"]
        else:
            lines += ["", "> 该运行的 arrival batch 全部已退休；服务次数分布没有尾部截尾。"]
    OUT.write_text("\n".join(lines) + "\n")
    print(f"wrote {OUT} ({len(analyses)} configs)")


if __name__ == "__main__":
    main()
