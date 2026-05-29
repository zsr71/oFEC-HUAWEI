#!/usr/bin/env python3
"""Offline dual-stream replay for shared-Chase validation.

This script consumes the invocation-level tile sample CSV emitted by
`apps/ofec_ber_window_probe.cpp` and performs the stage-2 validation
described in `doc/双虚拟流共享Chase方案验证分阶段实施建议.md`.

Default behavior:
- use `rows_need_siso_before_mux` as the demand metric
- split each `(ebn0_db, seed_index, tile_index)` time series into two virtual
  streams by invocation parity
- combine the two streams directly in timeline order
- compute peak / percentile / overflow statistics for the replayed demand
"""

from __future__ import annotations

import argparse
import csv
import math
from dataclasses import dataclass
from pathlib import Path
from statistics import mean
from typing import Dict, Iterable, List, Sequence, Tuple


GroupKey = Tuple[str, str, str]


@dataclass
class SampleRow:
    invocation: int
    rows_total: int
    rows_passed: int
    rows_hard_finish: int
    rows_need_siso_before_mux: int
    rows_unscheduled: int

    @property
    def need_decode_after_early_stop(self) -> int:
        return self.rows_total - self.rows_passed

    def demand_value(self, metric: str) -> int:
        if metric == "rows_need_siso_before_mux":
            return self.rows_need_siso_before_mux
        if metric == "need_decode_after_early_stop":
            return self.need_decode_after_early_stop
        if metric == "rows_unscheduled":
            return self.rows_unscheduled
        raise ValueError(f"Unsupported metric: {metric}")


@dataclass
class ReplayStats:
    split_mode: str
    count: int
    max_value: int
    mean_value: float
    p99: float
    p999: float
    p9999: float
    overflow_count: int
    overflow_ratio: float
    max_overflow_run: int
    avg_overflow_run: float


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Replay dual virtual streams with phase offsets from tile sample CSV."
    )
    parser.add_argument(
        "--input",
        type=Path,
        default=Path("data/ber_window_probe/per_seed_tile_early_stop_samples.csv"),
        help="Input CSV path from ofec_ber_window_probe.",
    )
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=Path("data/ber_window_probe/dual_stream_phase_replay"),
        help="Directory for summary CSV outputs.",
    )
    parser.add_argument(
        "--metric",
        choices=[
            "rows_need_siso_before_mux",
            "need_decode_after_early_stop",
            "rows_unscheduled",
        ],
        default="rows_need_siso_before_mux",
        help="Demand metric used as a(t).",
    )
    parser.add_argument(
        "--tile-index",
        type=int,
        default=None,
        help="Optional tile filter. When omitted, all tiles are analyzed independently.",
    )
    parser.add_argument(
        "--ebn0-db",
        type=str,
        default=None,
        help="Optional Eb/N0 filter using the exact CSV string value.",
    )
    parser.add_argument(
        "--seed-index",
        type=int,
        default=None,
        help="Optional seed filter.",
    )
    parser.add_argument(
        "--split-mode",
        choices=["parity", "chunk"],
        default="parity",
        help="How to construct the two virtual streams.",
    )
    parser.add_argument(
        "--chunk-size",
        type=int,
        default=8,
        help="Chunk size used when --split-mode=chunk.",
        )
    parser.add_argument(
        "--budget",
        type=int,
        default=None,
        help="Optional overflow threshold. When omitted, overflow metrics are skipped.",
    )
    return parser.parse_args()


def read_samples(
    csv_path: Path,
    tile_index: int | None,
    ebn0_db: str | None,
    seed_index: int | None,
) -> Dict[GroupKey, List[SampleRow]]:
    groups: Dict[GroupKey, List[SampleRow]] = {}
    with csv_path.open("r", newline="") as handle:
        reader = csv.DictReader(handle)
        for row in reader:
            row_tile = int(row["tile_index"])
            row_seed = int(row["seed_index"])
            row_ebn0 = row["ebn0_db"]
            if tile_index is not None and row_tile != tile_index:
                continue
            if ebn0_db is not None and row_ebn0 != ebn0_db:
                continue
            if seed_index is not None and row_seed != seed_index:
                continue

            key = (row_ebn0, row["seed_index"], row["tile_index"])
            groups.setdefault(key, []).append(
                SampleRow(
                    invocation=int(row["invocation"]),
                    rows_total=int(row["rows_total"]),
                    rows_passed=int(row["rows_passed"]),
                    rows_hard_finish=int(row["rows_hard_finish"]),
                    rows_need_siso_before_mux=int(row["rows_need_siso_before_mux"]),
                    rows_unscheduled=int(row["rows_unscheduled"]),
                )
            )
    for series in groups.values():
        series.sort(key=lambda item: item.invocation)
    return groups


def split_by_invocation_parity(
    samples: Sequence[SampleRow],
    metric: str,
) -> Tuple[List[int], List[int]]:
    stream_a: List[int] = []
    stream_b: List[int] = []
    for sample in samples:
        value = sample.demand_value(metric)
        if sample.invocation % 2 == 0:
            stream_a.append(value)
        else:
            stream_b.append(value)
    return stream_a, stream_b


def split_by_chunk(
    samples: Sequence[SampleRow],
    metric: str,
    chunk_size: int,
) -> Tuple[List[int], List[int]]:
    if chunk_size < 1:
        raise ValueError("chunk_size must be >= 1")
    stream_a: List[int] = []
    stream_b: List[int] = []
    for idx, sample in enumerate(samples):
        chunk_id = (idx // chunk_size) % 2
        value = sample.demand_value(metric)
        if chunk_id == 0:
            stream_a.append(value)
        else:
            stream_b.append(value)
    return stream_a, stream_b


def combine_streams(stream_a: Sequence[int], stream_b: Sequence[int]) -> List[int]:
    total_len = max(len(stream_a), len(stream_b))
    combined = [0] * total_len
    for idx, value in enumerate(stream_a):
        combined[idx] += value
    for idx, value in enumerate(stream_b):
        combined[idx] += value
    return combined


def percentile(values: Sequence[int], q: float) -> float:
    if not values:
        return 0.0
    if q <= 0.0:
        return float(min(values))
    if q >= 1.0:
        return float(max(values))
    ordered = sorted(values)
    pos = q * (len(ordered) - 1)
    lo = math.floor(pos)
    hi = math.ceil(pos)
    if lo == hi:
        return float(ordered[lo])
    frac = pos - lo
    return float(ordered[lo] * (1.0 - frac) + ordered[hi] * frac)


def overflow_runs(values: Sequence[int], budget: int | None) -> Tuple[int, int, float]:
    if budget is None:
        return 0, 0, 0.0
    overflow_count = 0
    run_lengths: List[int] = []
    current_run = 0
    for value in values:
        if value > budget:
            overflow_count += 1
            current_run += 1
        elif current_run > 0:
            run_lengths.append(current_run)
            current_run = 0
    if current_run > 0:
        run_lengths.append(current_run)
    if not run_lengths:
        return overflow_count, 0, 0.0
    return overflow_count, max(run_lengths), mean(run_lengths)


def summarize_replay(values: Sequence[int], split_mode: str, budget: int | None) -> ReplayStats:
    if not values:
        return ReplayStats(split_mode, 0, 0, 0.0, 0.0, 0.0, 0.0, 0, 0.0, 0, 0.0)
    overflow_count, max_run, avg_run = overflow_runs(values, budget)
    return ReplayStats(
        split_mode=split_mode,
        count=len(values),
        max_value=max(values),
        mean_value=mean(values),
        p99=percentile(values, 0.99),
        p999=percentile(values, 0.999),
        p9999=percentile(values, 0.9999),
        overflow_count=overflow_count,
        overflow_ratio=(overflow_count / len(values)) if values else 0.0,
        max_overflow_run=max_run,
        avg_overflow_run=avg_run,
    )


def write_per_group_csv(
    output_path: Path,
    rows: Iterable[Dict[str, object]],
) -> None:
    fieldnames = [
        "ebn0_db",
        "seed_index",
        "tile_index",
        "split_mode",
        "count",
        "max_value",
        "mean_value",
        "p99",
        "p999",
        "p9999",
        "overflow_count",
        "overflow_ratio",
        "max_overflow_run",
        "avg_overflow_run",
    ]
    with output_path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        for row in rows:
            writer.writerow(row)


def write_aggregated_csv(output_path: Path, stats: Sequence[ReplayStats]) -> None:
    fieldnames = [
        "split_mode",
        "count",
        "max_value",
        "mean_value",
        "p99",
        "p999",
        "p9999",
        "overflow_count",
        "overflow_ratio",
        "max_overflow_run",
        "avg_overflow_run",
    ]
    with output_path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        for item in stats:
            writer.writerow(
                {
                    "split_mode": item.split_mode,
                    "count": item.count,
                    "max_value": item.max_value,
                    "mean_value": f"{item.mean_value:.6f}",
                    "p99": f"{item.p99:.6f}",
                    "p999": f"{item.p999:.6f}",
                    "p9999": f"{item.p9999:.6f}",
                    "overflow_count": item.overflow_count,
                    "overflow_ratio": f"{item.overflow_ratio:.6f}",
                    "max_overflow_run": item.max_overflow_run,
                    "avg_overflow_run": f"{item.avg_overflow_run:.6f}",
                }
            )


def main() -> int:
    args = parse_args()
    groups = read_samples(args.input, args.tile_index, args.ebn0_db, args.seed_index)
    if not groups:
        raise SystemExit("No matching rows found in input CSV.")

    args.output_dir.mkdir(parents=True, exist_ok=True)

    per_group_rows: List[Dict[str, object]] = []
    pooled_values: List[int] = []
    for (ebn0_db, seed_str, tile_str), samples in sorted(groups.items()):
        if args.split_mode == "parity":
            stream_a, stream_b = split_by_invocation_parity(samples, args.metric)
        else:
            stream_a, stream_b = split_by_chunk(samples, args.metric, args.chunk_size)
        combined = combine_streams(stream_a, stream_b)
        pooled_values.extend(combined)
        stats = summarize_replay(combined, args.split_mode, args.budget)
        per_group_rows.append(
            {
                "ebn0_db": ebn0_db,
                "seed_index": seed_str,
                "tile_index": tile_str,
                "split_mode": args.split_mode,
                "count": stats.count,
                "max_value": stats.max_value,
                "mean_value": f"{stats.mean_value:.6f}",
                "p99": f"{stats.p99:.6f}",
                "p999": f"{stats.p999:.6f}",
                "p9999": f"{stats.p9999:.6f}",
                "overflow_count": stats.overflow_count,
                "overflow_ratio": f"{stats.overflow_ratio:.6f}",
                "max_overflow_run": stats.max_overflow_run,
                "avg_overflow_run": f"{stats.avg_overflow_run:.6f}",
            }
        )
    aggregated = [summarize_replay(pooled_values, args.split_mode, args.budget)]

    per_group_path = args.output_dir / "per_group_delta_stats.csv"
    aggregated_path = args.output_dir / "aggregated_delta_stats.csv"
    write_per_group_csv(per_group_path, per_group_rows)
    write_aggregated_csv(aggregated_path, aggregated)

    print(f"[INFO] input={args.input}")
    print(
        f"[INFO] groups={len(groups)} metric={args.metric}"
        f" split_mode={args.split_mode}"
        + (f" chunk_size={args.chunk_size}" if args.split_mode == "chunk" else "")
    )
    if args.budget is not None:
        print(f"[INFO] budget={args.budget}")
    print(f"[INFO] wrote {per_group_path}")
    print(f"[INFO] wrote {aggregated_path}")
    summary = aggregated[0]
    print(
        "[RESULT] replay summary:"
        f" max={summary.max_value}, p99={summary.p99:.3f},"
        f" p99.9={summary.p999:.3f}, p99.99={summary.p9999:.3f},"
        f" overflow_count={summary.overflow_count}"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
