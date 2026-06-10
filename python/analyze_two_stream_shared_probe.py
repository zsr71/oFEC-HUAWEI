#!/usr/bin/env python3
"""Analyze two-stream shared probe CSV outputs without extra dependencies.

This script reads the CSVs produced by:
`apps/ofec_two_stream_ber_window_probe.cpp`

Default input directory:
- data/two_stream_ber_window_probe

Default output directory:
- data/two_stream_ber_window_probe/analysis

It generates a few lightweight diagnosis tables:
- seed_ranking.csv
- suspicious_invocations.csv
- tile_index_summary.csv
- hybrid_class_summary.csv
- stream_class_summary.csv
"""

from __future__ import annotations

import argparse
import csv
from collections import defaultdict
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, Iterable, List, Sequence, Tuple


DEFAULT_INPUT_DIR = Path("data/two_stream_ber_window_probe")
DEFAULT_OUTPUT_SUBDIR = "analysis"


@dataclass(frozen=True)
class InvocationKey:
    run_id: str
    ebn0_db: str
    seed_index: int
    chunk_index: int
    invocation: int
    tile_index: int


@dataclass
class SeedSummary:
    run_id: str
    ebn0_db: str
    seed_index: int
    chunk_index: int
    bitgen_seed_a: int
    channel_seed_a: int
    bitgen_seed_b: int
    channel_seed_b: int
    pre_ber_a: float
    post_ber_a: float
    pre_ber_b: float
    post_ber_b: float
    tile_sample_rows: int
    hybrid_count_rows: int
    row_map_rows: int


@dataclass
class TileSample:
    key: InvocationKey
    rows_total: int
    rows_early_stop: int
    rows_not_early_stop: int
    rows_hard_finish: int
    rows_need_siso_before_mux: int
    rows_soft_scheduled: int
    rows_unscheduled: int
    produced_rows: int
    failed_rows: int
    produced_rows_a: int
    produced_rows_b: int
    failed_rows_a: int
    failed_rows_b: int
    stream_rows_a: int
    stream_rows_b: int


@dataclass
class HybridCounts:
    key: InvocationKey
    rows_seen_by_hybrid: int
    class_none_count: int
    class_bch_hard_decoded_count: int
    class_clean_count: int
    class_parity_only_count: int
    class_one_main_count: int
    class_one_main_plus_parity_count: int
    class_two_main_count: int
    class_suspicious_count: int
    class_hard_fail_count: int
    deferred_candidate_count: int
    deferred_priority_0_count: int
    deferred_priority_1_count: int
    deferred_priority_2_count: int
    deferred_reclaimed_to_hard_finish_count: int


@dataclass
class RowMap:
    key: InvocationKey
    merged_row: int
    stream_id: int
    stream_label: str
    source_local_row: int
    source_global_row: int
    early_stop_hit: bool
    hybrid_class: str
    final_tag: str
    scheduled_for_soft: bool
    produced_row: bool


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Analyze two-stream shared probe CSV outputs."
    )
    parser.add_argument(
        "--input-dir",
        type=Path,
        default=DEFAULT_INPUT_DIR,
        help="Directory containing probe CSV outputs.",
    )
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=None,
        help="Directory for analysis CSV outputs. Defaults to <input-dir>/analysis.",
    )
    parser.add_argument(
        "--run-id",
        type=str,
        default=None,
        help="Analyze only one run_id. Defaults to the latest run_id when multiple exist.",
    )
    parser.add_argument(
        "--all-runs",
        action="store_true",
        help="Analyze all run_id values together instead of selecting the latest one.",
    )
    parser.add_argument(
        "--ebn0-db",
        type=str,
        default=None,
        help="Optional exact Eb/N0 CSV string filter, for example 3.050000.",
    )
    parser.add_argument(
        "--seed-index",
        type=int,
        default=None,
        help="Optional seed_index filter.",
    )
    parser.add_argument(
        "--tile-index",
        type=int,
        default=None,
        help="Optional tile_index filter.",
    )
    parser.add_argument(
        "--top",
        type=int,
        default=20,
        help="How many rows to keep in the ranking-style outputs.",
    )
    return parser.parse_args()


def read_csv_rows(path: Path) -> List[dict]:
    with path.open("r", newline="") as handle:
        return list(csv.DictReader(handle))


def parse_bool_int(value: str) -> bool:
    return int(value) != 0


def invocation_key_from_row(row: dict) -> InvocationKey:
    return InvocationKey(
        run_id=row["run_id"],
        ebn0_db=row["ebn0_db"],
        seed_index=int(row["seed_index"]),
        chunk_index=int(row["chunk_index"]),
        invocation=int(row["invocation"]),
        tile_index=int(row["tile_index"]),
    )


def load_seed_summaries(rows: Sequence[dict]) -> List[SeedSummary]:
    out: List[SeedSummary] = []
    for row in rows:
        out.append(
            SeedSummary(
                run_id=row["run_id"],
                ebn0_db=row["ebn0_db"],
                seed_index=int(row["seed_index"]),
                chunk_index=int(row["chunk_index"]),
                bitgen_seed_a=int(row["bitgen_seed_a"]),
                channel_seed_a=int(row["channel_seed_a"]),
                bitgen_seed_b=int(row["bitgen_seed_b"]),
                channel_seed_b=int(row["channel_seed_b"]),
                pre_ber_a=float(row["pre_ber_a"]),
                post_ber_a=float(row["post_ber_a"]),
                pre_ber_b=float(row["pre_ber_b"]),
                post_ber_b=float(row["post_ber_b"]),
                tile_sample_rows=int(row["tile_sample_rows"]),
                hybrid_count_rows=int(row["hybrid_count_rows"]),
                row_map_rows=int(row["row_map_rows"]),
            )
        )
    return out


def load_tile_samples(rows: Sequence[dict]) -> List[TileSample]:
    out: List[TileSample] = []
    for row in rows:
        out.append(
            TileSample(
                key=invocation_key_from_row(row),
                rows_total=int(row["rows_total"]),
                rows_early_stop=int(row["rows_early_stop"]),
                rows_not_early_stop=int(row["rows_not_early_stop"]),
                rows_hard_finish=int(row["rows_hard_finish"]),
                rows_need_siso_before_mux=int(row["rows_need_siso_before_mux"]),
                rows_soft_scheduled=int(row["rows_soft_scheduled"]),
                rows_unscheduled=int(row["rows_unscheduled"]),
                produced_rows=int(row["produced_rows"]),
                failed_rows=int(row["failed_rows"]),
                produced_rows_a=int(row["produced_rows_a"]),
                produced_rows_b=int(row["produced_rows_b"]),
                failed_rows_a=int(row["failed_rows_a"]),
                failed_rows_b=int(row["failed_rows_b"]),
                stream_rows_a=int(row["stream_rows_a"]),
                stream_rows_b=int(row["stream_rows_b"]),
            )
        )
    return out


def load_hybrid_counts(rows: Sequence[dict]) -> List[HybridCounts]:
    out: List[HybridCounts] = []
    for row in rows:
        out.append(
            HybridCounts(
                key=invocation_key_from_row(row),
                rows_seen_by_hybrid=int(row["rows_seen_by_hybrid"]),
                class_none_count=int(row["class_none_count"]),
                class_bch_hard_decoded_count=int(row["class_bch_hard_decoded_count"]),
                class_clean_count=int(row["class_clean_count"]),
                class_parity_only_count=int(row["class_parity_only_count"]),
                class_one_main_count=int(row["class_one_main_count"]),
                class_one_main_plus_parity_count=int(
                    row["class_one_main_plus_parity_count"]
                ),
                class_two_main_count=int(row["class_two_main_count"]),
                class_suspicious_count=int(row["class_suspicious_count"]),
                class_hard_fail_count=int(row["class_hard_fail_count"]),
                deferred_candidate_count=int(row["deferred_candidate_count"]),
                deferred_priority_0_count=int(row["deferred_priority_0_count"]),
                deferred_priority_1_count=int(row["deferred_priority_1_count"]),
                deferred_priority_2_count=int(row["deferred_priority_2_count"]),
                deferred_reclaimed_to_hard_finish_count=int(
                    row["deferred_reclaimed_to_hard_finish_count"]
                ),
            )
        )
    return out


def load_row_maps(rows: Sequence[dict]) -> List[RowMap]:
    out: List[RowMap] = []
    for row in rows:
        out.append(
            RowMap(
                key=invocation_key_from_row(row),
                merged_row=int(row["merged_row"]),
                stream_id=int(row["stream_id"]),
                stream_label=row["stream_label"],
                source_local_row=int(row["source_local_row"]),
                source_global_row=int(row["source_global_row"]),
                early_stop_hit=parse_bool_int(row["early_stop_hit"]),
                hybrid_class=row["hybrid_class"],
                final_tag=row["final_tag"],
                scheduled_for_soft=parse_bool_int(row["scheduled_for_soft"]),
                produced_row=parse_bool_int(row["produced_row"]),
            )
        )
    return out


def available_run_ids(
    seed_summaries: Sequence[SeedSummary], tile_samples: Sequence[TileSample]
) -> List[str]:
    run_ids = {item.run_id for item in seed_summaries}
    run_ids.update(item.key.run_id for item in tile_samples)
    return sorted(run_ids)


def select_run_ids(
    run_ids: Sequence[str], requested_run_id: str | None, all_runs: bool
) -> List[str]:
    if requested_run_id is not None:
        if requested_run_id not in run_ids:
            raise ValueError(f"run_id not found: {requested_run_id}")
        return [requested_run_id]
    if all_runs or len(run_ids) <= 1:
        return list(run_ids)
    return [run_ids[-1]]


def filter_seed_summaries(
    items: Sequence[SeedSummary],
    run_ids: Sequence[str],
    ebn0_db: str | None,
    seed_index: int | None,
) -> List[SeedSummary]:
    selected = set(run_ids)
    return [
        item
        for item in items
        if item.run_id in selected
        and (ebn0_db is None or item.ebn0_db == ebn0_db)
        and (seed_index is None or item.seed_index == seed_index)
    ]


def filter_invocation_items(
    items: Sequence[TileSample] | Sequence[HybridCounts] | Sequence[RowMap],
    run_ids: Sequence[str],
    ebn0_db: str | None,
    seed_index: int | None,
    tile_index: int | None,
) -> list:
    selected = set(run_ids)
    out = []
    for item in items:
        key = item.key
        if key.run_id not in selected:
            continue
        if ebn0_db is not None and key.ebn0_db != ebn0_db:
            continue
        if seed_index is not None and key.seed_index != seed_index:
            continue
        if tile_index is not None and key.tile_index != tile_index:
            continue
        out.append(item)
    return out


def safe_div(numerator: float, denominator: float) -> float:
    if denominator == 0:
        return 0.0
    return numerator / denominator


def invocation_key_sort_value(key: InvocationKey) -> Tuple[str, float, int, int, int, int]:
    return (
        key.run_id,
        float(key.ebn0_db),
        key.seed_index,
        key.chunk_index,
        key.invocation,
        key.tile_index,
    )


def build_row_map_stats(row_maps: Sequence[RowMap]) -> Dict[InvocationKey, dict]:
    stats: Dict[InvocationKey, dict] = {}
    for row in row_maps:
        bucket = stats.setdefault(
            row.key,
            {
                "row_map_rows": 0,
                "early_stop_a": 0,
                "early_stop_b": 0,
                "hard_finish_a": 0,
                "hard_finish_b": 0,
                "softdecode_a": 0,
                "softdecode_b": 0,
                "unscheduled_a": 0,
                "unscheduled_b": 0,
                "produced_a": 0,
                "produced_b": 0,
            },
        )
        bucket["row_map_rows"] += 1
        stream_suffix = "a" if row.stream_label == "A" else "b"
        if row.early_stop_hit:
            bucket[f"early_stop_{stream_suffix}"] += 1
        if row.final_tag == "hard_finish":
            bucket[f"hard_finish_{stream_suffix}"] += 1
        elif row.final_tag == "soft_decode":
            bucket[f"softdecode_{stream_suffix}"] += 1
        elif row.final_tag == "unscheduled":
            bucket[f"unscheduled_{stream_suffix}"] += 1
        if row.produced_row:
            bucket[f"produced_{stream_suffix}"] += 1
    return stats


def build_seed_ranking(seed_summaries: Sequence[SeedSummary]) -> List[dict]:
    rows: List[dict] = []
    for item in seed_summaries:
        rows.append(
            {
                "run_id": item.run_id,
                "ebn0_db": item.ebn0_db,
                "seed_index": item.seed_index,
                "chunk_index": item.chunk_index,
                "bitgen_seed_a": item.bitgen_seed_a,
                "channel_seed_a": item.channel_seed_a,
                "bitgen_seed_b": item.bitgen_seed_b,
                "channel_seed_b": item.channel_seed_b,
                "pre_ber_a": item.pre_ber_a,
                "post_ber_a": item.post_ber_a,
                "pre_ber_b": item.pre_ber_b,
                "post_ber_b": item.post_ber_b,
                "max_post_ber": max(item.post_ber_a, item.post_ber_b),
                "mean_post_ber": 0.5 * (item.post_ber_a + item.post_ber_b),
                "ab_post_ber_gap": abs(item.post_ber_a - item.post_ber_b),
                "tile_sample_rows": item.tile_sample_rows,
                "hybrid_count_rows": item.hybrid_count_rows,
                "row_map_rows": item.row_map_rows,
            }
        )
    rows.sort(
        key=lambda row: (
            row["max_post_ber"],
            row["ab_post_ber_gap"],
            row["mean_post_ber"],
        ),
        reverse=True,
    )
    return rows


def build_suspicious_invocations(
    tile_samples: Sequence[TileSample],
    hybrid_counts: Sequence[HybridCounts],
    row_map_stats: Dict[InvocationKey, dict],
) -> List[dict]:
    hybrid_by_key = {item.key: item for item in hybrid_counts}
    rows: List[dict] = []
    for sample in tile_samples:
        key = sample.key
        hybrid = hybrid_by_key.get(key)
        row_stats = row_map_stats.get(key, {})
        unscheduled_a = int(row_stats.get("unscheduled_a", 0))
        unscheduled_b = int(row_stats.get("unscheduled_b", 0))
        softdecode_a = int(row_stats.get("softdecode_a", 0))
        softdecode_b = int(row_stats.get("softdecode_b", 0))
        early_stop_a = int(row_stats.get("early_stop_a", 0))
        early_stop_b = int(row_stats.get("early_stop_b", 0))

        row = {
            "run_id": key.run_id,
            "ebn0_db": key.ebn0_db,
            "seed_index": key.seed_index,
            "chunk_index": key.chunk_index,
            "invocation": key.invocation,
            "tile_index": key.tile_index,
            "rows_total": sample.rows_total,
            "rows_early_stop": sample.rows_early_stop,
            "rows_not_early_stop": sample.rows_not_early_stop,
            "rows_hard_finish": sample.rows_hard_finish,
            "rows_need_siso_before_mux": sample.rows_need_siso_before_mux,
            "rows_soft_scheduled": sample.rows_soft_scheduled,
            "rows_unscheduled": sample.rows_unscheduled,
            "produced_rows": sample.produced_rows,
            "failed_rows": sample.failed_rows,
            "produced_rows_a": sample.produced_rows_a,
            "produced_rows_b": sample.produced_rows_b,
            "failed_rows_a": sample.failed_rows_a,
            "failed_rows_b": sample.failed_rows_b,
            "stream_rows_a": sample.stream_rows_a,
            "stream_rows_b": sample.stream_rows_b,
            "early_stop_ratio": safe_div(sample.rows_early_stop, sample.rows_total),
            "hard_finish_ratio": safe_div(sample.rows_hard_finish, sample.rows_total),
            "soft_pressure": safe_div(
                sample.rows_need_siso_before_mux, sample.rows_total
            ),
            "mux_drop_ratio": safe_div(
                sample.rows_unscheduled, sample.rows_need_siso_before_mux
            ),
            "produced_gap": abs(sample.produced_rows_a - sample.produced_rows_b),
            "failed_gap": abs(sample.failed_rows_a - sample.failed_rows_b),
            "early_stop_a": early_stop_a,
            "early_stop_b": early_stop_b,
            "early_stop_gap": abs(early_stop_a - early_stop_b),
            "softdecode_a": softdecode_a,
            "softdecode_b": softdecode_b,
            "softdecode_gap": abs(softdecode_a - softdecode_b),
            "unscheduled_a": unscheduled_a,
            "unscheduled_b": unscheduled_b,
            "unscheduled_gap": abs(unscheduled_a - unscheduled_b),
            "row_map_rows": int(row_stats.get("row_map_rows", 0)),
            "rows_seen_by_hybrid": 0,
            "class_clean_count": 0,
            "class_parity_only_count": 0,
            "class_one_main_count": 0,
            "class_one_main_plus_parity_count": 0,
            "class_two_main_count": 0,
            "class_hard_fail_count": 0,
            "deferred_candidate_count": 0,
            "deferred_priority_0_count": 0,
            "deferred_priority_1_count": 0,
            "deferred_priority_2_count": 0,
            "deferred_reclaimed_to_hard_finish_count": 0,
        }
        if hybrid is not None:
            row["rows_seen_by_hybrid"] = hybrid.rows_seen_by_hybrid
            row["class_clean_count"] = hybrid.class_clean_count
            row["class_parity_only_count"] = hybrid.class_parity_only_count
            row["class_one_main_count"] = hybrid.class_one_main_count
            row["class_one_main_plus_parity_count"] = (
                hybrid.class_one_main_plus_parity_count
            )
            row["class_two_main_count"] = hybrid.class_two_main_count
            row["class_hard_fail_count"] = hybrid.class_hard_fail_count
            row["deferred_candidate_count"] = hybrid.deferred_candidate_count
            row["deferred_priority_0_count"] = hybrid.deferred_priority_0_count
            row["deferred_priority_1_count"] = hybrid.deferred_priority_1_count
            row["deferred_priority_2_count"] = hybrid.deferred_priority_2_count
            row["deferred_reclaimed_to_hard_finish_count"] = (
                hybrid.deferred_reclaimed_to_hard_finish_count
            )
        rows.append(row)

    rows.sort(
        key=lambda row: (
            row["rows_unscheduled"],
            row["mux_drop_ratio"],
            row["produced_gap"],
            row["failed_rows"],
            row["unscheduled_gap"],
            row["class_hard_fail_count"],
            row["class_two_main_count"],
        ),
        reverse=True,
    )
    return rows


def build_tile_index_summary(invocations: Sequence[dict]) -> List[dict]:
    grouped: Dict[int, dict] = {}
    for row in invocations:
        bucket = grouped.setdefault(
            int(row["tile_index"]),
            {
                "tile_index": int(row["tile_index"]),
                "invocation_count": 0,
                "rows_total_sum": 0,
                "rows_early_stop_sum": 0,
                "rows_hard_finish_sum": 0,
                "rows_need_siso_before_mux_sum": 0,
                "rows_unscheduled_sum": 0,
                "produced_gap_sum": 0,
                "unscheduled_gap_sum": 0,
                "failed_rows_sum": 0,
            },
        )
        bucket["invocation_count"] += 1
        bucket["rows_total_sum"] += int(row["rows_total"])
        bucket["rows_early_stop_sum"] += int(row["rows_early_stop"])
        bucket["rows_hard_finish_sum"] += int(row["rows_hard_finish"])
        bucket["rows_need_siso_before_mux_sum"] += int(row["rows_need_siso_before_mux"])
        bucket["rows_unscheduled_sum"] += int(row["rows_unscheduled"])
        bucket["produced_gap_sum"] += int(row["produced_gap"])
        bucket["unscheduled_gap_sum"] += int(row["unscheduled_gap"])
        bucket["failed_rows_sum"] += int(row["failed_rows"])

    out: List[dict] = []
    for tile_index, bucket in grouped.items():
        inv_count = bucket["invocation_count"]
        rows_total_sum = bucket["rows_total_sum"]
        rows_need_sum = bucket["rows_need_siso_before_mux_sum"]
        out.append(
            {
                "tile_index": tile_index,
                "invocation_count": inv_count,
                "avg_rows_total": safe_div(rows_total_sum, inv_count),
                "avg_early_stop_ratio": safe_div(
                    bucket["rows_early_stop_sum"], rows_total_sum
                ),
                "avg_hard_finish_ratio": safe_div(
                    bucket["rows_hard_finish_sum"], rows_total_sum
                ),
                "avg_soft_pressure": safe_div(rows_need_sum, rows_total_sum),
                "avg_mux_drop_ratio": safe_div(
                    bucket["rows_unscheduled_sum"], rows_need_sum
                ),
                "avg_produced_gap": safe_div(bucket["produced_gap_sum"], inv_count),
                "avg_unscheduled_gap": safe_div(
                    bucket["unscheduled_gap_sum"], inv_count
                ),
                "avg_failed_rows": safe_div(bucket["failed_rows_sum"], inv_count),
                "rows_unscheduled_sum": bucket["rows_unscheduled_sum"],
            }
        )
    out.sort(
        key=lambda row: (
            row["avg_mux_drop_ratio"],
            row["avg_unscheduled_gap"],
            row["avg_produced_gap"],
            row["rows_unscheduled_sum"],
        ),
        reverse=True,
    )
    return out


def build_hybrid_class_summary(row_maps: Sequence[RowMap]) -> Tuple[List[dict], List[dict]]:
    class_stats: Dict[str, dict] = {}
    stream_class_stats: Dict[Tuple[str, str], dict] = {}

    for row in row_maps:
        class_bucket = class_stats.setdefault(
            row.hybrid_class,
            {
                "hybrid_class": row.hybrid_class,
                "row_count": 0,
                "stream_a_rows": 0,
                "stream_b_rows": 0,
                "early_stop_rows": 0,
                "hard_finish_rows": 0,
                "soft_decode_rows": 0,
                "unscheduled_rows": 0,
                "scheduled_for_soft_rows": 0,
                "produced_rows": 0,
            },
        )
        stream_key = (row.stream_label, row.hybrid_class)
        stream_bucket = stream_class_stats.setdefault(
            stream_key,
            {
                "stream_label": row.stream_label,
                "hybrid_class": row.hybrid_class,
                "row_count": 0,
                "early_stop_rows": 0,
                "hard_finish_rows": 0,
                "soft_decode_rows": 0,
                "unscheduled_rows": 0,
                "scheduled_for_soft_rows": 0,
                "produced_rows": 0,
            },
        )

        for bucket in (class_bucket, stream_bucket):
            bucket["row_count"] += 1
            if row.early_stop_hit:
                bucket["early_stop_rows"] += 1
            if row.final_tag == "hard_finish":
                bucket["hard_finish_rows"] += 1
            elif row.final_tag == "soft_decode":
                bucket["soft_decode_rows"] += 1
            elif row.final_tag == "unscheduled":
                bucket["unscheduled_rows"] += 1
            if row.scheduled_for_soft:
                bucket["scheduled_for_soft_rows"] += 1
            if row.produced_row:
                bucket["produced_rows"] += 1

        if row.stream_label == "A":
            class_bucket["stream_a_rows"] += 1
        else:
            class_bucket["stream_b_rows"] += 1

    def finalize(rows: Iterable[dict]) -> List[dict]:
        out: List[dict] = []
        for bucket in rows:
            count = bucket["row_count"]
            finalized = dict(bucket)
            finalized["produced_rate"] = safe_div(bucket["produced_rows"], count)
            finalized["soft_decode_rate"] = safe_div(bucket["soft_decode_rows"], count)
            finalized["unscheduled_rate"] = safe_div(bucket["unscheduled_rows"], count)
            finalized["hard_finish_rate"] = safe_div(bucket["hard_finish_rows"], count)
            finalized["scheduled_for_soft_rate"] = safe_div(
                bucket["scheduled_for_soft_rows"], count
            )
            out.append(finalized)
        out.sort(
            key=lambda row: (
                row["unscheduled_rows"],
                row["soft_decode_rows"],
                row["row_count"],
            ),
            reverse=True,
        )
        return out

    return finalize(class_stats.values()), finalize(stream_class_stats.values())


def write_csv(path: Path, rows: Sequence[dict]) -> None:
    if not rows:
        path.parent.mkdir(parents=True, exist_ok=True)
        with path.open("w", newline="") as handle:
            handle.write("")
        return
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0].keys()))
        writer.writeheader()
        writer.writerows(rows)


def print_rows(title: str, rows: Sequence[dict], columns: Sequence[str], limit: int) -> None:
    print(f"\n[{title}]")
    if not rows:
        print("  (empty)")
        return
    selected = rows[:limit]
    widths = {}
    for col in columns:
        widths[col] = max(len(col), *(len(format_value(row.get(col))) for row in selected))
    header = "  " + "  ".join(col.ljust(widths[col]) for col in columns)
    print(header)
    print("  " + "  ".join("-" * widths[col] for col in columns))
    for row in selected:
        print(
            "  "
            + "  ".join(format_value(row.get(col)).ljust(widths[col]) for col in columns)
        )


def format_value(value: object) -> str:
    if isinstance(value, float):
        return f"{value:.6f}"
    return str(value)


def main() -> int:
    args = parse_args()
    input_dir: Path = args.input_dir
    output_dir = args.output_dir or (input_dir / DEFAULT_OUTPUT_SUBDIR)

    seed_rows = read_csv_rows(input_dir / "per_seed_summary.csv")
    tile_rows = read_csv_rows(input_dir / "per_invocation_shared_tile_samples.csv")
    hybrid_rows = read_csv_rows(
        input_dir / "per_invocation_shared_hybrid_class_counts.csv"
    )
    row_map_rows = read_csv_rows(input_dir / "per_invocation_shared_row_map.csv")

    seed_summaries = load_seed_summaries(seed_rows)
    tile_samples = load_tile_samples(tile_rows)
    hybrid_counts = load_hybrid_counts(hybrid_rows)
    row_maps = load_row_maps(row_map_rows)

    run_ids = available_run_ids(seed_summaries, tile_samples)
    selected_run_ids = select_run_ids(run_ids, args.run_id, args.all_runs)

    seed_summaries = filter_seed_summaries(
        seed_summaries, selected_run_ids, args.ebn0_db, args.seed_index
    )
    tile_samples = filter_invocation_items(
        tile_samples,
        selected_run_ids,
        args.ebn0_db,
        args.seed_index,
        args.tile_index,
    )
    hybrid_counts = filter_invocation_items(
        hybrid_counts,
        selected_run_ids,
        args.ebn0_db,
        args.seed_index,
        args.tile_index,
    )
    row_maps = filter_invocation_items(
        row_maps,
        selected_run_ids,
        args.ebn0_db,
        args.seed_index,
        args.tile_index,
    )

    row_map_stats = build_row_map_stats(row_maps)
    seed_ranking = build_seed_ranking(seed_summaries)
    suspicious_invocations = build_suspicious_invocations(
        tile_samples, hybrid_counts, row_map_stats
    )
    tile_index_summary = build_tile_index_summary(suspicious_invocations)
    hybrid_class_summary, stream_class_summary = build_hybrid_class_summary(row_maps)

    write_csv(output_dir / "seed_ranking.csv", seed_ranking)
    write_csv(output_dir / "suspicious_invocations.csv", suspicious_invocations[: args.top])
    write_csv(output_dir / "tile_index_summary.csv", tile_index_summary)
    write_csv(output_dir / "hybrid_class_summary.csv", hybrid_class_summary)
    write_csv(output_dir / "stream_class_summary.csv", stream_class_summary)

    print("two-stream shared probe analysis")
    print(f"input_dir: {input_dir}")
    print(f"output_dir: {output_dir}")
    print(f"selected_run_ids: {', '.join(selected_run_ids) if selected_run_ids else '(none)'}")
    print(
        "filtered rows:"
        f" seeds={len(seed_summaries)}"
        f" tile_samples={len(tile_samples)}"
        f" hybrid_counts={len(hybrid_counts)}"
        f" row_map={len(row_maps)}"
    )

    print_rows(
        "Worst Seeds",
        seed_ranking,
        [
            "run_id",
            "ebn0_db",
            "seed_index",
            "post_ber_a",
            "post_ber_b",
            "max_post_ber",
            "ab_post_ber_gap",
        ],
        min(args.top, 10),
    )
    print_rows(
        "Most Suspicious Invocations",
        suspicious_invocations,
        [
            "run_id",
            "ebn0_db",
            "seed_index",
            "invocation",
            "tile_index",
            "rows_unscheduled",
            "mux_drop_ratio",
            "produced_gap",
            "unscheduled_gap",
            "class_hard_fail_count",
            "class_two_main_count",
            "deferred_candidate_count",
        ],
        min(args.top, 10),
    )
    print_rows(
        "Tile Index Summary",
        tile_index_summary,
        [
            "tile_index",
            "invocation_count",
            "avg_early_stop_ratio",
            "avg_hard_finish_ratio",
            "avg_soft_pressure",
            "avg_mux_drop_ratio",
            "avg_produced_gap",
            "avg_unscheduled_gap",
        ],
        min(args.top, 10),
    )
    print_rows(
        "Hybrid Class Summary",
        hybrid_class_summary,
        [
            "hybrid_class",
            "row_count",
            "soft_decode_rows",
            "unscheduled_rows",
            "hard_finish_rows",
            "produced_rows",
            "produced_rate",
        ],
        min(args.top, 10),
    )

    print("\nwritten files:")
    for name in [
        "seed_ranking.csv",
        "suspicious_invocations.csv",
        "tile_index_summary.csv",
        "hybrid_class_summary.csv",
        "stream_class_summary.csv",
    ]:
        print(f"  {output_dir / name}")

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
