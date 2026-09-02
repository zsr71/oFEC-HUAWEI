#!/usr/bin/env python3
"""Estimate the scheduling opportunity for a two-batch FIFO policy.

This script reads existing single-ordinary-batch traces. It reports observed
entry fragmentation exactly, and separately runs an optimistic static replay.
The replay is not a decoder simulation and must not be treated as a BER or
FIFO-depth result for the proposed policy.
"""

from __future__ import annotations

import csv
from collections import defaultdict, deque
from dataclasses import dataclass
from pathlib import Path
from typing import Iterable


ROOT = Path(__file__).resolve().parent
DATA = ROOT / "data" / "level56_schedule"
ENTRY_BUDGET = 8
INVALID_BATCH_ID = "18446744073709551615"


@dataclass(frozen=True)
class Source:
    label: str
    rounds: str


PRIMARY_SOURCES = (
    Source(
        "3.065 dB long R16000 0/8",
        "ofec_single_ebn03.065_nbits60014592_fifo1_rbuf16000_drain1_"
        "hiso0_siso8_schedgroup4_level56_schedule_rounds.csv",
    ),
    Source(
        "3.065 dB long R16000 8/8",
        "ofec_single_ebn03.065_nbits60014592_fifo1_rbuf16000_drain1_"
        "hiso8_siso8_schedgroup4_level56_schedule_rounds.csv",
    ),
    Source(
        "3.09 dB short R1600 0/8",
        "ofec_single_ebn03.09_fifo1_rbuf1600_drain1_hiso0_siso8_"
        "schedgroup4_eqobs1_level56_schedule_rounds.csv",
    ),
    Source(
        "3.06 dB short R1600 0/8",
        "ofec_single_ebn03.06_fifo1_rbuf1600_drain1_hiso0_siso8_"
        "schedgroup4_eqobs1_level56_schedule_rounds.csv",
    ),
    Source(
        "3.06 dB short R3200 0/8",
        "ofec_single_ebn03.06_fifo1_rbuf3200_drain1_hiso0_siso8_"
        "schedgroup4_eqobs1_level56_schedule_rounds.csv",
    ),
    Source(
        "3.065 dB long R6400 0/8",
        "ofec_single_ebn03.065_nbits60014592_fifo1_rbuf6400_drain1_"
        "hiso0_siso8_schedgroup4_level56_schedule_rounds.csv",
    ),
)


def times_path(rounds_path: Path) -> Path:
    name = rounds_path.name.replace("ofec_single_", "debug_L6_", 1)
    name = name.replace(
        "_level56_schedule_rounds.csv", "_level56_buffered_times.csv"
    )
    return rounds_path.with_name(name)


def load_rounds(rounds_path: Path):
    """Return positive-entry services and each batch's old service chunks."""
    by_invocation: dict[int, tuple[int, int, int]] = {}
    with rounds_path.open(newline="") as stream:
        for row in csv.DictReader(stream):
            invocation = int(row["invocation"])
            service_time = int(row["buffered_t"])
            batch_id = int(row["buffered_batch"].removeprefix("B"))
            entries = int(row["total_group_entries"])
            previous = by_invocation.get(invocation)
            if previous is None or entries > previous[2]:
                by_invocation[invocation] = (service_time, batch_id, entries)

    services: list[tuple[int, int, int, int]] = []
    chunks: dict[int, list[int]] = defaultdict(list)
    for invocation, (service_time, batch_id, entries) in sorted(
        by_invocation.items(), key=lambda item: item[1][0]
    ):
        if entries <= 0:
            continue
        services.append((invocation, service_time, batch_id, entries))
        chunks[batch_id].append(entries)
    return services, chunks


def load_times(path: Path):
    depths_by_invocation: dict[int, int] = {}
    arrivals: list[int] = []
    input_time_count = 0
    old_peak = 0
    old_input_end_depth = 0
    rows = []
    with path.open(newline="") as stream:
        rows = list(csv.DictReader(stream))

    for row in rows:
        arrived = row["arrived_batch"]
        if arrived and INVALID_BATCH_ID not in arrived:
            arrivals.append(int(arrived.removeprefix("B")))
            input_time_count = max(input_time_count, int(row["t"]) + 1)
            old_peak = max(old_peak, int(row["fifo_depth_before"]))
            old_input_end_depth = int(row["fifo_depth_after"])
        if row["ordinary_service_used"] == "1":
            invocation = int(row["ordinary_schedule_invocation"])
            depths_by_invocation[invocation] = int(row["fifo_depth_before"])

    old_drain_times = len(rows) - input_time_count
    return (
        arrivals,
        input_time_count,
        depths_by_invocation,
        old_peak,
        old_input_end_depth,
        old_drain_times,
    )


def observed_metrics(
    services: Iterable[tuple[int, int, int, int]],
    input_time_count: int,
    depths_by_invocation: dict[int, int],
):
    input_services = [item for item in services if item[1] < input_time_count]
    count = len(input_services)
    entries = sum(item[3] for item in input_services)
    under_budget = sum(item[3] < ENTRY_BUDGET for item in input_services)
    candidate = sum(
        item[3] < ENTRY_BUDGET
        and depths_by_invocation.get(item[0], 0) >= 2
        for item in input_services
    )
    return {
        "services": count,
        "entries": entries,
        "mean_entries": entries / count,
        "utilization": entries / (ENTRY_BUDGET * count),
        "under_budget": under_budget / count,
        "candidate": candidate / count,
        "unused_entries": ENTRY_BUDGET * count - entries,
    }


def optimistic_replay(arrivals: list[int], old_chunks: dict[int, list[int]]):
    """Replay old work with at most two batches and splittable second work."""
    chunks = {batch: deque(values) for batch, values in old_chunks.items()}
    fifo: deque[int] = deque()
    peak = 0
    service_times = 0

    def serve_one_time(arrival: int | None):
        nonlocal peak, service_times
        if arrival is not None:
            fifo.append(arrival)
        peak = max(peak, len(fifo))

        remaining = ENTRY_BUDGET
        ordinary_batches = 0
        index = 0
        while index < len(fifo) and ordinary_batches < 2 and remaining > 0:
            batch = fifo[index]
            work = chunks.get(batch, deque())

            # No positive-entry chunk in the old trace is used as a proxy for
            # a zero-ordinary-service retirement (normally FullEarlyStop).
            if not work:
                del fifo[index]
                continue

            required = work[0]
            used = min(required, remaining)
            remaining -= used
            ordinary_batches += 1

            if used == required:
                work.popleft()
            else:
                work[0] = required - used
            chunks[batch] = work

            if not work:
                del fifo[index]
            else:
                index += 1

        service_times += 1

    for batch in arrivals:
        serve_one_time(batch)
    input_end_depth = len(fifo)
    input_service_times = service_times

    while fifo:
        before = (len(fifo), sum(sum(values) for values in chunks.values()))
        serve_one_time(None)
        after = (len(fifo), sum(sum(values) for values in chunks.values()))
        if after == before:
            raise RuntimeError("optimistic replay made no progress")

    return {
        "peak": peak,
        "input_end_depth": input_end_depth,
        "drain_times": service_times - input_service_times,
    }


def print_observed_table():
    print("Observed fragmentation from old single-batch traces")
    print(
        "label|ordinary services|mean entries|utilization|under 8|"
        "candidate|unused entries"
    )
    for source in PRIMARY_SOURCES:
        rounds_path = DATA / source.rounds
        services, _ = load_rounds(rounds_path)
        (
            _,
            input_time_count,
            depths,
            _,
            _,
            _,
        ) = load_times(times_path(rounds_path))
        metrics = observed_metrics(services, input_time_count, depths)
        print(
            f"{source.label}|{metrics['services']}|"
            f"{metrics['mean_entries']:.3f}|"
            f"{metrics['utilization']:.2%}|"
            f"{metrics['under_budget']:.2%}|"
            f"{metrics['candidate']:.2%}|"
            f"{metrics['unused_entries']}"
        )


def seed_files():
    pattern = (
        "ofec_single_ebn03.065_nbits60014592_bitseed*_fifo1_rbuf16000_"
        "drain1_hiso[08]_siso8_g4entries8_schedgroup4_"
        "level56_schedule_rounds.csv"
    )
    return sorted(DATA.glob(pattern))


def print_seed_replay_table():
    print("\nOptimistic replay; not a decoder simulation")
    print(
        "seed|mode|old peak|predicted peak|old input-end|"
        "predicted input-end|old drain|predicted drain"
    )
    for rounds_path in seed_files():
        name = rounds_path.name
        seed = name.split("bitseed", 1)[1].split("_", 1)[0]
        mode = name.split("_hiso", 1)[1].split("_", 1)[0] + "/8"
        _, chunks = load_rounds(rounds_path)
        (
            arrivals,
            _,
            _,
            old_peak,
            old_input_end,
            old_drain,
        ) = load_times(times_path(rounds_path))
        predicted = optimistic_replay(arrivals, chunks)
        print(
            f"{seed}|{mode}|{old_peak}|{predicted['peak']}|"
            f"{old_input_end}|{predicted['input_end_depth']}|"
            f"{old_drain}|{predicted['drain_times']}"
        )


def main():
    print_observed_table()
    print_seed_replay_table()


if __name__ == "__main__":
    main()
