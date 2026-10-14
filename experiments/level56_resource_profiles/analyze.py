#!/usr/bin/env python3
import csv
import json
import re
import sys
from pathlib import Path


def parse_run(folder):
    log = (folder / 'stdout.log').read_text(errors='replace')
    m = re.findall(r'Post-FEC BER=([0-9eE.+-]+)\s+\(errs=(\d+)\s*/\s*(\d+)', log)
    if not m:
        raise RuntimeError(f'no Post-FEC result in {folder}')
    _, errors, compared = m[-1]
    if '[DONE]' not in log or (folder / 'exit_code').read_text().strip() != '0':
        raise RuntimeError(f'run not complete: {folder}')
    parts = folder.name.split('_')
    profile = parts[-1]
    rbuf = int(parts[1][4:])
    sw, hw = {'1S1H': (1, 1), '2S1H': (2, 1), '1S2H': (1, 2)}[profile]
    time_files = list((folder / 'data/level56_schedule').glob('*dual_window_times.csv'))
    if len(time_files) != 1:
        raise RuntimeError(f'expected one time csv in {folder}')
    sums = dict(input_batches=0, output_batches=0, c5=0, c6=0,
                forced5=0, forced6=0, hiso=0, siso=0, max_entries=0,
                max_hiso=0, max_siso=0, times=0, entries=0,
                pending5_code_times=0, pending6_code_times=0,
                pending5_times=0, pending6_times=0, min_s5=rbuf+22, min_s6=rbuf,
                max_debt5=0, max_debt6=0, full_early_stop5=0, full_early_stop6=0)
    with time_files[0].open() as f:
        for row in csv.DictReader(f):
            sums['times'] += 1
            sums['input_batches'] += int(row['input_valid'])
            sums['output_batches'] += int(row['output_valid'])
            sums['c5'] += int(row['C5']); sums['c6'] += int(row['C6'])
            sums['forced5'] += int(row['forced_evicted5'])
            sums['forced6'] += int(row['forced_evicted6'])
            sums['hiso'] += int(row.get('hiso_calls', 0))
            sums['siso'] += int(row.get('siso_calls', 0))
            sums['max_entries'] = max(sums['max_entries'], int(row['entries_used']))
            sums['max_hiso'] = max(sums['max_hiso'], int(row.get('hiso_calls', 0)))
            sums['max_siso'] = max(sums['max_siso'], int(row.get('siso_calls', 0)))
            sums['entries'] += int(row['entries_used'])
            assert int(row['entries_used']) <= 8
            assert int(row['hiso_calls']) <= 8 * hw
            assert int(row['siso_calls']) <= 8 * sw
            for level in (5, 6):
                pending = int(row[f'pending{level}'])
                sums[f'pending{level}_code_times'] += pending
                sums[f'pending{level}_times'] += pending > 0
                sums[f'max_debt{level}'] = max(sums[f'max_debt{level}'],
                                             int(row[f'recovery_debt{level}']))
                sums[f'full_early_stop{level}'] += int(row[f'full_early_stop{level}'])
            for suffix in ('before', 'after'):
                s5, s6 = int(row[f'S5_{suffix}']), int(row[f'S6_{suffix}'])
                assert s5 % 2 == s6 % 2 == 0
                assert 0 <= s6 and s6 + 22 <= s5 and s5 + 22 <= rbuf + 44
                sums['min_s5'] = min(sums['min_s5'], s5)
                sums['min_s6'] = min(sums['min_s6'], s6)
    if sums['input_batches'] != sums['output_batches']:
        raise RuntimeError(f'input/output mismatch in {folder}: {sums}')
    if sums['c5'] + sums['forced5'] != sums['input_batches']:
        raise RuntimeError(f'L5 conservation failure in {folder}: {sums}')
    if sums['c6'] + sums['forced6'] != sums['input_batches']:
        raise RuntimeError(f'L6 conservation failure in {folder}: {sums}')
    steps_file = next((folder / 'data/level56_schedule').glob('*dual_window_steps.csv'))
    for row in csv.DictReader(steps_file.open()):
        s5, s6 = int(row['S5']), int(row['S6'])
        assert s5 % 2 == s6 % 2 == 0
        assert 0 <= s6 and s6 + 22 <= s5 and s5 + 22 <= rbuf + 44
        assert int(row['entries_after']) <= 8
    # Each scheduler invocation may have multiple round CSV rows. Count once.
    seen = set()
    planned_h = planned_s = consumed = 0
    rounds_file = next((folder / 'data/level56_schedule').glob('*schedule_rounds.csv'))
    for row in csv.DictReader(rounds_file.open()):
        if row['invocation'] in seen:
            continue
        seen.add(row['invocation'])
        planned_h += int(row['planned_hiso'])
        planned_s += int(row['planned_siso'])
        consumed += int(row['consumed_entries'])
    assert (planned_h, planned_s, consumed) == (sums['hiso'], sums['siso'], sums['entries'])
    timing = (folder / 'timing.log').read_text()
    assert re.search(r'Exit status: 0\s*$', timing)
    seconds = float(re.search(r'User time \(seconds\): ([\d.]+)', timing)[1])
    wall = re.search(r'Elapsed \(wall clock\) time.*: ([\d:.]+)', timing)[1]
    chunks = wall.split(':')
    wall_seconds = sum(float(v) * 60 ** i for i, v in enumerate(reversed(chunks)))
    return dict(folder=str(folder), ebn0=float(parts[0][4:]),
                rbuf=rbuf, profile=profile,
                post_fec_ber=int(errors) / int(compared), post_fec_errors=int(errors),
                compared_bits=int(compared), user_seconds=seconds,
                wall_seconds=wall_seconds,
                peak_rss_kib=int(re.search(r'Maximum resident set size \(kbytes\): (\d+)', timing)[1]),
                mean_hiso=sums['hiso']/sums['times'], mean_siso=sums['siso']/sums['times'],
                hiso_utilization=sums['hiso']/(sums['times']*8*hw),
                siso_utilization=sums['siso']/(sums['times']*8*sw), **sums)


def main():
    root = Path(sys.argv[1] if len(sys.argv) > 1 else '.')
    def complete(p):
        return (p / 'exit_code').exists() and (p / 'exit_code').read_text().strip() == '0'
    rows = [parse_run(p) for p in sorted(root.glob('ebn0*_rbuf*_1S1H')) if complete(p)]
    for profile in ('2S1H', '1S2H'):
        rows.extend(parse_run(p) for p in sorted(root.glob(f'ebn0*_rbuf*_{profile}')) if complete(p))
    rows.sort(key=lambda x: (x['ebn0'], x['rbuf'], x['profile']))
    (root / 'results.json').write_text(json.dumps(rows, indent=2) + '\n')
    if rows:
        with (root / 'results.csv').open('w', newline='') as f:
            w = csv.DictWriter(f, fieldnames=list(rows[0]))
            w.writeheader(); w.writerows(rows)
    print(json.dumps(rows, indent=2))


if __name__ == '__main__':
    main()
