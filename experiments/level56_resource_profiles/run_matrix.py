#!/usr/bin/env python3
"""Fixed-seed paired Dual resource comparison; never overwrite existing runs."""
import argparse
import concurrent.futures
import hashlib
import os
from pathlib import Path
import subprocess
import sys

REPO = Path(__file__).resolve().parents[2]
PROFILES = {'1S1H': (8, 8), '2S1H': (16, 8), '1S2H': (8, 16)}


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--output', type=Path, required=True)
    parser.add_argument('--binary', type=Path,
                        default=REPO / 'build_prior_fix_release/ofec_single')
    parser.add_argument('--short', action='store_true')
    parser.add_argument('--workers', type=int, choices=[1, 2], default=2)
    parser.add_argument('--resume', action='store_true')
    args = parser.parse_args()
    root = args.output.resolve()
    root.mkdir(parents=True, exist_ok=True)
    binary = args.binary.resolve()
    digest = hashlib.sha256(binary.read_bytes()).hexdigest()
    manifest = root / 'binary.sha256'
    if manifest.exists() and manifest.read_text().strip() != digest:
        raise RuntimeError('Output directory belongs to a different binary')
    manifest.write_text(digest + '\n')
    (root / 'commit.txt').write_text(subprocess.check_output(
        ['git', 'rev-parse', 'HEAD'], cwd=REPO, text=True))
    (root / 'source.patch').write_text(subprocess.check_output(
        ['git', 'diff'], cwd=REPO, text=True))
    bits = 2344320 if args.short else 60014592
    matrix = [(e, r, p) for e in (['3.06'] if args.short else ['3.06', '3.07'])
              for r in ([0, 8, 16, 32] if args.short else [0, 32, 64, 128])
              for p in PROFILES]

    def run(point):
        e, r, p = point
        folder = root / f'ebn0{e}_rbuf{r}_{p}'
        folder.mkdir(exist_ok=True)
        log = folder / 'stdout.log'
        if log.exists():
            if args.resume and '[DONE]' in log.read_text() and (
                    folder / 'exit_code').read_text().strip() == '0':
                print(f'SKIP completed {folder.name}', flush=True)
                return
            raise RuntimeError(f'Refusing to overwrite {log}')
        s, h = PROFILES[p]
        env = dict(os.environ, OFEC_EBN0_DB=e, OFEC_NUM_INFO_BITS=str(bits),
                   OFEC_BITGEN_SEED='20260319', OFEC_CHANNEL_SEED='3182026',
                   LEVEL56_BUFFER_ROWS=str(r), LEVEL56_BUFFERED_FIFO_ENABLE='0',
                   LEVEL56_DUAL_WINDOW_ENABLE='1',
                   LEVEL56_BUFFERED_FIFO_DRAIN_AT_FRAME_END='0',
                   LEVEL56_DUAL_WINDOW_DRAIN_AT_FRAME_END='1',
                   LEVEL56_SHARED_SISO_ACTIVE=str(s),
                   LEVEL56_SHARED_HISO_ACTIVE=str(h),
                   LEVEL56_GROUP4_MAX_ENTRIES='8', LEVEL56_SCHEDULE_MODE='group4',
                   POST_FEC_ERROR_POSITIONS_PATH=str(folder / 'errors.csv'))
        # Save every effective environment variable affecting this experiment.
        (folder / 'environment.txt').write_text('\n'.join(
            f'{key}={value}' for key, value in sorted(env.items())
            if key.startswith(('OFEC_', 'LEVEL56_', 'POST_FEC_'))) + '\n')
        print(f'START {folder.name}', flush=True)
        with log.open('x') as out, (folder / 'timing.log').open('x') as err:
            result = subprocess.run(['/usr/bin/time', '-v', str(binary)],
                                    cwd=folder, env=env, stdout=out, stderr=err)
        (folder / 'exit_code').write_text(str(result.returncode) + '\n')
        if result.returncode:
            raise RuntimeError(f'Failed {folder.name}: {result.returncode}')
        print(f'DONE {folder.name}', flush=True)

    with concurrent.futures.ThreadPoolExecutor(max_workers=args.workers) as pool:
        list(pool.map(run, matrix))
    subprocess.run([sys.executable, str(Path(__file__).with_name('analyze.py')),
                    str(root)], check=True)


if __name__ == '__main__':
    main()
