#!/usr/bin/env python3
"""执行冻结的 TwoMain 默认行为快速等价回归。"""

from __future__ import annotations

import argparse
import subprocess
import sys
from pathlib import Path


def main() -> int:
    parser = argparse.ArgumentParser(description="TwoMain 默认行为快速等价回归")
    parser.add_argument("--binary", type=Path, help="ofec_single 路径")
    parser.add_argument("--allow-dirty", action="store_true", help="允许在脏工作区运行")
    args = parser.parse_args()
    repo = Path(__file__).resolve().parents[1]
    runner = repo / "scripts" / "run_twomain_experiment.py"
    spec = repo / "experiments" / "twomain" / "specs" / "regression_default_quick.toml"
    command = [sys.executable, str(runner), "--spec", str(spec), "--verify-expected"]
    if args.binary:
        command.extend(["--binary", str(args.binary)])
    if args.allow_dirty:
        command.append("--allow-dirty")
    print("开始执行 TwoMain 默认行为快速等价回归")
    completed = subprocess.run(command, cwd=repo)
    if completed.returncode == 0:
        print("默认行为快速等价回归通过")
    else:
        print("默认行为快速等价回归失败", file=sys.stderr)
    return completed.returncode


if __name__ == "__main__":
    raise SystemExit(main())
