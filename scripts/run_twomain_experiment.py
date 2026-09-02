#!/usr/bin/env python3
"""运行一个受版本控制的 TwoMain 单点实验并保存可复现清单。"""

from __future__ import annotations

import argparse
import datetime as dt
import hashlib
import json
import os
import platform
import re
import shlex
import socket
import subprocess
import sys
import time
from pathlib import Path
from typing import Any


class SpecError(RuntimeError):
    pass


def parse_scalar(text: str) -> Any:
    text = text.strip()
    if not text:
        raise SpecError("配置值不能为空")
    if text.startswith('"'):
        try:
            return json.loads(text)
        except json.JSONDecodeError as exc:
            raise SpecError(f"字符串格式错误：{text}") from exc
    if text.startswith("["):
        if not text.endswith("]"):
            raise SpecError(f"数组缺少右方括号：{text}")
        body = text[1:-1].strip()
        if not body:
            return []
        values: list[Any] = []
        token = ""
        quoted = False
        escaped = False
        for char in body + ",":
            if escaped:
                token += char
                escaped = False
            elif char == "\\" and quoted:
                token += char
                escaped = True
            elif char == '"':
                token += char
                quoted = not quoted
            elif char == "," and not quoted:
                values.append(parse_scalar(token.strip()))
                token = ""
            else:
                token += char
        return values
    if text == "true":
        return True
    if text == "false":
        return False
    if re.fullmatch(r"[-+]?\d+", text):
        return int(text)
    if re.fullmatch(r"[-+]?(?:\d+\.\d*|\d*\.\d+)(?:[eE][-+]?\d+)?", text):
        return float(text)
    raise SpecError(f"不支持的配置值：{text}")


def load_simple_toml(path: Path) -> dict[str, Any]:
    """读取本项目使用的 TOML 子集，避免依赖系统外的 Python 包。"""
    root: dict[str, Any] = {}
    section: dict[str, Any] = root
    for line_number, raw_line in enumerate(path.read_text(encoding="utf-8").splitlines(), 1):
        line = raw_line.strip()
        if not line or line.startswith("#"):
            continue
        if line.startswith("[") and line.endswith("]"):
            section_name = line[1:-1].strip()
            if not section_name or "." in section_name:
                raise SpecError(f"{path}:{line_number}: 不支持的节名 {line}")
            if section_name in root:
                raise SpecError(f"{path}:{line_number}: 重复节 {section_name}")
            section = {}
            root[section_name] = section
            continue
        if "=" not in line:
            raise SpecError(f"{path}:{line_number}: 缺少等号")
        key, value_text = line.split("=", 1)
        key = key.strip()
        if not re.fullmatch(r"[A-Za-z_][A-Za-z0-9_]*", key):
            raise SpecError(f"{path}:{line_number}: 非法键名 {key}")
        if key in section:
            raise SpecError(f"{path}:{line_number}: 重复键 {key}")
        section[key] = parse_scalar(value_text)
    return root


def format_toml_scalar(value: Any) -> str:
    if isinstance(value, bool):
        return "true" if value else "false"
    if isinstance(value, str):
        return json.dumps(value, ensure_ascii=False)
    if isinstance(value, (int, float)):
        return str(value)
    if isinstance(value, list):
        return "[" + ", ".join(format_toml_scalar(item) for item in value) + "]"
    raise SpecError(f"无法写入TOML的配置值：{value!r}")


def dump_simple_toml(spec: dict[str, Any]) -> str:
    lines: list[str] = []
    for key, value in spec.items():
        if not isinstance(value, dict):
            lines.append(f"{key} = {format_toml_scalar(value)}")
    for section_name, section in spec.items():
        if not isinstance(section, dict):
            continue
        if lines:
            lines.append("")
        lines.append(f"[{section_name}]")
        for key, value in section.items():
            lines.append(f"{key} = {format_toml_scalar(value)}")
    return "\n".join(lines) + "\n"


def require(spec: dict[str, Any], section: str, key: str, value_type: type) -> Any:
    try:
        value = spec[section][key]
    except (KeyError, TypeError) as exc:
        raise SpecError(f"缺少必填配置 [{section}] {key}") from exc
    if value_type is float and isinstance(value, (int, float)) and not isinstance(value, bool):
        return float(value)
    if value_type is int and isinstance(value, int) and not isinstance(value, bool):
        return value
    if not isinstance(value, value_type):
        raise SpecError(f"[{section}] {key} 的类型应为 {value_type.__name__}")
    return value


def validate_spec(spec: dict[str, Any]) -> None:
    if require(spec, "identity", "experiment_id", str) == "":
        raise SpecError("experiment_id 不能为空")
    require(spec, "identity", "scheme_id", str)
    require(spec, "channel", "ebn0_db", float)
    require(spec, "channel", "bitgen_seed", int)
    require(spec, "channel", "channel_seed", int)
    nbits = require(spec, "frame", "num_info_bits", int)
    if nbits <= 0 or nbits % (32 * 111) != 0:
        raise SpecError("num_info_bits 必须是 3552 的正整数倍")
    mask = require(spec, "twomain", "hiso_class_mask", int)
    if not 0 <= mask <= 15:
        raise SpecError("hiso_class_mask 必须位于 0 到 15")
    output_policy = spec["twomain"].get("output_policy", "legacy_fixed")
    if output_policy not in {"legacy_fixed", "scheme6_parameterized"}:
        raise SpecError(
            "[twomain] output_policy 只能是 legacy_fixed 或 scheme6_parameterized")
    if output_policy == "scheme6_parameterized":
        m2 = require(spec, "twomain", "m2", float)
        rho_corr = require(spec, "twomain", "rho_corr", float)
        rho_keep = require(spec, "twomain", "rho_keep", float)
        if m2 < 0.0:
            raise SpecError("[twomain] m2 必须为非负数")
        if not 0.0 <= rho_corr <= 1.0 or not 0.0 <= rho_keep <= 1.0:
            raise SpecError("[twomain] rho_corr/rho_keep 必须位于 [0,1]")
    for key in ("buffer_rows", "hiso_active", "siso_active", "group4_max_entries"):
        if require(spec, "level56", key, int) < 0:
            raise SpecError(f"[level56] {key} 不能为负数")
    if require(spec, "level56", "schedule_mode", str) not in {"group4", "global"}:
        raise SpecError("schedule_mode 只能是 group4 或 global")
    require(spec, "level56", "buffered_fifo_enabled", bool)
    require(spec, "level56", "drain_at_frame_end", bool)
    require(spec, "outputs", "schedule_stats", bool)
    require(spec, "outputs", "post_fec_error_positions", bool)

    # 阶段二只把现有 ofec_single 已支持的环境变量交给运行时覆盖。
    # 其他公共字段必须与冻结基线默认值一致，防止配置声称采用某个值，
    # 实际程序却静默沿用另一个编译期常量。
    fixed_baseline = {
        ("channel", "generate_random_bits"): True,
        ("channel", "bits_per_symbol"): 1,
        ("frame", "ber_skip_prefix_windows"): 4,
        ("frame", "ber_skip_suffix_windows"): 2,
        ("decoder", "name"): "chase_baseline",
        ("decoder", "chase_l"): 6,
        ("decoder", "chase_n_test"): 64,
        ("decoder", "chase_topk_keep"): 24,
        ("decoder", "llr_bits"): 6,
        ("decoder", "quant_clip_ratio"): 0.5,
        ("decoder", "normalize_extrinsic"): False,
        ("decoder", "normalize_known_prefix_tail"): False,
        ("early_stop", "enabled"): True,
        ("early_stop", "per_tile_enabled"): [1, 1, 1, 1, 1, 1],
        ("early_stop", "condition_mode"): 1,
        ("early_stop", "action_mode"): 7,
        ("early_stop", "bind_group_size"): 1,
        ("early_stop", "condition_v1_require_bch"): True,
        ("early_stop", "condition_v1_require_overall"): True,
        ("early_stop", "action_beta"): [99.857143, 99.179301, 99.253626, 99.119585, 99.0, 99.0],
        ("level56", "shared_enabled"): True,
        ("level56", "temporal_lookahead_enabled"): False,
        ("level56", "priority_mode"): "level5_first",
        ("level56", "early_stop_group_update_mode"): "all_groups",
        ("level56", "unselected_early_stop_action_enabled"): True,
        ("outputs", "schedule_stats"): True,
        ("outputs", "target_trace"): False,
        ("outputs", "post_fec_error_positions"): True,
    }
    for (section, key), wanted in fixed_baseline.items():
        actual = spec.get(section, {}).get(key)
        if actual != wanted:
            raise SpecError(
                f"阶段二运行器尚不能覆盖 [{section}] {key}；"
                f"当前必须为冻结基线值 {wanted!r}，配置实际为 {actual!r}")


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def run_text(command: list[str], cwd: Path, *, check: bool = True) -> str:
    completed = subprocess.run(
        command, cwd=cwd, check=check, text=True,
        stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
    return completed.stdout.strip()


def git_metadata(repo: Path) -> dict[str, Any]:
    status = run_text(["git", "status", "--porcelain=v2"], repo)
    branch = run_text(["git", "branch", "--show-current"], repo)
    commit = run_text(["git", "rev-parse", "HEAD"], repo)
    tags_text = run_text(["git", "tag", "--points-at", "HEAD"], repo)
    return {
        "branch": branch,
        "commit_sha": commit,
        "tags_at_head": tags_text.splitlines() if tags_text else [],
        "worktree_clean": status == "",
        "status_porcelain_v2": status.splitlines() if status else [],
    }


def cmake_build_type(repo: Path, binary: Path) -> str:
    try:
        build_dir = binary.parent
        cache = (build_dir / "CMakeCache.txt").read_text(
            encoding="utf-8", errors="replace")
    except OSError:
        return "未知"
    match = re.search(r"^CMAKE_BUILD_TYPE:STRING=(.*)$", cache, re.MULTILINE)
    return match.group(1) if match and match.group(1) else "未指定"


def environment_from_spec(spec: dict[str, Any], run_dir: Path) -> dict[str, str]:
    channel = spec["channel"]
    frame = spec["frame"]
    level56 = spec["level56"]
    twomain = spec["twomain"]
    outputs = spec["outputs"]
    environment = {
        "OFEC_EBN0_DB": str(channel["ebn0_db"]),
        "OFEC_NUM_INFO_BITS": str(frame["num_info_bits"]),
        "OFEC_BITGEN_SEED": str(channel["bitgen_seed"]),
        "OFEC_CHANNEL_SEED": str(channel["channel_seed"]),
        "LEVEL56_BUFFERED_FIFO_ENABLE": "1" if level56["buffered_fifo_enabled"] else "0",
        "LEVEL56_BUFFER_ROWS": str(level56["buffer_rows"]),
        "LEVEL56_BUFFERED_FIFO_DRAIN_AT_FRAME_END": "1" if level56["drain_at_frame_end"] else "0",
        "LEVEL56_HISO_CLASS_MASK": str(twomain["hiso_class_mask"]),
        "LEVEL56_SHARED_HISO_ACTIVE": str(level56["hiso_active"]),
        "LEVEL56_SHARED_SISO_ACTIVE": str(level56["siso_active"]),
        "LEVEL56_GROUP4_MAX_ENTRIES": str(level56["group4_max_entries"]),
        "LEVEL56_SCHEDULE_MODE": str(level56["schedule_mode"]),
        "LEVEL56_EQUIVALENCE_OBSERVATION": "1" if outputs.get("equivalence_observation", False) else "0",
        "POST_FEC_ERROR_POSITIONS_PATH": str(run_dir / "post_fec_error_positions.txt"),
    }
    magnitudes = twomain.get("hybrid_hard_llr_mag_per_tile")
    if magnitudes is not None:
        environment["OFEC_HYBRID_HARD_LLR_MAG_LIST"] = ",".join(str(x) for x in magnitudes)
    if twomain.get("output_policy", "legacy_fixed") == "scheme6_parameterized":
        environment["OFEC_TWOMAIN_HISO_OUTPUT_MODE"] = "parameterized"
        environment["OFEC_TWOMAIN_HISO_M2"] = str(twomain["m2"])
        environment["OFEC_TWOMAIN_HISO_RHO_CORR"] = str(twomain["rho_corr"])
        environment["OFEC_TWOMAIN_HISO_RHO_KEEP"] = str(twomain["rho_keep"])
    else:
        environment["OFEC_TWOMAIN_HISO_OUTPUT_MODE"] = "legacy"
    return environment


RESULT_RE = re.compile(
    r"\[RESULT\] Pre-FEC BER=(?P<pre_ber>[^ ]+) \(errs=(?P<pre_errors>\d+)/(?P<pre_total>\d+)\)"
    r".*?Post-FEC BER=(?P<post_ber>[^ ]+) \(errs=(?P<post_errors>\d+)/(?P<post_total>\d+)\)")
FIFO_RE = re.compile(
    r"\[INFO\] Level56 buffered FIFO: t=(?P<times>\d+), completed=(?P<completed>\d+), "
    r"full_early_stop=(?P<full_early_stop>\d+), forced_evicted=(?P<forced>\d+), "
    r"max_fifo_depth=(?P<max_depth>\d+), min_S_t=(?P<min_s>\d+)")
ACTION_RE = re.compile(
    r"\[RESULT\] Level(?P<level>[56]) actions\(EarlyStop,HISO,SISO,Unscheduled\) = "
    r"\[(?P<values>\d+,\d+,\d+,\d+)\]")
CLASS_RE = re.compile(
    r"\[RESULT\] Level56 class (?P<name>[A-Za-z0-9]+) actions\(HISO,SISO,Unscheduled\) = "
    r"\[(?P<values>\d+,\d+,\d+)\]")
ENTRY_RE = re.compile(
    r"\[RESULT\] Level56 entry/HISO/SISO/idle-HISO/idle-SISO totals = "
    r"(?P<values>\d+/\d+/\d+/\d+/\d+)")


def parse_result(stdout: str, positions_path: Path, return_code: int) -> dict[str, Any]:
    result: dict[str, Any] = {"return_code": return_code}
    match = RESULT_RE.search(stdout)
    if match:
        result["pre_fec"] = {
            "ber": float(match.group("pre_ber")),
            "errors": int(match.group("pre_errors")),
            "total": int(match.group("pre_total")),
        }
        result["post_fec"] = {
            "ber": float(match.group("post_ber")),
            "errors": int(match.group("post_errors")),
            "total": int(match.group("post_total")),
        }
    fifo = FIFO_RE.search(stdout)
    if fifo:
        result["fifo"] = {
            "service_times": int(fifo.group("times")),
            "completed_batches": int(fifo.group("completed")),
            "full_early_stop_batches": int(fifo.group("full_early_stop")),
            "forced_evicted": int(fifo.group("forced")),
            "max_fifo_depth": int(fifo.group("max_depth")),
            "min_window_start": int(fifo.group("min_s")),
        }
    result["level_actions"] = {
        match.group("level"): [int(value) for value in match.group("values").split(",")]
        for match in ACTION_RE.finditer(stdout)
    }
    result["class_actions"] = {
        match.group("name"): [int(value) for value in match.group("values").split(",")]
        for match in CLASS_RE.finditer(stdout)
    }
    entry = ENTRY_RE.search(stdout)
    if entry:
        result["entry_hiso_siso_idle_totals"] = [
            int(value) for value in entry.group("values").split("/")]
    if positions_path.exists():
        result["post_fec_positions"] = {
            "count": sum(1 for line in positions_path.read_text(encoding="utf-8").splitlines() if line.strip()),
            "sha256": sha256_file(positions_path),
        }
    return result


def verify_expected(spec: dict[str, Any], result: dict[str, Any]) -> list[str]:
    expected = spec.get("expected", {})
    failures: list[str] = []
    comparisons = [
        ("post_fec_errors", result.get("post_fec", {}).get("errors")),
        ("post_fec_total", result.get("post_fec", {}).get("total")),
        ("post_fec_positions_sha256", result.get("post_fec_positions", {}).get("sha256")),
        ("forced_evicted", result.get("fifo", {}).get("forced_evicted")),
        ("completed_batches", result.get("fifo", {}).get("completed_batches")),
        ("max_fifo_depth", result.get("fifo", {}).get("max_fifo_depth")),
        ("level5_actions", result.get("level_actions", {}).get("5")),
        ("level6_actions", result.get("level_actions", {}).get("6")),
        ("entry_hiso_siso_idle_totals", result.get("entry_hiso_siso_idle_totals")),
        ("class_parityonly_actions", result.get("class_actions", {}).get("ParityOnly")),
        ("class_onemain_actions", result.get("class_actions", {}).get("OneMain")),
        ("class_onemainplusparity_actions", result.get("class_actions", {}).get("OneMainPlusParity")),
        ("class_twomain_actions", result.get("class_actions", {}).get("TwoMain")),
        ("class_hardfail_actions", result.get("class_actions", {}).get("HardFail")),
    ]
    for key, actual in comparisons:
        wanted = expected.get(key)
        if wanted is None or wanted == -1 or wanted == "待冻结":
            continue
        if isinstance(wanted, list) and wanted and all(value == -1 for value in wanted):
            continue
        if wanted != actual:
            failures.append(f"{key}: 期望 {wanted!r}，实际 {actual!r}")
    return failures


def main() -> int:
    parser = argparse.ArgumentParser(description="运行受版本控制的 TwoMain 实验")
    parser.add_argument("--spec", required=True, type=Path, help="实验 TOML 配置")
    parser.add_argument("--binary", type=Path, help="ofec_single 路径，默认 build/ofec_single")
    parser.add_argument("--runs-root", type=Path, help="原始运行根目录，默认 runs/twomain")
    parser.add_argument("--run-id", help="显式运行编号；默认使用时间和提交短号")
    parser.add_argument("--bitgen-seed", type=int,
                        help="覆盖配置中的信息比特随机种子；必须与--channel-seed同时使用")
    parser.add_argument("--channel-seed", type=int,
                        help="覆盖配置中的信道随机种子；必须与--bitgen-seed同时使用")
    parser.add_argument("--dry-run", action="store_true", help="只校验并打印，不执行")
    parser.add_argument("--verify-expected", action="store_true", help="核对配置内冻结签名")
    parser.add_argument("--allow-dirty", action="store_true", help="允许脏工作区调试运行")
    args = parser.parse_args()

    repo = Path(__file__).resolve().parents[1]
    spec_path = args.spec if args.spec.is_absolute() else (repo / args.spec)
    spec_path = spec_path.resolve()
    if not spec_path.is_file():
        raise SpecError(f"配置文件不存在：{spec_path}")
    spec = load_simple_toml(spec_path)
    if (args.bitgen_seed is None) != (args.channel_seed is None):
        raise SpecError("--bitgen-seed和--channel-seed必须成对使用")
    requested_overrides: dict[str, Any] = {}
    if args.bitgen_seed is not None and args.channel_seed is not None:
        spec["channel"]["bitgen_seed"] = args.bitgen_seed
        spec["channel"]["channel_seed"] = args.channel_seed
        requested_overrides = {
            "bitgen_seed": args.bitgen_seed,
            "channel_seed": args.channel_seed,
        }
    validate_spec(spec)
    git = git_metadata(repo)
    if not git["worktree_clean"] and not args.allow_dirty and not args.dry_run:
        raise SpecError("工作区存在未提交修改；正式运行前请提交，或调试时显式使用 --allow-dirty")

    binary = args.binary or (repo / "build" / "ofec_single")
    binary = binary if binary.is_absolute() else (repo / binary)
    binary = binary.resolve()
    runs_root = args.runs_root or (repo / "runs" / "twomain")
    runs_root = runs_root if runs_root.is_absolute() else (repo / runs_root)
    experiment_id = spec["identity"]["experiment_id"]
    timestamp = dt.datetime.now().astimezone().strftime("%Y%m%d_%H%M%S")
    run_id = args.run_id or f"{timestamp}_{git['commit_sha'][:8]}"
    run_dir = (runs_root / experiment_id / run_id).resolve()
    environment = environment_from_spec(spec, run_dir)

    if args.dry_run:
        print(f"配置：{spec_path}")
        print(f"方案：{spec['identity']['scheme_id']} / {spec['identity']['name']}")
        print(f"可执行文件：{binary}")
        print(f"运行目录：{run_dir}")
        print("环境变量：")
        for key in sorted(environment):
            print(f"  {key}={environment[key]}")
        print(f"工作区干净：{'是' if git['worktree_clean'] else '否'}")
        return 0

    if not binary.is_file() or not os.access(binary, os.X_OK):
        raise SpecError(f"找不到可执行文件：{binary}；请先构建 ofec_single")
    if run_dir.exists():
        raise SpecError(f"运行目录已经存在：{run_dir}")
    run_dir.mkdir(parents=True)
    (run_dir / "effective_config.toml").write_text(
        dump_simple_toml(spec), encoding="utf-8")

    manifest = {
        "schema_version": 1,
        "experiment_id": experiment_id,
        "scheme_id": spec["identity"]["scheme_id"],
        "run_id": run_id,
        "status": "running",
        "start_time": dt.datetime.now().astimezone().isoformat(),
        "repo_root": str(repo),
        "spec_path": str(spec_path.relative_to(repo)),
        "spec_sha256": sha256_file(spec_path),
        "requested_overrides": requested_overrides,
        "git": git,
        "host": {
            "hostname": socket.gethostname(),
            "platform": platform.platform(),
            "python": platform.python_version(),
        },
        "build_tools": {
            "build_type": cmake_build_type(repo, binary),
            "cmake": run_text(["cmake", "--version"], repo).splitlines()[0],
            "compiler": run_text(["c++", "--version"], repo).splitlines()[0],
        },
        "binary": {"path": str(binary), "sha256": sha256_file(binary)},
        "environment_overrides": environment,
        "effective_config": spec,
    }
    manifest_path = run_dir / "manifest.json"
    manifest_path.write_text(json.dumps(manifest, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")

    child_env = os.environ.copy()
    child_env.update(environment)
    command = [str(binary)]
    stdout_path = run_dir / "stdout.log"
    started = time.monotonic()
    with stdout_path.open("w", encoding="utf-8") as output:
        output.write("[实验运行器] 命令：" + shlex.join(command) + "\n")
        output.write("[实验运行器] 工作目录：" + str(run_dir) + "\n")
        output.flush()
        completed = subprocess.run(
            command, cwd=run_dir, env=child_env, text=True,
            stdout=output, stderr=subprocess.STDOUT)
    elapsed = time.monotonic() - started
    stdout = stdout_path.read_text(encoding="utf-8", errors="replace")
    result = parse_result(stdout, run_dir / "post_fec_error_positions.txt", completed.returncode)
    failures = verify_expected(spec, result) if args.verify_expected else []
    result["expected_verification"] = {
        "requested": args.verify_expected,
        "passed": not failures if args.verify_expected else None,
        "failures": failures,
    }
    result_path = run_dir / "result.json"
    result_path.write_text(json.dumps(result, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")
    manifest["status"] = "completed" if completed.returncode == 0 else "failed"
    manifest["end_time"] = dt.datetime.now().astimezone().isoformat()
    manifest["elapsed_seconds"] = elapsed
    manifest["return_code"] = completed.returncode
    manifest["expected_verification_passed"] = result["expected_verification"]["passed"]
    manifest_path.write_text(json.dumps(manifest, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")

    print(f"运行目录：{run_dir}")
    if "post_fec" in result:
        post = result["post_fec"]
        print(f"Post-FEC：{post['errors']} / {post['total']} = {post['ber']}")
    if "fifo" in result:
        fifo = result["fifo"]
        print(
            "FIFO：完成批次 {completed_batches}，强制顶出 {forced_evicted}，最大深度 {max_fifo_depth}".format(**fifo))
    if failures:
        print("冻结签名核对失败：", file=sys.stderr)
        for failure in failures:
            print(f"  - {failure}", file=sys.stderr)
        return 3
    if completed.returncode != 0:
        print(f"程序返回非零状态：{completed.returncode}", file=sys.stderr)
        return completed.returncode
    if args.verify_expected:
        print("冻结签名核对通过")
    return 0


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except SpecError as exc:
        print(f"[错误] {exc}", file=sys.stderr)
        raise SystemExit(2)
