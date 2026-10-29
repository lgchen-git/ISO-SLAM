#!/usr/bin/env python3
"""Run the bulk long-read SLAM-seq workflow branch."""

from __future__ import annotations

import argparse
import csv
import fcntl
import json
import shlex
import shutil
import subprocess
import sys
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path


SCRIPT_DIR = Path(__file__).resolve().parent
READ_FLOW_SCHEMA = "iso-slam-read-flow-v1"
REPORT_SCHEMA = "iso-slam-sample-report-v3"


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Run bulk long-read SLAM-seq analysis.")
    parser.add_argument("--config", required=True)
    parser.add_argument("--samples", required=True)
    parser.add_argument("--results-dir", required=True)
    parser.add_argument("--logs-dir", required=True)
    parser.add_argument("--tmp-dir", required=True)
    parser.add_argument("--threads", type=int, default=4)
    parser.add_argument("--run-mode", choices=["pilot", "full"], default="")
    parser.add_argument("--max-parallel-samples", type=int, default=0)
    return parser.parse_args()


def read_config(path: str) -> dict[str, str]:
    config: dict[str, str] = {}
    with open(path, "r", encoding="utf-8") as handle:
        for raw in handle:
            line = raw.strip()
            if not line or line.startswith("#") or raw[:1].isspace() or ":" not in raw:
                continue
            key, value = raw.split(":", 1)
            value = value.split("#", 1)[0].strip()
            if len(value) >= 2 and value[0] == value[-1] and value[0] in {"'", '"'}:
                value = value[1:-1]
            config[key.strip()] = value
    return config


def bool_config(config: dict[str, str], key: str, default: bool) -> bool:
    value = config.get(key, "")
    if not value:
        return default
    return value.lower() in {"1", "true", "yes", "on"}


def int_config(config: dict[str, str], key: str, default: int) -> int:
    value = config.get(key, "")
    return int(value) if value else default


def float_config(config: dict[str, str], key: str, default: float | None) -> float | None:
    value = config.get(key, "")
    return float(value) if value else default


def run_step(name: str, cmd: list[str], log_path: Path, lock_path: Path | None = None) -> None:
    log_path.parent.mkdir(parents=True, exist_ok=True)
    with log_path.open("w", encoding="utf-8") as log:
        log.write("$ " + " ".join(cmd) + "\n")
        log.flush()
        if lock_path is None:
            subprocess.run(cmd, check=True, stdout=log, stderr=subprocess.STDOUT)
            return
        lock_path.parent.mkdir(parents=True, exist_ok=True)
        with lock_path.open("w", encoding="utf-8") as lock:
            fcntl.flock(lock.fileno(), fcntl.LOCK_EX)
            log.write(f"# acquired global I/O lock: {lock_path}\n")
            log.flush()
            subprocess.run(cmd, check=True, stdout=log, stderr=subprocess.STDOUT)


def run_optional_step(cmd: list[str], log_path: Path) -> int:
    log_path.parent.mkdir(parents=True, exist_ok=True)
    with log_path.open("w", encoding="utf-8") as log:
        log.write("$ " + " ".join(cmd) + "\n")
        log.flush()
        completed = subprocess.run(cmd, check=False, stdout=log, stderr=subprocess.STDOUT)
        log.write(f"\nexit_code\t{completed.returncode}\n")
        return completed.returncode


def write_json(path: Path, payload: dict[str, object]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(payload, indent=2), encoding="utf-8")


def metadata_is_current(path: Path, schema_version: str) -> bool:
    if not path.is_file() or path.stat().st_size == 0:
        return False
    try:
        payload = json.loads(path.read_text(encoding="utf-8"))
    except (OSError, json.JSONDecodeError):
        return False
    return payload.get("ok", True) is not False and payload.get("schema_version") == schema_version


def any_source_newer(target: Path, sources: list[Path]) -> bool:
    if not target.is_file():
        return True
    return any(source.is_file() and source.stat().st_mtime > target.stat().st_mtime for source in sources)

def split_csv(value: str) -> list[str]:
    return [item.strip() for item in value.split(",") if item.strip()]


def resolve_tool_command(tool: str) -> tuple[list[str], str]:
    tool_path = Path(tool)
    resolved = tool_path if tool_path.is_file() else None
    if resolved is None and "/" not in tool:
        found = shutil.which(tool)
        resolved = Path(found) if found else None
    if resolved is None:
        return [tool], "unresolved"
    python = resolved.parent / "python"
    try:
        first_line = resolved.open("r", encoding="utf-8", errors="ignore").readline().strip()
    except OSError:
        first_line = ""
    if first_line.startswith("#!") and not Path(first_line[2:]).exists() and python.is_file():
        return [str(python), "-m", "nanoplot.NanoPlot"], "python_module_fallback"
    return [str(resolved)], "direct"


def nanoplot_is_current(status_path: Path, fastq: str, command: list[str], extra_args: list[str]) -> bool:
    if not status_path.is_file():
        return False
    try:
        status = json.loads(status_path.read_text(encoding="utf-8"))
    except (OSError, json.JSONDecodeError):
        return False
    if status.get("fastq") != fastq or status.get("command") != command or status.get("extra_args") != extra_args:
        return False
    if status.get("status") == "skipped":
        return True
    if status.get("status") != "complete":
        return False
    stats_path = Path(str(status.get("nanostats", "")))
    report_path = Path(str(status.get("report", "")))
    return stats_path.is_file() and report_path.is_file()


def run_nanoplot_qc(
    sample_id: str,
    label: str,
    fastq: str,
    qc_dir: Path,
    logs_dir: Path,
    base_command: list[str],
    command_mode: str,
    threads: int,
    extra_args: list[str],
) -> Path:
    outdir = qc_dir / "nanoplot" / label
    status_path = outdir / "nanoplot_status.json"
    if nanoplot_is_current(status_path, fastq, base_command, extra_args):
        return status_path
    stats_path = outdir / "NanoStats.txt"
    if not stats_path.is_file():
        stats_path = outdir / f"{sample_id}.{label}.NanoStats.txt"
    report_path = outdir / "NanoPlot-report.html"
    if not report_path.is_file():
        report_path = outdir / f"{sample_id}.{label}.NanoPlot-report.html"
    outdir.mkdir(parents=True, exist_ok=True)
    version_cmd = base_command + ["--version"]
    try:
        version = subprocess.run(version_cmd, check=False, text=True, capture_output=True, timeout=120)
    except (OSError, subprocess.TimeoutExpired) as exc:
        write_json(
            status_path,
            {
                "status": "skipped",
                "reason": str(exc),
                "sample_id": sample_id,
                "label": label,
                "fastq": fastq,
                "command": base_command,
                "command_mode": command_mode,
                "extra_args": extra_args,
            },
        )
        return status_path
    if version.returncode != 0:
        write_json(
            status_path,
            {
                "status": "skipped",
                "reason": (version.stderr or version.stdout).strip(),
                "sample_id": sample_id,
                "label": label,
                "fastq": fastq,
                "command": base_command,
                "command_mode": command_mode,
                "extra_args": extra_args,
            },
        )
        return status_path
    cmd = (
        base_command
        + [
            "--fastq",
            fastq,
            "--outdir",
            str(outdir),
            "--prefix",
            f"{sample_id}.{label}.",
            "--threads",
            str(threads),
        ]
        + extra_args
    )
    exit_code = run_optional_step(cmd, logs_dir / f"{sample_id}.00_nanoplot_{label}.log")
    write_json(
        status_path,
        {
            "status": "complete" if exit_code == 0 and stats_path.is_file() else "failed",
            "exit_code": exit_code,
            "sample_id": sample_id,
            "label": label,
            "fastq": fastq,
            "command": base_command,
            "command_mode": command_mode,
            "extra_args": extra_args,
            "version": (version.stdout or version.stderr).strip(),
            "outdir": str(outdir),
            "nanostats": str(stats_path),
            "report": str(report_path),
        },
    )
    return status_path


def values_match(observed: object, expected: object) -> bool:
    if expected is None:
        return observed is None or observed == ""
    if isinstance(expected, bool):
        if isinstance(observed, bool):
            return observed == expected
        return str(observed).lower() in {"1", "true", "yes", "on"} if expected else str(observed).lower() in {"0", "false", "no", "off", ""}
    if isinstance(expected, int):
        try:
            return int(observed) == expected
        except (TypeError, ValueError):
            return False
    if isinstance(expected, float):
        try:
            return abs(float(observed) - expected) < 1e-12
        except (TypeError, ValueError):
            return False
    return str(observed) == str(expected)


def slam_is_current(path: Path, expected: dict[str, object]) -> bool:
    if not path.is_file():
        return False
    try:
        qc = json.loads(path.read_text(encoding="utf-8"))
    except (OSError, json.JSONDecodeError):
        return False
    params = qc.get("calling_parameters")
    if not isinstance(params, dict):
        params = qc
    for key, expected_value in expected.items():
        observed = params.get(key, qc.get(key))
        if not values_match(observed, expected_value):
            return False
    return True


def file_nonempty(path: Path) -> bool:
    return path.is_file() and path.stat().st_size > 0


def alignment_is_current(sorted_bam: Path, primary_bam: Path, sample_id: str, min_mapq: int) -> bool:
    outdir = sorted_bam.parent
    required = [
        sorted_bam,
        sorted_bam.with_suffix(sorted_bam.suffix + ".bai"),
        outdir / f"{sample_id}.flagstat.txt",
        primary_bam,
        primary_bam.with_suffix(primary_bam.suffix + ".bai"),
        outdir / f"{sample_id}.primary.flagstat.txt",
    ]
    expected_primary_name = f"{sample_id}.primary.q{min_mapq}.bam"
    return primary_bam.name == expected_primary_name and all(file_nonempty(path) for path in required)


def link_full_reusable_stages(
    sample_dir: Path,
    reuse_results_root: Path | None,
    sample_id: str,
    min_mapq: int,
) -> list[dict[str, str]]:
    """Link validated immutable stages from the prior full run into a new result version."""
    if reuse_results_root is None:
        return []
    source_sample = reuse_results_root / sample_id
    source_trim = source_sample / "02_trim" / "04_tool_output" / "trim"
    source_align = source_sample / "03_alignment" / "04_tool_output" / "minimap2_samtools"
    source_qc = source_sample / "01_raw_qc" / "04_tool_output" / "longread_qc" / "fastq"
    qc_required = [
        source_qc / f"{sample_id}.raw_fastq_qc.json",
        source_qc / f"{sample_id}.insert_fastq_qc.json",
    ]
    trim_required = [
        source_trim / f"{sample_id}.insert.fastq.gz",
        source_trim / f"{sample_id}.trim_summary.json",
    ]
    source_sorted = source_align / f"{sample_id}.sorted.bam"
    source_primary = source_align / f"{sample_id}.primary.q{min_mapq}.bam"
    if not all(file_nonempty(path) for path in trim_required):
        raise RuntimeError(f"Reusable trim stage is incomplete for {sample_id}: {source_trim}")
    if not alignment_is_current(source_sorted, source_primary, sample_id, min_mapq):
        raise RuntimeError(f"Reusable alignment stage is incomplete for {sample_id}: {source_align}")

    stage_names = ("01_raw_qc", "02_trim", "03_alignment")
    records: list[dict[str, str]] = []
    sample_dir.mkdir(parents=True, exist_ok=True)
    for stage_name in stage_names:
        source = source_sample / stage_name
        target = sample_dir / stage_name
        if stage_name == "01_raw_qc" and not all(file_nonempty(path) for path in qc_required):
            # Do not link an incomplete QC stage: the V3 run must create a local copy instead.
            continue
        if not source.is_dir():
            if stage_name == "01_raw_qc":
                continue
            raise RuntimeError(f"Reusable stage is missing for {sample_id}: {source}")
        if target.is_symlink():
            if target.resolve() != source.resolve():
                raise RuntimeError(f"Existing reusable-stage link points elsewhere: {target}")
        elif target.exists():
            # Resume a local V3 stage without replacing or merging it.
            continue
        else:
            target.symlink_to(source, target_is_directory=True)
        if target.is_symlink():
            records.append(
                {
                    "stage": stage_name,
                    "source": str(source.resolve()),
                    "target": str(target),
                    "mode": "read_only_symlink",
                }
            )
    return records
