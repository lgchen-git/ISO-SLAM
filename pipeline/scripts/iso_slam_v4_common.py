#!/usr/bin/env python3
"""Small shared helpers for the independent-sample ISO-SLAM V4 workflow."""

from __future__ import annotations

import csv
import json
import shutil
import subprocess
from pathlib import Path


SUBDIRS = ("01_tables", "02_figures", "04_tool_output", "05_metadata")


def read_config(path: str | Path) -> dict[str, str]:
    values: dict[str, str] = {}
    with Path(path).open(encoding="utf-8") as handle:
        for raw in handle:
            if not raw.strip() or raw.lstrip().startswith("#") or raw[:1].isspace() or ":" not in raw:
                continue
            key, value = raw.split(":", 1)
            values[key.strip()] = value.split("#", 1)[0].strip().strip("'\"")
    return values


def as_bool(config: dict[str, str], key: str, default: bool = False) -> bool:
    value = config.get(key, "")
    return default if not value else value.lower() in {"1", "true", "yes", "on"}


def ensure_stage(path: Path) -> None:
    for name in SUBDIRS:
        (path / name).mkdir(parents=True, exist_ok=True)


def nonempty(path: Path) -> bool:
    return path.is_file() and path.stat().st_size > 0


def executable(value: str) -> str:
    path = Path(value)
    if path.is_file():
        return str(path.resolve())
    found = shutil.which(value)
    if not found:
        raise RuntimeError(f"Required executable is unavailable: {value}")
    return found


def run_logged(command: list[str], log_path: Path, cwd: Path | None = None) -> None:
    import os
    if os.environ.get("ISO_LOG_ROOT") and "03_logs" in log_path.parts:
        stage = log_path.parent.parent
        log_path = Path(os.environ["ISO_LOG_ROOT"]) / stage.parent.name / stage.name / log_path.name
    log_path.parent.mkdir(parents=True, exist_ok=True)
    with log_path.open("w", encoding="utf-8") as log:
        log.write("$ " + " ".join(command) + "\n")
        log.flush()
        subprocess.run(command, check=True, stdout=log, stderr=subprocess.STDOUT, cwd=cwd)


def write_json(path: Path, payload: dict[str, object]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(payload, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")


def write_tsv(path: Path, rows: list[dict[str, object]], fields: list[str] | None = None) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    if fields is None:
        fields = list(rows[0]) if rows else []
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def find_one(root: Path, patterns: tuple[str, ...], label: str) -> Path:
    matches: list[Path] = []
    for pattern in patterns:
        matches.extend(path for path in root.rglob(pattern) if nonempty(path))
    unique = sorted(set(matches))
    if len(unique) != 1:
        raise RuntimeError(f"Expected one {label} under {root}, found {len(unique)}: {unique[:5]}")
    return unique[0]
