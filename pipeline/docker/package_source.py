#!/usr/bin/env python3
"""Create a deterministic, data-free ISO-SLAM source snapshot."""

from __future__ import annotations

import csv
import gzip
import hashlib
import os
import shutil
import tarfile
import tempfile
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]
ARCHIVE = ROOT / "outputs/docker/iso-slam_20260930_v1.0.2_source.tar.gz"
EXCLUDED_ROOTS = {
    "inputs", "reference", "log", "outputs", "tmp", ".envs", ".mamba", ".codex", ".aws", ".agents",
    "results", "raw_data", "data", ".snakemake", "__pycache__", ".pytest_cache", ".ruff_cache", ".mypy_cache",
}
EXCLUDED_NAMES = {"__pycache__", ".pytest_cache", ".ruff_cache", ".mypy_cache", ".snakemake"}
EXCLUDED_SUFFIXES = {".pyc", ".pyo", ".bam", ".bai", ".fastq", ".fq"}


def ignored(directory: Path, names: list[str]) -> set[str]:
    directory = Path(directory)
    if directory == ROOT:
        return set(names) & EXCLUDED_ROOTS
    result = set(names) & EXCLUDED_NAMES
    if directory.name == ".git":
        result |= {"hooks", "logs"} & set(names)
    return result | {
        name for name in names
        if Path(name).suffix.lower() in EXCLUDED_SUFFIXES or name.lower().endswith((".fastq.gz", ".fq.gz"))
    }


def digest(path: Path) -> str:
    value = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            value.update(block)
    return value.hexdigest()


def write_manifest(stage: Path) -> None:
    manifest = stage / "docs/source_snapshot_manifest.tsv"
    rows = []
    for path in sorted(stage.rglob("*")):
        if not path.is_file() or ".git" in path.relative_to(stage).parts or path == manifest:
            continue
        rows.append((path.relative_to(stage).as_posix(), digest(path)))
    manifest.parent.mkdir(parents=True, exist_ok=True)
    with manifest.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(("path", "sha256"))
        writer.writerows(rows)


def package(stage: Path) -> int:
    with ARCHIVE.open("wb") as raw:
        with gzip.GzipFile(fileobj=raw, mode="wb", mtime=0, filename="") as compressed:
            with tarfile.open(fileobj=compressed, mode="w", format=tarfile.PAX_FORMAT) as archive:
                for path in [stage, *sorted(stage.rglob("*"))]:
                    relative = Path("ISO-SLAM") if path == stage else Path("ISO-SLAM") / path.relative_to(stage)
                    info = archive.gettarinfo(str(path), arcname=relative.as_posix())
                    info.uid = info.gid = 0
                    info.uname = info.gname = ""
                    info.mtime = 0
                    if info.issym() or info.islnk():
                        raise RuntimeError(f"Source snapshot cannot contain symbolic links: {path}")
                    if info.isfile():
                        with path.open("rb") as handle:
                            archive.addfile(info, handle)
                    else:
                        archive.addfile(info)
    return len([p for p in stage.rglob("*") if p.is_file()])


def main() -> int:
    if ARCHIVE.exists():
        raise SystemExit(f"Refusing to overwrite existing snapshot: {ARCHIVE}")
    ARCHIVE.parent.mkdir(parents=True, exist_ok=True)
    with tempfile.TemporaryDirectory(prefix="source_package_", dir=ROOT / "tmp") as temporary:
        stage = Path(temporary) / "ISO-SLAM"
        shutil.copytree(ROOT, stage, ignore=ignored, symlinks=True)
        git_config = stage / ".git/config"
        if git_config.is_file():
            git_config.write_text(
                "[core]\n\trepositoryformatversion = 0\n\tfilemode = true\n\tbare = false\n\tlogallrefupdates = true\n",
                encoding="utf-8",
            )
        for name in ("config", "inputs", "log", "outputs", "pipeline", "reference", "tmp"):
            (stage / name).mkdir(parents=True, exist_ok=True)
        for name, content in {
            "inputs/README.md": "Place the FASTQ files at the paths recorded in the configured sample sheet.\n",
            "reference/README.md": "Mount the matching GRCh38-2024-A reference files here; do not modify them during a run.\n",
        }.items():
            (stage / name).write_text(content, encoding="utf-8")
        write_manifest(stage)
        shutil.copy2(stage / "docs/source_snapshot_manifest.tsv", ROOT / "docs/source_snapshot_manifest.tsv")
        count = package(stage)
    print(f"Source snapshot created: {ARCHIVE} ({count} files; excludes inputs, reference data, logs, outputs, and tmp)")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
