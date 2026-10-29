#!/usr/bin/env python3
"""Archive the exact duplicate pre-v2 sample report after final validation."""

from __future__ import annotations

import argparse
import csv
import hashlib
from pathlib import Path


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--sample-root", required=True, type=Path)
    parser.add_argument("--sample-id", required=True)
    parser.add_argument("--archive-root", required=True, type=Path)
    args = parser.parse_args()

    source = args.sample_root / "05_sample_report" / "04_tool_output" / "delivery" / f"{args.sample_id}.bulk_slam_report.html"
    if not source.is_file():
        print(f"legacy report absent: {source}")
        return 0
    digest = sha256(source)
    sample_archive = args.archive_root / args.sample_id
    sample_archive.mkdir(parents=True, exist_ok=True)
    target = sample_archive / source.name
    if target.exists() and sha256(target) != digest:
        target = sample_archive / f"{source.stem}.{digest[:12]}{source.suffix}"
    if target.exists():
        source.unlink()
        status = "removed_duplicate_after_existing_archive_match"
    else:
        source.replace(target)
        status = "archived"

    manifest = args.archive_root / "archive_manifest.tsv"
    rows: list[dict[str, str]] = []
    if manifest.is_file():
        with manifest.open(encoding="utf-8", newline="") as handle:
            rows = list(csv.DictReader(handle, delimiter="\t"))
    row = {
        "sample_id": args.sample_id,
        "source": str(source.resolve()),
        "archive": str(target.resolve()),
        "sha256": digest,
        "status": status,
    }
    rows = [item for item in rows if item.get("sample_id") != args.sample_id]
    rows.append(row)
    with manifest.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=["sample_id", "source", "archive", "sha256", "status"],
            delimiter="\t",
            lineterminator="\n",
        )
        writer.writeheader()
        writer.writerows(sorted(rows, key=lambda item: item["sample_id"]))
    print(f"{status}: {source} -> {target}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
