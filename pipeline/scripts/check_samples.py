#!/usr/bin/env python3
"""Validate a simple TSV sample sheet for the template."""

from __future__ import annotations

import argparse
import csv
from pathlib import Path
import sys


REQUIRED_COLUMNS = ("sample_id", "group", "input_1", "input_2")
LR_SCSLAM_BD_COLUMNS = (
    "fastq",
    "species",
    "chemistry",
    "genome_fa",
    "gtf",
    "bd_cell_label_py",
)
BULK_SLAM_COLUMNS = (
    "fastq",
    "species",
    "genome_fa",
    "gtf",
)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Validate sample sheet TSV.")
    parser.add_argument("--samples", required=True, help="Path to TSV sample sheet")
    parser.add_argument(
        "--analysis-type",
        default="generic",
        help="Optional analysis type for extended sample contracts",
    )
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    sample_path = Path(args.samples)

    if not sample_path.is_file():
        print(f"Sample sheet not found: {sample_path}", file=sys.stderr)
        return 1

    with sample_path.open("r", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        header = reader.fieldnames or []

        required = list(REQUIRED_COLUMNS)
        if args.analysis_type == "lr_scslam_bd":
            required.extend(LR_SCSLAM_BD_COLUMNS)
        if args.analysis_type == "lr_bulk_slam":
            required.extend(BULK_SLAM_COLUMNS)
        missing_columns = [name for name in required if name not in header]
        if missing_columns:
            print(
                "Missing required columns: " + ", ".join(missing_columns),
                file=sys.stderr,
            )
            return 1

        seen_ids: set[str] = set()
        row_count = 0

        for row_number, row in enumerate(reader, start=2):
            row_count += 1
            sample_id = (row.get("sample_id") or "").strip()
            group = (row.get("group") or "").strip()
            input_1 = (row.get("input_1") or "").strip()
            input_2 = (row.get("input_2") or "").strip()

            if not sample_id:
                print(f"Row {row_number}: sample_id is empty", file=sys.stderr)
                return 1
            if sample_id in seen_ids:
                print(f"Row {row_number}: duplicate sample_id {sample_id}", file=sys.stderr)
                return 1
            seen_ids.add(sample_id)

            if not group:
                print(f"Row {row_number}: group is empty", file=sys.stderr)
                return 1
            if not input_1:
                print(f"Row {row_number}: input_1 is empty", file=sys.stderr)
                return 1
            if not input_2:
                print(f"Row {row_number}: input_2 is empty", file=sys.stderr)
                return 1

            for input_name in (input_1, input_2):
                if not Path(input_name).exists():
                    print(
                        f"Row {row_number}: input path does not exist: {input_name}",
                        file=sys.stderr,
                    )
                    return 1

            if args.analysis_type == "lr_scslam_bd":
                for column in LR_SCSLAM_BD_COLUMNS:
                    value = (row.get(column) or "").strip()
                    if not value:
                        print(f"Row {row_number}: {column} is empty", file=sys.stderr)
                        return 1
                for column in ("fastq", "genome_fa", "gtf", "bd_cell_label_py"):
                    value = (row.get(column) or "").strip()
                    if not Path(value).exists():
                        print(
                            f"Row {row_number}: {column} path does not exist: {value}",
                            file=sys.stderr,
                        )
                        return 1
            if args.analysis_type == "lr_bulk_slam":
                for column in BULK_SLAM_COLUMNS:
                    value = (row.get(column) or "").strip()
                    if not value:
                        print(f"Row {row_number}: {column} is empty", file=sys.stderr)
                        return 1
                for column in ("fastq", "genome_fa", "gtf"):
                    value = (row.get(column) or "").strip()
                    if not Path(value).exists():
                        print(
                            f"Row {row_number}: {column} path does not exist: {value}",
                            file=sys.stderr,
                        )
                        return 1

    if row_count == 0:
        print("Sample sheet has no data rows", file=sys.stderr)
        return 1

    print(f"Sample sheet OK: {sample_path} ({row_count} samples)")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
