#!/usr/bin/env python3
"""Build a comparable LR-bulkSLAM/LR-scSLAM summary table."""

from __future__ import annotations

import argparse
import csv
from pathlib import Path


FIELDS = [
    "sample_id",
    "analysis_type",
    "unit_type",
    "total_units",
    "assigned_units",
    "new_units",
    "old_units",
    "ambiguous_units",
    "new_fraction",
    "new_fraction_excluding_ambiguous",
    "total_t_covered",
    "total_tc_count",
    "global_tc_rate",
    "bg_tc_rate",
    "bg_rate_source",
    "assignment_mode",
]


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Summarize LR-SLAM outputs for cross-workflow comparison.")
    parser.add_argument("--results-dir", required=True)
    parser.add_argument("--analysis-type", choices=["lr_bulk_slam", "lr_scslam_bd"], default="")
    parser.add_argument("--out", default="")
    return parser.parse_args()


def read_first_row(path: Path) -> dict[str, str]:
    with path.open("r", encoding="utf-8", newline="") as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    return rows[0] if rows else {}


def find_summaries(results_dir: Path, analysis_type: str) -> list[tuple[str, Path]]:
    sample_root = results_dir / "samples"
    if not sample_root.is_dir():
        return []
    out: list[tuple[str, Path]] = []
    for sample_dir in sorted(path for path in sample_root.iterdir() if path.is_dir()):
        sample_id = sample_dir.name
        if analysis_type in {"", "lr_bulk_slam"}:
            path = sample_dir / "03_bulk_slam" / f"{sample_id}.bulk_slam_summary.tsv"
            if path.is_file():
                out.append(("lr_bulk_slam", path))
        if analysis_type in {"", "lr_scslam_bd"}:
            path = sample_dir / "03_slam_molecules" / f"{sample_id}.scslam_summary.tsv"
            if path.is_file():
                out.append(("lr_scslam_bd", path))
    return out


def comparable_row(analysis_type: str, row: dict[str, str]) -> dict[str, str]:
    sample_id = row.get("sample", "")
    if analysis_type == "lr_bulk_slam":
        unit_type = "read"
        total_units = row.get("total_reads", "0")
        assigned_units = row.get("assigned_reads", "0")
        new_units = row.get("new_read_count", "0")
        old_units = row.get("old_read_count", "0")
        ambiguous_units = row.get("ambiguous_read_count", "0")
    else:
        unit_type = "umi_molecule"
        total_units = row.get("molecules_total", "0")
        assigned_units = row.get("molecules_total", "0")
        new_units = row.get("new_molecule_count", "0")
        old_units = row.get("old_molecule_count", "0")
        ambiguous_units = row.get("ambiguous_molecule_count", "0")
    return {
        "sample_id": sample_id,
        "analysis_type": analysis_type,
        "unit_type": unit_type,
        "total_units": total_units,
        "assigned_units": assigned_units,
        "new_units": new_units,
        "old_units": old_units,
        "ambiguous_units": ambiguous_units,
        "new_fraction": row.get("new_fraction", "0"),
        "new_fraction_excluding_ambiguous": row.get("new_fraction_excluding_ambiguous", "0"),
        "total_t_covered": row.get("total_t_covered", "0"),
        "total_tc_count": row.get("total_tc_count", "0"),
        "global_tc_rate": row.get("global_tc_rate", "0"),
        "bg_tc_rate": row.get("bg_tc_rate", ""),
        "bg_rate_source": row.get("bg_rate_source", ""),
        "assignment_mode": row.get("assignment_mode", ""),
    }


def main() -> None:
    args = parse_args()
    results_dir = Path(args.results_dir)
    out_path = Path(args.out) if args.out else results_dir / "tables" / "lr_slam_comparable_summary.tsv"
    out_path.parent.mkdir(parents=True, exist_ok=True)

    rows = []
    for analysis_type, path in find_summaries(results_dir, args.analysis_type):
        summary = read_first_row(path)
        if summary:
            rows.append(comparable_row(analysis_type, summary))

    with out_path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, delimiter="\t", lineterminator="\n", fieldnames=FIELDS)
        writer.writeheader()
        writer.writerows(rows)


if __name__ == "__main__":
    main()
