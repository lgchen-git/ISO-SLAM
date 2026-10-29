#!/usr/bin/env python3
"""Aggregate only technical QC for the independent-sample V4 project report."""

from __future__ import annotations

import argparse
import csv
import json
from pathlib import Path

from iso_slam_v4_common import ensure_stage, write_json, write_tsv


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--results-root", required=True, type=Path)
    parser.add_argument("--samples", required=True, type=Path)
    parser.add_argument("--baseline-results-root", type=Path)
    args = parser.parse_args()
    stage = args.results_root / "_project/08_project_qc"
    ensure_stage(stage)
    with args.samples.open(encoding="utf-8", newline="") as handle:
        samples = list(csv.DictReader(handle, delimiter="\t"))
    rows: list[dict[str, object]] = []
    for sample_row in samples:
        sample = sample_row["sample_id"]
        root = Path(sample_row.get("source_results_root") or args.results_root) / sample
        trim = json.loads((root / f"02_trim/04_tool_output/trim/{sample}.trim_summary.json").read_text(encoding="utf-8"))
        catalog = json.loads((root / "04_isoform_catalog/05_metadata/catalog_validation.json").read_text(encoding="utf-8"))
        grand = json.loads((root / "05_slam_model/05_metadata/grand3_validation.json").read_text(encoding="utf-8"))
        dynamics = json.loads((root / "06_isoform_dynamics/05_metadata/dynamics_validation.json").read_text(encoding="utf-8"))
        calls = dynamics["posterior_call_counts"]
        baseline_no_t_sites = ""
        if args.baseline_results_root:
            baseline = json.loads((args.baseline_results_root / sample / "06_isoform_dynamics/05_metadata/dynamics_validation.json").read_text(encoding="utf-8"))
            baseline_no_t_sites = int(baseline["posterior_call_counts"].get("No_T_sites", 0))
        rows.append(
            {
                "sample_id": sample,
                "group_label_metadata_only": sample_row.get("group", ""),
                "source_results_root": str(root.parent.resolve()),
                "raw_reads": int(trim["total_reads"]),
                "insert_reads": int(trim["trimmed_reads"]),
                "insert_rate": int(trim["trimmed_reads"]) / int(trim["total_reads"]) if int(trim["total_reads"]) else 0,
                "isoform_assigned_reads": dynamics["joined_reads"],
                "dynamics_eligible_reads": dynamics.get("dynamics_eligible_reads", dynamics["joined_reads"]),
                "no_t_sites_reads": int(calls.get("No_T_sites", 0)),
                "assignment_conflict_reads": int(calls.get("Assignment_conflict", 0)),
                "assignment_conflict_rate": int(calls.get("Assignment_conflict", 0)) / int(dynamics["joined_reads"]) if int(dynamics["joined_reads"]) else 0,
                "baseline_no_t_sites_reads": baseline_no_t_sites,
                "no_t_sites_reclassified_reads": (baseline_no_t_sites - int(calls.get("No_T_sites", 0))) if baseline_no_t_sites != "" else "",
                "catalog_isoforms": catalog["catalog_isoforms"],
                "filtered_isoforms": catalog["filtered_isoforms"],
                "grand3_global_binom_ntr": grand["global_parameters"]["binom_ntr"],
                "posterior_mean": dynamics["posterior_mean"],
                "posterior_grand3_absolute_difference": dynamics["posterior_grand3_absolute_difference"],
                "switch_candidates": dynamics["switch_candidates"],
                "sample_report": str((root / f"07_sample_report/{sample}_report.html").resolve()),
            }
        )
    write_tsv(stage / "01_tables/sample_qc_summary.tsv", rows)
    comparison_rows = [
        {
            "sample_id": row["sample_id"],
            "baseline_no_t_sites_reads": row["baseline_no_t_sites_reads"],
            "revised_no_t_sites_reads": row["no_t_sites_reads"],
            "assignment_conflict_reads": row["assignment_conflict_reads"],
            "no_t_sites_reclassified_reads": row["no_t_sites_reclassified_reads"],
        }
        for row in rows
    ]
    write_tsv(stage / "01_tables/no_t_sites_reclassification_summary.tsv", comparison_rows)
    write_json(stage / "05_metadata/project_qc_validation.json", {"schema_version": "iso-slam-project-qc-v6", "status": "PASS", "sample_count": len(rows), "scope": "technical_qc_only", "group_inference": False, "baseline_results_root": str(args.baseline_results_root) if args.baseline_results_root else "", "assignment_conflict_reads": sum(int(row["assignment_conflict_reads"]) for row in rows)})
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
