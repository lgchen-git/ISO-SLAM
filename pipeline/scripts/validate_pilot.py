#!/usr/bin/env python3
"""Apply fixed data, strand, report-input, and method readiness pilot gates."""

from __future__ import annotations

import argparse
import csv
import json
from pathlib import Path


def read_config(path: str) -> dict[str, str]:
    if not path:
        return {}
    config: dict[str, str] = {}
    with Path(path).open(encoding="utf-8") as handle:
        for raw in handle:
            if not raw.strip() or raw.lstrip().startswith("#") or raw[:1].isspace() or ":" not in raw:
                continue
            key, value = raw.split(":", 1)
            config[key.strip()] = value.split("#", 1)[0].strip().strip("'\"")
    return config


def as_bool(value: str, default: bool = False) -> bool:
    if not value:
        return default
    return value.lower() in {"1", "true", "yes", "on"}


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--results-root", required=True)
    parser.add_argument("--samples", required=True)
    parser.add_argument("--config", default="")
    args = parser.parse_args()
    root = Path(args.results_root)
    config = read_config(args.config)
    require_grand_slam = as_bool(config.get("grand_slam_enabled", ""), False)
    grand_pilot_ok = True
    if require_grand_slam:
        grand_pilot_path = root / "_project/06_cohort_qc/05_metadata/grand3_pilot_validation.json"
        grand_pilot_ok = grand_pilot_path.is_file() and json.loads(
            grand_pilot_path.read_text(encoding="utf-8")
        ).get("status") == "PASS"
    with Path(args.samples).open(encoding="utf-8", newline="") as handle:
        samples = list(csv.DictReader(handle, delimiter="\t"))

    rows = []
    passed = True
    observed_modes: set[str] = set()
    for meta in samples:
        sample_id = meta["sample_id"]
        sample_root = root / sample_id
        trim = json.loads(
            (sample_root / "02_trim/04_tool_output/trim" / f"{sample_id}.trim_summary.json").read_text()
        )
        align = sample_root / "03_alignment/04_tool_output/minimap2_samtools" / f"{sample_id}.primary.flagstat.txt"
        primary = int(align.read_text().splitlines()[0].split()[0])
        caller = sample_root / "04_slam_calling/04_tool_output/caller"
        qc = json.loads((caller / f"{sample_id}.bulk_slam_qc.json").read_text())
        strand = json.loads(
            (
                sample_root
                / "04_slam_calling/05_metadata/strand_orientation"
                / f"{sample_id}.strand_orientation_audit.json"
            ).read_text()
        )
        grand = json.loads((caller / f"{sample_id}.grand_slam_compatibility.json").read_text())
        raw = int(trim["total_reads"])
        trimmed = int(trim["trimmed_reads"])
        assigned = int(qc.get("assigned_reads", 0))
        trim_rate = trimmed / raw if raw else 0.0
        map_rate = primary / trimmed if trimmed else 0.0
        assign_rate = assigned / primary if primary else 0.0
        strand_ok = strand.get("status") == "PASS"
        selected_mode = str(strand.get("selected_mode", ""))
        if strand_ok:
            observed_modes.add(selected_mode)
        grand_ok = grand.get("status") == "ready" if require_grand_slam else True
        exact_inputs_ok = all(
            (caller / f"{sample_id}.{suffix}").is_file()
            for suffix in ("assignment_summary.tsv", "substitution_summary.tsv", "ntr_sufficient_stats.tsv")
        )
        ok = (
            trim_rate >= 0.85
            and map_rate >= 0.70
            and assign_rate >= 0.80
            and int(qc.get("t_covered", 0)) > 0
            and strand_ok
            and exact_inputs_ok
            and grand_ok
        )
        passed &= ok
        rows.append(
            [
                sample_id,
                raw,
                trimmed,
                trim_rate,
                primary,
                map_rate,
                assigned,
                assign_rate,
                qc.get("t_covered", 0),
                selected_mode,
                strand.get("selected_same_strand_fraction", 0),
                strand.get("predicted_antisense_fraction", 0),
                grand.get("status", "missing"),
                "PASS" if exact_inputs_ok else "FAIL",
                "PASS" if ok else "FAIL",
            ]
        )
    if len(observed_modes) != 1:
        passed = False
    passed &= grand_pilot_ok

    out = root / "_project/06_cohort_qc/01_tables"
    out.mkdir(parents=True, exist_ok=True)
    with (out / "pilot_acceptance.tsv").open("w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(
            [
                "sample_id",
                "raw_reads",
                "trimmed_reads",
                "trim_rate",
                "primary_reads",
                "mapping_rate",
                "assigned_reads",
                "assignment_rate",
                "t_covered",
                "strand_mode",
                "strand_consistency",
                "antisense_fraction",
                "grand_slam_status",
                "exact_report_inputs",
                "status",
            ]
        )
        writer.writerows(rows)
    (out / "PILOT_PASS").write_text("PASS\n" if passed else "FAIL\n", encoding="utf-8")
    raise SystemExit(0 if passed else 1)


if __name__ == "__main__":
    main()
