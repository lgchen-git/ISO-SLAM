#!/usr/bin/env python3
"""Assemble cohort matrices and sample QC from completed ISO-SLAM samples."""

from __future__ import annotations

import argparse
import csv
import json
from pathlib import Path


def args_parse():
    p = argparse.ArgumentParser()
    p.add_argument("--results-root", required=True)
    p.add_argument("--samples", required=True)
    return p.parse_args()


def read_tsv(path: Path):
    with path.open(encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def write_matrix(path: Path, genes: list[str], names: dict[str, str], samples: list[str], values):
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(["gene_id", "gene_name", *samples])
        for gene in genes:
            writer.writerow([gene, names.get(gene, gene), *[values[s].get(gene, 0) for s in samples]])


def flagstat_total(path: Path) -> int:
    if not path.is_file():
        return 0
    line = path.read_text(encoding="utf-8", errors="ignore").splitlines()[0]
    try:
        return int(line.split(" ", 1)[0])
    except (IndexError, ValueError):
        return 0


def main():
    args = args_parse()
    root = Path(args.results_root)
    sheet = read_tsv(Path(args.samples))
    sample_ids = [row["sample_id"] for row in sheet]
    groups = {row["sample_id"]: row["group"] for row in sheet}
    out_stage = root / "_project" / "06_cohort_qc"
    for child in ("01_tables", "02_figures", "03_logs", "04_tool_output", "05_metadata"):
        (out_stage / child).mkdir(parents=True, exist_ok=True)
    tables = out_stage / "01_tables"
    total = {s: {} for s in sample_ids}
    new = {s: {} for s in sample_ids}
    frac = {s: {} for s in sample_ids}
    names = {}
    metric_rows = []
    for sample in sample_ids:
        sample_root = root / sample
        trim_path = sample_root / "02_trim" / "04_tool_output" / "trim" / f"{sample}.trim_summary.json"
        align_dir = sample_root / "03_alignment" / "04_tool_output" / "minimap2_samtools"
        slam_dir = sample_root / "04_slam_calling" / "04_tool_output" / "caller"
        gene_path = slam_dir / f"{sample}.bulk_gene_summary.tsv"
        qc_path = slam_dir / f"{sample}.bulk_slam_qc.json"
        report = sample_root / "05_sample_report" / f"{sample}.report.html"
        required = [trim_path, gene_path, qc_path, report]
        missing = [str(p) for p in required if not p.is_file() or p.stat().st_size == 0]
        if missing:
            raise SystemExit(f"{sample}: missing completed outputs: {', '.join(missing)}")
        trim = json.loads(trim_path.read_text(encoding="utf-8"))
        qc = json.loads(qc_path.read_text(encoding="utf-8"))
        for row in read_tsv(gene_path):
            gene = row["gene_id"]
            names[gene] = row.get("gene_name") or gene
            total[sample][gene] = int(row["total_reads"])
            new[sample][gene] = int(row["new_reads"])
            frac[sample][gene] = float(row.get("new_fraction_excluding_ambiguous", 0) or 0)
        raw = int(trim.get("total_reads", 0))
        trimmed = int(trim.get("trimmed_reads", 0))
        primary = flagstat_total(align_dir / f"{sample}.primary.flagstat.txt")
        assigned = int(qc.get("assigned_reads", 0))
        metric_rows.append({
            "sample_id": sample, "group": groups[sample], "raw_reads": raw,
            "trimmed_reads": trimmed, "trimmed_rate": trimmed / raw if raw else 0,
            "primary_reads": primary, "primary_rate": primary / trimmed if trimmed else 0,
            "assigned_reads": assigned, "assigned_rate": assigned / primary if primary else 0,
            "new_reads": int(qc.get("new_reads", 0)), "old_reads": int(qc.get("old_reads", 0)),
            "ambiguous_reads": int(qc.get("ambiguous_reads", 0)),
            "tc_rate": float(qc.get("t_to_c_rate", 0) or 0), "genes_detected": int(qc.get("genes_detected", 0)),
            "sample_report": str(report),
        })
    genes = sorted(set().union(*(set(total[s]) for s in sample_ids)))
    write_matrix(tables / "total_rna_counts.tsv", genes, names, sample_ids, total)
    write_matrix(tables / "new_rna_counts.tsv", genes, names, sample_ids, new)
    write_matrix(tables / "new_fraction.tsv", genes, names, sample_ids, frac)
    with (tables / "sample_metrics.tsv").open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, delimiter="\t", lineterminator="\n", fieldnames=list(metric_rows[0]))
        writer.writeheader(); writer.writerows(metric_rows)


if __name__ == "__main__":
    main()
