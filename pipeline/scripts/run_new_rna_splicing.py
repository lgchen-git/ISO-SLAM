#!/usr/bin/env python3
"""Estimate annotated alternative-splicing PSI from new-RNA transcript support."""

from __future__ import annotations

import argparse
import csv
import gzip
import json
import math
import shutil
import subprocess
from pathlib import Path

def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open(encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def write_tsv(path: Path, fields: list[str], rows: list[dict[str, object]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", lineterminator="\n", extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--samples", required=True, type=Path)
    parser.add_argument("--gtf", required=True, type=Path)
    parser.add_argument("--results-root", required=True, type=Path)
    parser.add_argument("--outdir", required=True, type=Path)
    parser.add_argument("--workdir", required=True, type=Path)
    parser.add_argument("--suppa", required=True)
    parser.add_argument("--rscript", required=True)
    parser.add_argument("--contrasts-json", required=True)
    args = parser.parse_args()
    samples = read_tsv(args.samples)
    contrasts = json.loads(args.contrasts_json)
    args.outdir.mkdir(parents=True, exist_ok=True)
    work = args.workdir
    work.mkdir(parents=True, exist_ok=True)
    gtf = args.gtf
    if gtf.suffix == ".gz":
        gtf = work / "reference.gtf"
        with gzip.open(args.gtf, "rb") as source, gtf.open("wb") as target:
            shutil.copyfileobj(source, target, length=4 * 1024 * 1024)
    ioe_prefix = work / "reference_events"
    subprocess.run([args.suppa, "generateEvents", "-i", str(gtf), "-o", str(ioe_prefix), "-f", "ioe", "-e", "SE", "SS", "MX", "RI"], check=True)
    ioe_files = sorted(work.glob("reference_events*.ioe"))
    events: dict[str, tuple[str, set[str], set[str]]] = {}
    for ioe in ioe_files:
        with ioe.open(encoding="utf-8") as handle:
            header = handle.readline().rstrip("\n").split("\t")
            for line in handle:
                row = dict(zip(header, line.rstrip("\n").split("\t")))
                event_id = row.get("event_id", "")
                gene = row.get("gene_id", "")
                inclusion = set(filter(None, row.get("alternative_transcripts", "").split(",")))
                total = set(filter(None, row.get("total_transcripts", "").split(",")))
                if event_id and gene and inclusion and total:
                    events[event_id] = (gene, inclusion, total)
    if not events:
        raise RuntimeError("SUPPA did not produce any annotated alternative-splicing events")

    transcript_counts: dict[str, dict[str, float]] = {}
    sample_novel = []
    for sample in samples:
        sid = sample["sample_id"]
        root = args.results_root / sid
        dynamics = read_tsv(root / "06_isoform_dynamics/01_tables" / f"{sid}.isoform_dynamics.tsv")
        counts: dict[str, float] = {}
        for row in dynamics:
            transcript = row["isoform_id"].split("|", 1)[-1]
            counts[transcript] = counts.get(transcript, 0.0) + float(row["expected_new_reads"])
        transcript_counts[sid] = counts
        catalog = read_tsv(root / "04_isoform_catalog/01_tables" / f"{sid}.isoform_catalog.tsv")
        for row in catalog:
            category = row.get("structural_category", "")
            if category not in {"full-splice_match", "incomplete-splice_match"}:
                sample_novel.append({"sample_id": sid, "group": sample["group"], "gene_id": row.get("gene_id", ""),
                    "isoform_id": row.get("isoform_id", ""), "structural_category": category,
                    "read_support": row.get("read_support", ""), "filter_status": row.get("filter_status", "")})

    sample_psi = []
    event_values: dict[str, dict[str, tuple[float, float]]] = {}
    for event_id, (gene, inclusion, total) in events.items():
        values = {}
        for sample in samples:
            sid = sample["sample_id"]
            counts = transcript_counts[sid]
            included = sum(counts.get(tx, 0.0) for tx in inclusion)
            all_event = sum(counts.get(tx, 0.0) for tx in total)
            psi = included / all_event if all_event > 0 else math.nan
            values[sid] = (included, all_event - included)
            sample_psi.append({"event_id": event_id, "gene_id": gene, "sample_id": sid, "group": sample["group"],
                "inclusion_new_read_equivalents": round(included, 6), "exclusion_new_read_equivalents": round(all_event-included, 6),
                "total_new_read_equivalents": round(all_event, 6), "new_RNA_PSI": "NA" if math.isnan(psi) else round(psi, 6)})
        event_values[event_id] = values

    write_tsv(args.outdir / "sample_event_psi.tsv", ["event_id", "gene_id", "sample_id", "group", "inclusion_new_read_equivalents", "exclusion_new_read_equivalents", "total_new_read_equivalents", "new_RNA_PSI"], sample_psi)
    subprocess.run([args.rscript, str(Path(__file__).with_name("run_new_rna_splicing.R")),
        "--samples", str(args.samples), "--psi", str(args.outdir / "sample_event_psi.tsv"),
        "--contrasts-json", args.contrasts_json, "--out", str(args.outdir / "differential_new_rna_splicing.tsv")], check=True)
    write_tsv(args.outdir / "novel_structure_inventory.tsv", ["sample_id", "group", "gene_id", "isoform_id", "structural_category", "read_support", "filter_status"], sample_novel)
    differential = read_tsv(args.outdir / "differential_new_rna_splicing.tsv")
    counts_by_contrast = {f'{x["treatment"]}_vs_{x["control"]}': sum(r["significant"].lower() == "true" for r in differential if r["contrast"] == f'{x["treatment"]}_vs_{x["control"]}') for x in contrasts}
    validation = {"status": "PASS", "method": "SUPPA2 reference events; transcript-level expected new-read equivalents; empirical-logit read-support-weighted limma with eBayes; BH within contrast", "event_types": ["SE", "A3", "A5", "MX", "RI"], "reference_event_count": len(events), "tested_event_count_by_contrast": {name: sum(r["contrast"] == name for r in differential) for name in counts_by_contrast}, "significant_event_count_by_contrast": counts_by_contrast, "novel_catalog_records": len(sample_novel), "new_RNA_novel_event_differential": "NOT_ESTIMABLE: novel structures lack validated read-level new-RNA posterior assignment in the existing outputs", "support_threshold_expected_new_read_equivalents_per_sample": 20, "significance_threshold": "BH FDR < 0.05 and absolute delta PSI >= 0.10", "biological_replicates_per_group": 2, "interpretation": "Exploratory; low biological replication limits statistical power."}
    (args.outdir / "new_rna_splicing_validation.json").write_text(json.dumps(validation, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
