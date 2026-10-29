#!/usr/bin/env python3
"""Run IsoQuant/SQANTI3 and normalize one sample-specific isoform catalog."""

from __future__ import annotations

import argparse
import csv
import fcntl
import gzip
import os
import re
import shutil
from collections import Counter, defaultdict
from pathlib import Path

from iso_slam_v4_common import as_bool, ensure_stage, executable, find_one, nonempty, read_config, run_logged, write_json, write_tsv


CATALOG_SCHEMA = "iso-slam-isoform-catalog-v5"


def arguments() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--config", required=True)
    parser.add_argument("--sample-id", required=True)
    parser.add_argument("--bam", required=True, type=Path)
    parser.add_argument("--fastq", required=True, type=Path)
    parser.add_argument("--fasta", required=True, type=Path)
    parser.add_argument("--gtf", required=True, type=Path)
    parser.add_argument("--outdir", required=True, type=Path)
    parser.add_argument("--threads", type=int, default=8)
    parser.add_argument("--tmp-dir", type=Path)
    parser.add_argument("--skip-tools", action="store_true")
    parser.add_argument("--isoquant-dir", type=Path)
    parser.add_argument("--sqanti-classification", type=Path)
    return parser.parse_args()


def open_text(path: Path):
    return gzip.open(path, "rt", encoding="utf-8", newline="") if path.suffix == ".gz" else path.open(encoding="utf-8", newline="")


def iter_rows(path: Path):
    with open_text(path) as handle:
        lines = (line for line in handle if line.strip() and not line.startswith("#"))
        yield from csv.DictReader(lines, delimiter="\t")


def rows(path: Path) -> list[dict[str, str]]:
    return list(iter_rows(path))


def pick(row: dict[str, str], aliases: tuple[str, ...], default: str = "") -> str:
    normalized = {re.sub(r"[^a-z0-9]", "", key.lower()): value for key, value in row.items()}
    for alias in aliases:
        key = re.sub(r"[^a-z0-9]", "", alias.lower())
        if key in normalized and normalized[key] not in {None, ""}:
            return str(normalized[key])
    return default


def parse_assignments(path: Path, sample_id: str) -> tuple[list[dict[str, object]], Counter, dict[str, set[str]]]:
    read_states: dict[str, dict[str, object] | None | object] = {}
    ambiguous = object()
    counts: Counter = Counter()
    assignment_type_rows: Counter = Counter()
    isoform_reads: dict[str, set[str]] = defaultdict(set)
    for row in iter_rows(path):
        read_id = pick(row, ("read_id", "read", "readname", "read_name"))
        transcript_id = pick(row, ("transcript_id", "isoform_id", "isoform", "transcript"), ".")
        gene_id = pick(row, ("gene_id", "gene", "geneid"), ".")
        classification = pick(row, ("assignment_type", "classification", "assignment", "status"), "unknown")
        polya_value = pick(row, ("polya_position", "polyA_position", "polya", "polyA", "polya_info"), "")
        if not read_id:
            continue
        counts["assignment_rows"] += 1
        if read_id not in read_states:
            read_states[read_id] = None
        if classification.strip().lower() not in {"unique", "unique_minor_difference"}:
            continue
        if transcript_id in {"", ".", "NA", "None", "unassigned"}:
            continue
        candidate = {
            "sample_id": sample_id,
            "read_id": read_id,
            "isoform_id": transcript_id,
            "gene_id": gene_id,
            "isoquant_assignment": classification,
            "isoform_assigned": 1,
            "polya_value": polya_value,
        }
        previous = read_states[read_id]
        if previous is None:
            read_states[read_id] = candidate
        elif previous is not ambiguous and (
            previous["isoform_id"] != transcript_id or previous["gene_id"] != gene_id
        ):
            read_states[read_id] = ambiguous

    normalized: list[dict[str, object]] = []
    for read_id, state in read_states.items():
        if state is ambiguous:
            counts["multiple_confident_assignments"] += 1
            continue
        if state is None:
            continue
        normalized.append(state)
        isoform_reads[str(state["isoform_id"])].add(read_id)
        if state["polya_value"] not in {"", ".", "NA", "None", "False", "false", "0"}:
            counts["polya_evaluable_reads"] += 1
    counts["total_reads"] = len(read_states)
    counts["assigned_reads"] = len(normalized)
    counts["unassigned_reads"] = counts["total_reads"] - counts["assigned_reads"]
    counts["duplicate_assignment_rows"] = counts["assignment_rows"] - counts["total_reads"]
    counts["assignment_type_rows"] = dict(assignment_type_rows)
    if not normalized:
        raise RuntimeError(f"No readable assignment records in {path}")
    return normalized, counts, isoform_reads


def support_from_counts(path: Path | None) -> dict[str, int]:
    support: dict[str, int] = {}
    if path is None:
        return support
    for row in rows(path):
        isoform = pick(row, ("transcript_id", "isoform_id", "transcript", "feature_id"))
        raw_count = pick(row, ("count", "read_count", "reads", "estimated_count"))
        if not isoform or not raw_count:
            continue
        try:
            support[isoform] = int(round(float(raw_count)))
        except ValueError:
            continue
    return support


def truth(value: str) -> bool:
    return value.strip().lower() in {"1", "true", "yes", "y", "t"}


def filter_reason(row: dict[str, str], support: int, config: dict[str, str]) -> str:
    category = pick(row, ("structural_category", "structuralcategory", "category"), "unknown").lower()
    exons_raw = pick(row, ("exons", "num_exons", "exon_count"), "0")
    try:
        exons = int(float(exons_raw))
    except ValueError:
        exons = 0
    if category in {"full-splice_match", "fsm"}:
        minimum = int(config.get("sqanti_min_reads_fsm", "1"))
    elif category in {"novel_in_catalog", "nic"}:
        minimum = int(config.get("sqanti_min_reads_nic", "2"))
    else:
        minimum = int(config.get("sqanti_min_reads_other", "3"))
    if exons <= 1 and category not in {"full-splice_match", "fsm"}:
        minimum = max(minimum, int(config.get("sqanti_min_reads_novel_single_exon", "5")))
    if support < minimum:
        return f"support_lt_{minimum}"
    if as_bool(config, "sqanti_reject_rtswitching", True) and truth(pick(row, ("RTS_stage", "RT_switching", "rts"))):
        return "rt_switching"
    if as_bool(config, "sqanti_reject_intrapriming", True):
        intra = pick(row, ("perc_A_downstream_TTS", "intrapriming", "intra_priming"), "0")
        try:
            if float(intra) >= 60:
                return "intrapriming"
        except ValueError:
            pass
    if as_bool(config, "sqanti_reject_noncanonical", True):
        canonical = pick(row, ("canonical", "all_canonical", "canonical_splice_sites"), "")
        if canonical and canonical.strip().lower() in {"0", "false", "no", "n"}:
            return "noncanonical_splice"
    return "PASS"


def write_filtered_gtf(source: Path, target: Path, keep: set[str]) -> None:
    opener = gzip.open if source.suffix == ".gz" else open
    target.parent.mkdir(parents=True, exist_ok=True)
    with opener(source, "rt", encoding="utf-8") as src, target.open("w", encoding="utf-8") as dst:
        for line in src:
            if line.startswith("#"):
                dst.write(line)
                continue
            match = re.search(r'transcript_id[ ="]+([^"; ]+)', line)
            if match and match.group(1) in keep:
                dst.write(line)


def plain_gtf(source: Path, tmp_dir: Path) -> Path:
    if source.suffix != ".gz":
        return source
    target = tmp_dir / source.with_suffix("").name
    lock_path = tmp_dir / (target.name + ".lock")
    tmp_dir.mkdir(parents=True, exist_ok=True)
    with lock_path.open("w", encoding="utf-8") as lock:
        fcntl.flock(lock.fileno(), fcntl.LOCK_EX)
        if nonempty(target):
            return target
        partial = tmp_dir / f"{target.name}.part.{os.getpid()}"
        with gzip.open(source, "rb") as src, partial.open("wb") as dst:
            shutil.copyfileobj(src, dst, length=1024 * 1024)
        partial.replace(target)
    return target


def main() -> int:
    args = arguments()
    config = read_config(args.config)
    ensure_stage(args.outdir)
    isoquant_root = args.isoquant_dir or args.outdir / "04_tool_output/isoquant"
    sqanti_root = args.outdir / "04_tool_output/sqanti3"
    isoquant_root.mkdir(parents=True, exist_ok=True)
    sqanti_root.mkdir(parents=True, exist_ok=True)

    if not args.skip_tools:
        isoquant = executable(config.get("isoquant_bin", "isoquant"))
        sqanti = executable(config.get("sqanti3_bin", "sqanti3_qc.py"))
        reference_gtf = plain_gtf(args.gtf, args.tmp_dir or args.outdir / "04_tool_output/reference")
        marker = args.outdir / "05_metadata/isoquant.complete.json"
        if not nonempty(marker):
            command = [
                isoquant, "--reference", str(args.fasta), "--genedb", str(reference_gtf),
                "--bam", str(args.bam), "--data_type", config.get("isoquant_data_type", "nanopore"),
                "--prefix", args.sample_id, "--output", str(isoquant_root), "--threads", str(args.threads),
            ]
            if as_bool(config, "isoquant_complete_genedb", True):
                command.append("--complete_genedb")
            if as_bool(config, "isoquant_full_length", True):
                command.append("--fl_data")
            run_logged(command, args.outdir / f"03_logs/{args.sample_id}.isoquant.log")
            write_json(marker, {"status": "complete", "command": command})
        model_gtf = find_one(isoquant_root, ("*.transcript_models.gtf", "*transcript_models.gtf"), "IsoQuant transcript GTF")
        classification = args.sqanti_classification
        if classification is None or not nonempty(classification):
            command = [
                sqanti, "--isoforms", str(model_gtf), "--refGTF", str(reference_gtf), "--refFasta", str(args.fasta),
                "--output", args.sample_id, "--dir", str(sqanti_root), "--cpus", str(args.threads), "--report", "skip",
            ]
            run_logged(command, args.outdir / f"03_logs/{args.sample_id}.sqanti3.log")
            classification = find_one(sqanti_root, ("*classification.txt", "*classification.tsv"), "SQANTI3 classification")
    else:
        model_gtf = find_one(isoquant_root, ("*.transcript_models.gtf", "*transcript_models.gtf"), "IsoQuant transcript GTF")
        classification = args.sqanti_classification or find_one(
            sqanti_root, ("*classification.txt", "*classification.tsv"), "SQANTI3 classification"
        )

    assignments_path = find_one(isoquant_root, ("*read_assignments.tsv.gz", "*read_assignments.tsv"), "IsoQuant read assignments")
    count_matches = sorted(isoquant_root.rglob("*transcript_model_counts.tsv"))
    count_path = count_matches[0] if len(count_matches) == 1 else None
    normalized, assignment_counts, isoform_reads = parse_assignments(assignments_path, args.sample_id)
    count_support = support_from_counts(count_path)
    assignment_table = args.outdir / f"01_tables/{args.sample_id}.read_isoform_assignments.tsv.gz"
    with gzip.open(assignment_table, "wt", encoding="utf-8", newline="") as handle:
        fields = list(normalized[0])
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(normalized)

    catalog: list[dict[str, object]] = []
    keep: set[str] = set()
    class_counts: Counter = Counter()
    for row in rows(classification):
        isoform = pick(row, ("isoform", "isoform_id", "transcript_id"))
        if not isoform:
            continue
        support = max(count_support.get(isoform, 0), len(isoform_reads.get(isoform, set())))
        reason = filter_reason(row, support, config)
        category = pick(row, ("structural_category", "category"), "unknown")
        if reason == "PASS":
            keep.add(isoform)
        class_counts[category] += 1
        catalog.append(
            {
                "sample_id": args.sample_id,
                "isoform_id": f"{args.sample_id}|{isoform}",
                "source_isoform_id": isoform,
                "gene_id": pick(row, ("associated_gene", "gene_id", "gene"), "."),
                "associated_transcript": pick(row, ("associated_transcript", "transcript_id"), "."),
                "structural_category": category,
                "length": pick(row, ("length", "isoform_length"), "0"),
                "exons": pick(row, ("exons", "num_exons", "exon_count"), "0"),
                "read_support": support,
                "filter_status": "PASS" if reason == "PASS" else "FAIL",
                "filter_reason": reason,
            }
        )
    if not catalog:
        raise RuntimeError(f"No isoforms parsed from {classification}")
    write_tsv(args.outdir / f"01_tables/{args.sample_id}.isoform_catalog.tsv", catalog)
    filtered_gtf = args.outdir / f"04_tool_output/filtered/{args.sample_id}.filtered_isoforms.gtf"
    write_filtered_gtf(model_gtf, filtered_gtf, keep)
    validation = {
        "schema_version": CATALOG_SCHEMA,
        "status": "PASS" if keep and assignment_counts["assigned_reads"] else "FAIL",
        "sample_id": args.sample_id,
        "catalog_isoforms": len(catalog),
        "filtered_isoforms": len(keep),
        "assignment_counts": dict(assignment_counts),
        "polya_evaluable_fraction": assignment_counts["polya_evaluable_reads"] / assignment_counts["total_reads"] if assignment_counts["total_reads"] else 0.0,
        "structural_categories": dict(class_counts),
        "structural_category_closure_error": abs(sum(class_counts.values()) - len(catalog)),
        "assignment_partition_closure_error": abs(assignment_counts["total_reads"] - assignment_counts["assigned_reads"] - assignment_counts["unassigned_reads"]),
        "isoquant_read_assignments": str(assignments_path),
        "sqanti_classification": str(classification),
        "filtered_gtf": str(filtered_gtf),
        "reference_gtf": str(reference_gtf) if not args.skip_tools else str(args.gtf),
        "sample_prefixed_ids": True,
        "transcript_assignment_policy": ["unique", "unique_minor_difference"],
    }
    write_json(args.outdir / "05_metadata/catalog_validation.json", validation)
    if validation["status"] != "PASS" or validation["structural_category_closure_error"] or validation["assignment_partition_closure_error"]:
        raise RuntimeError("Isoform catalog validation failed")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
