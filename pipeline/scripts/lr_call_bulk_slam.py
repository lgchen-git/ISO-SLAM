#!/usr/bin/env python3
"""Call bulk long-read SLAM-seq read, site, and gene summaries."""

from __future__ import annotations

import argparse
import csv
import gzip
import json
import logging
import math
import os
import re
import shutil
import subprocess
from collections import Counter, defaultdict
from dataclasses import dataclass
from pathlib import Path
from typing import Iterable

import pysam

from lr_gene_assignment import (
    AssignmentResult,
    assign_region,
    load_exon_annotations,
    load_gene_annotations,
    load_transcript_models,
)


try:
    from scipy.stats import binom as scipy_binom
except ImportError:  # pragma: no cover - exercised only when scipy is absent
    scipy_binom = None


@dataclass
class ReadMeta:
    read_id: str
    chrom: str
    start: int
    end: int
    mapq: int
    read_length: int
    gene_id: str
    gene_name: str
    transcript_id: str
    strand: str
    region_class: str
    assignment_class: str = ""
    assignment_level: str = ""
    assignment_support_tx_count: int = 0
    assignment_support_gene_count: int = 0
    assignment_resolution_method: str = ""


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Call bulk SLAM-seq old/new reads.")
    parser.add_argument("--bam", required=True)
    parser.add_argument("--fasta", required=True)
    parser.add_argument("--gtf", required=True)
    parser.add_argument("--sample-id", required=True)
    parser.add_argument("--outdir", required=True)
    parser.add_argument("--new-threshold", type=int, default=1)
    parser.add_argument("--min-tc-baseq", type=int, default=15)
    parser.add_argument("--min-mapq", type=int, default=20)
    parser.add_argument("--min-t-covered", type=int, default=10)
    parser.add_argument("--bg-tc-rate", type=float, default=None)
    parser.add_argument("--alpha", type=float, default=0.01)
    parser.add_argument("--main-call", choices=["site_binomial", "fixed_threshold", "long_tuc_threshold"], default="site_binomial")
    parser.add_argument("--long-tuc-formal-threshold", type=int, default=20)
    parser.add_argument("--long-tuc-report-thresholds", default="6,20,30")
    parser.add_argument("--min-site-depth", type=int, default=10)
    parser.add_argument("--min-site-tc-count", type=int, default=3)
    parser.add_argument("--min-site-tc-rate", type=float, default=0.2)
    parser.add_argument("--max-site-other-mismatch-rate", type=float, default=0.05)
    parser.add_argument("--max-read-other-mismatch-rate", type=float, default=0.1)
    parser.add_argument("--min-high-conf-tc-sites-hit", type=int, default=1)
    parser.add_argument("--site-evidence-mode", choices=["supportive", "required"], default="supportive")
    parser.add_argument("--assignment-mode", choices=["gene_span", "exon_gene", "transcript_aware_gene"], default="gene_span")
    parser.add_argument(
        "--library-strand-mode",
        choices=["legacy", "forward", "reverse", "unstranded"],
        default="legacy",
        help="Transcript orientation model. legacy preserves --assignment-require-same-strand behavior.",
    )
    parser.add_argument("--assignment-min-exon-overlap-bp", type=int, default=20)
    parser.add_argument("--assignment-require-same-strand", action="store_true")
    parser.add_argument("--assignment-end-distance-tiebreak", action="store_true")
    parser.add_argument("--assignment-allow-gene-level-salvage", action="store_true")
    parser.add_argument("--assignment-count-shared-resolved-reads", action="store_true")
    parser.add_argument("--noisy-site-mode", choices=["exclude", "ambiguous"], default="exclude")
    parser.add_argument("--low-t-zero-tc-as-old", action="store_true")
    parser.add_argument("--low-t-high-conf-tc-as-new", action="store_true")
    parser.add_argument("--bg-from-background-sites", action="store_true",
                        help="Compute bg_tc_rate from background_site/low_depth_site only, excluding high_conf_tc_site signal.")
    parser.add_argument("--grand-slam-bin", default="", help="Optional pinned gedi/GRAND-SLAM executable for readiness audit.")
    parser.add_argument("--grand-slam-java-bin", default="", help="Optional directory containing the pinned Java runtime.")
    return parser.parse_args()


def setup_logging(outdir: Path, sample_id: str) -> None:
    import os
    log_path = Path(os.environ.get("ISO_LOG_ROOT", str(outdir))) / f"{sample_id}.bulk_slam_calling.log"
    log_path.parent.mkdir(parents=True, exist_ok=True)
    logging.basicConfig(
        level=logging.INFO,
        format="%(asctime)s %(levelname)s %(message)s",
        handlers=[logging.FileHandler(log_path, mode="w", encoding="utf-8"), logging.StreamHandler()],
    )


def fetch_ref_base(fasta: pysam.FastaFile, chrom: str, rpos: int) -> str:
    return fasta.fetch(chrom, rpos, rpos + 1).upper()


def is_transcript_t(ref_base: str, strand: str) -> bool:
    return (strand == "+" and ref_base == "T") or (strand == "-" and ref_base == "A")


def classify_t_base(ref_base: str, query_base: str, strand: str) -> str | None:
    ref_base = ref_base.upper()
    query_base = query_base.upper()
    if strand == "+" and ref_base == "T":
        return { "T": "match_t_count", "C": "tc_count", "A": "ta_count", "G": "tg_count" }.get(query_base)
    if strand == "-" and ref_base == "A":
        return { "A": "match_t_count", "G": "tc_count", "T": "ta_count", "C": "tg_count" }.get(query_base)
    return None


def check_base_qualities(read: pysam.AlignedSegment) -> None:
    if read.query_qualities is None:
        raise SystemExit(
            f"Read {read.query_name} has no query base qualities. "
            "Formal site_binomial mode requires base qualities for --min-tc-baseq filtering."
        )


def iter_t_events(
    read: pysam.AlignedSegment,
    fasta: pysam.FastaFile,
    strand: str,
    min_baseq: int,
) -> Iterable[tuple[str, int, str]]:
    check_base_qualities(read)
    if read.has_tag("MD"):
        yield from iter_t_events_from_cigar_md(read, strand, min_baseq)
        return
    query_sequence = (read.query_sequence or "").upper()
    qualities = read.query_qualities
    block_cache = None
    for qpos, rpos in read.get_aligned_pairs(matches_only=False):
        if rpos is None:
            continue
        if block_cache is None:
            block_cache = [
                (start, end, fasta.fetch(read.reference_name, start, end).upper())
                for start, end in read.get_blocks()
            ]
        ref_base = ref_base_from_blocks(block_cache, rpos)
        if not ref_base:
            continue
        ref_base = ref_base.upper()
        if not is_transcript_t(ref_base, strand):
            continue
        if qpos is None:
            yield read.reference_name, rpos, "deletion_count"
            continue
        if qualities[qpos] < min_baseq:
            continue
        event = classify_t_base(ref_base, query_sequence[qpos], strand)
        if event:
            yield read.reference_name, rpos, event


def iter_t_events_from_cigar_md(
    read: pysam.AlignedSegment,
    strand: str,
    min_baseq: int,
) -> Iterable[tuple[str, int, str]]:
    query_sequence = (read.query_sequence or "").upper()
    qualities = read.query_qualities
    ref_sequence = read.get_reference_sequence().upper()
    qpos = 0
    rpos = read.reference_start
    ref_idx = 0
    for op, length in read.cigartuples or []:
        if op in {0, 7, 8}:  # M, =, X
            for offset in range(length):
                ref_base = ref_sequence[ref_idx + offset]
                if is_transcript_t(ref_base, strand) and qualities[qpos + offset] >= min_baseq:
                    event = classify_t_base(ref_base, query_sequence[qpos + offset], strand)
                    if event:
                        yield read.reference_name, rpos + offset, event
            qpos += length
            rpos += length
            ref_idx += length
        elif op == 2:  # deletion
            for offset in range(length):
                ref_base = ref_sequence[ref_idx + offset]
                if is_transcript_t(ref_base, strand):
                    yield read.reference_name, rpos + offset, "deletion_count"
            rpos += length
            ref_idx += length
        elif op == 3:  # skipped reference region
            rpos += length
        elif op in {1, 4}:  # insertion or soft clip
            qpos += length
        elif op in {5, 6}:  # hard clip or padding
            continue


def ref_base_from_blocks(block_cache: list[tuple[int, int, str]], rpos: int) -> str:
    for start, end, ref_seq in block_cache:
        if start <= rpos < end:
            return ref_seq[rpos - start]
    return ""


def empty_base_counter() -> Counter:
    return Counter(
        {
            "match_t_count": 0,
            "tc_count": 0,
            "ta_count": 0,
            "tg_count": 0,
            "deletion_count": 0,
        }
    )


def read_t_counts(read: pysam.AlignedSegment, fasta: pysam.FastaFile, strand: str, min_baseq: int) -> Counter:
    counts = empty_base_counter()
    for _chrom, _rpos, event in iter_t_events(read, fasta, strand, min_baseq):
        counts[event] += 1
    counts["t_covered"] = counts["match_t_count"] + counts["tc_count"] + counts["ta_count"] + counts["tg_count"]
    return counts


def read_t_counts_and_site_hits(
    read: pysam.AlignedSegment,
    fasta: pysam.FastaFile,
    strand: str,
    min_baseq: int,
    site_calls: dict[tuple[str, int], str],
    noisy_site_mode: str,
) -> tuple[Counter, int, int]:
    counts = empty_base_counter()
    high_conf_hits = 0
    noisy_hits = 0
    for chrom, rpos, event in iter_t_events(read, fasta, strand, min_baseq):
        if event == "deletion_count":
            counts[event] += 1
            continue
        site_call = site_calls.get((chrom, rpos), "")
        if site_call == "noisy_site":
            noisy_hits += int(event in {"tc_count", "ta_count", "tg_count"})
            if noisy_site_mode == "exclude":
                continue
        counts[event] += 1
        high_conf_hits += int(site_call == "high_conf_tc_site")
    counts["t_covered"] = counts["match_t_count"] + counts["tc_count"] + counts["ta_count"] + counts["tg_count"]
    return counts, high_conf_hits, noisy_hits


def site_rates(counts: Counter) -> dict[str, float | int | str]:
    t_depth = counts["match_t_count"] + counts["tc_count"] + counts["ta_count"] + counts["tg_count"]
    return {
        "t_depth": t_depth,
        "tc_rate": counts["tc_count"] / t_depth if t_depth else 0.0,
        "ta_rate": counts["ta_count"] / t_depth if t_depth else 0.0,
        "tg_rate": counts["tg_count"] / t_depth if t_depth else 0.0,
        "other_mismatch_rate": (counts["ta_count"] + counts["tg_count"]) / t_depth if t_depth else 0.0,
    }


def classify_site(counts: Counter, args: argparse.Namespace) -> str:
    rates = site_rates(counts)
    if rates["t_depth"] < args.min_site_depth:
        return "low_depth_site"
    if rates["other_mismatch_rate"] > args.max_site_other_mismatch_rate:
        return "noisy_site"
    if (
        counts["tc_count"] >= args.min_site_tc_count
        and rates["tc_rate"] >= args.min_site_tc_rate
        and rates["other_mismatch_rate"] <= args.max_site_other_mismatch_rate
    ):
        return "high_conf_tc_site"
    return "background_site"


def binomial_tail(k: int, n: int, p: float) -> float:
    if n <= 0:
        return math.nan
    if k <= 0:
        return 1.0
    p = min(max(p, 0.0), 1.0)
    if scipy_binom is not None:
        return float(scipy_binom.sf(k - 1, n, p))
    if p <= 0.0:
        return 1.0 if k <= 0 else 0.0
    if p >= 1.0:
        return 1.0
    if k <= n * p:
        return 1.0
    log_pmf = (
        math.lgamma(n + 1)
        - math.lgamma(k + 1)
        - math.lgamma(n - k + 1)
        + k * math.log(p)
        + (n - k) * math.log1p(-p)
    )
    term = math.exp(log_pmf) if log_pmf > -745 else 0.0
    tail = 0.0
    x = k
    while x <= n and term > 0.0:
        tail += term
        if tail >= 1.0:
            return 1.0
        x += 1
        if x > n:
            break
        term *= ((n - x + 1) / x) * (p / (1.0 - p))
        if term < 1e-320:
            break
    return min(max(tail, 0.0), 1.0)


def format_pvalue(value: float) -> str:
    if math.isnan(value):
        return "NA"
    return f"{value:.6g}"


def collect_site_counts(
    bam_path: str,
    fasta: pysam.FastaFile,
    exon_annotations,
    gene_annotations,
    model_lists,
    start_lists,
    args: argparse.Namespace,
) -> tuple[dict[tuple[str, int], Counter], dict[str, ReadMeta], Counter]:
    site_counts: dict[tuple[str, int], Counter] = defaultdict(empty_base_counter)
    read_meta: dict[str, ReadMeta] = {}
    qc = Counter()
    with pysam.AlignmentFile(bam_path, "rb") as bam:
        for read in bam.fetch(until_eof=True):
            qc["total_reads"] += 1
            if qc["total_reads"] % 100000 == 0:
                logging.info("site pass scanned %s reads", qc["total_reads"])
            if read.is_unmapped or read.is_secondary or read.is_supplementary:
                qc["skipped_non_primary"] += 1
                continue
            assignment = assign_region(
                read,
                exon_annotations,
                gene_annotations,
                model_lists,
                start_lists,
                args,
            )
            qc[f"region_{assignment.region_class}_reads"] += 1
            qc[f"assignment_{assignment.assignment_class}_reads"] += 1
            meta = ReadMeta(
                read_id=read.query_name,
                chrom=read.reference_name,
                start=read.reference_start,
                end=read.reference_end,
                mapq=read.mapping_quality,
                read_length=read.query_length or len(read.query_sequence or ""),
                gene_id=assignment.gene_id,
                gene_name=assignment.gene_name,
                transcript_id=assignment.transcript_id,
                strand=assignment.strand,
                region_class=assignment.region_class,
                assignment_class=assignment.assignment_class,
                assignment_level=assignment.assignment_level,
                assignment_support_tx_count=assignment.assignment_support_tx_count,
                assignment_support_gene_count=assignment.assignment_support_gene_count,
                assignment_resolution_method=assignment.assignment_resolution_method,
            )
            read_meta[read.query_name] = meta
            if read.mapping_quality < args.min_mapq:
                qc["low_mapq_reads"] += 1
                continue
            if meta.assignment_level == "unresolved":
                qc["unassigned_gene"] += 1
                continue
            for chrom, rpos, event in iter_t_events(read, fasta, meta.strand, args.min_tc_baseq):
                site_counts[(chrom, rpos)][event] += 1
    return site_counts, read_meta, qc


def write_site_table(
    outdir: Path,
    sample_id: str,
    site_counts: dict[tuple[str, int], Counter],
    fasta: pysam.FastaFile,
    args: argparse.Namespace,
) -> tuple[dict[tuple[str, int], str], Counter, Counter, Counter]:
    site_calls: dict[tuple[str, int], str] = {}
    site_call_counts = Counter()
    totals = Counter()
    bg_totals = Counter()
    path = outdir / f"{sample_id}.bulk_t_sites.tsv.gz"
    with gzip.open(path, "wt") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(
            [
                "chrom",
                "pos",
                "ref_base",
                "t_depth",
                "match_t_count",
                "tc_count",
                "ta_count",
                "tg_count",
                "deletion_count",
                "tc_rate",
                "ta_rate",
                "tg_rate",
                "other_mismatch_rate",
                "site_call",
            ]
        )
        for (chrom, rpos), counts in sorted(site_counts.items()):
            rates = site_rates(counts)
            site_call = classify_site(counts, args)
            site_calls[(chrom, rpos)] = site_call
            site_call_counts[site_call] += 1
            for key in ("match_t_count", "tc_count", "ta_count", "tg_count", "deletion_count"):
                totals[key] += counts[key]
            totals["t_depth"] += rates["t_depth"]
            if site_call in {"background_site", "low_depth_site"}:
                for key in ("match_t_count", "tc_count", "ta_count", "tg_count", "deletion_count"):
                    bg_totals[key] += counts[key]
                bg_totals["t_depth"] += rates["t_depth"]
            writer.writerow(
                [
                    chrom,
                    rpos + 1,
                    fetch_ref_base(fasta, chrom, rpos),
                    rates["t_depth"],
                    counts["match_t_count"],
                    counts["tc_count"],
                    counts["ta_count"],
                    counts["tg_count"],
                    counts["deletion_count"],
                    f"{rates['tc_rate']:.8f}",
                    f"{rates['ta_rate']:.8f}",
                    f"{rates['tg_rate']:.8f}",
                    f"{rates['other_mismatch_rate']:.8f}",
                    site_call,
                ]
            )
    return site_calls, site_call_counts, totals, bg_totals


def choose_background_rate(args: argparse.Namespace, totals: Counter, bg_totals: Counter | None = None) -> tuple[float, str]:
    if args.bg_tc_rate is not None:
        if args.bg_tc_rate < 0 or args.bg_tc_rate > 1:
            raise SystemExit("--bg-tc-rate must be between 0 and 1")
        return args.bg_tc_rate, "provided_bg_tc_rate"
    # When bg-from-background-sites is enabled, use only background/low_depth sites
    # to avoid inflating bg_rate with true T>C signal from high_conf_tc_sites.
    source_counts = bg_totals if (args.bg_from_background_sites and bg_totals and bg_totals.get("t_depth", 0) > 0) else totals
    source_label = "background_site_filtered" if source_counts is bg_totals else "global"
    t_depth = source_counts["t_depth"]
    ta_rate = source_counts["ta_count"] / t_depth if t_depth else 0.0
    tg_rate = source_counts["tg_count"] / t_depth if t_depth else 0.0
    return max(ta_rate, tg_rate, 1e-5), f"max_ta_tg_rate_{source_label}"


def call_for_read(
    meta: ReadMeta,
    counts: Counter,
    high_conf_hits: int,
    noisy_hits: int,
    binom_p: float,
    args: argparse.Namespace,
) -> tuple[str, str, str]:
    t_covered = counts["t_covered"]
    other_rate = (counts["ta_count"] + counts["tg_count"]) / t_covered if t_covered else 0.0
    reasons = []
    warnings = []
    if meta.mapq < args.min_mapq:
        reasons.append("low_mapq")
    if meta.assignment_class == "intergenic":
        reasons.append("no_gene")
    if meta.assignment_class == "shared_gene_ambiguous":
        reasons.append("gene_assignment_tie")
    if meta.assignment_class == "antisense_only":
        reasons.append("antisense_only")
    if meta.assignment_class == "shared_antisense_ambiguous":
        reasons.append("shared_antisense_ambiguous")
    if other_rate > args.max_read_other_mismatch_rate:
        reasons.append("high_other_mismatch")
    if args.noisy_site_mode == "ambiguous" and noisy_hits > 0:
        reasons.append("noisy_site_mismatch")
    if reasons:
        return "ambiguous", ",".join(reasons), "technical_ambiguous"
    if args.main_call == "long_tuc_threshold":
        if counts["tc_count"] > args.long_tuc_formal_threshold:
            return "new", f"long_tuc_gt{args.long_tuc_formal_threshold}", "high_conf_new"
        return "old", f"long_tuc_le{args.long_tuc_formal_threshold}", "high_conf_old"
    if args.main_call == "fixed_threshold":
        if t_covered < args.min_t_covered and counts["tc_count"] == 0 and args.low_t_zero_tc_as_old:
            return "old", "low_t_zero_tc", "low_conf_old"
        if t_covered < args.min_t_covered:
            return "ambiguous", "low_t_covered", "low_t_ambiguous"
        return ("new", "fixed_threshold", "high_conf_new") if counts["tc_count"] >= args.new_threshold else ("old", "no_tc_threshold", "high_conf_old")
    if t_covered < args.min_t_covered:
        if counts["tc_count"] == 0 and args.low_t_zero_tc_as_old:
            return "old", "low_t_zero_tc", "low_conf_old"
        low_t_new = (
            args.low_t_high_conf_tc_as_new
            and counts["tc_count"] > 0
            and high_conf_hits >= args.min_high_conf_tc_sites_hit
            and not math.isnan(binom_p)
            and binom_p < args.alpha
        )
        if low_t_new:
            return "new", "low_t_high_conf_tc", "low_conf_new"
        return "ambiguous", "low_t_covered", "low_t_ambiguous"
    is_new = binom_p < args.alpha
    if args.site_evidence_mode == "required":
        is_new = is_new and high_conf_hits >= args.min_high_conf_tc_sites_hit
    if is_new:
        return "new", "binomial_tc_signal", "high_conf_new"
    return "old", "not_significant", "high_conf_old"


def write_read_and_gene_tables(
    bam_path: str,
    fasta: pysam.FastaFile,
    read_meta: dict[str, ReadMeta],
    site_calls: dict[tuple[str, int], str],
    bg_rate: float,
    args: argparse.Namespace,
) -> tuple[Counter, dict[str, Counter], dict[str, str]]:
    qc = Counter()
    gene_summary: dict[str, Counter] = defaultdict(Counter)
    gene_names: dict[str, str] = {}
    read_path = Path(args.outdir) / f"{args.sample_id}.bulk_read_calls.tsv.gz"
    long_tuc_thresholds = parse_thresholds(args.long_tuc_report_thresholds)
    with pysam.AlignmentFile(bam_path, "rb") as bam, gzip.open(read_path, "wt") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        header = [
            "read_id",
            "gene_id",
            "transcript_id",
            "chrom",
            "start",
            "end",
            "strand",
            "mapq",
            "read_length",
            "region_class",
            "assignment_class",
            "assignment_level",
            "assignment_support_tx_count",
            "assignment_support_gene_count",
            "assignment_resolution_method",
            "t_covered",
            "tc_count",
            "ta_count",
            "tg_count",
            "tc_rate",
            "other_mismatch_rate",
            "high_conf_tc_sites_hit",
            "noisy_sites_hit",
            "binomial_p",
            "final_call",
            "ambiguity_reason",
            "call_confidence",
            "gene_name",
            "tc_events",
            "t_to_a",
            "t_to_c",
            "t_to_g",
            "long_tuc_formal_threshold",
        ]
        header.extend([f"long_tuc_gt{threshold}" for threshold in long_tuc_thresholds])
        header.append("assignment")
        writer.writerow(header)
        for read in bam.fetch(until_eof=True):
            if read.is_unmapped or read.is_secondary or read.is_supplementary:
                continue
            qc["read_pass_reads"] += 1
            if qc["read_pass_reads"] % 100000 == 0:
                logging.info("read pass processed %s reads", qc["read_pass_reads"])
            meta = read_meta.get(read.query_name)
            if meta is None:
                continue
            if meta.region_class in {"intergenic", "ambiguous"}:
                counts = empty_base_counter()
                counts["t_covered"] = 0
                high_conf_hits = 0
                noisy_hits = 0
            else:
                counts, high_conf_hits, noisy_hits = read_t_counts_and_site_hits(
                    read,
                    fasta,
                    meta.strand,
                    args.min_tc_baseq,
                    site_calls,
                    args.noisy_site_mode,
                )
            t_covered = counts["t_covered"]
            tc_rate = counts["tc_count"] / t_covered if t_covered else 0.0
            other_rate = (counts["ta_count"] + counts["tg_count"]) / t_covered if t_covered else 0.0
            binom_p = (
                binomial_tail(counts["tc_count"], t_covered, bg_rate)
                if t_covered > 0
                else math.nan
            )
            final_call, ambiguity_reason, call_confidence = call_for_read(
                meta,
                counts,
                high_conf_hits,
                noisy_hits,
                binom_p,
                args,
            )
            long_tuc_flags = {threshold: int(counts["tc_count"] > threshold) for threshold in long_tuc_thresholds}
            writer.writerow(
                [
                    meta.read_id,
                    meta.gene_id,
                    meta.transcript_id,
                    meta.chrom,
                    meta.start,
                    meta.end,
                    meta.strand,
                    meta.mapq,
            meta.read_length,
            meta.region_class,
            meta.assignment_class,
            meta.assignment_level,
            meta.assignment_support_tx_count,
            meta.assignment_support_gene_count,
            meta.assignment_resolution_method,
            t_covered,
                    counts["tc_count"],
                    counts["ta_count"],
                    counts["tg_count"],
                    f"{tc_rate:.8f}",
                    f"{other_rate:.8f}",
                    high_conf_hits,
                    noisy_hits,
                    format_pvalue(binom_p),
                    final_call,
                    ambiguity_reason,
                    call_confidence,
                    meta.gene_name,
                    counts["tc_count"],
                    counts["ta_count"],
                    counts["tc_count"],
                    counts["tg_count"],
                    args.long_tuc_formal_threshold,
                    *[long_tuc_flags[threshold] for threshold in long_tuc_thresholds],
                    final_call,
                ]
            )
            qc[f"{final_call}_reads"] += 1
            qc[f"{call_confidence}_reads"] += 1
            for threshold in long_tuc_thresholds:
                if long_tuc_flags[threshold]:
                    qc[f"long_tuc_gt{threshold}_reads"] += 1
            if ambiguity_reason:
                for reason in ambiguity_reason.split(","):
                    qc[f"reason_{reason}_reads"] += 1
            qc["high_quality_reads"] += int(
                meta.mapq >= args.min_mapq
                and t_covered >= args.min_t_covered
                and other_rate <= args.max_read_other_mismatch_rate
            )
            qc["total_t_covered"] += t_covered
            qc["total_tc_count"] += counts["tc_count"]
            qc["total_ta_count"] += counts["ta_count"]
            qc["total_tg_count"] += counts["tg_count"]
            qc["assigned_reads"] += int(meta.assignment_level != "unresolved")
            qc[f"assignment_level_{meta.assignment_level}_reads"] += 1
            if meta.assignment_level == "gene" and meta.assignment_class in {"gene_unique_multi_tx", "antisense_gene_salvaged"}:
                qc["gene_level_salvaged_reads"] += 1
            if meta.assignment_level == "gene" and meta.assignment_class == "shared_gene_resolved_by_ends":
                qc["shared_gene_resolved_reads"] += 1
            if meta.assignment_level == "unresolved":
                continue
            summary = gene_summary[meta.gene_id]
            summary["total_reads"] += 1
            summary[f"{final_call}_reads"] += 1
            summary[f"{call_confidence}_reads"] += 1
            summary["tc_events"] += counts["tc_count"]
            summary["t_covered"] += t_covered
            summary["t_to_a"] += counts["ta_count"]
            summary["t_to_c"] += counts["tc_count"]
            summary["t_to_g"] += counts["tg_count"]
            for threshold in long_tuc_thresholds:
                if long_tuc_flags[threshold]:
                    summary[f"long_tuc_gt{threshold}_reads"] += 1
            summary[f"{meta.region_class}_reads"] += 1
            summary[f"{meta.assignment_class}_reads"] += 1
            gene_names[meta.gene_id] = meta.gene_name
    return qc, gene_summary, gene_names


def parse_thresholds(raw: str) -> list[int]:
    thresholds = []
    for item in str(raw or "").split(","):
        item = item.strip()
        if not item:
            continue
        threshold = int(item)
        if threshold < 0:
            raise SystemExit("--long-tuc-report-thresholds must be non-negative integers")
        thresholds.append(threshold)
    return sorted(set(thresholds))


def write_gene_summary(
    outdir: Path,
    sample_id: str,
    gene_summary: dict[str, Counter],
    gene_names: dict[str, str],
    long_tuc_thresholds: list[int],
) -> None:
    gene_path = outdir / f"{sample_id}.bulk_gene_summary.tsv"
    with gene_path.open("w", encoding="utf-8") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        header = [
            "gene_id",
            "gene_name",
            "total_reads",
            "new_reads",
            "old_reads",
            "ambiguous_reads",
            "new_fraction",
            "new_fraction_excluding_ambiguous",
            "gene_reads",
            "exon_reads",
            "intron_reads",
            "tx_unique_reads",
            "gene_unique_multi_tx_reads",
            "shared_gene_resolved_by_ends_reads",
            "intron_only_unique_gene_reads",
            "tc_events",
            "t_covered",
            "tc_rate",
            "t_to_a",
            "t_to_g",
        ]
        for threshold in long_tuc_thresholds:
            header.extend([f"long_tuc_gt{threshold}_reads", f"long_tuc_gt{threshold}_fraction"])
        writer.writerow(header)
        for gene_id, summary in sorted(gene_summary.items(), key=lambda item: item[1]["total_reads"], reverse=True):
            total_reads = summary["total_reads"]
            called_reads = summary["new_reads"] + summary["old_reads"]
            values = [
                gene_id,
                gene_names.get(gene_id, gene_id),
                total_reads,
                summary["new_reads"],
                summary["old_reads"],
                summary["ambiguous_reads"],
                summary["new_reads"] / total_reads if total_reads else 0.0,
                summary["new_reads"] / called_reads if called_reads else 0.0,
                summary["gene_reads"],
                summary["exon_reads"],
                summary["intron_reads"],
                summary["tx_unique_reads"],
                summary["gene_unique_multi_tx_reads"],
                summary["shared_gene_resolved_by_ends_reads"],
                summary["intron_only_unique_gene_reads"],
                summary["tc_events"],
                summary["t_covered"],
                summary["t_to_c"] / summary["t_covered"] if summary["t_covered"] else 0.0,
                summary["t_to_a"],
                summary["t_to_g"],
            ]
            for threshold in long_tuc_thresholds:
                threshold_reads = summary[f"long_tuc_gt{threshold}_reads"]
                values.extend([threshold_reads, threshold_reads / total_reads if total_reads else 0.0])
            writer.writerow(values)


def write_ntr_sufficient_stats(
    outdir: Path,
    sample_id: str,
    gene_summary: dict[str, Counter],
    gene_names: dict[str, str],
    bg_rate: float,
    bg_rate_source: str,
) -> None:
    """Write the gene-level k/n statistics required by mixture-model pilots."""
    path = outdir / f"{sample_id}.ntr_sufficient_stats.tsv"
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(
            [
                "gene_id",
                "gene_name",
                "informative_reads",
                "t_covered",
                "tc_count",
                "ta_count",
                "tg_count",
                "diagnostic_background_rate",
                "background_source",
            ]
        )
        for gene_id, summary in sorted(gene_summary.items(), key=lambda item: item[1]["total_reads"], reverse=True):
            writer.writerow(
                [
                    gene_id,
                    gene_names.get(gene_id, gene_id),
                    summary["total_reads"],
                    summary["t_covered"],
                    summary["t_to_c"],
                    summary["t_to_a"],
                    summary["t_to_g"],
                    f"{bg_rate:.10f}",
                    bg_rate_source,
                ]
            )


def write_report_summaries(
    outdir: Path,
    sample_id: str,
    qc: Counter,
    site_call_counts: Counter,
    bg_rate: float,
    bg_rate_source: str,
) -> None:
    """Write exact, compact report inputs without sampling large detail tables."""
    total_reads = int(qc.get("read_pass_reads", qc.get("total_reads", 0)))
    raw_classes = [
        "tx_unique",
        "gene_unique_multi_tx",
        "shared_gene_resolved_by_ends",
        "intron_only_unique_gene",
        "gene_overlap_unresolved",
        "antisense_gene_salvaged",
        "shared_gene_ambiguous",
        "shared_antisense_ambiguous",
        "antisense_only",
        "intergenic",
    ]
    grouped = [
        (
            "transcript_supported",
            qc.get("assignment_tx_unique_reads", 0) + qc.get("assignment_gene_unique_multi_tx_reads", 0),
        ),
        (
            "same_strand_gene_only",
            qc.get("assignment_shared_gene_resolved_by_ends_reads", 0)
            + qc.get("assignment_intron_only_unique_gene_reads", 0)
            + qc.get("assignment_gene_overlap_unresolved_reads", 0),
        ),
        ("antisense_rescue", qc.get("assignment_antisense_gene_salvaged_reads", 0)),
        (
            "ambiguous_shared",
            qc.get("assignment_shared_gene_ambiguous_reads", 0)
            + qc.get("assignment_shared_antisense_ambiguous_reads", 0)
            + qc.get("assignment_antisense_only_reads", 0),
        ),
        ("unassigned_intergenic", qc.get("assignment_intergenic_reads", 0)),
    ]
    assignment_path = outdir / f"{sample_id}.assignment_summary.tsv"
    with assignment_path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(["level", "category", "reads", "denominator", "fraction"])
        for category, count in grouped:
            writer.writerow(["report_group", category, count, total_reads, count / total_reads if total_reads else 0.0])
        for category in raw_classes:
            count = int(qc.get(f"assignment_{category}_reads", 0))
            writer.writerow(["raw_class", category, count, total_reads, count / total_reads if total_reads else 0.0])

    total_t = int(qc.get("total_t_covered", 0))
    substitution_path = outdir / f"{sample_id}.substitution_summary.tsv"
    with substitution_path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(["substitution", "count", "denominator", "rate", "role", "source"])
        for label, key in (("T>A", "total_ta_count"), ("T>C", "total_tc_count"), ("T>G", "total_tg_count")):
            count = int(qc.get(key, 0))
            writer.writerow([label, count, total_t, count / total_t if total_t else 0.0, "observed", "all_assigned_T_bases"])
        writer.writerow(["diagnostic_background", "", "", bg_rate, "diagnostic_only", bg_rate_source])

    site_path = outdir / f"{sample_id}.site_call_summary.tsv"
    with site_path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(["site_call", "sites"])
        for label in ("high_conf_tc_site", "background_site", "low_depth_site", "noisy_site"):
            writer.writerow([label, int(site_call_counts.get(label, 0))])


def write_grand_slam_compatibility(
    outdir: Path,
    sample_id: str,
    bam_path: str,
    library_strand_mode: str,
    executable_hint: str = "",
    java_bin: str = "",
) -> None:
    """Record BAM compatibility and verify that the pinned Grand3 command starts."""
    examined = 0
    md_present = 0
    nh_present = 0
    spliced = 0
    with pysam.AlignmentFile(bam_path, "rb") as bam:
        for read in bam.fetch(until_eof=True):
            if read.is_unmapped or read.is_secondary or read.is_supplementary:
                continue
            examined += 1
            md_present += int(read.has_tag("MD"))
            nh_present += int(read.has_tag("NH"))
            spliced += int(any(op == 3 for op, _length in (read.cigartuples or [])))
            if examined >= 10000:
                break
    hinted = Path(executable_hint) if executable_hint else None
    executable = (
        str(hinted.resolve())
        if hinted is not None and hinted.is_file()
        else shutil.which("gedi") or shutil.which("GRAND-SLAM") or shutil.which("grand-slam")
    )
    reasons = []
    warnings = []
    if not executable:
        reasons.append("grand_slam_executable_not_found")
    if examined and md_present != examined:
        reasons.append("md_tag_incomplete")
    if examined and nh_present != examined:
        warnings.append("nh_tag_incomplete_primary_q20_bam")
    command = []
    command_exit_code = None
    command_output = ""
    grand3_version = ""
    if executable:
        gedi_home = outdir / ".gedi_home"
        (gedi_home / ".gedi").mkdir(parents=True, exist_ok=True)
        command = [
            executable,
            "-mem",
            "2G",
            f"-Duser.home={gedi_home}",
            "-e",
            "Grand3",
            "-h",
        ]
        command_env = os.environ.copy()
        if java_bin:
            command_env["PATH"] = java_bin + os.pathsep + command_env.get("PATH", "")
        try:
            completed = subprocess.run(
                command,
                check=False,
                capture_output=True,
                text=True,
                timeout=120,
                env=command_env,
            )
            command_exit_code = completed.returncode
            command_output = (completed.stdout + "\n" + completed.stderr).strip()
            match = re.search(r"GRAND3 version\s+([^\s]+)", command_output, flags=re.IGNORECASE)
            grand3_version = match.group(1) if match else ""
            if completed.returncode != 0 or not grand3_version:
                reasons.append("grand3_command_self_test_failed")
        except (OSError, subprocess.TimeoutExpired) as exc:
            command_output = str(exc)
            reasons.append("grand3_command_self_test_failed")
    command_output = re.sub(r"\x1b\[[0-9;]*m", "", command_output)
    payload = {
        "schema_version": "iso-slam-grand-slam-compatibility-v1",
        "sample_id": sample_id,
        "status": "ready" if not reasons else "blocked",
        "library_strand_mode": library_strand_mode,
        "examined_primary_reads": examined,
        "md_tag_fraction": md_present / examined if examined else 0.0,
        "nh_tag_fraction": nh_present / examined if examined else 0.0,
        "spliced_read_fraction": spliced / examined if examined else 0.0,
        "executable": executable or "",
        "java_bin": java_bin,
        "command": command,
        "command_exit_code": command_exit_code,
        "grand3_version": grand3_version,
        "command_output_tail": command_output[-2000:],
        "reasons": reasons,
        "warnings": warnings,
        "note": "Compatibility audit only; no per-read new/old certainty is inferred here.",
    }
    (outdir / f"{sample_id}.grand_slam_compatibility.json").write_text(
        json.dumps(payload, indent=2) + "\n", encoding="utf-8"
    )


def write_sample_summary(
    outdir: Path,
    sample_id: str,
    qc: Counter,
    site_call_counts: Counter,
    bg_rate: float,
    bg_rate_source: str,
    long_tuc_thresholds: list[int],
    assignment_mode: str,
    library_strand_mode: str,
) -> dict[str, object]:
    total_called = qc["new_reads"] + qc["old_reads"] + qc["ambiguous_reads"]
    total_t = qc["total_t_covered"]
    summary = {
        "sample": sample_id,
        "total_reads": qc["total_reads"],
        "assigned_reads": qc["assigned_reads"],
        "assignment_mode": assignment_mode,
        "library_strand_mode": library_strand_mode,
        "high_quality_reads": qc["high_quality_reads"],
        "total_t_covered": total_t,
        "total_tc_count": qc["total_tc_count"],
        "total_ta_count": qc["total_ta_count"],
        "total_tg_count": qc["total_tg_count"],
        "global_tc_rate": qc["total_tc_count"] / total_t if total_t else 0.0,
        "global_ta_rate": qc["total_ta_count"] / total_t if total_t else 0.0,
        "global_tg_rate": qc["total_tg_count"] / total_t if total_t else 0.0,
        "global_other_mismatch_rate": (qc["total_ta_count"] + qc["total_tg_count"]) / total_t if total_t else 0.0,
        "bg_tc_rate": bg_rate,
        "bg_rate_source": bg_rate_source,
        "high_conf_tc_site_count": site_call_counts["high_conf_tc_site"],
        "noisy_site_count": site_call_counts["noisy_site"],
        "low_depth_site_count": site_call_counts["low_depth_site"],
        "background_site_count": site_call_counts["background_site"],
        "new_read_count": qc["new_reads"],
        "old_read_count": qc["old_reads"],
        "ambiguous_read_count": qc["ambiguous_reads"],
        "gene_reads": qc.get("region_gene_reads", 0),
        "exon_reads": qc.get("region_exon_reads", 0),
        "intron_reads": qc.get("region_intron_reads", 0),
        "intergenic_reads": qc.get("region_intergenic_reads", 0),
        "region_ambiguous_reads": qc.get("region_ambiguous_reads", 0),
        "tx_unique_reads": qc.get("assignment_tx_unique_reads", 0),
        "gene_unique_multi_tx_reads": qc.get("assignment_gene_unique_multi_tx_reads", 0),
        "shared_gene_resolved_reads": qc.get("assignment_shared_gene_resolved_by_ends_reads", 0),
        "shared_gene_ambiguous_reads": qc.get("assignment_shared_gene_ambiguous_reads", 0),
        "gene_overlap_unresolved_reads": qc.get("assignment_gene_overlap_unresolved_reads", 0),
        "antisense_only_reads": qc.get("assignment_antisense_only_reads", 0),
        "antisense_gene_salvaged_reads": qc.get("assignment_antisense_gene_salvaged_reads", 0),
        "shared_antisense_ambiguous_reads": qc.get("assignment_shared_antisense_ambiguous_reads", 0),
        "intron_only_unique_gene_reads": qc.get("assignment_intron_only_unique_gene_reads", 0),
        "gene_level_salvaged_reads": qc.get("gene_level_salvaged_reads", 0),
        "low_conf_old_read_count": qc["low_conf_old_reads"],
        "low_conf_new_read_count": qc["low_conf_new_reads"],
        "technical_ambiguous_read_count": qc["technical_ambiguous_reads"],
        "low_t_ambiguous_read_count": qc["low_t_ambiguous_reads"],
        "low_t_zero_tc_read_count": qc["reason_low_t_zero_tc_reads"],
        "low_t_high_conf_tc_read_count": qc["reason_low_t_high_conf_tc_reads"],
        "no_gene_read_count": qc["reason_no_gene_reads"],
        "gene_assignment_tie_read_count": qc["reason_gene_assignment_tie_reads"],
        "gene_overlap_unresolved_read_count": qc.get("assignment_gene_overlap_unresolved_reads", 0),
        "antisense_only_read_count": qc["reason_antisense_only_reads"],
        "high_other_mismatch_read_count": qc["reason_high_other_mismatch_reads"],
        "noisy_site_mismatch_read_count": qc["reason_noisy_site_mismatch_reads"],
        "new_fraction": qc["new_reads"] / total_called if total_called else 0.0,
        "old_fraction": qc["old_reads"] / total_called if total_called else 0.0,
        "ambiguous_fraction": qc["ambiguous_reads"] / total_called if total_called else 0.0,
        "new_fraction_excluding_ambiguous": qc["new_reads"] / (qc["new_reads"] + qc["old_reads"])
        if qc["new_reads"] + qc["old_reads"]
        else 0.0,
    }
    for threshold in long_tuc_thresholds:
        count = qc[f"long_tuc_gt{threshold}_reads"]
        summary[f"long_tuc_gt{threshold}_read_count"] = count
        summary[f"long_tuc_gt{threshold}_fraction"] = count / total_called if total_called else 0.0
    path = outdir / f"{sample_id}.bulk_slam_summary.tsv"
    with path.open("w", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, delimiter="\t", lineterminator="\n", fieldnames=list(summary.keys()))
        writer.writeheader()
        writer.writerow(summary)
    return summary


def main() -> None:
    args = parse_args()
    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    setup_logging(outdir, args.sample_id)
    logging.info("bulk SLAM caller parameters: %s", json.dumps(vars(args), sort_keys=True))
    if args.main_call == "site_binomial":
        logging.info("--new-threshold is retained for compatibility and ignored by site_binomial main calling")
    if args.main_call == "long_tuc_threshold":
        logging.info("long-TUC formal new call uses T>C events > %s", args.long_tuc_formal_threshold)
    logging.info("loading annotations assignment_mode=%s", args.assignment_mode)
    long_tuc_thresholds = parse_thresholds(args.long_tuc_report_thresholds)
    exon_annotations = load_exon_annotations(args.gtf) if args.assignment_mode == "exon_gene" else None
    gene_annotations = load_gene_annotations(args.gtf)
    _transcript_models, model_lists, start_lists = load_transcript_models(args.gtf) if args.assignment_mode == "transcript_aware_gene" else ({}, {}, {})
    logging.info("annotations loaded")
    fasta = pysam.FastaFile(args.fasta)

    site_counts, read_meta, first_qc = collect_site_counts(args.bam, fasta, exon_annotations, gene_annotations, model_lists, start_lists, args)
    site_calls, site_call_counts, site_totals, bg_site_totals = write_site_table(outdir, args.sample_id, site_counts, fasta, args)
    bg_rate, bg_rate_source = choose_background_rate(args, site_totals, bg_site_totals)
    logging.info("bg_tc_rate=%s bg_rate_source=%s", bg_rate, bg_rate_source)
    read_qc, gene_summary, gene_names = write_read_and_gene_tables(args.bam, fasta, read_meta, site_calls, bg_rate, args)
    write_gene_summary(outdir, args.sample_id, gene_summary, gene_names, long_tuc_thresholds)

    read_qc.update(first_qc)
    read_qc["genes_detected"] = len(gene_summary)
    read_qc["assigned_gene_reads"] = read_qc["assigned_reads"]
    read_qc["gene_reads"] = read_qc["region_gene_reads"]
    read_qc["exon_reads"] = read_qc["region_exon_reads"]
    read_qc["intron_reads"] = read_qc["region_intron_reads"]
    read_qc["intergenic_reads"] = read_qc["region_intergenic_reads"]
    read_qc["region_ambiguous_reads"] = (
        read_qc.get("assignment_shared_gene_ambiguous_reads", 0)
        + read_qc.get("assignment_antisense_only_reads", 0)
        + read_qc.get("assignment_shared_antisense_ambiguous_reads", 0)
    )
    read_qc["total_t_covered"] = site_totals["t_depth"]
    read_qc["total_tc_count"] = site_totals["tc_count"]
    read_qc["total_ta_count"] = site_totals["ta_count"]
    read_qc["total_tg_count"] = site_totals["tg_count"]
    read_qc["t_covered"] = read_qc["total_t_covered"]
    read_qc["t_to_a"] = read_qc["total_ta_count"]
    read_qc["t_to_c"] = read_qc["total_tc_count"]
    read_qc["t_to_g"] = read_qc["total_tg_count"]
    for key in ("t_to_a", "t_to_c", "t_to_g"):
        read_qc[f"{key}_rate"] = read_qc[key] / read_qc["t_covered"] if read_qc["t_covered"] else 0.0
    read_qc["main_call"] = args.main_call
    read_qc["long_tuc_formal_threshold"] = args.long_tuc_formal_threshold
    read_qc["long_tuc_report_thresholds"] = long_tuc_thresholds
    read_qc["site_evidence_mode"] = args.site_evidence_mode
    read_qc["bg_tc_rate"] = bg_rate
    read_qc["bg_rate_source"] = bg_rate_source
    read_qc["calling_parameters"] = vars(args)
    sample_summary = write_sample_summary(
        outdir,
        args.sample_id,
        read_qc,
        site_call_counts,
        bg_rate,
        bg_rate_source,
        long_tuc_thresholds,
        args.assignment_mode,
        args.library_strand_mode,
    )
    for key, value in sample_summary.items():
        read_qc[key] = value
    write_ntr_sufficient_stats(outdir, args.sample_id, gene_summary, gene_names, bg_rate, bg_rate_source)
    write_report_summaries(outdir, args.sample_id, read_qc, site_call_counts, bg_rate, bg_rate_source)
    write_grand_slam_compatibility(
        outdir,
        args.sample_id,
        args.bam,
        args.library_strand_mode,
        args.grand_slam_bin,
        args.grand_slam_java_bin,
    )
    qc_path = outdir / f"{args.sample_id}.bulk_slam_qc.json"
    qc_path.write_text(json.dumps(read_qc, indent=2), encoding="utf-8")


if __name__ == "__main__":
    main()
