#!/usr/bin/env python3
"""Summarize all high-quality primary-alignment base substitutions."""

from __future__ import annotations

import argparse
import csv
from collections import Counter
from pathlib import Path

import pysam


BASES = ("A", "C", "G", "T")


def arguments() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--sample-id", required=True)
    parser.add_argument("--bam", required=True, type=Path)
    parser.add_argument("--fasta", required=True, type=Path)
    parser.add_argument("--outdir", required=True, type=Path)
    parser.add_argument("--min-baseq", type=int, default=15)
    return parser.parse_args()


def main() -> int:
    args = arguments()
    counts: Counter[tuple[str, str]] = Counter()
    denominators: Counter[str] = Counter()
    with pysam.AlignmentFile(args.bam, "rb") as bam, pysam.FastaFile(args.fasta) as fasta:
        for read in bam.fetch(until_eof=True):
            if read.is_unmapped or read.is_secondary or read.is_supplementary or read.query_sequence is None or read.query_qualities is None:
                continue
            sequence = read.query_sequence.upper()
            qualities = read.query_qualities
            blocks = [(start, end, fasta.fetch(read.reference_name, start, end).upper()) for start, end in read.get_blocks()]
            for qpos, rpos in read.get_aligned_pairs(matches_only=False):
                if qpos is None or rpos is None or qualities[qpos] < args.min_baseq:
                    continue
                ref = next((seq[rpos - start] for start, end, seq in blocks if start <= rpos < end), "")
                query = sequence[qpos]
                if ref not in BASES or query not in BASES:
                    continue
                denominators[ref] += 1
                if query != ref:
                    counts[(ref, query)] += 1
    args.outdir.mkdir(parents=True, exist_ok=True)
    output = args.outdir / f"{args.sample_id}.base_mismatch_spectrum.tsv"
    with output.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(["reference_base", "observed_base", "count", "denominator", "probability", "scope"])
        for ref in BASES:
            for query in BASES:
                if query == ref:
                    continue
                total = denominators[ref]
                writer.writerow([ref, query, counts[(ref, query)], total, f"{counts[(ref, query)] / total if total else 0.0:.12g}", "primary_q20_aligned_bases"])
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
