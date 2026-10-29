#!/usr/bin/env python3
"""Compute long-read FASTQ QC metrics."""

from __future__ import annotations

import argparse
import csv
import gzip
import json
import math
import statistics
from collections import Counter
from pathlib import Path


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Compute long-read FASTQ QC metrics.")
    parser.add_argument("--fastq", required=True)
    parser.add_argument("--sample-id", required=True)
    parser.add_argument("--outdir", required=True)
    parser.add_argument("--label", choices=("raw", "insert"), default="raw")
    parser.add_argument("--max-reads", type=int, default=0)
    parser.add_argument("--write-read-metrics", action="store_true")
    return parser.parse_args()


def open_text(path: str):
    if path.endswith(".gz"):
        return gzip.open(path, "rt")
    return open(path, "r", encoding="utf-8")


def iter_fastq(path: str):
    with open_text(path) as handle:
        while True:
            header = handle.readline()
            if not header:
                break
            seq = handle.readline().strip()
            _plus = handle.readline()
            qual = handle.readline().strip()
            if not qual:
                break
            yield header.rstrip("\n")[1:].split()[0], seq, qual


def n50(length_counts: Counter, total_bases: int) -> int:
    if not length_counts or total_bases <= 0:
        return 0
    half = total_bases / 2
    running = 0
    for length in sorted(length_counts, reverse=True):
        running += length * length_counts[length]
        if running >= half:
            return length
    return 0


def mean(total: float, count: int) -> float:
    return total / count if count else 0.0


def median_counted(counts: Counter, total_count: int) -> float:
    if not counts or total_count <= 0:
        return 0.0
    target = (total_count + 1) / 2
    running = 0
    for value in sorted(counts):
        running += counts[value]
        if running >= target:
            return float(value)
    return 0.0


def length_bins(counts: Counter) -> dict[str, int]:
    bins = {
        "lt_200": 0,
        "200_499": 0,
        "500_999": 0,
        "1000_1499": 0,
        "1500_1999": 0,
        "ge_2000": 0,
    }
    for length, count in counts.items():
        if length < 200:
            bins["lt_200"] += count
        elif length < 500:
            bins["200_499"] += count
        elif length < 1000:
            bins["500_999"] += count
        elif length < 1500:
            bins["1000_1499"] += count
        elif length < 2000:
            bins["1500_1999"] += count
        else:
            bins["ge_2000"] += count
    return bins


def q_bins_for_read(q: float) -> str:
    if q < 7:
        return "lt_7"
    if q < 10:
        return "7_9"
    if q < 15:
        return "10_14"
    if q < 20:
        return "15_19"
    if q < 30:
        return "20_29"
    return "ge_30"


def gc_bin_for_read(gc: float) -> str:
    if gc < 0.2:
        return "0_20"
    if gc < 0.4:
        return "20_40"
    if gc < 0.6:
        return "40_60"
    if gc < 0.8:
        return "60_80"
    return "80_100"


def main() -> None:
    args = parse_args()
    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    metrics_path = outdir / f"{args.sample_id}.{args.label}_read_metrics.tsv.gz"
    summary_path = outdir / f"{args.sample_id}.{args.label}_fastq_qc.json"

    read_count = 0
    base_count = 0
    q_counts = Counter()
    length_counts = Counter()
    mean_q_bins = Counter()
    gc_bin_counts = Counter()
    n_count = 0
    gc_count = 0
    total_q = 0.0
    length_sum = 0

    metrics_handle = gzip.open(metrics_path, "wt", encoding="utf-8", newline="") if args.write_read_metrics else None
    try:
        writer = csv.writer(metrics_handle, delimiter="\t", lineterminator="\n") if metrics_handle else None
        if writer:
            writer.writerow(["read_id", "length", "mean_q", "median_q", "min_q", "max_q", "gc_fraction", "n_count"])
        for read_number, (read_id, seq, qual) in enumerate(iter_fastq(args.fastq), start=1):
            if args.max_reads and read_number > args.max_reads:
                break
            q_values = [ord(char) - 33 for char in qual]
            seq_upper = seq.upper()
            length = len(seq_upper)
            read_gc = seq_upper.count("G") + seq_upper.count("C")
            read_n = seq_upper.count("N")
            read_count += 1
            base_count += length
            gc_count += read_gc
            n_count += read_n
            length_sum += length
            length_counts[length] += 1

            read_q_total = 0.0
            read_min_q = math.inf
            read_max_q = -math.inf
            for q in q_values:
                q_counts[q] += 1
                total_q += q
                read_q_total += q
                read_min_q = min(read_min_q, q)
                read_max_q = max(read_max_q, q)
            read_mean_q = read_q_total / len(q_values) if q_values else 0.0
            read_median_q = statistics.median(q_values) if q_values else 0.0
            mean_q_bins[q_bins_for_read(read_mean_q)] += 1
            gc_fraction = read_gc / length if length else 0.0
            gc_bin_counts[gc_bin_for_read(gc_fraction)] += 1
            if writer:
                writer.writerow(
                    [
                        read_id,
                        length,
                        f"{read_mean_q:.4f}",
                        f"{read_median_q:.4f}",
                        read_min_q if read_min_q is not math.inf else 0,
                        read_max_q if read_max_q is not -math.inf else 0,
                        f"{gc_fraction:.6f}",
                        read_n,
                    ]
                )
    finally:
        if metrics_handle:
            metrics_handle.close()

    q_total = sum(q_counts.values())
    summary = {
        "sample_id": args.sample_id,
        "label": args.label,
        "fastq": args.fastq,
        "max_reads": args.max_reads,
        "read_count": read_count,
        "base_count": base_count,
        "mean_length": mean(float(length_sum), read_count),
        "median_length": median_counted(length_counts, read_count),
        "min_length": min(length_counts) if read_count else 0,
        "max_length": max(length_counts) if read_count else 0,
        "n50": n50(length_counts, length_sum),
        "mean_read_q": total_q / q_total if q_total else 0.0,
        "median_read_q": median_counted(q_counts, q_total),
        "q10_base_fraction": sum(count for q, count in q_counts.items() if q >= 10) / q_total if q_total else 0.0,
        "q20_base_fraction": sum(count for q, count in q_counts.items() if q >= 20) / q_total if q_total else 0.0,
        "q30_base_fraction": sum(count for q, count in q_counts.items() if q >= 30) / q_total if q_total else 0.0,
        "gc_fraction": gc_count / base_count if base_count else 0.0,
        "n_fraction": n_count / base_count if base_count else 0.0,
        "length_bins": length_bins(length_counts),
        "mean_q_bins": {
            "lt_7": mean_q_bins.get("lt_7", 0),
            "7_9": mean_q_bins.get("7_9", 0),
            "10_14": mean_q_bins.get("10_14", 0),
            "15_19": mean_q_bins.get("15_19", 0),
            "20_29": mean_q_bins.get("20_29", 0),
            "ge_30": mean_q_bins.get("ge_30", 0),
        },
        "gc_bins": {
            "0_20": gc_bin_counts.get("0_20", 0),
            "20_40": gc_bin_counts.get("20_40", 0),
            "40_60": gc_bin_counts.get("40_60", 0),
            "60_80": gc_bin_counts.get("60_80", 0),
            "80_100": gc_bin_counts.get("80_100", 0),
        },
        "read_metrics": str(metrics_path) if args.write_read_metrics else "",
    }
    summary_path.write_text(json.dumps(summary, indent=2), encoding="utf-8")


if __name__ == "__main__":
    main()
