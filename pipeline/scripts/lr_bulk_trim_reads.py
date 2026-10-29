#!/usr/bin/env python3
"""Extract inserts from bulk long-read SLAM-seq FASTQ."""

from __future__ import annotations

import argparse
import gzip
import json
import math
from collections import Counter
from pathlib import Path


RC_TABLE = str.maketrans("ACGTNacgtn", "TGCANtgcan")


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Trim bulk long-read SLAM-seq reads.")
    parser.add_argument("--fastq", required=True)
    parser.add_argument("--sample-id", required=True)
    parser.add_argument("--outdir", required=True)
    parser.add_argument("--primer-5p", required=True)
    parser.add_argument("--adapter-3p", required=True)
    parser.add_argument(
        "--trim-mode",
        choices=("primer_adapter", "ont_cdna_polyt"),
        default="primer_adapter",
    )
    parser.add_argument("--poly-t-min", type=int, default=12)
    parser.add_argument("--max-reads", type=int, default=0)
    parser.add_argument("--min-insert-len", type=int, default=50)
    parser.add_argument("--seed-len", type=int, default=18)
    parser.add_argument("--search-window", type=int, default=160)
    parser.add_argument("--adapter-max-error-rate", type=float, default=1.0)
    parser.add_argument("--adapter-min-overlap", type=int, default=12)
    parser.add_argument("--poly-t-max-non-t", type=int, default=0)
    return parser.parse_args()


def rc(text: str) -> str:
    return text.translate(RC_TABLE)[::-1]


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
            plus = handle.readline()
            qual = handle.readline().strip()
            if not qual:
                break
            yield header.rstrip("\n")[1:].split()[0], seq, qual


def seeds(adapter: str, seed_len: int) -> list[tuple[str, int]]:
    if len(adapter) <= seed_len:
        return [(adapter, 0)]
    starts = sorted({0, max(0, len(adapter) // 2 - seed_len // 2), len(adapter) - seed_len})
    return [(adapter[start : start + seed_len], start) for start in starts]


def adapter_error_rate(seq: str, adapter: str, candidate_start: int) -> tuple[int, float]:
    adapter_start = max(0, -candidate_start)
    seq_start = max(0, candidate_start)
    overlap = min(len(adapter) - adapter_start, len(seq) - seq_start)
    if overlap <= 0:
        return 0, 1.0
    mismatches = 0
    for idx in range(overlap):
        mismatches += int(seq[seq_start + idx] != adapter[adapter_start + idx])
    return overlap, mismatches / overlap


def find_adapter(
    seq: str,
    adapter: str,
    seed_len: int,
    start: int,
    end: int,
    max_error_rate: float = 1.0,
    min_overlap: int = 12,
    prefer_right: bool = False,
) -> int:
    seq_region = seq[start:end]
    if max_error_rate >= 1.0:
        raw_hits = []
        for seed, offset in seeds(adapter, seed_len):
            pos = seq_region.rfind(seed) if prefer_right else seq_region.find(seed)
            if pos >= 0:
                raw_hits.append(start + pos - offset)
        if not raw_hits:
            return -1
        return max(0, max(raw_hits) if prefer_right else min(raw_hits))
    hits: list[tuple[float, int]] = []
    for seed, offset in seeds(adapter, seed_len):
        search_from = 0
        checked = 0
        while True:
            pos = seq_region.find(seed, search_from)
            if pos < 0:
                break
            candidate = start + pos - offset
            overlap, error_rate = adapter_error_rate(seq, adapter, candidate)
            if overlap >= min_overlap and error_rate <= max_error_rate:
                hits.append((error_rate, candidate))
            search_from = pos + 1
            checked += 1
            if checked >= 80:
                break
    if not hits:
        return -1
    best_error = min(error_rate for error_rate, _candidate in hits)
    best = [candidate for error_rate, candidate in hits if error_rate == best_error]
    return max(0, max(best) if prefer_right else min(best))


def find_left_adapter(
    seq: str,
    adapter: str,
    seed_len: int,
    search_window: int,
    max_error_rate: float,
    min_overlap: int,
) -> int:
    return find_adapter(seq, adapter, seed_len, 0, min(len(seq), search_window), max_error_rate, min_overlap)


def find_right_adapter(
    seq: str,
    adapter: str,
    seed_len: int,
    search_window: int,
    max_error_rate: float,
    min_overlap: int,
) -> int:
    start = max(0, len(seq) - search_window)
    return find_adapter(seq, adapter, seed_len, start, len(seq), max_error_rate, min_overlap, prefer_right=True)


def trim_read(
    seq: str,
    qual: str,
    primer_5p: str,
    adapter_3p: str,
    seed_len: int,
    search_window: int,
    min_insert_len: int,
    max_error_rate: float,
    min_overlap: int,
):
    forward_start = find_left_adapter(seq, primer_5p, seed_len, search_window, max_error_rate, min_overlap)
    reverse_end = find_right_adapter(seq, adapter_3p, seed_len, search_window, max_error_rate, min_overlap)
    if forward_start >= 0 and reverse_end > forward_start:
        insert_start = forward_start + len(primer_5p)
        insert_end = reverse_end
        if insert_end - insert_start >= min_insert_len:
            return "forward", seq[insert_start:insert_end], qual[insert_start:insert_end]

    rc_adapter = rc(adapter_3p)
    rc_primer = rc(primer_5p)
    reverse_start = find_left_adapter(seq, rc_adapter, seed_len, search_window, max_error_rate, min_overlap)
    forward_end = find_right_adapter(seq, rc_primer, seed_len, search_window, max_error_rate, min_overlap)
    if reverse_start >= 0 and forward_end > reverse_start:
        insert_start = reverse_start + len(rc_adapter)
        insert_end = forward_end
        if insert_end - insert_start >= min_insert_len:
            return "reverse", rc(seq[insert_start:insert_end]), qual[insert_start:insert_end][::-1]

    return "unassigned", "", ""


def find_homopolymer(
    seq: str,
    base: str,
    min_len: int,
    start: int,
    end: int,
    max_non_base: int = 0,
) -> tuple[int, int]:
    if max_non_base > 0:
        limit = min(end, len(seq))
        window_start = start
        non_base = 0
        base_count = 0
        for pos in range(start, limit):
            if seq[pos] == base:
                base_count += 1
            else:
                non_base += 1
            while non_base > max_non_base and window_start <= pos:
                if seq[window_start] == base:
                    base_count -= 1
                else:
                    non_base -= 1
                window_start += 1
            if base_count >= min_len:
                return window_start, pos + 1
        return -1, -1
    run_start = -1
    run_len = 0
    for pos in range(start, min(end, len(seq))):
        if seq[pos] == base:
            if run_start < 0:
                run_start = pos
            run_len += 1
            continue
        if run_len >= min_len:
            return run_start, pos
        run_start = -1
        run_len = 0
    if run_len >= min_len:
        return run_start, min(end, len(seq))
    return -1, -1


def trim_ont_cdna_polyt(
    seq: str,
    qual: str,
    left_adapter: str,
    right_adapter: str,
    seed_len: int,
    search_window: int,
    min_insert_len: int,
    poly_t_min: int,
    max_error_rate: float,
    min_overlap: int,
    poly_t_max_non_t: int,
):
    def extract_after_polyt(oriented_seq: str, oriented_qual: str):
        poly_start, poly_end = find_homopolymer(
            oriented_seq,
            "T",
            poly_t_min,
            0,
            min(len(oriented_seq), search_window),
            poly_t_max_non_t,
        )
        if poly_start < 0:
            return "", ""
        insert_seq = oriented_seq[poly_end:]
        insert_qual = oriented_qual[poly_end:]
        if len(insert_seq) < min_insert_len:
            return "", ""
        return insert_seq, insert_qual

    left_start = find_left_adapter(seq, left_adapter, seed_len, search_window, max_error_rate, min_overlap)
    right_end = find_right_adapter(seq, right_adapter, seed_len, search_window, max_error_rate, min_overlap)
    if left_start >= 0 and right_end > left_start:
        core_start = left_start + len(left_adapter)
        core_end = right_end
        insert_seq, insert_qual = extract_after_polyt(seq[core_start:core_end], qual[core_start:core_end])
        if insert_seq:
            return "forward", insert_seq, insert_qual

    rc_right = rc(right_adapter)
    rc_left = rc(left_adapter)
    reverse_start = find_left_adapter(seq, rc_right, seed_len, search_window, max_error_rate, min_overlap)
    forward_end = find_right_adapter(seq, rc_left, seed_len, search_window, max_error_rate, min_overlap)
    if reverse_start >= 0 and forward_end > reverse_start:
        core_start = reverse_start + len(rc_right)
        core_end = forward_end
        oriented_seq = rc(seq[core_start:core_end])
        oriented_qual = qual[core_start:core_end][::-1]
        insert_seq, insert_qual = extract_after_polyt(oriented_seq, oriented_qual)
        if insert_seq:
            return "reverse", insert_seq, insert_qual

    return "unassigned", "", ""


def main() -> None:
    args = parse_args()
    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    out_fastq = outdir / f"{args.sample_id}.insert.fastq.gz"
    summary_path = outdir / f"{args.sample_id}.trim_summary.json"
    length_path = outdir / f"{args.sample_id}.insert_lengths.tsv"
    counts = Counter()
    assigned_count = 0
    insert_length_sum = 0
    insert_length_min = math.inf
    insert_length_max = 0

    with gzip.open(out_fastq, "wt") as out_handle, length_path.open("w", encoding="utf-8") as length_handle:
        length_handle.write("read_id\torientation\tinsert_len\n")
        for read_number, (read_id, seq, qual) in enumerate(iter_fastq(args.fastq), start=1):
            if args.max_reads and read_number > args.max_reads:
                break
            counts["total_reads"] += 1
            upper_seq = seq.upper()
            if args.trim_mode == "ont_cdna_polyt":
                orientation, insert_seq, insert_qual = trim_ont_cdna_polyt(
                    upper_seq,
                    qual,
                    args.primer_5p.upper(),
                    args.adapter_3p.upper(),
                    args.seed_len,
                    args.search_window,
                    args.min_insert_len,
                    args.poly_t_min,
                    args.adapter_max_error_rate,
                    args.adapter_min_overlap,
                    args.poly_t_max_non_t,
                )
            else:
                orientation, insert_seq, insert_qual = trim_read(
                    upper_seq,
                    qual,
                    args.primer_5p.upper(),
                    args.adapter_3p.upper(),
                    args.seed_len,
                    args.search_window,
                    args.min_insert_len,
                    args.adapter_max_error_rate,
                    args.adapter_min_overlap,
                )
            counts[orientation] += 1
            if orientation == "unassigned":
                continue
            insert_len = len(insert_seq)
            assigned_count += 1
            insert_length_sum += insert_len
            insert_length_min = min(insert_length_min, insert_len)
            insert_length_max = max(insert_length_max, insert_len)
            length_handle.write(f"{read_id}\t{orientation}\t{len(insert_seq)}\n")
            out_handle.write(f"@{read_id}\n{insert_seq}\n+\n{insert_qual}\n")

    assigned = counts["forward"] + counts["reverse"]
    total = counts["total_reads"]
    mean_insert_len = insert_length_sum / assigned_count if assigned_count else 0.0
    if insert_length_min == math.inf:
        insert_length_min = 0
    summary = {
        "total_reads": total,
        "trimmed_reads": assigned,
        "trimmed_rate": assigned / total if total else 0.0,
        "forward_reads": counts["forward"],
        "reverse_reads": counts["reverse"],
        "unassigned_reads": counts["unassigned"],
        "mean_insert_len": mean_insert_len,
        "min_insert_len": insert_length_min,
        "max_insert_len": insert_length_max,
    }
    summary_path.write_text(json.dumps(summary, indent=2), encoding="utf-8")


if __name__ == "__main__":
    main()
