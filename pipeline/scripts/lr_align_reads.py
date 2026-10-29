#!/usr/bin/env python3
"""Align long-read cDNA FASTQ with minimap2 and samtools."""

from __future__ import annotations

import argparse
import subprocess
from pathlib import Path


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Align cDNA reads.")
    parser.add_argument("--fastq", required=True)
    parser.add_argument("--reference", required=True)
    parser.add_argument("--outdir", required=True)
    parser.add_argument("--sample-id", required=True)
    parser.add_argument("--threads", type=int, default=4)
    parser.add_argument("--min-mapq", type=int, default=20)
    parser.add_argument("--minimap2-bin", default="minimap2")
    parser.add_argument("--samtools-bin", default="samtools")
    parser.add_argument("--junc-bed", default="")
    return parser.parse_args()


def require_file(path: str) -> None:
    if not Path(path).is_file():
        raise SystemExit(f"Required file not found: {path}")


def main() -> None:
    args = parse_args()
    require_file(args.fastq)
    require_file(args.reference)
    if args.junc_bed:
        require_file(args.junc_bed)
    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    sorted_bam = outdir / f"{args.sample_id}.sorted.bam"
    primary_bam = outdir / f"{args.sample_id}.primary.q{args.min_mapq}.bam"

    mm2_cmd = [
        args.minimap2_bin,
        "-ax",
        "splice",
        "-uf",
        "--secondary=no",
        "--MD",
        "--cs=long",
        "-t",
        str(max(1, args.threads - 3)),
    ]
    if args.junc_bed:
        mm2_cmd.extend(["--junc-bed", args.junc_bed])
    mm2_cmd.extend([args.reference, args.fastq])
    sort_cmd = [args.samtools_bin, "sort", "-@", "1", "-m", "2G", "-T", str(Path(__import__("os").environ["TMPDIR"]) / args.sample_id), "-o", str(sorted_bam), "-"]

    mm2 = subprocess.Popen(mm2_cmd, stdout=subprocess.PIPE)
    assert mm2.stdout is not None
    sort_proc = subprocess.Popen(sort_cmd, stdin=mm2.stdout)
    mm2.stdout.close()
    if sort_proc.wait() != 0:
        raise SystemExit("samtools sort failed")
    if mm2.wait() != 0:
        raise SystemExit("minimap2 failed")

    subprocess.run([args.samtools_bin, "index", "-@", str(max(0, args.threads - 1)), str(sorted_bam)], check=True)
    with (outdir / f"{args.sample_id}.flagstat.txt").open("w", encoding="utf-8") as handle:
        subprocess.run([args.samtools_bin, "flagstat", str(sorted_bam)], check=True, stdout=handle)
    with primary_bam.open("wb") as handle:
        subprocess.run(
            [args.samtools_bin, "view", "-b", "-q", str(args.min_mapq), "-F", "2304", str(sorted_bam)],
            check=True,
            stdout=handle,
        )
    subprocess.run([args.samtools_bin, "index", "-@", str(max(0, args.threads - 1)), str(primary_bam)], check=True)
    with (outdir / f"{args.sample_id}.primary.flagstat.txt").open("w", encoding="utf-8") as handle:
        subprocess.run([args.samtools_bin, "flagstat", str(primary_bam)], check=True, stdout=handle)


if __name__ == "__main__":
    main()
