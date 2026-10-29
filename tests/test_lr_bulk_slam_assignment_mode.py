#!/usr/bin/env python3
"""Synthetic test for transcript-aware bulk SLAM gene assignment."""

from __future__ import annotations

import csv
import gzip
import subprocess
import sys
import tempfile
from pathlib import Path

import pysam


ROOT = Path(__file__).resolve().parents[1]


def write_fasta(path: Path) -> None:
    path.write_text(">chr1\n" + "T" * 500 + "\n", encoding="utf-8")
    pysam.faidx(str(path))


def write_gtf(path: Path) -> None:
    path.write_text(
        "\n".join(
            [
                'chr1\ttest\tgene\t1\t140\t.\t+\t.\tgene_id "GENE1"; gene_name "Gene1";',
                'chr1\ttest\ttranscript\t1\t120\t.\t+\t.\tgene_id "GENE1"; transcript_id "TX1"; gene_name "Gene1";',
                'chr1\ttest\texon\t1\t40\t.\t+\t.\tgene_id "GENE1"; transcript_id "TX1"; gene_name "Gene1";',
                'chr1\ttest\texon\t81\t120\t.\t+\t.\tgene_id "GENE1"; transcript_id "TX1"; gene_name "Gene1";',
                'chr1\ttest\ttranscript\t1\t115\t.\t+\t.\tgene_id "GENE1"; transcript_id "TX1B"; gene_name "Gene1";',
                'chr1\ttest\texon\t1\t35\t.\t+\t.\tgene_id "GENE1"; transcript_id "TX1B"; gene_name "Gene1";',
                'chr1\ttest\texon\t81\t115\t.\t+\t.\tgene_id "GENE1"; transcript_id "TX1B"; gene_name "Gene1";',
                'chr1\ttest\ttranscript\t1\t115\t.\t+\t.\tgene_id "GENE1"; transcript_id "TX1C"; gene_name "Gene1";',
                'chr1\ttest\texon\t1\t35\t.\t+\t.\tgene_id "GENE1"; transcript_id "TX1C"; gene_name "Gene1";',
                'chr1\ttest\texon\t81\t115\t.\t+\t.\tgene_id "GENE1"; transcript_id "TX1C"; gene_name "Gene1";',
                'chr1\ttest\tgene\t141\t265\t.\t+\t.\tgene_id "GENE2"; gene_name "Gene2";',
                'chr1\ttest\ttranscript\t141\t260\t.\t+\t.\tgene_id "GENE2"; transcript_id "TX2"; gene_name "Gene2";',
                'chr1\ttest\texon\t141\t180\t.\t+\t.\tgene_id "GENE2"; transcript_id "TX2"; gene_name "Gene2";',
                'chr1\ttest\texon\t221\t260\t.\t+\t.\tgene_id "GENE2"; transcript_id "TX2"; gene_name "Gene2";',
                'chr1\ttest\tgene\t141\t265\t.\t+\t.\tgene_id "GENE3"; gene_name "Gene3";',
                'chr1\ttest\ttranscript\t141\t265\t.\t+\t.\tgene_id "GENE3"; transcript_id "TX3"; gene_name "Gene3";',
                'chr1\ttest\texon\t141\t180\t.\t+\t.\tgene_id "GENE3"; transcript_id "TX3"; gene_name "Gene3";',
                'chr1\ttest\texon\t221\t265\t.\t+\t.\tgene_id "GENE3"; transcript_id "TX3"; gene_name "Gene3";',
                'chr1\ttest\tgene\t281\t400\t.\t+\t.\tgene_id "GENE4"; gene_name "Gene4";',
                'chr1\ttest\ttranscript\t281\t400\t.\t+\t.\tgene_id "GENE4"; transcript_id "TX4"; gene_name "Gene4";',
                'chr1\ttest\texon\t281\t320\t.\t+\t.\tgene_id "GENE4"; transcript_id "TX4"; gene_name "Gene4";',
                'chr1\ttest\texon\t361\t400\t.\t+\t.\tgene_id "GENE4"; transcript_id "TX4"; gene_name "Gene4";',
                'chr1\ttest\tgene\t281\t400\t.\t+\t.\tgene_id "GENE5"; gene_name "Gene5";',
                'chr1\ttest\ttranscript\t281\t400\t.\t+\t.\tgene_id "GENE5"; transcript_id "TX5"; gene_name "Gene5";',
                'chr1\ttest\texon\t281\t320\t.\t+\t.\tgene_id "GENE5"; transcript_id "TX5"; gene_name "Gene5";',
                'chr1\ttest\texon\t361\t400\t.\t+\t.\tgene_id "GENE5"; transcript_id "TX5"; gene_name "Gene5";',
                "",
            ]
        ),
        encoding="utf-8",
    )


def make_read(name: str, start: int, cigar: list[tuple[int, int]], query_len: int) -> pysam.AlignedSegment:
    read = pysam.AlignedSegment()
    read.query_name = name
    read.query_sequence = "T" * query_len
    read.flag = 0
    read.reference_id = 0
    read.reference_start = start
    read.mapping_quality = 60
    read.cigar = cigar
    read.query_qualities = [35] * query_len
    return read


def write_bam(path: Path) -> None:
    header = {"HD": {"VN": "1.6"}, "SQ": [{"SN": "chr1", "LN": 500}]}
    with pysam.AlignmentFile(path, "wb", header=header) as bam:
        bam.write(make_read("tx_unique", 0, [(0, 40), (3, 40), (0, 40)], 80))
        bam.write(make_read("gene_multi", 0, [(0, 35), (3, 45), (0, 35)], 70))
        bam.write(make_read("overlap_unresolved", 29, [(0, 31)], 31))
        bam.write(make_read("shared_resolved", 140, [(0, 40), (3, 40), (0, 40)], 80))
        bam.write(make_read("intergenic", 269, [(0, 10)], 10))
        bam.write(make_read("shared_ambiguous", 280, [(0, 40), (3, 40), (0, 40)], 80))
    pysam.index(str(path))


def read_rows(path: Path) -> list[dict[str, str]]:
    opener = gzip.open if path.suffix == ".gz" else open
    with opener(path, "rt", encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def main() -> None:
    with tempfile.TemporaryDirectory(prefix="lr_bulk_slam_assignment_") as tmp_raw:
        tmp = Path(tmp_raw)
        fasta = tmp / "ref.fa"
        gtf = tmp / "genes.gtf"
        bam = tmp / "reads.bam"
        outdir = tmp / "out"
        write_fasta(fasta)
        write_gtf(gtf)
        write_bam(bam)
        subprocess.run(
            [
                sys.executable,
                str(ROOT / "pipeline/scripts/lr_call_bulk_slam.py"),
                "--bam",
                str(bam),
                "--fasta",
                str(fasta),
                "--gtf",
                str(gtf),
                "--sample-id",
                "synthetic",
                "--outdir",
                str(outdir),
                "--assignment-mode",
                "transcript_aware_gene",
                "--assignment-min-exon-overlap-bp",
                "20",
                "--assignment-require-same-strand",
                "--assignment-end-distance-tiebreak",
                "--assignment-allow-gene-level-salvage",
                "--assignment-count-shared-resolved-reads",
                "--min-t-covered",
                "1",
                "--min-site-depth",
                "1",
                "--min-site-tc-count",
                "1",
                "--min-site-tc-rate",
                "0.1",
                "--max-site-other-mismatch-rate",
                "0.5",
            ],
            check=True,
            cwd=ROOT,
        )
        read_calls = {row["read_id"]: row for row in read_rows(outdir / "synthetic.bulk_read_calls.tsv.gz")}
        summary = read_rows(outdir / "synthetic.bulk_slam_summary.tsv")[0]
        assert read_calls["tx_unique"]["assignment_class"] == "tx_unique"
        assert read_calls["gene_multi"]["assignment_class"] == "gene_unique_multi_tx"
        assert read_calls["shared_resolved"]["assignment_class"] == "shared_gene_resolved_by_ends"
        assert read_calls["shared_ambiguous"]["assignment_class"] == "shared_gene_ambiguous"
        assert read_calls["overlap_unresolved"]["assignment_class"] == "gene_overlap_unresolved"
        assert read_calls["intergenic"]["assignment_class"] == "intergenic"
        assert summary["tx_unique_reads"] == "1"
        assert summary["gene_unique_multi_tx_reads"] == "1"
        assert summary["shared_gene_resolved_reads"] == "1"
        assert summary["shared_gene_ambiguous_reads"] == "1"
        assert summary["gene_overlap_unresolved_reads"] == "1"
        assert summary["gene_level_salvaged_reads"] == "1"


if __name__ == "__main__":
    main()
