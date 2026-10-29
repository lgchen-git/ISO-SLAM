#!/usr/bin/env python3
"""Minimal synthetic test for formal bulk SLAM calling."""

from __future__ import annotations

import csv
import gzip
import json
import subprocess
import sys
import tempfile
from pathlib import Path

import pysam


ROOT = Path(__file__).resolve().parents[1]


def write_fasta(path: Path) -> None:
    path.write_text(">chr1\n" + "T" * 80 + "\n", encoding="utf-8")
    pysam.faidx(str(path))


def write_gtf(path: Path) -> None:
    path.write_text(
        'chr1\ttest\tgene\t1\t80\t.\t+\t.\tgene_id "GENE1"; gene_name "Gene1";\n'
        'chr1\ttest\texon\t1\t80\t.\t+\t.\tgene_id "GENE1"; gene_name "Gene1";\n',
        encoding="utf-8",
    )


def make_read(name: str, seq: str, mapq: int = 60, qualities: list[int] | None = None) -> pysam.AlignedSegment:
    read = pysam.AlignedSegment()
    read.query_name = name
    read.query_sequence = seq
    read.flag = 0
    read.reference_id = 0
    read.reference_start = 0
    read.mapping_quality = mapq
    read.cigar = [(0, len(seq))]
    read.query_qualities = qualities or [35] * len(seq)
    return read


def write_bam(path: Path) -> None:
    header = {"HD": {"VN": "1.6"}, "SQ": [{"SN": "chr1", "LN": 80}]}
    with pysam.AlignmentFile(path, "wb", header=header) as bam:
        bam.write(make_read("new_1", "CCCTTTTTTTTTTTTTTTTT"))
        bam.write(make_read("new_2", "CCCTTTTTTTTTTTTTTTTT"))
        bam.write(make_read("new_3", "CCCTTTTTTTTTTTTTTTTT"))
        bam.write(make_read("long_tuc_20", "C" * 20 + "T" * 5))
        bam.write(make_read("long_tuc_21", "C" * 21 + "T" * 4))
        bam.write(make_read("old_1", "TTTTTTTTTTTTTTTTTTTT"))
        bam.write(make_read("low_t", "C" * 5, qualities=[35] * 5))
        bam.write(make_read("low_baseq_tc", "CCCTTTTTTTTTTTTTTTTT", qualities=[5, 5, 5] + [35] * 17))
        bam.write(make_read("noisy", "AAGGTTTTTTTTTTTTTTTT"))
    pysam.index(str(path))


def read_rows(path: Path) -> list[dict[str, str]]:
    opener = gzip.open if path.suffix == ".gz" else open
    with opener(path, "rt", encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def main() -> None:
    with tempfile.TemporaryDirectory(prefix="lr_bulk_slam_formal_") as tmp_raw:
        tmp = Path(tmp_raw)
        fasta = tmp / "ref.fa"
        gtf = tmp / "genes.gtf"
        bam = tmp / "reads.bam"
        results_dir = tmp / "results"
        sample_root = results_dir / "synthetic"
        outdir = sample_root / "04_slam_calling" / "04_tool_output" / "caller"
        trim_dir = sample_root / "02_trim" / "04_tool_output" / "trim"
        align_dir = sample_root / "03_alignment" / "04_tool_output" / "minimap2_samtools"
        trim_dir.mkdir(parents=True)
        align_dir.mkdir(parents=True)
        (trim_dir / "synthetic.trim_summary.json").write_text(
            json.dumps({"total_reads": 9, "trimmed_reads": 9, "forward_reads": 9, "reverse_reads": 0}),
            encoding="utf-8",
        )
        (trim_dir / "synthetic.insert_lengths.tsv").write_text(
            "read_id\torientation\tinsert_length\nsynthetic\tforward\t20\n",
            encoding="utf-8",
        )
        (align_dir / "synthetic.primary.flagstat.txt").write_text("9 + 0 in total (QC-passed reads + QC-failed reads)\n", encoding="utf-8")
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
                "--min-t-covered",
                "6",
                "--min-site-depth",
                "3",
                "--min-site-tc-count",
                "2",
                "--min-site-tc-rate",
                "0.2",
                "--max-site-other-mismatch-rate",
                "0.2",
            ],
            check=True,
            cwd=ROOT,
        )
        read_calls = {row["read_id"]: row for row in read_rows(outdir / "synthetic.bulk_read_calls.tsv.gz")}
        site_rows = read_rows(outdir / "synthetic.bulk_t_sites.tsv.gz")
        summary = read_rows(outdir / "synthetic.bulk_slam_summary.tsv")[0]
        assert read_calls["new_1"]["final_call"] == "new"
        assert read_calls["old_1"]["final_call"] == "old"
        assert read_calls["low_t"]["final_call"] == "ambiguous"
        assert read_calls["low_baseq_tc"]["tc_count"] == "0"
        assert summary["bg_rate_source"] == "max_ta_tg_rate_global"
        assert any(row["site_call"] == "high_conf_tc_site" for row in site_rows)
        assert int(summary["new_read_count"]) >= 3
        config = tmp / "config.yaml"
        samples = tmp / "samples.tsv"
        config.write_text("main_call: site_binomial\nnew_threshold: 1\n", encoding="utf-8")
        samples.write_text("sample_id\nsynthetic\n", encoding="utf-8")
        subprocess.run(
            [
                sys.executable,
                str(ROOT / "pipeline/scripts/build_lr_bulk_slam_report.py"),
                "--results-dir",
                str(results_dir),
                "--sample-id",
                "synthetic",
                "--project-name",
                "synthetic",
                "--report-title",
                "Synthetic",
                "--config-path",
                str(config),
                "--samples-path",
                str(samples),
            ],
            check=True,
            cwd=ROOT,
        )
        report_path = sample_root / "05_sample_report/synthetic.bulk_slam_report.html"
        html = report_path.read_text(encoding="utf-8")
        assert "NanoPlot QC 明细" in html
        assert "long_tuc_gt20_reads" not in html
        subprocess.run(
            [
                sys.executable,
                str(ROOT / "pipeline/scripts/00_shared/validate_delivery_report.py"),
                "--report",
                str(report_path),
                "--require-section",
                "1. 样本概览",
            ],
            check=True,
            cwd=ROOT,
        )
        long_tuc_outdir = sample_root / "03_bulk_slam_long_tuc"
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
                str(long_tuc_outdir),
                "--main-call",
                "long_tuc_threshold",
                "--long-tuc-formal-threshold",
                "20",
                "--long-tuc-report-thresholds",
                "6,20,30",
                "--min-t-covered",
                "6",
                "--min-site-depth",
                "3",
                "--min-site-tc-count",
                "2",
                "--min-site-tc-rate",
                "0.2",
                "--max-site-other-mismatch-rate",
                "0.2",
            ],
            check=True,
            cwd=ROOT,
        )
        long_tuc_calls = {row["read_id"]: row for row in read_rows(long_tuc_outdir / "synthetic.bulk_read_calls.tsv.gz")}
        long_tuc_summary = read_rows(long_tuc_outdir / "synthetic.bulk_slam_summary.tsv")[0]
        long_tuc_genes = read_rows(long_tuc_outdir / "synthetic.bulk_gene_summary.tsv")
        assert long_tuc_calls["long_tuc_20"]["final_call"] == "old"
        assert long_tuc_calls["long_tuc_20"]["long_tuc_gt20"] == "0"
        assert long_tuc_calls["long_tuc_21"]["final_call"] == "new"
        assert long_tuc_calls["long_tuc_21"]["long_tuc_gt20"] == "1"
        assert long_tuc_summary["long_tuc_gt20_read_count"] == "1"
        assert "long_tuc_gt20_reads" in long_tuc_genes[0]


if __name__ == "__main__":
    main()
