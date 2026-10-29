#!/usr/bin/env python3
"""Regression tests for lr_bulk_slam runner helper functions."""

from __future__ import annotations

import importlib.util
import json
import tempfile
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
SPEC = importlib.util.spec_from_file_location("run_lr_bulk_slam", ROOT / "pipeline/scripts/lr_bulk_slam_helpers.py")
assert SPEC and SPEC.loader
run_lr_bulk_slam = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(run_lr_bulk_slam)


def test_resolve_broken_nanoplot_shebang() -> None:
    with tempfile.TemporaryDirectory() as tmp:
        bindir = Path(tmp)
        python = bindir / "python"
        nanoplot = bindir / "NanoPlot"
        python.write_text("#!/usr/bin/env sh\nexit 0\n", encoding="utf-8")
        nanoplot.write_text("#!/missing/python\n", encoding="utf-8")
        python.chmod(0o755)
        nanoplot.chmod(0o755)

        command, mode = run_lr_bulk_slam.resolve_tool_command(str(nanoplot))
        assert command == [str(python), "-m", "nanoplot.NanoPlot"]
        assert mode == "python_module_fallback"


def test_nanoplot_status_current() -> None:
    with tempfile.TemporaryDirectory() as tmp:
        tmpdir = Path(tmp)
        stats = tmpdir / "sample.raw.NanoStats.txt"
        report = tmpdir / "sample.raw.NanoPlot-report.html"
        status = tmpdir / "nanoplot_status.json"
        stats.write_text("Metric\tValue\n", encoding="utf-8")
        report.write_text("<html></html>\n", encoding="utf-8")
        payload = {
            "status": "complete",
            "fastq": "reads.fastq.gz",
            "command": ["python", "-m", "nanoplot.NanoPlot"],
            "extra_args": ["--no_static", "--tsv_stats"],
            "nanostats": str(stats),
            "report": str(report),
        }
        status.write_text(json.dumps(payload), encoding="utf-8")

        assert run_lr_bulk_slam.nanoplot_is_current(
            status,
            "reads.fastq.gz",
            ["python", "-m", "nanoplot.NanoPlot"],
            ["--no_static", "--tsv_stats"],
        )
        assert not run_lr_bulk_slam.nanoplot_is_current(
            status,
            "other.fastq.gz",
            ["python", "-m", "nanoplot.NanoPlot"],
            ["--no_static", "--tsv_stats"],
        )


def test_alignment_is_current_requires_nonempty_auxiliary_files() -> None:
    with tempfile.TemporaryDirectory() as tmp:
        outdir = Path(tmp)
        sorted_bam = outdir / "BC80.sorted.bam"
        primary_bam = outdir / "BC80.primary.q20.bam"
        required = [
            sorted_bam,
            outdir / "BC80.sorted.bam.bai",
            outdir / "BC80.flagstat.txt",
            primary_bam,
            outdir / "BC80.primary.q20.bam.bai",
            outdir / "BC80.primary.flagstat.txt",
        ]
        for path in required:
            path.write_text("ok\n", encoding="utf-8")

        assert run_lr_bulk_slam.alignment_is_current(sorted_bam, primary_bam, "BC80", 20)

        (outdir / "BC80.primary.flagstat.txt").write_text("", encoding="utf-8")
        assert not run_lr_bulk_slam.alignment_is_current(sorted_bam, primary_bam, "BC80", 20)


def test_link_full_reusable_stages() -> None:
    with tempfile.TemporaryDirectory() as tmp:
        root = Path(tmp)
        source_sample = root / "old" / "S1"
        qc = source_sample / "01_raw_qc"
        trim = source_sample / "02_trim" / "04_tool_output" / "trim"
        align = source_sample / "03_alignment" / "04_tool_output" / "minimap2_samtools"
        qc.mkdir(parents=True)
        trim.mkdir(parents=True)
        align.mkdir(parents=True)
        (qc / "source.txt").write_text("old\n", encoding="utf-8")
        qc_fastq = qc / "04_tool_output" / "longread_qc" / "fastq"
        qc_fastq.mkdir(parents=True)
        (qc_fastq / "S1.raw_fastq_qc.json").write_text("{}\n", encoding="utf-8")
        (qc_fastq / "S1.insert_fastq_qc.json").write_text("{}\n", encoding="utf-8")
        for path in (
            trim / "S1.insert.fastq.gz",
            trim / "S1.trim_summary.json",
            align / "S1.sorted.bam",
            align / "S1.sorted.bam.bai",
            align / "S1.flagstat.txt",
            align / "S1.primary.q20.bam",
            align / "S1.primary.q20.bam.bai",
            align / "S1.primary.flagstat.txt",
        ):
            path.write_text("ok\n", encoding="utf-8")

        target_sample = root / "new" / "S1"
        records = run_lr_bulk_slam.link_full_reusable_stages(target_sample, root / "old", "S1", 20)
        assert [record["stage"] for record in records] == ["01_raw_qc", "02_trim", "03_alignment"]
        assert all((target_sample / stage).is_symlink() for stage in ("01_raw_qc", "02_trim", "03_alignment"))
        assert (target_sample / "03_alignment/04_tool_output/minimap2_samtools/S1.primary.q20.bam").read_text() == "ok\n"


def main() -> None:
    test_resolve_broken_nanoplot_shebang()
    test_nanoplot_status_current()
    test_alignment_is_current_requires_nonempty_auxiliary_files()
    test_link_full_reusable_stages()


if __name__ == "__main__":
    main()
