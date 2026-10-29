#!/usr/bin/env python3
"""Regression tests for NanoPlot-driven bulk SLAM report helpers."""

from __future__ import annotations

import importlib.util
import json
import tempfile
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
SPEC = importlib.util.spec_from_file_location("build_lr_bulk_slam_report_v3", ROOT / "pipeline/scripts/build_lr_bulk_slam_report_v2.py")
assert SPEC and SPEC.loader
report = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(report)


def test_parse_nanostats() -> None:
    with tempfile.TemporaryDirectory() as tmp:
        path = Path(tmp) / "sample.insert.NanoStats.txt"
        path.write_text(
            "\n".join(
                [
                    "Metrics\tdataset",
                    "number_of_reads\t6209919",
                    "number_of_bases\t3933491723.0",
                    "median_read_length\t549.0",
                    "mean_read_length\t633.4",
                    "read_length_stdev\t342.2",
                    "n50\t701.0",
                    "mean_qual\t17.1",
                    "median_qual\t17.7",
                    "Reads >Q20:\t892507 (14.4%) 455.3Mb",
                    "Reads >Q30:\t1809 (0.0%) 0.4Mb",
                ]
            )
            + "\n",
            encoding="utf-8",
        )
        metrics = report.parse_nanostats(path)
        assert metrics["number_of_reads"] == 6209919
        assert metrics["n50"] == 701.0
        assert metrics["mean_read_length"] == 633.4
        assert report.threshold_metric(metrics, 20, "count") == 892507.0
        assert abs(report.threshold_metric(metrics, 20, "fraction") - 0.144) < 1e-9
        assert report.threshold_metric(metrics, 25, "count") == 0.0


def test_read_nanoplot_rows_failed_status_keeps_reason() -> None:
    with tempfile.TemporaryDirectory() as tmp:
        qc_dir = Path(tmp)
        insert_dir = qc_dir / "nanoplot" / "insert"
        insert_dir.mkdir(parents=True)
        status = {
            "status": "failed",
            "sample_id": "sampleA",
            "label": "insert",
            "version": "NanoPlot 1.46.2",
            "reason": "Fatal: No reads found in input.",
            "report": str(insert_dir / "sampleA.insert.NanoPlot-report.html"),
            "nanostats": str(insert_dir / "sampleA.insert.NanoStats.txt"),
        }
        (insert_dir / "nanoplot_status.json").write_text(json.dumps(status), encoding="utf-8")

        rows = report.read_nanoplot_rows(qc_dir, "sampleA")
        assert len(rows) == 1
        assert rows[0]["status"] == "failed"
        assert rows[0]["reason"] == "Fatal: No reads found in input."
        assert report.available_nanoplot_rows(rows) == []


def test_internal_fastq_qc_fills_table_when_nanoplot_skipped() -> None:
    with tempfile.TemporaryDirectory() as tmp:
        qc_dir = Path(tmp)
        internal = qc_dir / "fastq" / "sampleA.raw_fastq_qc.json"
        status_path = qc_dir / "nanoplot" / "raw" / "nanoplot_status.json"
        internal.parent.mkdir(parents=True)
        status_path.parent.mkdir(parents=True)
        internal.write_text(json.dumps({
            "read_count": 200000,
            "base_count": 272999168,
            "median_length": 1277,
            "n50": 1508,
            "median_read_q": 33,
            "q20_base_fraction": 0.8631646,
        }), encoding="utf-8")
        status_path.write_text(json.dumps({"status": "skipped", "reason": "timed out"}), encoding="utf-8")

        rows = report.read_nanoplot_rows(qc_dir, "sampleA")
        assert len(rows) == 1
        assert rows[0]["metric_source"] == "internal_qc"
        assert rows[0]["nanoplot_status"] == "skipped"
        table = report.nanoplot_table(rows)
        assert "200,000" in table
        assert "86.32% (bases)" in table


def main() -> None:
    test_parse_nanostats()
    test_read_nanoplot_rows_failed_status_keeps_reason()
    test_internal_fastq_qc_fills_table_when_nanoplot_skipped()


if __name__ == "__main__":
    main()
