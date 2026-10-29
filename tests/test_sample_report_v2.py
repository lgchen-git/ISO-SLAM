#!/usr/bin/env python3
from __future__ import annotations

import csv
import gzip
import json
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]


def write_flagstat(path: Path, total: int, primary: int, primary_mapped: int, supplementary: int = 0) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(
        "\n".join(
            [
                f"{total} + 0 in total (QC-passed reads + QC-failed reads)",
                f"{primary} + 0 primary",
                "0 + 0 secondary",
                f"{supplementary} + 0 supplementary",
                f"{primary_mapped + supplementary} + 0 mapped (0.00% : N/A)",
                f"{primary_mapped} + 0 primary mapped (0.00% : N/A)",
            ]
        )
        + "\n",
        encoding="utf-8",
    )


class SampleReportV2Tests(unittest.TestCase):
    def test_self_contained_paginated_report_with_read_flow(self) -> None:
        with tempfile.TemporaryDirectory(prefix="iso-slam-report-v2-") as raw:
            root = Path(raw)
            results = root / "results"
            sample = "S-1"
            sample_root = results / sample
            trim_dir = sample_root / "02_trim/04_tool_output/trim"
            align_dir = sample_root / "03_alignment/04_tool_output/minimap2_samtools"
            slam_dir = sample_root / "04_slam_calling/04_tool_output/caller"
            trim_dir.mkdir(parents=True)
            slam_dir.mkdir(parents=True)
            (trim_dir / f"{sample}.trim_summary.json").write_text(
                json.dumps(
                    {
                        "total_reads": 100,
                        "trimmed_reads": 80,
                        "unassigned_reads": 20,
                        "forward_reads": 80,
                        "reverse_reads": 0,
                    }
                ),
                encoding="utf-8",
            )
            write_flagstat(align_dir / f"{sample}.flagstat.txt", 85, 80, 70, 5)
            write_flagstat(align_dir / f"{sample}.primary.flagstat.txt", 60, 60, 60)
            qc = {
                "total_reads": 60,
                "assigned_reads": 50,
                "assigned_gene_reads": 50,
                "unassigned_gene": 10,
                "assignment_intergenic_reads": 4,
                "assignment_shared_antisense_ambiguous_reads": 3,
                "assignment_shared_gene_ambiguous_reads": 3,
                "assignment_antisense_only_reads": 0,
                "new_reads": 10,
                "old_reads": 40,
                "ambiguous_reads": 10,
                "genes_detected": 25,
                "t_covered": 1000,
                "t_to_a": 1,
                "t_to_c": 20,
                "t_to_g": 2,
                "global_tc_rate": 0.02,
                "bg_tc_rate": 0.002,
                "bg_rate_source": "synthetic",
                "main_call": "site_binomial",
            }
            (slam_dir / f"{sample}.bulk_slam_qc.json").write_text(json.dumps(qc), encoding="utf-8")
            with (slam_dir / f"{sample}.bulk_gene_summary.tsv").open("w", encoding="utf-8", newline="") as handle:
                writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
                writer.writerow(
                    ["gene_id", "gene_name", "total_reads", "new_reads", "old_reads", "ambiguous_reads", "new_fraction", "gene_reads", "tc_events", "tc_rate"]
                )
                for index in range(25):
                    writer.writerow([f"G{index}", f"Gene{index}", 10, 2, 7, 1, 0.2, 10, 3, 0.03])
            (slam_dir / f"{sample}.bulk_slam_summary.tsv").write_text("sample_id\tnew_reads\nS-1\t10\n", encoding="utf-8")
            (slam_dir / f"{sample}.assignment_summary.tsv").write_text(
                "level\tcategory\treads\tdenominator\tfraction\n"
                "report_group\ttranscript_supported\t30\t60\t0.5\n"
                "report_group\tsame_strand_gene_only\t20\t60\t0.333333\n"
                "report_group\tantisense_rescue\t0\t60\t0\n"
                "report_group\tambiguous_shared\t6\t60\t0.1\n"
                "report_group\tunassigned_intergenic\t4\t60\t0.066667\n",
                encoding="utf-8",
            )
            (slam_dir / f"{sample}.substitution_summary.tsv").write_text(
                "substitution\tcount\tdenominator\trate\trole\tsource\n"
                "T>A\t1\t1000\t0.001\tobserved\tall_assigned_T_bases\n"
                "T>C\t20\t1000\t0.02\tobserved\tall_assigned_T_bases\n"
                "T>G\t2\t1000\t0.002\tobserved\tall_assigned_T_bases\n"
                "diagnostic_background\t\t\t0.002\tdiagnostic_only\tsynthetic\n",
                encoding="utf-8",
            )
            (slam_dir / f"{sample}.ntr_sufficient_stats.tsv").write_text(
                "gene_id\tgene_name\tinformative_reads\tt_covered\ttc_count\tta_count\ttg_count\tdiagnostic_background_rate\tbackground_source\n"
                "G1\tGene1\t10\t100\t2\t0\t0\t0.002\tsynthetic\n",
                encoding="utf-8",
            )
            (slam_dir / f"{sample}.grand_slam_compatibility.json").write_text(
                json.dumps({"schema_version": "iso-slam-grand-slam-compatibility-v1", "status": "blocked"}),
                encoding="utf-8",
            )
            with gzip.open(slam_dir / f"{sample}.bulk_read_calls.tsv.gz", "wt", encoding="utf-8", newline="") as handle:
                writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
                writer.writerow(["read_id", "read_length", "tc_events", "t_covered", "other_mismatch_rate"])
                writer.writerow(["r1", 1000, 2, 20, 0.01])
            with gzip.open(slam_dir / f"{sample}.bulk_t_sites.tsv.gz", "wt", encoding="utf-8", newline="") as handle:
                writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
                writer.writerow(["tc_rate", "other_mismatch_rate", "site_call"])
                writer.writerow([0.2, 0.01, "high_conf_tc_site"])
            config = root / "config.yaml"
            samples = root / "samples.tsv"
            config.write_text("main_call: site_binomial\nassignment_mode: transcript_aware_gene\n", encoding="utf-8")
            samples.write_text("sample_id\tgroup\nS-1\tA\n", encoding="utf-8")
            flow = sample_root / "05_sample_report/04_tool_output/read_flow"
            subprocess.run(
                [
                    sys.executable,
                    str(ROOT / "pipeline/scripts/build_lr_bulk_read_flow.py"),
                    "--sample-id", sample,
                    "--sample-root", str(sample_root),
                    "--output-dir", str(flow),
                ],
                cwd=ROOT,
                check=True,
            )
            legacy_report = sample_root / "05_sample_report/report.html"
            legacy_report.write_text("legacy unprefixed report", encoding="utf-8")
            subprocess.run(
                [
                    sys.executable,
                    str(ROOT / "pipeline/scripts/build_lr_bulk_slam_report_v2.py"),
                    "--results-dir", str(results),
                    "--sample-id", sample,
                    "--project-name", "Synthetic",
                    "--report-title", "Synthetic",
                    "--config-path", str(config),
                    "--samples-path", str(samples),
                    "--read-flow-dir", str(flow),
                ],
                cwd=ROOT,
                check=True,
            )
            stage = sample_root / "05_sample_report"
            report = stage / f"{sample}.report.html"
            manifest = stage / "05_metadata/delivery_manifest.tsv"
            checked = subprocess.run(
                [
                    sys.executable,
                    str(ROOT / "pipeline/scripts/00_shared/validate_delivery_report_v2.py"),
                    "--report", str(report),
                    "--manifest", str(manifest),
                    "--require-section", "1. 样本概览",
                    "--require-section", "2. Reads处理流与长读长质量",
                    "--forbid-text", "tab-panel",
                    "--forbid-text", "数据保留漏斗",
                    "--json",
                ],
                cwd=ROOT,
                text=True,
                capture_output=True,
                check=False,
            )
            self.assertEqual(checked.returncode, 0, checked.stdout + checked.stderr)
            validation = json.loads(checked.stdout)
            self.assertTrue(validation["ok"])
            self.assertEqual(validation["external_hrefs"], 0)
            self.assertEqual(validation["external_srcs"], 0)
            text = report.read_text(encoding="utf-8")
            self.assertIn("Plotly.react", text)
            self.assertIn('class="sidebar"', text)
            self.assertIn('id="genePageInfo"', text)
            self.assertEqual(text.count("class='kpi'"), 6)
            self.assertIn("T&gt;C signal versus diagnostic background", text)
            self.assertIn("Transcript-aware gene assignment", text)
            self.assertNotIn("NanoPlot reads / bases summary", text)
            self.assertNotIn("read length vs tc_count", text)
            self.assertNotIn("tab-panel", text)
            with manifest.open(encoding="utf-8", newline="") as handle:
                rows = list(csv.DictReader(handle, delimiter="\t"))
            self.assertEqual(rows[-1]["artifact"], "final_sample_report")
            self.assertEqual(sum(Path(row["target"]).resolve() == report.resolve() for row in rows), 1)
            self.assertFalse(legacy_report.exists())


if __name__ == "__main__":
    unittest.main(verbosity=2)
