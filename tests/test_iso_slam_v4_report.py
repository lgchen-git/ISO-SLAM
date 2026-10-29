#!/usr/bin/env python3
"""Synthetic offline report and conserved V4 read-flow regression."""

from __future__ import annotations

import base64
import csv
import re
import json
import subprocess
import sys
import tempfile
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]


def flagstat(path: Path, total: int, mapped: int) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(
        f"{total} + 0 in total (QC-passed reads + QC-failed reads)\n{total} + 0 primary\n0 + 0 secondary\n0 + 0 supplementary\n{mapped} + 0 mapped (0.00% : N/A)\n{mapped} + 0 primary mapped (0.00% : N/A)\n",
        encoding="utf-8",
    )


def write_tsv(path: Path, fields: list[str], rows: list[list[object]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(fields)
        writer.writerows(rows)


def main() -> None:
    test_tmp = ROOT / "tmp" / "test_runs"
    test_tmp.mkdir(exist_ok=True)
    with tempfile.TemporaryDirectory(prefix="v4-report-", dir=test_tmp) as raw:
        root = Path(raw) / "results/S1"
        trim = root / "02_trim/04_tool_output/trim/S1.trim_summary.json"
        trim.parent.mkdir(parents=True)
        trim.write_text(json.dumps({"total_reads": 100, "trimmed_reads": 80}), encoding="utf-8")
        align = root / "03_alignment/04_tool_output/minimap2_samtools"
        flagstat(align / "S1.flagstat.txt", 80, 70)
        flagstat(align / "S1.primary.flagstat.txt", 60, 60)
        catalog_meta = root / "04_isoform_catalog/05_metadata/catalog_validation.json"
        catalog_meta.parent.mkdir(parents=True)
        catalog_meta.write_text(json.dumps({"schema_version": "iso-slam-isoform-catalog-v4", "status": "PASS", "catalog_isoforms": 3, "filtered_isoforms": 2}), encoding="utf-8")
        grand_meta = root / "05_slam_model/05_metadata/grand3_validation.json"
        grand_meta.parent.mkdir(parents=True)
        grand_meta.write_text(json.dumps({"schema_version": "iso-slam-grand3-sample-v4", "status": "PASS", "global_parameters": {"binom_ntr": 0.2}}), encoding="utf-8")
        write_tsv(root / "05_slam_model/01_tables/S1.grand3_gene_ntr.tsv", ["gene_id", "read_count", "grand3_binom_ntr_map"], [["g1", 50, 0.2]])
        dyn_meta = root / "06_isoform_dynamics/05_metadata/dynamics_validation.json"
        dyn_meta.parent.mkdir(parents=True)
        dyn_meta.write_text(json.dumps({"schema_version": "iso-slam-isoform-dynamics-v6", "status": "PASS", "joined_reads": 50, "dynamics_eligible_reads": 49, "posterior_call_counts": {"New": 10, "Old": 29, "Ambiguous": 5, "No_T_sites": 5, "Assignment_conflict": 1}, "assignment_conflict_subtype_counts": {"Caller_intergenic": 1}, "posterior_mean": 0.2, "posterior_grand3_absolute_difference": 0.0, "switch_candidates": 1}), encoding="utf-8")
        caller = Path(raw) / "v3_caller"
        caller.mkdir()
        (caller / "S1.bulk_slam_qc.json").write_text(json.dumps({"total_reads": 60, "t_to_c_rate": 0.02, "bg_tc_rate": 0.002}), encoding="utf-8")
        write_tsv(caller / "S1.substitution_summary.tsv", ["substitution", "count", "denominator", "rate", "role", "source"], [["T>C", 1, 50, 0.02, "observed", "site_calls"], ["diagnostic_background", "", "", 0.002, "diagnostic_only", "max_ta_tg_rate_background_site_filtered"]])
        write_tsv(root / "04_isoform_catalog/01_tables/S1.isoform_catalog.tsv", ["isoform_id", "gene_id", "structural_category", "read_support", "filter_status"], [["S1|tx1", "g1", "full-splice_match", 30, "PASS"], ["S1|tx2", "g1", "novel_in_catalog", 20, "PASS"], ["S1|tx3", "g2", "novel_not_in_catalog", 1, "FAIL"]])
        fields = ["isoform_id", "gene_id", "total_reads", "expected_new_reads", "posterior_new_fraction", "technical_bootstrap_ci_low", "technical_bootstrap_ci_high", "total_isoform_usage", "new_isoform_usage", "delta_isoform_usage"]
        write_tsv(root / "06_isoform_dynamics/01_tables/S1.isoform_dynamics.tsv", fields, [["S1|tx1", "g1", 30, 9, 0.3, 0.2, 0.4, 0.6, 0.9, 0.3], ["S1|tx2", "g1", 20, 1, 0.05, 0.01, 0.1, 0.4, 0.1, -0.3]])
        write_tsv(root / "06_isoform_dynamics/01_tables/S1.isoform_switch_candidates.tsv", fields, [["S1|tx1", "g1", 30, 9, 0.3, 0.2, 0.4, 0.6, 0.9, 0.3]])
        report_stage = root / "07_sample_report"
        flow = report_stage / "04_tool_output/read_flow"
        subprocess.run([sys.executable, "pipeline/scripts/build_lr_iso_slam_read_flow.py", "--sample-id", "S1", "--sample-root", str(root), "--outdir", str(flow)], cwd=ROOT, check=True)
        subprocess.run([sys.executable, "pipeline/scripts/build_lr_iso_slam_report_v4.py", "--sample-id", "S1", "--sample-root", str(root), "--v3-caller", str(caller), "--report-title", "Synthetic result report"], cwd=ROOT, check=True)
        report = report_stage / "S1_report.html"
        checked = subprocess.run([sys.executable, "pipeline/scripts/00_shared/validate_delivery_report_v2.py", "--report", str(report), "--manifest", str(report_stage / "05_metadata/delivery_manifest.tsv"), "--require-section", "1. 样本概览", "--require-section", "3. 转录本结构概览", "--require-section", "4. RNA 标记信号", "--forbid-text", "Assignment conflict", "--json"], cwd=ROOT, text=True, capture_output=True, check=False)
        assert checked.returncode == 0, checked.stdout + checked.stderr
        validation = json.loads(checked.stdout)
        assert validation["ok"]
        production = subprocess.run([sys.executable, "pipeline/scripts/00_shared/validate_bio_report.py", "--report", str(report), "--manifest", str(report_stage / "05_metadata/delivery_manifest.tsv"), "--require-expandable-figures", "--require-print-css", "--json"], cwd=ROOT, text=True, capture_output=True)
        assert production.returncode == 0, production.stdout + production.stderr
        text = report.read_text(encoding="utf-8")
        library = re.search(r"eval\(atob\((.*?)\)\)", text)
        assert library is not None
        decoded = base64.b64decode(json.loads(library.group(1)), validate=True)
        assert decoded == (ROOT / "pipeline/assets/plotly-3.5.1.min.js").read_bytes()
        assert text.count("class='kpi'") == 6
        assert "T&gt;C conversion and diagnostic background" in text
        assert "Diagnostic background" in text
        assert "诊断背景" not in text
        assert 'text-anchor="end"' in text
        assert "全碱基错配概率谱" not in text
        assert "Excluded after assignment screen" in text
        assert "No usable T sites" in text
        assert "Assignment conflict" not in text
        assert "ISO-SLAM V4" not in text
        assert "id='bio-flow-view'" not in text
        assert "2. 数据处理与长读长质量" in text
        assert "SQANTI3 未运行 ORF 预测" not in text
        assert "方法与限制" not in text
        payload = json.loads(dyn_meta.read_text(encoding="utf-8"))
        payload["posterior_call_counts"]["No_T_sites"] = 0
        dyn_meta.write_text(json.dumps(payload), encoding="utf-8")
        subprocess.run([sys.executable, "pipeline/scripts/build_lr_iso_slam_report_v4.py", "--sample-id", "S1", "--sample-root", str(root), "--v3-caller", str(caller), "--report-title", "Synthetic result report"], cwd=ROOT, check=True)
        assert "No usable T sites" not in report.read_text(encoding="utf-8")


if __name__ == "__main__":
    main()
