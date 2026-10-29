#!/usr/bin/env python3
"""Regression test for copy-only customer delivery packaging."""

from __future__ import annotations

import csv
import json
import shutil
import subprocess
import sys
import tempfile
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "pipeline/scripts/10_delivery"))
from build_iso_slam_delivery_package import build_analysis_summary, ensure_delivery_desensitized, group_delivery_target  # noqa: E402


def write_tsv(path: Path, fields: list[str], rows: list[list[object]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(fields)
        writer.writerows(rows)


def main() -> None:
    assert str(group_delivery_target(Path("tables/differential/Dtag-short_vs_DMSO/gene_abundance.tsv"))) == "04_Differential_Analysis/differential/Dtag-short_vs_DMSO/gene_abundance.tsv"
    groups = ["DMSO", "Dtag-short", "Dtag-long", "Dtag-FL", "Dtag-PR"]
    contrasts = ["Dtag-short_vs_DMSO", "Dtag-long_vs_DMSO", "Dtag-FL_vs_DMSO", "Dtag-PR_vs_Dtag-FL", "Dtag-long_vs_Dtag-short"]
    validation = {"groups": groups, "contrasts": contrasts, "biological_replicates_per_group": {group: 2 for group in groups}}
    sample_metrics = [{"raw_reads": 100, "insert_rate": 0.8, "isoform_assigned_reads": 50,
                       "filtered_isoforms": 3, "global_new_rna_fraction": 0.2,
                       "transcript_assignment_screen_rate": 0.02}]
    contrast_rows = [{"contrast": contrast, "significant_features": "0"} for contrast in contrasts]
    summary = build_analysis_summary(sample_metrics, [], [], contrast_rows, validation, False)
    assert "5 个独立组" in summary and "5 项预设比较" in summary
    assert "四个独立组" not in summary and "sample_reports/" not in summary

    test_tmp = ROOT / "tmp" / "test_runs"
    test_tmp.mkdir(exist_ok=True)
    with tempfile.TemporaryDirectory(prefix="delivery-package-", dir=test_tmp) as raw:
        base = Path(raw)
        results, scratch = base / "results", base / "tmp"
        sample = results / "S1"
        (results / "_project").mkdir(parents=True)
        (results / "_project/RUN_COMPLETE.tsv").write_text("status\ncompleted\n", encoding="utf-8")
        project = results / "_project/09_project_report/report.html"
        project.parent.mkdir(parents=True)
        project.write_text("<html><body>长读长 RNA 标记分析 逐样本审核指标 FASTQ 碱基 Q20 MAPQ≥20 检测基因数</body></html>", encoding="utf-8")
        review_fields = ["sample_id", "raw_reads", "raw_reads_status", "fastq_q20_base_fraction", "fastq_q20_status", "primary_mapq20_fraction", "primary_mapq20_status", "genes_detected", "genes_detected_status", "insert_retention_rate", "insert_mapping_rate", "overall_review_status"]
        write_tsv(results / "_project/09_project_report/01_tables/sample_qc_review.tsv", review_fields,
                  [["S1", 100, "FAIL", 0.85, "PASS", 0.95, "PASS", 12000, "PASS", 0.8, 0.9, "FAIL"]])
        report = sample / "07_sample_report/S1_report.html"
        report.parent.mkdir(parents=True)
        report.write_text("<html><body>样本报告</body></html>", encoding="utf-8")
        qc_fields = ["sample_id", "raw_reads", "insert_reads", "insert_rate", "isoform_assigned_reads", "dynamics_eligible_reads", "no_t_sites_reads", "assignment_conflict_reads", "assignment_conflict_rate", "catalog_isoforms", "filtered_isoforms", "grand3_global_binom_ntr", "posterior_mean", "posterior_grand3_absolute_difference", "switch_candidates", "source_results_root"]
        write_tsv(results / "_project/08_project_qc/01_tables/sample_qc_summary.tsv", qc_fields, [["S1", 100, 80, 0.8, 50, 45, 0, 5, 0.02, 3, 2, 0.2, 0.21, 0.01, 1, "/internal/private"]])
        group_root = results / "_project/08_project_qc/02_group_analysis"
        (group_root / "metadata").mkdir(parents=True)
        (group_root / "metadata/group_analysis_validation.json").write_text(json.dumps({
            "status": "PASS", "sample_count": 1, "groups": ["A"],
            "biological_replicates_per_group": {"A": 1}, "contrasts": ["A_vs_A"],
        }), encoding="utf-8")
        write_tsv(group_root / "tables/contrast_summary.tsv", ["contrast", "significant_features"], [["A_vs_A", 0]])
        as_root = group_root / "tables/differential/new_rna_splicing"
        (as_root / "new_rna_splicing_validation.json").parent.mkdir(parents=True, exist_ok=True)
        (as_root / "new_rna_splicing_validation.json").write_text(json.dumps({
            "status": "PASS", "method": "test", "event_types": ["SE"], "reference_event_count": 1,
            "tested_event_count_by_contrast": {"A_vs_A": 1},
            "significant_event_count_by_contrast": {"A_vs_A": 1},
            "significance_threshold": "FDR < 0.05", "biological_replicates_per_group": 1,
            "interpretation": "Exploratory", "internal_path": "/internal/private/results",
        }), encoding="utf-8")
        write_tsv(as_root / "differential_new_rna_splicing.tsv",
                  ["contrast", "event_id", "gene_id", "group_1", "group_2", "mean_PSI_group_1", "mean_PSI_group_2", "delta_PSI", "p_value", "FDR", "significant", "support_group_1", "support_group_2"],
                  [["A_vs_A", "ENSG000001;SE:chr1:1-2", "ENSG000001", "A", "A", 0.7, 0.3, 0.4, 0.01, 0.02, "TRUE", 20, 20]])
        write_tsv(as_root / "sample_event_psi.tsv",
                  ["event_id", "gene_id", "sample_id", "group", "inclusion_new_read_equivalents", "exclusion_new_read_equivalents", "total_new_read_equivalents", "new_RNA_PSI"],
                  [["ENSG000001;SE:chr1:1-2", "ENSG000001", "S1", "A", 20, 10, 30, 0.667]])
        write_tsv(as_root / "novel_structure_inventory.tsv",
                  ["sample_id", "group", "gene_id", "isoform_id", "structural_category", "read_support", "filter_status"],
                  [["S1", "A", "ENSG000001", "tx1", "novel", 5, "PASS"]])
        validation = sample / "06_isoform_dynamics/05_metadata/dynamics_validation.json"
        validation.parent.mkdir(parents=True)
        validation.write_text(json.dumps({"posterior_call_counts": {"New": 10, "Old": 30, "Ambiguous": 5, "No_T_sites": 0, "Assignment_conflict": 5}}), encoding="utf-8")
        fields = ["sample_id", "gene_id", "isoform_id", "total_reads", "expected_new_reads", "expected_old_reads", "posterior_new_fraction", "technical_bootstrap_ci_low", "technical_bootstrap_ci_high", "total_isoform_usage", "new_isoform_usage", "delta_isoform_usage", "switch_status"]
        row = ["S1", "g1", "S1|tx1", 50, 10, 40, 0.2, 0.1, 0.3, 1.0, 1.0, 0.0, "candidate"]
        write_tsv(sample / "06_isoform_dynamics/01_tables/S1.isoform_dynamics.tsv", fields, [row])
        write_tsv(sample / "06_isoform_dynamics/01_tables/S1.isoform_switch_candidates.tsv", fields, [row])
        write_tsv(sample / "04_isoform_catalog/01_tables/S1.isoform_catalog.tsv", ["sample_id", "isoform_id"], [["S1", "tx1"]])
        write_tsv(sample / "05_slam_model/01_tables/S1.grand3_gene_ntr.tsv", ["gene_id", "grand3_binom_ntr_map"], [["g1", 0.2]])
        write_tsv(sample / "05_slam_model/01_tables/S1.grand3_global_parameters.tsv", ["sample_id", "binom_ntr"], [["S1", 0.2]])
        write_tsv(sample / "06_isoform_dynamics/01_tables/S1.gene_dynamics.tsv", ["sample_id", "gene_id"], [["S1", "g1"]])
        samples = base / "samples.tsv"
        write_tsv(samples, ["sample_id", "group"], [["S1", "A"]])
        delivery = results / "_project/10_delivery/ISO_SLAM_Delivery"
        config = base / "config.yaml"
        rscript = shutil.which("Rscript")
        if not rscript:
            raise RuntimeError("Rscript is required for the integrated AS delivery test")
        config.write_text(f"results_root: {results}\ntmp_root: {scratch}\ndelivery_root: {delivery}\ndelivery_name: Test delivery\ninclude_sample_reports: false\nsamples_tsv: {samples}\nrscript_bin: {rscript}\n", encoding="utf-8")
        subprocess.run([sys.executable, "pipeline/scripts/10_delivery/build_iso_slam_delivery_package.py", "--config", str(config)], cwd=ROOT, check=True)
        assert (delivery / "00_Report_and_Figures/ISO_SLAM_project_report.html").is_file()
        assert (delivery / "RESULTS_ANALYSIS_CN.md").is_file()
        assert not (delivery / "sample_reports").exists()
        assert (delivery / "01_QC_and_Preprocessing/sample_overview.tsv").is_file()
        assert (delivery / "01_QC_and_Preprocessing/project_qc_detailed.tsv").is_file()
        assert (delivery / "01_QC_and_Preprocessing/sample_qc_review.tsv").is_file()
        assert (delivery / "01_QC_and_Preprocessing/REVIEW_CRITERIA.md").is_file()
        assert (delivery / "03_Expression_Quantification_and_Matrices/S1/isoform_dynamics.tsv").is_file()
        assert "gene_name" in (delivery / "03_Expression_Quantification_and_Matrices/S1/isoform_dynamics.tsv").read_text(encoding="utf-8").splitlines()[0]
        assert (delivery / "TABLE_GUIDE.md").is_file()
        as_delivery = delivery / "04_Differential_Analysis/differential/new_rna_splicing"
        assert (as_delivery / "AS_结果摘要.md").is_file()
        assert "AS_结果摘要.md" in (delivery / "RESULTS_ANALYSIS_CN.md").read_text(encoding="utf-8")
        assert "AS_结果摘要.md" in (delivery / "README.md").read_text(encoding="utf-8")
        for table in ("differential_new_rna_splicing.tsv", "sample_event_psi.tsv", "novel_structure_inventory.tsv"):
            assert "gene_name" in (as_delivery / table).read_text(encoding="utf-8").splitlines()[0]
        assert (as_delivery / "figures/AS_contrast_overview.png").is_file()
        assert (as_delivery / "figures/AS_significant_delta_psi.pdf").is_file()
        delivered_as_validation = json.loads((as_delivery / "new_rna_splicing_validation.json").read_text(encoding="utf-8"))
        assert "internal_path" not in delivered_as_validation
        delivery_validation = json.loads((delivery / "DELIVERY_VALIDATION.json").read_text(encoding="utf-8"))
        assert delivery_validation["alternative_splicing"]["distinct_significant_event_count"] == 1
        assert "reads_without_usable_t_sites" not in (delivery / "01_QC_and_Preprocessing/read_classification_summary.tsv").read_text(encoding="utf-8")
        assert all(not path.is_symlink() for path in delivery.rglob("*"))
        manifest = (delivery / "delivery_manifest.tsv").read_text(encoding="utf-8")
        assert "sample_reports/S1.html" not in manifest
        assert "03_Expression_Quantification_and_Matrices/S1/gene_rna_labeling.tsv" in manifest
        assert "04_Differential_Analysis/differential/new_rna_splicing/AS_结果摘要.md" in manifest
        assert "/internal/private" not in (delivery / "01_QC_and_Preprocessing/project_qc_detailed.tsv").read_text(encoding="utf-8")
        marker_text = (results / "_project/RUN_COMPLETE.tsv").read_text(encoding="utf-8")
        assert marker_text == "status\ncompleted\n"
        leak = delivery / "privacy-probe.html"
        leak.write_text('<script>const localPath="/internal/private/results";</script>', encoding="utf-8")
        try:
            ensure_delivery_desensitized(delivery)
            raise AssertionError("Expected hidden HTML paths to fail de-identification scanning")
        except RuntimeError:
            pass
        leak.unlink()
        binary_leak = delivery / "privacy-probe.png"
        binary_leak.write_bytes(b"PNG\x00/work/private/results\x00")
        try:
            ensure_delivery_desensitized(delivery)
            raise AssertionError("Expected binary metadata paths to fail de-identification scanning")
        except RuntimeError:
            pass
        binary_leak.unlink()


if __name__ == "__main__":
    main()
