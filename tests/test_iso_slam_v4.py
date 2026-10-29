#!/usr/bin/env python3
"""Focused V4 posterior and isoform-filter regression tests."""

from __future__ import annotations

import csv
import gzip
import importlib.util
import sys
import tempfile
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "pipeline/scripts"))


def load(name: str, path: Path):
    spec = importlib.util.spec_from_file_location(name, path)
    assert spec and spec.loader
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def main() -> None:
    dynamics = load("build_isoform_dynamics", ROOT / "pipeline/scripts/build_isoform_dynamics.py")
    assert dynamics.read_posterior(0, 100, 0.2, 0.001, 0.25) < 0.1
    assert dynamics.read_posterior(30, 100, 0.2, 0.001, 0.25) > 0.9
    underflow_safe = dynamics.read_posterior(500, 2000, 0.2, 0.001, 0.25)
    assert 0 <= underflow_safe <= 1
    assert abs(dynamics.read_posterior(0, 0, 0.2, 0.001, 0.25) - 0.2) < 1e-12

    catalog = load("run_isoform_catalog", ROOT / "pipeline/scripts/run_isoform_catalog.py")
    config = {
        "sqanti_min_reads_fsm": "1",
        "sqanti_min_reads_nic": "2",
        "sqanti_min_reads_other": "3",
        "sqanti_min_reads_novel_single_exon": "5",
        "sqanti_reject_rtswitching": "true",
        "sqanti_reject_intrapriming": "true",
        "sqanti_reject_noncanonical": "true",
    }
    assert catalog.filter_reason({"structural_category": "full-splice_match", "exons": "2"}, 1, config) == "PASS"
    assert catalog.filter_reason({"structural_category": "novel_in_catalog", "exons": "2"}, 1, config) == "support_lt_2"
    assert catalog.filter_reason({"structural_category": "novel_not_in_catalog", "exons": "1"}, 3, config) == "support_lt_5"
    assert catalog.filter_reason({"structural_category": "novel_not_in_catalog", "exons": "2", "RTS_stage": "TRUE"}, 5, config) == "rt_switching"

    test_tmp = ROOT / "tmp" / "test_runs"
    test_tmp.mkdir(exist_ok=True)
    with tempfile.TemporaryDirectory(prefix="iso-slam-v4-", dir=test_tmp) as raw:
        root = Path(raw)
        read_calls = root / "calls.tsv.gz"
        assignments = root / "assignments.tsv.gz"
        catalog_path = root / "catalog.tsv"
        global_path = root / "global.tsv"
        genes_path = root / "genes.tsv"
        with gzip.open(read_calls, "wt", encoding="utf-8", newline="") as handle:
            writer = csv.DictWriter(handle, fieldnames=["read_id", "gene_id", "region_class", "assignment_level", "t_covered", "tc_count"], delimiter="\t", lineterminator="\n")
            writer.writeheader()
            for index in range(11):
                conflict = index == 10
                writer.writerow({"read_id": f"r{index}", "gene_id": "." if conflict else "g1", "region_class": "intergenic" if conflict else "gene", "assignment_level": "unresolved" if conflict else "gene", "t_covered": 0 if index >= 9 else 100, "tc_count": 30 if index < 2 else 0})
        with gzip.open(assignments, "wt", encoding="utf-8", newline="") as handle:
            writer = csv.DictWriter(handle, fieldnames=["read_id", "isoform_id", "gene_id"], delimiter="\t", lineterminator="\n")
            writer.writeheader()
            for index in range(11):
                writer.writerow({"read_id": f"r{index}", "isoform_id": "tx1", "gene_id": "g1"})
        catalog_path.write_text("source_isoform_id\tfilter_status\tgene_id\tstructural_category\ntx1\tPASS\tg1\tfull-splice_match\n", encoding="utf-8")
        global_path.write_text("binom_ntr\tbinom_p_err\tbinom_p_conv\n0.2\t0.001\t0.25\n", encoding="utf-8")
        genes_path.write_text("gene_id\tgrand3_binom_ntr_map\ng1\t0.2\n", encoding="utf-8")
        output = root / "out"
        import subprocess, sys
        subprocess.run([sys.executable, str(ROOT / "pipeline/scripts/build_isoform_dynamics.py"), "--sample-id", "S1", "--read-calls", str(read_calls), "--read-assignments", str(assignments), "--catalog", str(catalog_path), "--grand3-global", str(global_path), "--grand3-genes", str(genes_path), "--outdir", str(output), "--bootstrap-iterations", "20"], cwd=ROOT, check=True)
        validation = __import__("json").loads((output / "05_metadata/dynamics_validation.json").read_text(encoding="utf-8"))
        assert validation["posterior_concordance_status"] == "PASS"
        assert validation["weighted_closure_error"] <= 1e-6
        assert sum(validation["posterior_call_counts"].values()) == 11
        assert validation["posterior_call_counts"]["New"] == 2
        assert validation["posterior_call_counts"]["Old"] == 7
        assert validation["posterior_call_counts"]["No_T_sites"] == 1
        assert validation["posterior_call_counts"]["Ambiguous"] == 0
        assert validation["posterior_call_counts"]["Assignment_conflict"] == 1
        assert validation["assignment_conflict_subtype_counts"]["Caller_intergenic"] == 1
        assert validation["dynamics_eligible_reads"] == 10

        isoquant = root / "isoquant"
        isoquant.mkdir()
        (isoquant / "S1.transcript_models.gtf").write_text(
            'chr1\tIsoQuant\ttranscript\t1\t100\t.\t+\t.\tgene_id "g1"; transcript_id "tx1";\n'
            'chr1\tIsoQuant\texon\t1\t100\t.\t+\t.\tgene_id "g1"; transcript_id "tx1";\n'
            'chr1\tIsoQuant\ttranscript\t201\t300\t.\t+\t.\tgene_id "g1"; transcript_id "tx2";\n',
            encoding="utf-8",
        )
        with gzip.open(isoquant / "S1.read_assignments.tsv.gz", "wt", encoding="utf-8", newline="") as handle:
            writer = csv.DictWriter(handle, fieldnames=["read_id", "transcript_id", "gene_id", "assignment_type", "polya_position"], delimiter="\t", lineterminator="\n")
            writer.writeheader()
            writer.writerow({"read_id": "r1", "transcript_id": "tx1", "gene_id": "g1", "assignment_type": "unique", "polya_position": 100})
            writer.writerow({"read_id": "r2", "transcript_id": "tx2", "gene_id": "g1", "assignment_type": "unique", "polya_position": ""})
        (isoquant / "S1.transcript_model_counts.tsv").write_text("transcript_id\tcount\ntx1\t1\ntx2\t1\n", encoding="utf-8")
        classification = root / "classification.txt"
        classification.write_text("isoform\tassociated_gene\tstructural_category\tlength\texons\tRTS_stage\n"
                                  "tx1\tg1\tfull-splice_match\t100\t1\tFALSE\n"
                                  "tx2\tg1\tnovel_in_catalog\t100\t2\tFALSE\n", encoding="utf-8")
        config_path = root / "config.yaml"
        config_path.write_text("sqanti_min_reads_fsm: 1\nsqanti_min_reads_nic: 2\nsqanti_min_reads_other: 3\nsqanti_min_reads_novel_single_exon: 5\n", encoding="utf-8")
        catalog_out = root / "catalog_out"
        subprocess.run([sys.executable, str(ROOT / "pipeline/scripts/run_isoform_catalog.py"), "--config", str(config_path), "--sample-id", "S1", "--bam", str(root / "unused.bam"), "--fastq", str(root / "unused.fastq"), "--fasta", str(root / "unused.fa"), "--gtf", str(root / "unused.gtf"), "--outdir", str(catalog_out), "--skip-tools", "--isoquant-dir", str(isoquant), "--sqanti-classification", str(classification)], cwd=ROOT, check=True)
        catalog_validation = __import__("json").loads((catalog_out / "05_metadata/catalog_validation.json").read_text(encoding="utf-8"))
        assert catalog_validation["filtered_isoforms"] == 1
        assert catalog_validation["polya_evaluable_fraction"] == 0.5
        filtered = (catalog_out / "04_tool_output/filtered/S1.filtered_isoforms.gtf").read_text(encoding="utf-8")
        assert "tx1" in filtered and "tx2" not in filtered


if __name__ == "__main__":
    main()
