#!/usr/bin/env python3
from __future__ import annotations

import csv
import json
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
SCRIPT = ROOT / "pipeline/scripts" / "build_lr_bulk_read_flow.py"


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


class ReadFlowAdapterTests(unittest.TestCase):
    def setUp(self) -> None:
        self.temp = tempfile.TemporaryDirectory(prefix="iso-slam-flow-")
        self.root = Path(self.temp.name)
        self.sample = "S-1"
        self.sample_root = self.root / self.sample
        trim = self.sample_root / "02_trim/04_tool_output/trim/S-1.trim_summary.json"
        align = self.sample_root / "03_alignment/04_tool_output/minimap2_samtools"
        qc = self.sample_root / "04_slam_calling/04_tool_output/caller/S-1.bulk_slam_qc.json"
        trim.parent.mkdir(parents=True, exist_ok=True)
        qc.parent.mkdir(parents=True, exist_ok=True)
        trim.write_text(json.dumps({"total_reads": 100, "trimmed_reads": 80, "unassigned_reads": 20}), encoding="utf-8")
        write_flagstat(align / "S-1.flagstat.txt", total=85, primary=80, primary_mapped=70, supplementary=5)
        write_flagstat(align / "S-1.primary.flagstat.txt", total=60, primary=60, primary_mapped=60)
        qc.write_text(
            json.dumps(
                {
                    "total_reads": 60,
                    "assigned_reads": 50,
                    "unassigned_gene": 10,
                    "assignment_intergenic_reads": 4,
                    "assignment_shared_antisense_ambiguous_reads": 3,
                    "assignment_shared_gene_ambiguous_reads": 3,
                    "assignment_antisense_only_reads": 0,
                }
            ),
            encoding="utf-8",
        )
        self.qc = qc

    def tearDown(self) -> None:
        self.temp.cleanup()

    def run_adapter(self) -> subprocess.CompletedProcess[str]:
        return subprocess.run(
            [
                sys.executable,
                str(SCRIPT),
                "--sample-id", self.sample,
                "--sample-root", str(self.sample_root),
                "--output-dir", str(self.root / "out"),
            ],
            cwd=ROOT,
            text=True,
            capture_output=True,
            check=False,
        )

    def test_builds_exact_conserved_flow_and_preserves_zero_branch(self) -> None:
        result = self.run_adapter()
        self.assertEqual(result.returncode, 0, result.stderr)
        out = self.root / "out"
        with (out / "tables/read_flow_nodes.tsv").open(encoding="utf-8", newline="") as handle:
            nodes = {row["node_id"]: row for row in csv.DictReader(handle, delimiter="\t")}
        self.assertEqual(nodes["raw"]["value"], "100")
        self.assertEqual(nodes["primary_mapped"]["value"], "70")
        self.assertEqual(nodes["primary_q20"]["value"], "60")
        self.assertEqual(nodes["gene_assigned"]["value"], "50")
        self.assertEqual(nodes["antisense_only"]["value"], "0")
        validation = json.loads((out / "reports/read_flow_adapter_validation.json").read_text(encoding="utf-8"))
        self.assertTrue(validation["ok"])
        self.assertEqual(validation["transformations"], 0)
        self.assertEqual(validation["unresolved_residuals"], 0)
        self.assertTrue((out / "reports/read_flow_sankey_section.html").is_file())
        self.assertTrue((out / "tables/read_flow_source_manifest.tsv").is_file())
        fragment = (out / "reports/read_flow_sankey_section.html").read_text(encoding="utf-8")
        self.assertIn("本视图不包含跨单位转换", fragment)
        self.assertNotIn("Read-pair", fragment)
        self.assertNotIn("Bin50", fragment)
        self.assertIn("flowShortLabels", fragment)

    def test_rejects_broken_gene_partition(self) -> None:
        values = json.loads(self.qc.read_text(encoding="utf-8"))
        values["assignment_intergenic_reads"] = 5
        self.qc.write_text(json.dumps(values), encoding="utf-8")
        result = self.run_adapter()
        self.assertEqual(result.returncode, 2)
        self.assertIn("unassigned gene fates", result.stderr)

    def test_represents_reads_without_primary_alignment_record(self) -> None:
        align = self.sample_root / "03_alignment/04_tool_output/minimap2_samtools/S-1.flagstat.txt"
        write_flagstat(align, total=84, primary=79, primary_mapped=70, supplementary=5)
        result = self.run_adapter()
        self.assertEqual(result.returncode, 0, result.stderr)
        with (self.root / "out/tables/read_flow_nodes.tsv").open(encoding="utf-8", newline="") as handle:
            nodes = {row["node_id"]: row for row in csv.DictReader(handle, delimiter="\t")}
        self.assertEqual(nodes["primary_unrepresented"]["value"], "1")

    def test_rejects_primary_records_above_trimmed_reads(self) -> None:
        align = self.sample_root / "03_alignment/04_tool_output/minimap2_samtools/S-1.flagstat.txt"
        write_flagstat(align, total=86, primary=81, primary_mapped=70, supplementary=5)
        result = self.run_adapter()
        self.assertEqual(result.returncode, 2)
        self.assertIn("Alignment-derived loss counts must be non-negative", result.stderr)


if __name__ == "__main__":
    unittest.main(verbosity=2)
