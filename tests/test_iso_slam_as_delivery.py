#!/usr/bin/env python3
"""Smoke test for alternative-splicing delivery figures and privacy scanning."""

from __future__ import annotations

import csv
import shutil
import subprocess
import sys
import tempfile
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "pipeline/scripts/10_delivery"))
from build_iso_slam_as_supplement import ensure_customer_package_desensitized  # noqa: E402


def write_tsv(path: Path, fields: list[str], rows: list[list[object]]) -> None:
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(fields)
        writer.writerows(rows)


def main() -> None:
    rscript = shutil.which("Rscript")
    if not rscript:
        raise RuntimeError("Rscript is required for the AS delivery plot smoke test")

    temp_root = ROOT / "tmp/test_runs"
    temp_root.mkdir(parents=True, exist_ok=True)
    with tempfile.TemporaryDirectory(prefix="as-delivery-", dir=temp_root) as raw:
        base = Path(raw)
        differential = base / "differential.tsv"
        sample_psi = base / "sample_psi.tsv"
        figure_dir = base / "figures"
        rows = [
            ["Dtag-short_vs_DMSO", "ENSG000001;SE:chr1:1-2", "GENE1", -0.30, 0.01, "TRUE"],
            ["Dtag-short_vs_DMSO", "ENSG000002;RI:chr1:3-4", "GENE2", 0.20, 0.04, "TRUE"],
            ["Dtag-long_vs_DMSO", "ENSG000001;SE:chr1:1-2", "GENE1", -0.12, 0.02, "TRUE"],
            ["Dtag-long_vs_DMSO", "ENSG000003;A5:chr1:5-6", "GENE3", 0.41, 0.03, "TRUE"],
            ["Dtag-PR_vs_Dtag-FL", "ENSG000004;MX:chr1:7-8", "GENE4", 0.01, 0.80, "FALSE"],
        ]
        write_tsv(differential, ["contrast", "event_id", "gene_name", "delta_PSI", "FDR", "significant"], rows)
        groups = ["DMSO", "Dtag-short", "Dtag-long", "Dtag-FL", "Dtag-PR"]
        psi_rows = []
        for event_id, gene_name, baseline in [
            ("ENSG000001;SE:chr1:1-2", "GENE1", 0.35),
            ("ENSG000002;RI:chr1:3-4", "GENE2", 0.55),
            ("ENSG000003;A5:chr1:5-6", "GENE3", 0.70),
        ]:
            for group_index, group in enumerate(groups):
                for replicate in (1, 2):
                    value = "NA" if event_id.startswith("ENSG000002") and group == "Dtag-PR" and replicate == 2 else f"{min(0.95, baseline + 0.02 * group_index + 0.01 * replicate):.3f}"
                    psi_rows.append([event_id, gene_name, f"{group}-{replicate}", group, value])
        write_tsv(sample_psi, ["event_id", "gene_name", "sample_id", "group", "new_RNA_PSI"], psi_rows)

        subprocess.run(
            [rscript, str(ROOT / "pipeline/scripts/10_delivery/plot_iso_slam_as_supplement.R"), str(differential), str(sample_psi), str(figure_dir)],
            cwd=ROOT,
            check=True,
        )
        expected = [
            "AS_contrast_overview.png", "AS_contrast_overview.pdf",
            "AS_significant_delta_psi.png", "AS_significant_delta_psi.pdf",
            "AS_candidate_sample_psi_page1.png", "AS_candidate_sample_psi.pdf",
        ]
        for filename in expected:
            output = figure_dir / filename
            assert output.is_file() and output.stat().st_size > 1000, filename
        assert (figure_dir / "AS_contrast_overview.png").read_bytes().startswith(b"\x89PNG\r\n\x1a\n")
        assert (figure_dir / "AS_significant_delta_psi.pdf").read_bytes().startswith(b"%PDF")

        ensure_customer_package_desensitized(figure_dir)
        privacy_probe = figure_dir / "privacy-probe.png"
        privacy_probe.write_bytes(b"PNG\x00/internal/internal/results\x00")
        try:
            ensure_customer_package_desensitized(figure_dir)
            raise AssertionError("Expected an absolute path in image metadata to fail the privacy scan")
        except RuntimeError:
            pass


if __name__ == "__main__":
    main()
