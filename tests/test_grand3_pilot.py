#!/usr/bin/env python3
"""Synthetic regression test for Grand3 sample-level concordance parsing."""

from __future__ import annotations

import csv
import gzip
import importlib.util
import tempfile
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
SPEC = importlib.util.spec_from_file_location("run_grand3_pilot", ROOT / "pipeline/scripts/run_grand3_pilot.py")
assert SPEC and SPEC.loader
run_grand3_pilot = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(run_grand3_pilot)


def main() -> None:
    with tempfile.TemporaryDirectory(prefix="iso-slam-grand3-") as raw:
        root = Path(raw)
        results = root / "results"
        for sample_id, fraction in (("S1", 0.2), ("S2", 0.4)):
            caller = results / sample_id / "04_slam_calling/04_tool_output/caller"
            caller.mkdir(parents=True)
            with (caller / f"{sample_id}.bulk_slam_summary.tsv").open(
                "w", encoding="utf-8", newline=""
            ) as handle:
                writer = csv.DictWriter(
                    handle,
                    fieldnames=["sample", "new_fraction_excluding_ambiguous"],
                    delimiter="\t",
                    lineterminator="\n",
                )
                writer.writeheader()
                writer.writerow({"sample": sample_id, "new_fraction_excluding_ambiguous": fraction})

        grand3 = root / "data.tsv.gz"
        models = root / "model.parameters.tsv"
        fields = [
            "Gene",
            "S1.4sU.bam Read count",
            "S2.4sU.bam Read count",
            "S1.4sU.bam 4sU Binom NTR MAP",
            "S2.4sU.bam 4sU Binom NTR MAP",
        ]
        with gzip.open(grand3, "wt", encoding="utf-8", newline="") as handle:
            writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", lineterminator="\n")
            writer.writeheader()
            writer.writerow(
                {
                    "Gene": "G1",
                    "S1.4sU.bam Read count": 10,
                    "S2.4sU.bam Read count": 20,
                    "S1.4sU.bam 4sU Binom NTR MAP": 0.2,
                    "S2.4sU.bam 4sU Binom NTR MAP": 0.4,
                }
            )
        with models.open("w", encoding="utf-8", newline="") as handle:
            writer = csv.DictWriter(
                handle,
                fieldnames=["Estimator", "Condition", "Subread", "Label", "Binom ntr"],
                delimiter="\t",
                lineterminator="\n",
            )
            writer.writeheader()
            writer.writerow({"Estimator": "Separate", "Condition": "/tmp/S1.4sU", "Subread": "Antisense", "Label": "4sU", "Binom ntr": 0.2})
            writer.writerow({"Estimator": "Separate", "Condition": "/tmp/S2.4sU", "Subread": "Antisense", "Label": "4sU", "Binom ntr": 0.4})
        output = root / "comparison.tsv"
        passed, rows = run_grand3_pilot.compare_sample_ntr(
            grand3, models, results, ["S1", "S2"], "Antisense", 0.05, output
        )
        assert passed
        assert all(float(row["absolute_error"]) == 0.0 for row in rows)
        assert all(float(row["grand3_global_binomial_ntr"]) in {0.2, 0.4} for row in rows)
        assert output.is_file()


if __name__ == "__main__":
    main()
