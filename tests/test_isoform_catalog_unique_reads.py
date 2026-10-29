#!/usr/bin/env python3
"""IsoQuant candidate rows must collapse to one confident assignment per read."""

from __future__ import annotations

import csv
import gzip
import sys
import tempfile
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "pipeline/scripts"))
from run_isoform_catalog import parse_assignments


def main() -> None:
    with tempfile.TemporaryDirectory(prefix="isoform-assignments-", dir=ROOT / "tmp/test_runs") as raw:
        path = Path(raw) / "read_assignments.tsv.gz"
        fields = ["read_id", "chr", "strand", "isoform_id", "gene_id", "assignment_type", "polya_position"]
        records = [
            ["R1", "chr1", "+", "TX1", "G1", "unique", ""],
            ["R1", "chr1", "+", ".", ".", "noninformative", ""],
            ["R2", "chr1", "+", "TX2", "G1", "inconsistent_ambiguous", ""],
            ["R2", "chr1", "+", "TX3", "G1", "inconsistent_ambiguous", ""],
            ["R3", "chr2", "-", "TX4", "G2", "unique_minor_difference", ""],
            ["R4", "chr2", "+", "TX5", "G3", "unique", ""],
            ["R4", "chr2", "+", "TX6", "G3", "unique_minor_difference", ""],
        ]
        with gzip.open(path, "wt", encoding="utf-8", newline="") as handle:
            handle.write("# IsoQuant fixture\n# version: test\n")
            writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
            writer.writerow(fields)
            writer.writerows(records)

        normalized, counts, support = parse_assignments(path, "S1")
        assert counts["assignment_rows"] == 7
        assert counts["total_reads"] == 4
        assert counts["duplicate_assignment_rows"] == 3
        assert counts["assigned_reads"] == 2
        assert counts["unassigned_reads"] == 2
        assert counts["multiple_confident_assignments"] == 1
        assert [row["read_id"] for row in normalized] == ["R1", "R3"]
        assert support == {"TX1": {"R1"}, "TX4": {"R3"}}
    print("IsoQuant unique read assignment accounting PASS")


if __name__ == "__main__":
    main()
