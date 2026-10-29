#!/usr/bin/env python3
"""Select one project-wide strand mode from passed sample audits."""

from __future__ import annotations

import argparse
import csv
import json
from pathlib import Path


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--samples", required=True)
    parser.add_argument("--audit-root", required=True)
    parser.add_argument("--output", required=True)
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    with Path(args.samples).open(encoding="utf-8", newline="") as handle:
        sample_ids = [row["sample_id"] for row in csv.DictReader(handle, delimiter="\t")]
    audits = []
    for sample_id in sample_ids:
        path = (
            Path(args.audit_root)
            / sample_id
            / "04_slam_calling"
            / "05_metadata"
            / "strand_orientation"
            / f"{sample_id}.strand_orientation_audit.json"
        )
        if not path.is_file():
            raise SystemExit(f"Missing strand audit: {path}")
        audits.append(json.loads(path.read_text(encoding="utf-8")))
    modes = {str(item.get("selected_mode", "")) for item in audits}
    failed = [str(item.get("sample_id", "")) for item in audits if item.get("status") != "PASS"]
    passed = not failed and len(modes) == 1 and modes <= {"forward", "reverse"}
    payload = {
        "schema_version": "iso-slam-project-strand-selection-v1",
        "status": "PASS" if passed else "FAIL",
        "selected_mode": next(iter(modes)) if len(modes) == 1 else "",
        "sample_count": len(audits),
        "samples": [item.get("sample_id", "") for item in audits],
        "failed_samples": failed,
        "observed_modes": sorted(modes),
    }
    output = Path(args.output)
    output.parent.mkdir(parents=True, exist_ok=True)
    output.write_text(json.dumps(payload, indent=2) + "\n", encoding="utf-8")
    print(output)
    raise SystemExit(0 if passed else 1)


if __name__ == "__main__":
    main()
