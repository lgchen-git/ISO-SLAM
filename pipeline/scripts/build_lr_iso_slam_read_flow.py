#!/usr/bin/env python3
"""Adapt V4 counts to one conserved integer-read Sankey."""

from __future__ import annotations

import argparse
import csv
import json
import re
import subprocess
import sys
from pathlib import Path

from iso_slam_v4_common import nonempty, write_json, write_tsv


NODE_FIELDS = ["view_id", "node_id", "stage", "label", "value", "unit", "scope", "node_kind"]
LINK_FIELDS = ["view_id", "link_id", "source", "target", "value", "reason"]


def flagstat(path: Path) -> dict[str, int]:
    values: dict[str, int] = {}
    for line in path.read_text(encoding="utf-8").splitlines():
        match = re.match(r"^(\d+) \+ \d+ (.+)$", line)
        if not match:
            continue
        if match.group(2).startswith("in total"):
            values["total"] = int(match.group(1))
        elif match.group(2) == "primary":
            values["primary"] = int(match.group(1))
        elif match.group(2).startswith("primary mapped"):
            values["mapped"] = int(match.group(1))
    return values


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--sample-id", required=True)
    parser.add_argument("--sample-root", required=True, type=Path)
    parser.add_argument("--outdir", required=True, type=Path)
    parser.add_argument("--plotly-js", type=Path, default=Path("pipeline/assets/plotly-3.5.1.min.js"))
    args = parser.parse_args()
    sample = args.sample_id
    trim_path = args.sample_root / f"02_trim/04_tool_output/trim/{sample}.trim_summary.json"
    align_path = args.sample_root / f"03_alignment/04_tool_output/minimap2_samtools/{sample}.flagstat.txt"
    q20_path = args.sample_root / f"03_alignment/04_tool_output/minimap2_samtools/{sample}.primary.flagstat.txt"
    dyn_path = args.sample_root / "06_isoform_dynamics/05_metadata/dynamics_validation.json"
    for path in (trim_path, align_path, q20_path, dyn_path):
        if not nonempty(path):
            raise RuntimeError(f"Missing Sankey source: {path}")
    trim = json.loads(trim_path.read_text(encoding="utf-8"))
    aligned = flagstat(align_path)
    q20 = flagstat(q20_path)
    dynamics = json.loads(dyn_path.read_text(encoding="utf-8"))
    raw = int(trim["total_reads"])
    insert = int(trim["trimmed_reads"])
    structure_loss = raw - insert
    mapped = int(aligned["mapped"])
    unmapped = insert - mapped
    mapq20 = int(q20["total"])
    low_mapq = mapped - mapq20
    assigned = int(dynamics["joined_reads"])
    isoform_unassigned = mapq20 - assigned
    calls = dynamics["posterior_call_counts"]
    new = int(calls.get("New", 0)); old = int(calls.get("Old", 0))
    ambiguous = int(calls.get("Ambiguous", 0)); no_t_sites = int(calls.get("No_T_sites", 0))
    assignment_conflict = int(calls.get("Assignment_conflict", 0))
    checks = [raw == insert + structure_loss, insert == mapped + unmapped, mapped == mapq20 + low_mapq, mapq20 == assigned + isoform_unassigned, assigned == new + old + ambiguous + no_t_sites + assignment_conflict]
    if min(structure_loss, unmapped, low_mapq, isoform_unassigned, new, old, ambiguous, no_t_sites, assignment_conflict) < 0 or not all(checks):
        raise RuntimeError("V4 read-flow counts are not conserved")
    values = {"raw": raw, "insert": insert, "structure_loss": structure_loss, "mapped": mapped, "unmapped": unmapped, "mapq20": mapq20, "low_mapq": low_mapq, "assigned": assigned, "isoform_unassigned": isoform_unassigned, "new": new, "old": old, "ambiguous": ambiguous, "no_t_sites": no_t_sites, "assignment_conflict": assignment_conflict}
    stages = {"raw": 0, "insert": 1, "structure_loss": 1, "mapped": 2, "unmapped": 2, "mapq20": 3, "low_mapq": 3, "assigned": 4, "isoform_unassigned": 4, "new": 5, "old": 5, "ambiguous": 5, "no_t_sites": 5, "assignment_conflict": 5}
    labels = {"raw": "Raw reads", "insert": "Insert reads", "structure_loss": "Structure loss", "mapped": "Primary mapped", "unmapped": "Unmapped", "mapq20": "Primary MAPQ≥20", "low_mapq": "MAPQ<20", "assigned": "Filtered isoform assigned", "isoform_unassigned": "No filtered isoform", "new": "New (P≥0.50)", "old": "Old (P≤0.30)", "ambiguous": "Ambiguous (0.30<P<0.50)", "no_t_sites": "No T sites", "assignment_conflict": "Assignment conflict"}
    retained = {"raw", "insert", "mapped", "mapq20", "assigned"}
    nodes = [{"view_id": "sample_reads", "node_id": key, "stage": stages[key], "label": labels[key], "value": value, "unit": "reads", "scope": "sample_library", "node_kind": "retained" if key in retained else "branch" if key in {"new", "old", "ambiguous", "no_t_sites", "assignment_conflict"} else "loss"} for key, value in values.items()]
    specs = [("raw_insert", "raw", "insert"), ("raw_loss", "raw", "structure_loss"), ("insert_mapped", "insert", "mapped"), ("insert_unmapped", "insert", "unmapped"), ("mapped_q20", "mapped", "mapq20"), ("mapped_low", "mapped", "low_mapq"), ("q20_assigned", "mapq20", "assigned"), ("q20_unassigned", "mapq20", "isoform_unassigned"), ("assigned_new", "assigned", "new"), ("assigned_old", "assigned", "old"), ("assigned_ambiguous", "assigned", "ambiguous"), ("assigned_no_t_sites", "assigned", "no_t_sites"), ("assigned_assignment_conflict", "assigned", "assignment_conflict")]
    links = [{"view_id": "sample_reads", "link_id": link, "source": source, "target": target, "value": values[target], "reason": labels[target]} for link, source, target in specs]
    inputs = args.outdir / "inputs"
    write_tsv(inputs / "read_flow_nodes.input.tsv", nodes, NODE_FIELDS)
    write_tsv(inputs / "read_flow_links.input.tsv", links, LINK_FIELDS)
    write_tsv(inputs / "read_flow_views.tsv", [{"view_id": "sample_reads", "display_name": f"{sample} reads flow", "order": 0, "is_default": 1}], ["view_id", "display_name", "order", "is_default"])
    subprocess.run([sys.executable, "pipeline/scripts/build_read_flow_sankey.py", "--nodes", str(inputs / "read_flow_nodes.input.tsv"), "--links", str(inputs / "read_flow_links.input.tsv"), "--views", str(inputs / "read_flow_views.tsv"), "--default-mode", "view", "--node-layout", "main_top_losses", "--plotly-js", str(args.plotly_js), "--output-dir", str(args.outdir), "--title", f"{sample} conserved read flow"], check=True)
    subprocess.run([sys.executable, "pipeline/scripts/validate_read_flow_outputs.py", "--output-dir", str(args.outdir)], check=True)
    write_json(args.outdir / "reports/read_flow_adapter_validation.json", {"ok": True, "status": "PASS", "schema_version": "iso-slam-read-flow-v6", "sample_id": sample, "counts": values, "unit": "reads", "fractional_expected_counts_in_sankey": False})
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
