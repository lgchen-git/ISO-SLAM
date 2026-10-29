#!/usr/bin/env python3
"""Build one conserved long-read ISO-SLAM processing-flow module."""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import re
import subprocess
import sys
from pathlib import Path


NODE_FIELDS = ["view_id", "node_id", "stage", "label", "value", "unit", "scope", "node_kind"]
LINK_FIELDS = ["view_id", "link_id", "source", "target", "value", "reason"]
VIEW_FIELDS = ["view_id", "display_name", "order", "is_default"]
SCHEMA_VERSION = "iso-slam-read-flow-v1"


class FlowInputError(RuntimeError):
    """Raised when source metrics do not define a conserved processing flow."""


def arguments() -> argparse.Namespace:
    script_dir = Path(__file__).resolve().parent
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--sample-id", required=True)
    parser.add_argument("--sample-root", required=True, type=Path)
    parser.add_argument("--output-dir", required=True, type=Path)
    parser.add_argument("--renderer", type=Path, default=script_dir / "build_read_flow_sankey.py")
    parser.add_argument("--validator", type=Path, default=script_dir / "validate_read_flow_outputs.py")
    parser.add_argument("--plotly-js", type=Path, default=script_dir.parent / "assets" / "plotly-3.5.1.min.js")
    return parser.parse_args()


def read_json(path: Path) -> dict[str, object]:
    if not path.is_file() or path.stat().st_size == 0:
        raise FlowInputError(f"Missing or empty source file: {path}")
    value = json.loads(path.read_text(encoding="utf-8"))
    if not isinstance(value, dict):
        raise FlowInputError(f"Expected a JSON object: {path}")
    return value


def integer(source: dict[str, object], key: str, source_name: str) -> int:
    value = source.get(key)
    if isinstance(value, bool):
        raise FlowInputError(f"{source_name}: {key} must be a non-negative integer")
    try:
        number = int(value)  # type: ignore[arg-type]
    except (TypeError, ValueError) as exc:
        raise FlowInputError(f"{source_name}: missing or invalid integer {key}") from exc
    if number < 0 or str(value).strip() != str(number):
        raise FlowInputError(f"{source_name}: {key} must be a non-negative integer")
    return number


def parse_flagstat(path: Path) -> dict[str, int]:
    if not path.is_file() or path.stat().st_size == 0:
        raise FlowInputError(f"Missing or empty source file: {path}")
    metrics: dict[str, int] = {}
    for line in path.read_text(encoding="utf-8", errors="strict").splitlines():
        match = re.match(r"^(\d+) \+ \d+ (.+)$", line.strip())
        if not match:
            continue
        count = int(match.group(1))
        label = match.group(2)
        if label.startswith("in total"):
            metrics["total"] = count
        elif label == "primary":
            metrics["primary"] = count
        elif label.startswith("primary mapped"):
            metrics["primary_mapped"] = count
    missing = [key for key in ("total", "primary", "primary_mapped") if key not in metrics]
    if missing:
        raise FlowInputError(f"{path}: missing flagstat metrics: {', '.join(missing)}")
    return metrics


def require_equal(label: str, left: int, right: int) -> None:
    if left != right:
        raise FlowInputError(f"Conservation failure for {label}: {left} != {right}")


def write_tsv(path: Path, fields: list[str], rows: list[dict[str, object]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def source_paths(sample_root: Path, sample_id: str) -> dict[str, Path]:
    return {
        "trim_summary": sample_root / "02_trim" / "04_tool_output" / "trim" / f"{sample_id}.trim_summary.json",
        "alignment_flagstat": sample_root / "03_alignment" / "04_tool_output" / "minimap2_samtools" / f"{sample_id}.flagstat.txt",
        "primary_q20_flagstat": sample_root / "03_alignment" / "04_tool_output" / "minimap2_samtools" / f"{sample_id}.primary.flagstat.txt",
        "slam_qc": sample_root / "04_slam_calling" / "04_tool_output" / "caller" / f"{sample_id}.bulk_slam_qc.json",
    }


def contract(sample_id: str, paths: dict[str, Path]) -> tuple[list[dict[str, object]], list[dict[str, object]], list[dict[str, object]], dict[str, int]]:
    trim = read_json(paths["trim_summary"])
    qc = read_json(paths["slam_qc"])
    aligned = parse_flagstat(paths["alignment_flagstat"])
    q20 = parse_flagstat(paths["primary_q20_flagstat"])

    raw = integer(trim, "total_reads", "trim summary")
    insert = integer(trim, "trimmed_reads", "trim summary")
    structure_loss = integer(trim, "unassigned_reads", "trim summary")
    primary_mapped = aligned["primary_mapped"]
    unmapped = aligned["primary"] - primary_mapped
    primary_unrepresented = insert - aligned["primary"]
    primary_q20 = q20["total"]
    low_mapq = primary_mapped - primary_q20
    assigned = integer(qc, "assigned_reads", "SLAM QC")
    intergenic = integer(qc, "assignment_intergenic_reads", "SLAM QC")
    shared_antisense = integer(qc, "assignment_shared_antisense_ambiguous_reads", "SLAM QC")
    shared_gene = integer(qc, "assignment_shared_gene_ambiguous_reads", "SLAM QC")
    antisense_only = integer(qc, "assignment_antisense_only_reads", "SLAM QC")
    unassigned_gene = integer(qc, "unassigned_gene", "SLAM QC")
    caller_total = integer(qc, "total_reads", "SLAM QC")

    if unmapped < 0 or low_mapq < 0 or primary_unrepresented < 0:
        raise FlowInputError("Alignment-derived loss counts must be non-negative")
    require_equal("raw read partition", raw, insert + structure_loss)
    require_equal("primary mapping partition", aligned["primary"], primary_mapped + unmapped)
    require_equal("MAPQ partition", primary_mapped, primary_q20 + low_mapq)
    require_equal("Q20 BAM total versus primary", q20["total"], q20["primary"])
    require_equal("Q20 BAM total versus primary mapped", q20["total"], q20["primary_mapped"])
    require_equal("caller denominator versus Q20 reads", caller_total, primary_q20)
    require_equal("unassigned gene fates", unassigned_gene, intergenic + shared_antisense + shared_gene + antisense_only)
    require_equal("gene assignment partition", primary_q20, assigned + unassigned_gene)

    values = {
        "raw": raw,
        "insert": insert,
        "structure_loss": structure_loss,
        "primary_unrepresented": primary_unrepresented,
        "primary_mapped": primary_mapped,
        "unmapped": unmapped,
        "primary_q20": primary_q20,
        "low_mapq": low_mapq,
        "gene_assigned": assigned,
        "intergenic": intergenic,
        "shared_antisense_ambiguous": shared_antisense,
        "shared_gene_ambiguous": shared_gene,
        "antisense_only": antisense_only,
    }
    labels = {
        "raw": "Raw reads",
        "insert": "Insert reads",
        "structure_loss": "文库结构未识别",
        "primary_unrepresented": "无 primary 对齐记录",
        "primary_mapped": "Primary mapped reads",
        "unmapped": "Unmapped reads",
        "primary_q20": "Primary MAPQ≥20 reads",
        "low_mapq": "Primary MAPQ<20 reads",
        "gene_assigned": "基因分配 reads",
        "intergenic": "Intergenic",
        "shared_antisense_ambiguous": "Shared antisense ambiguous",
        "shared_gene_ambiguous": "Shared gene ambiguous",
        "antisense_only": "Antisense only",
    }
    stages = {
        "raw": 0,
        "insert": 1,
        "structure_loss": 1,
        "primary_unrepresented": 2,
        "primary_mapped": 2,
        "unmapped": 2,
        "primary_q20": 3,
        "low_mapq": 3,
        "gene_assigned": 4,
        "intergenic": 4,
        "shared_antisense_ambiguous": 4,
        "shared_gene_ambiguous": 4,
        "antisense_only": 4,
    }
    retained = {"raw", "insert", "primary_mapped", "primary_q20", "gene_assigned"}
    nodes = [
        {
            "view_id": "sample_processing",
            "node_id": node_id,
            "stage": stages[node_id],
            "label": labels[node_id],
            "value": value,
            "unit": "reads",
            "scope": "sample_library",
            "node_kind": "retained" if node_id in retained else "loss",
        }
        for node_id, value in values.items()
    ]
    link_specs = [
        ("raw_to_insert", "raw", "insert", "成功识别文库结构并提取 insert"),
        ("raw_to_structure_loss", "raw", "structure_loss", "未识别完整文库结构或 insert 不满足要求"),
        ("insert_to_primary_unrepresented", "insert", "primary_unrepresented", "未产生 primary 对齐记录"),
        ("insert_to_primary_mapped", "insert", "primary_mapped", "minimap2 primary alignment mapped"),
        ("insert_to_unmapped", "insert", "unmapped", "无 primary mapped alignment"),
        ("mapped_to_q20", "primary_mapped", "primary_q20", "primary alignment MAPQ≥20"),
        ("mapped_to_low_mapq", "primary_mapped", "low_mapq", "primary alignment MAPQ<20"),
        ("q20_to_assigned", "primary_q20", "gene_assigned", "分配到唯一基因或通过 gene-level salvage"),
        ("q20_to_intergenic", "primary_q20", "intergenic", "无基因重叠"),
        ("q20_to_shared_antisense", "primary_q20", "shared_antisense_ambiguous", "多个反义候选基因且无法唯一消歧"),
        ("q20_to_shared_gene", "primary_q20", "shared_gene_ambiguous", "多个候选基因且无法唯一消歧"),
        ("q20_to_antisense_only", "primary_q20", "antisense_only", "仅反义重叠且未被 salvage"),
    ]
    links = [
        {
            "view_id": "sample_processing",
            "link_id": link_id,
            "source": source,
            "target": target,
            "value": values[target],
            "reason": reason,
        }
        for link_id, source, target, reason in link_specs
    ]
    views = [{"view_id": "sample_processing", "display_name": f"{sample_id} reads处理流", "order": 0, "is_default": 1}]
    return nodes, links, views, values


def main() -> int:
    cli = arguments()
    try:
        paths = source_paths(cli.sample_root.resolve(), cli.sample_id)
        nodes, links, views, values = contract(cli.sample_id, paths)
        inputs = cli.output_dir / "inputs"
        nodes_path = inputs / "read_flow_nodes.input.tsv"
        links_path = inputs / "read_flow_links.input.tsv"
        views_path = inputs / "read_flow_views.tsv"
        write_tsv(nodes_path, NODE_FIELDS, nodes)
        write_tsv(links_path, LINK_FIELDS, links)
        write_tsv(views_path, VIEW_FIELDS, views)

        subprocess.run(
            [
                sys.executable,
                str(cli.renderer),
                "--nodes", str(nodes_path),
                "--links", str(links_path),
                "--views", str(views_path),
                "--default-mode", "view",
                "--node-layout", "main_top_losses",
                "--plotly-js", str(cli.plotly_js),
                "--output-dir", str(cli.output_dir),
                "--title", f"{cli.sample_id} reads处理流",
            ],
            check=True,
        )
        subprocess.run([sys.executable, str(cli.validator), "--output-dir", str(cli.output_dir)], check=True)

        tables = cli.output_dir / "tables"
        reports = cli.output_dir / "reports"
        source_rows = [
            {
                "source": name,
                "path": str(path.resolve()),
                "bytes": path.stat().st_size,
                "sha256": sha256(path),
            }
            for name, path in paths.items()
        ]
        write_tsv(tables / "read_flow_source_manifest.tsv", ["source", "path", "bytes", "sha256"], source_rows)
        validation = {
            "ok": True,
            "schema_version": SCHEMA_VERSION,
            "sample_id": cli.sample_id,
            "view_order": ["sample_processing"],
            "unit": "reads",
            "scope": "sample_library",
            "transformations": 0,
            "unresolved_residuals": 0,
            "counts": values,
            "source_manifest": str((tables / "read_flow_source_manifest.tsv").resolve()),
        }
        reports.mkdir(parents=True, exist_ok=True)
        (reports / "read_flow_adapter_validation.json").write_text(
            json.dumps(validation, ensure_ascii=False, indent=2) + "\n", encoding="utf-8"
        )
        print(json.dumps(validation, ensure_ascii=False))
        return 0
    except (FlowInputError, OSError, ValueError, json.JSONDecodeError, subprocess.CalledProcessError) as exc:
        print(json.dumps({"ok": False, "error": str(exc)}, ensure_ascii=False), file=sys.stderr)
        return 2


if __name__ == "__main__":
    raise SystemExit(main())
