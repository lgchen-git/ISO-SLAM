#!/usr/bin/env python3
"""Add configured GTF gene names to selected tabular results."""

from __future__ import annotations

import argparse
import csv
import json
import os
import tempfile
from pathlib import Path

from gene_annotation import fields_with_gene_name, load_gtf_gene_map


def annotate(path: Path, mapping: dict[str, str]) -> tuple[int, int]:
    with path.open(encoding="utf-8", newline="") as source:
        reader = csv.DictReader(source, delimiter="\t")
        fields = reader.fieldnames
        if not fields or "gene_id" not in fields:
            raise ValueError(f"{path.name} has no gene_id column")
        output_fields = fields_with_gene_name(fields)
        rows = list(reader)
    mapped = 0
    unmapped = 0
    for row in rows:
        gene_id = row.get("gene_id", "")
        if gene_id:
            row["gene_name"] = mapping.get(gene_id, gene_id)
            mapped += gene_id in mapping
            unmapped += gene_id not in mapping
    fd, temporary = tempfile.mkstemp(prefix=f".{path.name}.", suffix=".tmp", dir=path.parent)
    try:
        with os.fdopen(fd, "w", encoding="utf-8", newline="") as target:
            writer = csv.DictWriter(target, fieldnames=output_fields, delimiter="\t", lineterminator="\n", extrasaction="ignore")
            writer.writeheader()
            writer.writerows(rows)
            target.flush()
            os.fsync(target.fileno())
        os.replace(temporary, path)
    finally:
        if os.path.exists(temporary):
            os.unlink(temporary)
    return len(rows), mapped, unmapped


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--gtf", required=True, type=Path)
    parser.add_argument("--tables", required=True, nargs="+", type=Path)
    parser.add_argument("--validation", required=True, type=Path)
    args = parser.parse_args()
    mapping = load_gtf_gene_map(args.gtf)
    summaries = []
    for path in args.tables:
        rows, mapped, unmapped = annotate(path, mapping)
        summaries.append({"file": path.name, "rows": rows, "mapped_gene_rows": mapped, "unmapped_gene_rows": unmapped})
    result = {"status": "PASS", "table_count": len(summaries), "row_count": sum(x["rows"] for x in summaries), "tables": summaries}
    args.validation.parent.mkdir(parents=True, exist_ok=True)
    temporary = args.validation.with_name(f".{args.validation.name}.tmp")
    temporary.write_text(json.dumps(result, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")
    os.replace(temporary, args.validation)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
