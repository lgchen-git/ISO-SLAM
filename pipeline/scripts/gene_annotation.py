#!/usr/bin/env python3
"""Helpers for adding GRCh38 gene names to reports and delivery tables."""

from __future__ import annotations

import gzip
import re
from pathlib import Path


ATTRIBUTE_RE = re.compile(r'(\w+) "([^"]+)"')


def load_gtf_gene_map(path: Path) -> dict[str, str]:
    mapping: dict[str, str] = {}
    with gzip.open(path, "rt", encoding="utf-8") as handle:
        for line in handle:
            if line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) != 9 or fields[2] != "gene":
                continue
            attrs = dict(ATTRIBUTE_RE.findall(fields[8]))
            gene_id = attrs.get("gene_id", "")
            if gene_id:
                mapping.setdefault(gene_id, attrs.get("gene_name") or gene_id)
    if not mapping:
        raise RuntimeError(f"No gene annotations were read from {path}")
    return mapping


def add_gene_names(rows: list[dict[str, object]], mapping: dict[str, str]) -> None:
    for row in rows:
        gene_id = str(row.get("gene_id", ""))
        if gene_id:
            row["gene_name"] = mapping.get(gene_id, gene_id)


def fields_with_gene_name(fields: list[str]) -> list[str]:
    if "gene_id" not in fields:
        return fields
    updated = list(fields)
    if "gene_name" in updated:
        updated.remove("gene_name")
    updated.insert(updated.index("gene_id") + 1, "gene_name")
    return updated
