#!/usr/bin/env python3
"""Validate an offline bioinformatics HTML report and optional delivery manifest."""

from __future__ import annotations

import argparse
import base64
import csv
import json
import re
from pathlib import Path


def add_error(result: dict, message: str) -> None:
    result["errors"].append(message)


def read_manifest(path: Path) -> list[dict[str, str]]:
    with path.open(encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def validate_data_uris(text: str, result: dict) -> None:
    uris = re.findall(r'\bsrc="(data:[^"]+)"', text)
    result["data_uri_count"] = len(uris)
    for index, uri in enumerate(uris, start=1):
        if ";base64," not in uri:
            add_error(result, f"data URI {index} is not base64 encoded")
            continue
        payload = uri.split(";base64,", 1)[1]
        try:
            decoded = base64.b64decode(payload, validate=True)
        except ValueError:
            add_error(result, f"data URI {index} has invalid base64")
            continue
        if not decoded:
            add_error(result, f"data URI {index} is empty")


def validate_tables(text: str, max_rows: int, result: dict) -> None:
    tables = re.findall(r"<table\b([^>]*)>(.*?)</table>", text, flags=re.IGNORECASE | re.DOTALL)
    result["table_count"] = len(tables)
    for index, (attrs, body) in enumerate(tables, start=1):
        rows = len(re.findall(r"<tr\b", body, flags=re.IGNORECASE)) - 1
        if rows > max_rows and "data-page-size=" not in attrs:
            add_error(result, f"table {index} has {rows} rows without pagination")


def validate_manifest(report: Path, manifest: Path, result: dict) -> None:
    if not manifest.is_file() or manifest.stat().st_size == 0:
        add_error(result, f"manifest missing or empty: {manifest}")
        return
    rows = read_manifest(manifest)
    result["manifest_rows"] = len(rows)
    for row in rows:
        if row.get("status") not in {"copied", "up_to_date"}:
            continue
        target = row.get("target", "")
        if not target or not Path(target).is_file() or Path(target).stat().st_size == 0:
            add_error(result, f"manifest target missing for {row.get('description', 'unnamed item')}: {target}")
    report_rows = [row for row in rows if Path(row.get("target", "")) == report]
    if len(report_rows) != 1:
        add_error(result, f"manifest must contain exactly one final-report target; found {len(report_rows)}")


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--report", required=True, type=Path)
    parser.add_argument("--manifest", type=Path)
    parser.add_argument("--require-section", action="append", default=[])
    parser.add_argument("--forbid-text", action="append", default=[])
    parser.add_argument("--max-unpaginated-rows", type=int, default=12)
    parser.add_argument("--json", action="store_true", dest="as_json")
    args = parser.parse_args()

    report = args.report.resolve()
    result: dict[str, object] = {"report": str(report), "errors": []}
    if not report.is_file() or report.stat().st_size == 0:
        add_error(result, f"report missing or empty: {report}")
    else:
        text = report.read_text(encoding="utf-8", errors="replace")
        hrefs = re.findall(r'\bhref="([^"]+)"', text)
        srcs = re.findall(r'\bsrc="([^"]+)"', text)
        result["external_hrefs"] = sum(not href.startswith("#") for href in hrefs)
        result["external_srcs"] = sum(not src.startswith("data:") for src in srcs)
        if result["external_hrefs"]:
            add_error(result, "report contains external or file href values")
        if result["external_srcs"]:
            add_error(result, "report contains non-data src values")
        ids = set(re.findall(r'\bid="([^"]+)"', text))
        stale = [href for href in hrefs if href.startswith("#") and href[1:] not in ids]
        result["stale_anchor_count"] = len(stale)
        if stale:
            add_error(result, f"report contains stale anchors: {', '.join(stale[:5])}")
        if 'class="sidebar"' not in text:
            add_error(result, "report sidebar is missing")
        if 'class="footer"' not in text:
            add_error(result, "report footer is missing")
        validate_data_uris(text, result)
        validate_tables(text, args.max_unpaginated_rows, result)
        for section in args.require_section:
            if section not in text:
                add_error(result, f"required section is missing: {section}")
        for term in args.forbid_text:
            if term and term in text:
                add_error(result, f"forbidden text found in report: {term}")
    if args.manifest:
        validate_manifest(report, args.manifest.resolve(), result)
        if args.manifest.is_file():
            manifest_text = args.manifest.read_text(encoding="utf-8", errors="replace")
            for term in args.forbid_text:
                if term and term in manifest_text:
                    add_error(result, f"forbidden text found in manifest: {term}")
    result["ok"] = not result["errors"]
    if args.as_json:
        print(json.dumps(result, ensure_ascii=False, indent=2))
    else:
        print(f"report: {result['report']}")
        print(f"ok: {result['ok']}")
        for error in result["errors"]:
            print(f"ERROR: {error}")
    return 0 if result["ok"] else 1


if __name__ == "__main__":
    raise SystemExit(main())
