#!/usr/bin/env python3
"""Validate an offline report by parsing HTML tags rather than script text."""

from __future__ import annotations

import argparse
import base64
import csv
import json
from html.parser import HTMLParser
from pathlib import Path


class DeliveryHTMLParser(HTMLParser):
    def __init__(self) -> None:
        super().__init__(convert_charrefs=True)
        self.hrefs: list[str] = []
        self.srcs: list[str] = []
        self.ids: list[str] = []
        self.classes: set[str] = set()
        self.tables: list[dict[str, object]] = []
        self._table: dict[str, object] | None = None

    def handle_starttag(self, tag: str, attrs: list[tuple[str, str | None]]) -> None:
        values = {key: value or "" for key, value in attrs}
        if values.get("href"):
            self.hrefs.append(values["href"])
        if values.get("src"):
            self.srcs.append(values["src"])
        if values.get("id"):
            self.ids.append(values["id"])
        if values.get("class"):
            self.classes.update(values["class"].split())
        if tag.lower() == "table":
            self._table = {"rows": 0, "paginated": bool(values.get("data-page-size"))}
        elif tag.lower() == "tr" and self._table is not None:
            self._table["rows"] = int(self._table["rows"]) + 1

    def handle_endtag(self, tag: str) -> None:
        if tag.lower() == "table" and self._table is not None:
            self.tables.append(self._table)
            self._table = None


def add_error(result: dict[str, object], message: str) -> None:
    errors = result.setdefault("errors", [])
    assert isinstance(errors, list)
    errors.append(message)


def read_manifest(path: Path) -> list[dict[str, str]]:
    with path.open(encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def validate_data_uris(srcs: list[str], result: dict[str, object]) -> None:
    uris = [src for src in srcs if src.startswith("data:")]
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


def validate_manifest(report: Path, manifest: Path, result: dict[str, object]) -> None:
    if not manifest.is_file() or manifest.stat().st_size == 0:
        add_error(result, f"manifest missing or empty: {manifest}")
        return
    rows = read_manifest(manifest)
    result["manifest_rows"] = len(rows)
    for row in rows:
        if row.get("status") not in {"copied", "up_to_date"}:
            continue
        target = Path(row.get("target", ""))
        if not target.is_absolute():
            target = manifest.parent / target
        if not target.is_file() or target.stat().st_size == 0:
            add_error(result, f"manifest target missing for {row.get('description', 'unnamed item')}: {target}")
    report_rows = []
    for row in rows:
        target = Path(row.get("target", ""))
        if not target.is_absolute():
            target = manifest.parent / target
        if target.resolve() == report:
            report_rows.append(row)
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
    text = ""
    if not report.is_file() or report.stat().st_size == 0:
        add_error(result, f"report missing or empty: {report}")
    else:
        text = report.read_text(encoding="utf-8", errors="strict")
        html_parser = DeliveryHTMLParser()
        html_parser.feed(text)
        external_hrefs = [href for href in html_parser.hrefs if not href.startswith(("#", "data:"))]
        external_srcs = [src for src in html_parser.srcs if not src.startswith("data:")]
        result["external_hrefs"] = len(external_hrefs)
        result["external_srcs"] = len(external_srcs)
        if external_hrefs:
            add_error(result, f"report contains external or file href values: {external_hrefs[:3]}")
        if external_srcs:
            add_error(result, f"report contains non-data src values: {external_srcs[:3]}")
        ids = set(html_parser.ids)
        duplicates = sorted({value for value in html_parser.ids if html_parser.ids.count(value) > 1})
        if duplicates:
            add_error(result, f"report contains duplicate HTML IDs: {duplicates[:5]}")
        stale = [href for href in html_parser.hrefs if href.startswith("#") and href[1:] not in ids]
        result["stale_anchor_count"] = len(stale)
        if stale:
            add_error(result, f"report contains stale anchors: {', '.join(stale[:5])}")
        if "sidebar" not in html_parser.classes:
            add_error(result, "report sidebar is missing")
        if "footer" not in html_parser.classes:
            add_error(result, "report footer is missing")
        validate_data_uris(html_parser.srcs, result)
        result["table_count"] = len(html_parser.tables)
        for index, table in enumerate(html_parser.tables, start=1):
            data_rows = max(0, int(table["rows"]) - 1)
            if data_rows > args.max_unpaginated_rows and not bool(table["paginated"]):
                add_error(result, f"table {index} has {data_rows} rows without pagination")
        for section in args.require_section:
            if section not in text:
                add_error(result, f"required section is missing: {section}")
        for term in args.forbid_text:
            if term and term in text:
                add_error(result, f"forbidden text found in report: {term}")
    if args.manifest:
        manifest = args.manifest.resolve()
        validate_manifest(report, manifest, result)
        if manifest.is_file():
            manifest_text = manifest.read_text(encoding="utf-8", errors="strict")
            for term in args.forbid_text:
                if term and term in manifest_text:
                    add_error(result, f"forbidden text found in manifest: {term}")
    result["ok"] = not result["errors"]
    if args.as_json:
        print(json.dumps(result, ensure_ascii=False, indent=2))
    else:
        print(f"report: {result['report']}")
        print(f"ok: {result['ok']}")
        for error in result["errors"]:  # type: ignore[index]
            print(f"ERROR: {error}")
    return 0 if result["ok"] else 1


if __name__ == "__main__":
    raise SystemExit(main())
