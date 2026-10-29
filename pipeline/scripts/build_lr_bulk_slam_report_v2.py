#!/usr/bin/env python3
"""Build a Chinese self-contained bulk long-read SLAM-seq HTML report."""

from __future__ import annotations

import argparse
import csv
import gzip
import html
import json
import math
import re
from pathlib import Path


COLORS = ["#E64B35", "#4DBBD5", "#00A087", "#3C5488", "#F39B7F", "#8491B4"]
NPG_RED = "#E64B35"
NPG_BLUE = "#3C5488"
NPG_GREEN = "#00A087"
NPG_CYAN = "#4DBBD5"
NEUTRAL = "#94A3B8"
QUALITY_THRESHOLDS = [10, 15, 20, 25, 30]


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Build bulk SLAM-seq report.")
    parser.add_argument("--results-dir", required=True)
    parser.add_argument("--sample-id", required=True)
    parser.add_argument("--project-name", required=True)
    parser.add_argument("--report-title", required=True)
    parser.add_argument("--config-path", required=True)
    parser.add_argument("--samples-path", required=True)
    parser.add_argument("--report-stage", type=Path)
    parser.add_argument("--read-flow-dir", type=Path)
    parser.add_argument(
        "--plotly-js",
        type=Path,
        default=Path(__file__).resolve().parents[1] / "assets" / "plotly-3.5.1.min.js",
    )
    return parser.parse_args()


def read_json(path: Path) -> dict[str, object]:
    if not path.is_file():
        return {}
    return json.loads(path.read_text(encoding="utf-8"))


def read_config_scalar(path: Path, key: str, default: str = "") -> str:
    if not path.is_file():
        return default
    with path.open("r", encoding="utf-8") as handle:
        for raw_line in handle:
            line = raw_line.rstrip("\n")
            stripped = line.strip()
            if not stripped or stripped.startswith("#"):
                continue
            if raw_line[:1].isspace() or ":" not in line:
                continue
            name, raw_value = line.split(":", 1)
            if name.strip() != key:
                continue
            value = raw_value.split("#", 1)[0].strip()
            if len(value) >= 2 and value[0] == value[-1] and value[0] in {"'", '"'}:
                value = value[1:-1]
            return value or default
    return default


def read_config_int(path: Path, key: str, default: int) -> int:
    try:
        return int(read_config_scalar(path, key, str(default)))
    except ValueError:
        return default


def read_config_list_int(path: Path, key: str, default: list[int]) -> list[int]:
    raw = read_config_scalar(path, key, "")
    if not raw:
        return default
    values = []
    for item in raw.split(","):
        item = item.strip()
        if item:
            values.append(int(item))
    return values or default


def read_top_genes(path: Path, limit: int = 30) -> list[dict[str, str]]:
    if not path.is_file():
        return []
    with path.open("r", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        return [row for _, row in zip(range(limit), reader)]


def read_tsv_rows(path: Path, limit: int | None = None) -> list[dict[str, str]]:
    if not path.is_file():
        return []
    opener = gzip.open if path.suffix == ".gz" else open
    rows = []
    with opener(path, "rt", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        for row in reader:
            rows.append(row)
            if limit is not None and len(rows) >= limit:
                break
    return rows


def read_first_tsv_row(path: Path) -> dict[str, str]:
    rows = read_tsv_rows(path, limit=1)
    return rows[0] if rows else {}


def read_sample_group(path: Path, sample_id: str) -> str:
    for row in read_tsv_rows(path):
        if row.get("sample_id") == sample_id:
            return row.get("group", "")
    return ""


def read_flagstat_total(path: Path) -> int:
    if not path.is_file():
        return 0
    for line in path.read_text(encoding="utf-8").splitlines():
        if " in total " in line:
            return int(line.split()[0])
    return 0


def fmt_int(value: object) -> str:
    try:
        return f"{int(float(value)):,}"
    except (TypeError, ValueError):
        return "0"


def fmt_float(value: object, digits: int = 1) -> str:
    try:
        return f"{float(value):,.{digits}f}"
    except (TypeError, ValueError):
        return "0.0"


def fmt_pct(value: object) -> str:
    try:
        return f"{float(value) * 100:.2f}%"
    except (TypeError, ValueError):
        return "0.00%"


def fmt_chart_value(value: object) -> str:
    try:
        number = float(value)
    except (TypeError, ValueError):
        return "0"
    if not math.isfinite(number):
        return "0"
    if abs(number) >= 1000:
        return fmt_int(number)
    if number.is_integer():
        return str(int(number))
    return f"{number:.2f}"


def esc(value: object) -> str:
    return html.escape(str(value))


def float_or_zero(value: object) -> float:
    try:
        return float(value)
    except (TypeError, ValueError):
        return 0.0


def int_or_zero(value: object) -> int:
    try:
        return int(float(value))
    except (TypeError, ValueError):
        return 0


def svg_header(width: int = 900, height: int = 360) -> str:
    return f'<svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 {width} {height}" role="img">'


def bar_chart(title: str, data: list[tuple[str, float]], y_label: str = "") -> str:
    width, height = 900, 380
    left, right, top, bottom = 72, 28, 58, 92
    plot_w = width - left - right
    plot_h = height - top - bottom
    max_value = max([max(float(value), 0.0) for _label, value in data] + [1.0])
    bar_w = plot_w / max(len(data), 1) * 0.58
    gap = plot_w / max(len(data), 1)
    parts = [
        svg_header(width, height),
        '<rect width="100%" height="100%" fill="#ffffff"/>',
        f'<text x="{width / 2}" y="30" text-anchor="middle" class="svg-title">{esc(title)}</text>',
        f'<line x1="{left}" y1="{top + plot_h}" x2="{width - right}" y2="{top + plot_h}" class="axis"/>',
        f'<line x1="{left}" y1="{top}" x2="{left}" y2="{top + plot_h}" class="axis"/>',
    ]
    if y_label:
        parts.append(
            f'<text x="20" y="{top + plot_h / 2}" transform="rotate(-90 20 {top + plot_h / 2})" '
            f'text-anchor="middle" class="svg-small">{esc(y_label)}</text>'
        )
    for idx, (label, value) in enumerate(data):
        numeric = max(float(value), 0.0)
        x = left + idx * gap + (gap - bar_w) / 2
        bar_h = numeric / max_value * plot_h
        y = top + plot_h - bar_h
        color = COLORS[idx % len(COLORS)]
        parts.append(f'<rect x="{x:.1f}" y="{y:.1f}" width="{bar_w:.1f}" height="{bar_h:.1f}" fill="{color}" rx="4"/>')
        parts.append(f'<text x="{x + bar_w / 2:.1f}" y="{max(y - 8, 48):.1f}" text-anchor="middle" class="svg-small">{fmt_chart_value(numeric)}</text>')
        parts.append(
            f'<text x="{x + bar_w / 2:.1f}" y="{top + plot_h + 26}" text-anchor="middle" '
            f'transform="rotate(35 {x + bar_w / 2:.1f} {top + plot_h + 26})" class="svg-small">{esc(label)}</text>'
        )
    parts.append("</svg>")
    return "\n".join(parts)


def grouped_metric_chart(title: str, metrics: list[tuple[str, str, float]]) -> str:
    return bar_chart(title, [(f"{label} {metric}", value) for label, metric, value in metrics], "")


def horizontal_ratio_chart(title: str, rows: list[dict[str, str]], limit: int = 20) -> str:
    width, height = 900, 600
    left, right, top, row_h = 160, 42, 56, 24
    rows = rows[:limit]
    parts = [
        svg_header(width, height),
        '<rect width="100%" height="100%" fill="#ffffff"/>',
        f'<text x="{width / 2}" y="30" text-anchor="middle" class="svg-title">{esc(title)}</text>',
    ]
    for idx, row in enumerate(rows):
        y = top + idx * row_h
        gene = row.get("gene_name") or row.get("gene_id") or "."
        total = int_or_zero(row.get("total_reads"))
        frac = float_or_zero(row.get("new_fraction"))
        old_frac = max(0.0, 1.0 - frac)
        bar_w = width - left - right
        parts.append(f'<text x="{left - 10}" y="{y + 14}" text-anchor="end" class="svg-small">{esc(gene[:18])}</text>')
        parts.append(f'<rect x="{left}" y="{y}" width="{bar_w * old_frac:.1f}" height="16" fill="#cbd5e1" rx="3"/>')
        parts.append(f'<rect x="{left + bar_w * old_frac:.1f}" y="{y}" width="{bar_w * frac:.1f}" height="16" fill="#2563eb" rx="3"/>')
        parts.append(f'<text x="{width - right + 6}" y="{y + 13}" class="svg-small">{fmt_pct(frac)} / {fmt_int(total)}</text>')
    parts.append("</svg>")
    return "\n".join(parts)


def funnel_chart(values: list[tuple[str, int]], title: str) -> str:
    width, height = 900, 390
    left, top = 220, 78
    max_value = max([value for _label, value in values] + [1])
    parts = [
        svg_header(width, height),
        '<rect width="100%" height="100%" fill="#ffffff"/>',
        f'<text x="{width / 2}" y="32" text-anchor="middle" class="svg-title">{esc(title)}</text>',
    ]
    for idx, (label, value) in enumerate(values):
        y = top + idx * 64
        w = 580 * value / max_value if max_value else 0
        pct = value / values[0][1] if values and values[0][1] else 0.0
        parts.append(f'<text x="{left - 16}" y="{y + 21}" text-anchor="end" class="svg-label">{esc(label)}</text>')
        parts.append(f'<rect x="{left}" y="{y}" width="580" height="30" fill="#e2e8f0" rx="5"/>')
        parts.append(f'<rect x="{left}" y="{y}" width="{w:.1f}" height="30" fill="{COLORS[idx % len(COLORS)]}" rx="5"/>')
        parts.append(f'<text x="{left + 12}" y="{y + 21}" class="svg-inbar">{fmt_int(value)}</text>')
        parts.append(f'<text x="820" y="{y + 21}" class="svg-label">{fmt_pct(pct)}</text>')
    parts.append("</svg>")
    return "\n".join(parts)


def composition_three_svg(new_reads: int, old_reads: int, ambiguous_reads: int) -> str:
    width, height = 900, 260
    total = max(new_reads + old_reads + ambiguous_reads, 1)
    bar_w = 700
    new_w = bar_w * new_reads / total
    old_w = bar_w * old_reads / total
    ambiguous_w = bar_w - new_w - old_w
    parts = [
        svg_header(width, height),
        '<rect width="100%" height="100%" fill="#ffffff"/>',
        '<text x="450" y="32" text-anchor="middle" class="svg-title">new / old / ambiguous 构成</text>',
        '<text x="92" y="98" class="svg-label">formal read calls</text>',
        f'<rect x="92" y="120" width="{new_w:.1f}" height="42" fill="{NPG_GREEN}" rx="6"/>',
        f'<rect x="{92 + new_w:.1f}" y="120" width="{old_w:.1f}" height="42" fill="{NEUTRAL}" rx="0"/>',
        f'<rect x="{92 + new_w + old_w:.1f}" y="120" width="{ambiguous_w:.1f}" height="42" fill="#8491B4" rx="6"/>',
        f'<text x="92" y="190" class="svg-label">new {fmt_int(new_reads)} ({fmt_pct(new_reads / total)})</text>',
        f'<text x="342" y="190" class="svg-label">old {fmt_int(old_reads)} ({fmt_pct(old_reads / total)})</text>',
        f'<text x="592" y="190" class="svg-label">ambiguous {fmt_int(ambiguous_reads)} ({fmt_pct(ambiguous_reads / total)})</text>',
    ]
    parts.append("</svg>")
    return "\n".join(parts)


def assignment_composition_svg(rows: list[dict[str, str]]) -> str:
    labels = {
        "transcript_supported": "Transcript-supported",
        "same_strand_gene_only": "Same-strand gene-only",
        "antisense_rescue": "Antisense rescue",
        "ambiguous_shared": "Ambiguous/shared",
        "unassigned_intergenic": "Unassigned/intergenic",
    }
    colors = {
        "transcript_supported": NPG_GREEN,
        "same_strand_gene_only": NPG_CYAN,
        "antisense_rescue": NPG_RED,
        "ambiguous_shared": "#8491B4",
        "unassigned_intergenic": "#CBD5E1",
    }
    grouped = [row for row in rows if row.get("level") == "report_group"]
    total = max([int_or_zero(row.get("denominator")) for row in grouped] + [0])
    width, height = 900, 310
    left, bar_y, bar_w, bar_h = 70, 82, 760, 44
    parts = [
        svg_header(width, height),
        '<rect width="100%" height="100%" fill="#ffffff"/>',
        '<text x="450" y="32" text-anchor="middle" class="svg-title">Transcript-aware gene assignment</text>',
        f'<rect x="{left}" y="{bar_y}" width="{bar_w}" height="{bar_h}" fill="#F1F5F9" rx="6"/>',
    ]
    x = float(left)
    for row in grouped:
        category = row.get("category", "")
        count = int_or_zero(row.get("reads"))
        fraction = count / total if total else 0.0
        segment_w = bar_w * fraction
        if segment_w > 0:
            parts.append(
                f'<rect x="{x:.1f}" y="{bar_y}" width="{segment_w:.1f}" height="{bar_h}" '
                f'fill="{colors.get(category, NEUTRAL)}"/>'
            )
        x += segment_w
    for idx, row in enumerate(grouped):
        category = row.get("category", "")
        count = int_or_zero(row.get("reads"))
        fraction = count / total if total else 0.0
        col = idx % 2
        row_idx = idx // 2
        lx = 76 + col * 410
        ly = 168 + row_idx * 42
        parts.append(f'<rect x="{lx}" y="{ly - 12}" width="14" height="14" rx="2" fill="{colors.get(category, NEUTRAL)}"/>')
        parts.append(
            f'<text x="{lx + 22}" y="{ly}" class="svg-label">{esc(labels.get(category, category))}: '
            f'{fmt_int(count)} ({fmt_pct(fraction)})</text>'
        )
    parts.append("</svg>")
    return "\n".join(parts)


def t_enrichment_svg(slam_qc: dict[str, object]) -> str:
    tc_rate = float_or_zero(slam_qc.get("global_tc_rate"))
    bg_rate = float_or_zero(slam_qc.get("bg_tc_rate"))
    fold = tc_rate / bg_rate if bg_rate else math.nan
    maximum = max(tc_rate, bg_rate, 1e-8)
    width, height = 900, 360
    left, right, top, bottom = 110, 70, 68, 78
    plot_h = height - top - bottom
    baseline = top + plot_h
    bar_w = 150
    parts = [
        svg_header(width, height),
        '<rect width="100%" height="100%" fill="#ffffff"/>',
        '<text x="450" y="32" text-anchor="middle" class="svg-title">T&gt;C signal versus diagnostic background</text>',
        f'<line x1="{left}" y1="{baseline}" x2="{width-right}" y2="{baseline}" class="axis"/>',
        f'<line x1="{left}" y1="{top}" x2="{left}" y2="{baseline}" class="axis"/>',
        f'<text x="24" y="{top + plot_h/2}" transform="rotate(-90 24 {top + plot_h/2})" text-anchor="middle" class="svg-label">% of covered T</text>',
    ]
    for idx, (label, value, color) in enumerate(
        (("Diagnostic background", bg_rate, NEUTRAL), ("Observed T>C", tc_rate, NPG_RED))
    ):
        x = 250 + idx * 300
        h = value / maximum * plot_h if maximum else 0.0
        y = baseline - h
        parts.append(f'<rect x="{x}" y="{y:.1f}" width="{bar_w}" height="{h:.1f}" fill="{color}" rx="5"/>')
        parts.append(f'<text x="{x + bar_w/2}" y="{max(y-10, 54):.1f}" text-anchor="middle" class="svg-label">{value*100:.4f}%</text>')
        parts.append(f'<text x="{x + bar_w/2}" y="{baseline+28}" text-anchor="middle" class="svg-label">{esc(label)}</text>')
    fold_label = f"{fold:.2f}× enrichment" if math.isfinite(fold) else "fold unavailable"
    parts.append(f'<text x="450" y="{height-18}" text-anchor="middle" class="svg-title" fill="{NPG_RED}">{esc(fold_label)}</text>')
    parts.append("</svg>")
    return "\n".join(parts)


def substitution_table(rows: list[dict[str, str]]) -> str:
    observed = [row for row in rows if row.get("role") == "observed"]
    if not observed:
        return "<p class='note'>替换频率汇总不可用。</p>"
    out = ["<table><thead><tr><th>替换</th><th class='num'>事件数</th><th class='num'>Covered T</th><th class='num'>频率</th></tr></thead><tbody>"]
    for row in observed:
        out.append(
            "<tr>"
            f"<td>{esc(row.get('substitution', ''))}</td>"
            f"<td class='num'>{fmt_int(row.get('count'))}</td>"
            f"<td class='num'>{fmt_int(row.get('denominator'))}</td>"
            f"<td class='num'>{fmt_pct(row.get('rate'))}</td>"
            "</tr>"
        )
    out.append("</tbody></table>")
    return "\n".join(out)


def histogram_chart(title: str, values: list[float], bins: list[tuple[str, float, float]], y_label: str = "reads") -> str:
    counts = []
    for label, start, end in bins:
        if label.startswith(">="):
            count = sum(1 for value in values if value >= start)
        else:
            count = sum(1 for value in values if start <= value < end)
        counts.append((label, count))
    return bar_chart(title, counts, y_label)


def scatter_chart(title: str, points: list[tuple[float, float]], x_label: str, y_label: str) -> str:
    width, height = 900, 420
    left, right, top, bottom = 76, 34, 58, 78
    plot_w = width - left - right
    plot_h = height - top - bottom
    if not points:
        points = [(0.0, 0.0)]
    max_x = max([x for x, _y in points] + [1.0])
    max_y = max([y for _x, y in points] + [1.0])
    sampled = points[:5000]
    parts = [
        svg_header(width, height),
        '<rect width="100%" height="100%" fill="#ffffff"/>',
        f'<text x="{width / 2}" y="30" text-anchor="middle" class="svg-title">{esc(title)}</text>',
        f'<line x1="{left}" y1="{top + plot_h}" x2="{width - right}" y2="{top + plot_h}" class="axis"/>',
        f'<line x1="{left}" y1="{top}" x2="{left}" y2="{top + plot_h}" class="axis"/>',
        f'<text x="{left + plot_w / 2}" y="{height - 24}" text-anchor="middle" class="svg-label">{esc(x_label)}</text>',
        f'<text x="22" y="{top + plot_h / 2}" transform="rotate(-90 22 {top + plot_h / 2})" text-anchor="middle" class="svg-label">{esc(y_label)}</text>',
    ]
    for x, y in sampled:
        px = left + (x / max_x) * plot_w if max_x else left
        py = top + plot_h - (y / max_y) * plot_h if max_y else top + plot_h
        parts.append(f'<circle cx="{px:.1f}" cy="{py:.1f}" r="2.2" fill="#2563eb" opacity="0.28"/>')
    parts.append(f'<text x="{left}" y="{top + plot_h + 20}" class="svg-small">0</text>')
    parts.append(f'<text x="{width - right}" y="{top + plot_h + 20}" text-anchor="end" class="svg-small">{fmt_chart_value(max_x)}</text>')
    parts.append(f'<text x="{left - 8}" y="{top + 4}" text-anchor="end" class="svg-small">{fmt_chart_value(max_y)}</text>')
    parts.append("</svg>")
    return "\n".join(parts)


def save_svg(figures_dir: Path, name: str, svg: str) -> None:
    figures_dir.mkdir(parents=True, exist_ok=True)
    (figures_dir / name).write_text(svg, encoding="utf-8")


def metric_card(label: str, value: str, note: str = "") -> str:
    return (
        "<div class='kpi'>"
        f"<div class='kpi-val'>{esc(value)}</div>"
        f"<div class='kpi-label'>{esc(label)}</div>"
        f"<div class='kpi-note'>{esc(note)}</div>"
        "</div>"
    )


def table_rows(genes: list[dict[str, str]]) -> str:
    rows = []
    for row in genes:
        rows.append(
            "<tr>"
            f"<td>{esc(row.get('gene_name', ''))}</td>"
            f"<td class='mono'>{esc(row.get('gene_id', ''))}</td>"
            f"<td class='num'>{fmt_int(row.get('total_reads'))}</td>"
            f"<td class='num'>{fmt_int(row.get('new_reads'))}</td>"
            f"<td class='num'>{fmt_int(row.get('old_reads'))}</td>"
            f"<td class='num'>{fmt_int(row.get('ambiguous_reads'))}</td>"
            f"<td class='num'>{fmt_pct(row.get('new_fraction'))}</td>"
            f"<td class='num'>{fmt_int(row.get('gene_reads'))}</td>"
            f"<td class='num'>{fmt_int(row.get('tc_events'))}</td>"
            f"<td class='num'>{fmt_pct(row.get('tc_rate'))}</td>"
            "</tr>"
        )
    return "\n".join(rows)


def parse_metric_number(text: str) -> float:
    clean = text.replace(",", "").strip()
    match = re.search(r"[-+]?\d+(?:\.\d+)?", clean)
    if not match:
        return 0.0
    return float(match.group(0))


def parse_threshold_line(value: str) -> dict[str, float]:
    count_match = re.search(r"^\s*([\d,]+)", value)
    frac_match = re.search(r"\(([\d.]+)%\)", value)
    bases_match = re.search(r"([\d.]+)\s*Mb", value)
    return {
        "count": float(count_match.group(1).replace(",", "")) if count_match else 0.0,
        "fraction": float(frac_match.group(1)) / 100.0 if frac_match else 0.0,
        "bases_mb": float(bases_match.group(1)) if bases_match else 0.0,
    }


def parse_nanostats(path: Path) -> dict[str, object]:
    metrics: dict[str, object] = {"thresholds": {}, "raw_rows": []}
    if not path.is_file():
        return metrics
    with path.open("r", encoding="utf-8", errors="ignore") as handle:
        for line in handle:
            stripped = line.strip()
            if not stripped:
                continue
            if "\t" in stripped:
                key, value = stripped.split("\t", 1)
            elif ":" in stripped:
                key, value = stripped.split(":", 1)
            else:
                continue
            key = key.strip()
            value = value.strip()
            metrics["raw_rows"].append((key, value))
            if key == "number_of_reads":
                metrics["number_of_reads"] = int(parse_metric_number(value))
            elif key == "number_of_bases":
                metrics["number_of_bases"] = parse_metric_number(value)
            elif key == "median_read_length":
                metrics["median_read_length"] = parse_metric_number(value)
            elif key == "mean_read_length":
                metrics["mean_read_length"] = parse_metric_number(value)
            elif key == "read_length_stdev":
                metrics["read_length_stdev"] = parse_metric_number(value)
            elif key == "n50":
                metrics["n50"] = parse_metric_number(value)
            elif key == "mean_qual":
                metrics["mean_qual"] = parse_metric_number(value)
            elif key == "median_qual":
                metrics["median_qual"] = parse_metric_number(value)
            elif key.startswith("Reads >Q"):
                threshold_match = re.search(r"Q(\d+)", key)
                if threshold_match:
                    threshold = threshold_match.group(1)
                    metrics["thresholds"][threshold] = parse_threshold_line(value)
    return metrics


def threshold_metric(metrics: dict[str, object], threshold: int, field: str) -> float:
    thresholds = metrics.get("thresholds", {})
    if not isinstance(thresholds, dict):
        return 0.0
    item = thresholds.get(str(threshold), {})
    if not isinstance(item, dict):
        return 0.0
    return float_or_zero(item.get(field))


def parse_internal_fastq_qc(path: Path) -> dict[str, object]:
    """Normalize the workflow's deterministic FASTQ QC to the report metric schema."""
    payload = read_json(path)
    if not payload:
        return {"thresholds": {}, "raw_rows": []}
    return {
        "number_of_reads": int_or_zero(payload.get("read_count")),
        "number_of_bases": int_or_zero(payload.get("base_count")),
        "median_read_length": float_or_zero(payload.get("median_length")),
        "mean_read_length": float_or_zero(payload.get("mean_length")),
        "n50": float_or_zero(payload.get("n50")),
        "mean_qual": float_or_zero(payload.get("mean_read_q")),
        "median_qual": float_or_zero(payload.get("median_read_q")),
        "thresholds": {
            "20": {
                "fraction": float_or_zero(payload.get("q20_base_fraction")),
                "basis": "bases",
            }
        },
        "raw_rows": [],
    }


def read_nanoplot_rows(qc_dir: Path, sample_id: str) -> list[dict[str, object]]:
    rows: list[dict[str, object]] = []
    for label in ("raw", "insert"):
        status_path = qc_dir / "nanoplot" / label / "nanoplot_status.json"
        status = read_json(status_path)
        internal_path = qc_dir / "fastq" / f"{sample_id}.{label}_fastq_qc.json"
        internal_metrics = parse_internal_fastq_qc(internal_path)
        stats_path = Path(str(status.get("nanostats", ""))) if status else Path()
        report_path = Path(str(status.get("report", ""))) if status else Path()
        if internal_path.is_file() and internal_metrics.get("number_of_reads", 0):
            metrics = internal_metrics
            display_status = "内置 FASTQ QC"
            metric_source = "internal_qc"
        elif status.get("status") == "complete":
            metrics = parse_nanostats(stats_path)
            display_status = "NanoPlot"
            metric_source = "nanoplot"
        elif status:
            metrics = {"thresholds": {}, "raw_rows": []}
            display_status = str(status.get("status", "missing"))
            metric_source = "unavailable"
        else:
            continue
        rows.append(
            {
                "label": label,
                "status": display_status,
                "metric_source": metric_source,
                "nanoplot_status": status.get("status", "not_run"),
                "version": status.get("version", ""),
                "reason": status.get("reason", ""),
                "metrics": metrics,
                "stats_path": stats_path,
                "report_path": report_path,
                "report_available": bool(report_path.name and report_path.is_file()),
            }
        )
    return rows


def nanoplot_row_map(rows: list[dict[str, object]]) -> dict[str, dict[str, object]]:
    return {str(row.get("label", "")): row for row in rows}


def available_nanoplot_rows(rows: list[dict[str, object]]) -> list[dict[str, object]]:
    return [row for row in rows if row.get("metric_source") in {"nanoplot", "internal_qc"}]


def metric_from_row(row: dict[str, object], key: str) -> float:
    metrics = row.get("metrics", {})
    if not isinstance(metrics, dict):
        return 0.0
    return float_or_zero(metrics.get(key))


def nanoplot_table(rows: list[dict[str, object]]) -> str:
    if not rows:
        return "<p class='note'>长读长 QC 尚未生成。</p>"
    out = [
        "<table>",
        "<thead><tr><th>FASTQ</th><th>指标来源</th><th class='num'>Reads</th><th class='num'>Median length</th><th class='num'>N50</th><th class='num'>Median Q</th><th class='num'>Q≥20 保留率</th></tr></thead><tbody>",
    ]
    for row in rows:
        metrics = row.get("metrics", {})
        metrics = metrics if isinstance(metrics, dict) else {}
        q20 = threshold_metric(metrics, 20, "fraction")
        available = row.get("metric_source") in {"nanoplot", "internal_qc"}
        q20_basis = "bases" if row.get("metric_source") == "internal_qc" else "reads"
        out.append(
            "<tr>"
            f"<td>{esc(row.get('label', ''))}</td>"
            f"<td>{esc(row.get('status', ''))}</td>"
            f"<td class='num'>{fmt_int(metrics.get('number_of_reads', 0)) if available else ''}</td>"
            f"<td class='num'>{fmt_int(metrics.get('median_read_length', 0)) if available else ''}</td>"
            f"<td class='num'>{fmt_int(metrics.get('n50', 0)) if available else ''}</td>"
            f"<td class='num'>{fmt_float(metrics.get('median_qual', 0), 1) if available else ''}</td>"
            f"<td class='num'>{fmt_pct(q20) + ' (' + q20_basis + ')' if available else ''}</td>"
            "</tr>"
        )
    out.append("</tbody></table>")
    return "\n".join(out)


def available_qc_figure_note(available_rows: list[dict[str, object]]) -> str:
    if not available_rows:
        return "当前没有成功完成的 NanoPlot 结果，QC 图无法绘制。"
    labels = ", ".join(str(row.get("label", "")) for row in available_rows)
    return f"以下 QC 图仅使用 NanoPlot 成功完成的 FASTQ 结果绘制：{labels}。"


def qc_count_chart(rows: list[dict[str, object]]) -> str:
    data = []
    for row in rows:
        label = str(row.get("label", ""))
        data.append((f"{label} reads", metric_from_row(row, "number_of_reads")))
        data.append((f"{label} bases", metric_from_row(row, "number_of_bases")))
    return grouped_metric_chart("NanoPlot reads / bases summary", [(a.split()[0], a.split()[1], b) for a, b in data])


def qc_length_chart(rows: list[dict[str, object]]) -> str:
    metrics: list[tuple[str, str, float]] = []
    for row in rows:
        label = str(row.get("label", ""))
        metrics.extend(
            [
                (label, "mean", metric_from_row(row, "mean_read_length")),
                (label, "median", metric_from_row(row, "median_read_length")),
                (label, "N50", metric_from_row(row, "n50")),
            ]
        )
    return grouped_metric_chart("NanoPlot read length summary", metrics)


def qc_quality_chart(rows: list[dict[str, object]]) -> str:
    metrics: list[tuple[str, str, float]] = []
    for row in rows:
        label = str(row.get("label", ""))
        metrics.extend(
            [
                (label, "meanQ", metric_from_row(row, "mean_qual")),
                (label, "medianQ", metric_from_row(row, "median_qual")),
            ]
        )
    return grouped_metric_chart("NanoPlot read quality summary", metrics)


def qc_threshold_chart(rows: list[dict[str, object]]) -> str:
    metrics: list[tuple[str, str, float]] = []
    for row in rows:
        label = str(row.get("label", ""))
        row_metrics = row.get("metrics", {})
        row_metrics = row_metrics if isinstance(row_metrics, dict) else {}
        for threshold in QUALITY_THRESHOLDS:
            metrics.append((label, f">Q{threshold}", threshold_metric(row_metrics, threshold, "fraction") * 100.0))
    return grouped_metric_chart("NanoPlot quality threshold retention (%)", metrics)


def read_tc_event_bins(path: Path, new_threshold: int, extra_levels: int, strict_gt: bool = False) -> list[tuple[str, int]]:
    threshold = max(1, new_threshold)
    last_explicit = threshold + max(0, extra_levels) - 1
    counts = {str(value): 0 for value in range(last_explicit + 1)}
    tail_label = f">{last_explicit}" if strict_gt else f">={last_explicit + 1}"
    counts[tail_label] = 0
    if path.is_file():
        with gzip.open(path, "rt", encoding="utf-8") as handle:
            reader = csv.DictReader(handle, delimiter="\t")
            for row in reader:
                tc_events = int_or_zero(row.get("tc_events"))
                if tc_events > last_explicit:
                    counts[tail_label] += 1
                else:
                    counts[str(tc_events)] += 1
    return list(counts.items())


def main() -> None:
    args = parse_args()
    results_dir = Path(args.results_dir)
    sample_root = results_dir / args.sample_id
    qc_dir = sample_root / "01_raw_qc" / "04_tool_output" / "longread_qc"
    trim_dir = sample_root / "02_trim" / "04_tool_output" / "trim"
    align_dir = sample_root / "03_alignment" / "04_tool_output" / "minimap2_samtools"
    slam_dir = sample_root / "04_slam_calling" / "04_tool_output" / "caller"
    trim_qc = read_json(trim_dir / f"{args.sample_id}.trim_summary.json")
    nanoplot_rows = read_nanoplot_rows(qc_dir, args.sample_id)
    nanoplot_by_label = nanoplot_row_map(nanoplot_rows)
    available_qc_rows = available_nanoplot_rows(nanoplot_rows)
    raw_nanoplot = nanoplot_by_label.get("raw", {})
    insert_nanoplot = nanoplot_by_label.get("insert", {})
    raw_metrics = raw_nanoplot.get("metrics", {}) if raw_nanoplot else {}
    insert_metrics = insert_nanoplot.get("metrics", {}) if insert_nanoplot else {}
    raw_metrics = raw_metrics if isinstance(raw_metrics, dict) else {}
    insert_metrics = insert_metrics if isinstance(insert_metrics, dict) else {}
    slam_qc = read_json(slam_dir / f"{args.sample_id}.bulk_slam_qc.json")
    genes = read_top_genes(slam_dir / f"{args.sample_id}.bulk_gene_summary.tsv", limit=30)
    sample_summary = read_first_tsv_row(slam_dir / f"{args.sample_id}.bulk_slam_summary.tsv")
    assignment_rows = read_tsv_rows(slam_dir / f"{args.sample_id}.assignment_summary.tsv")
    substitution_rows = read_tsv_rows(slam_dir / f"{args.sample_id}.substitution_summary.tsv")
    grand_slam_compatibility = read_json(slam_dir / f"{args.sample_id}.grand_slam_compatibility.json")
    sample_group = read_sample_group(Path(args.samples_path), args.sample_id)
    primary_reads = read_flagstat_total(align_dir / f"{args.sample_id}.primary.flagstat.txt")

    new_threshold = read_config_int(Path(args.config_path), "new_threshold", 1)
    long_tuc_formal_threshold = read_config_int(Path(args.config_path), "long_tuc_formal_threshold", int_or_zero(slam_qc.get("long_tuc_formal_threshold")) or 20)
    long_tuc_report_thresholds = read_config_list_int(Path(args.config_path), "long_tuc_report_thresholds", [6, 20, 30])
    main_call = read_config_scalar(Path(args.config_path), "main_call", str(slam_qc.get("main_call", "site_binomial")))
    assignment_mode = read_config_scalar(Path(args.config_path), "assignment_mode", "gene_span")
    library_strand_mode = str(slam_qc.get("library_strand_mode", read_config_scalar(Path(args.config_path), "library_strand_mode", "legacy")))

    if main_call == "site_binomial":
        threshold_text = "read-level binomial + site evidence"
        threshold_rule_text = "binomial p-value 和 site-level T>C evidence"
        read_length_tc_note = "read length 与 T>C 事件数的关系；当前正式判定使用 site_binomial，T>C event burden 仅作为辅助观察指标。"
    elif main_call == "long_tuc_threshold":
        threshold_text = f"Long-TUC 文献阈值：T>C 事件数 > {long_tuc_formal_threshold}"
        threshold_rule_text = f"T>C 事件数 > {long_tuc_formal_threshold}"
        read_length_tc_note = "read length 与 T>C 事件数的关系；Long-TUC 对齐模式使用每条 read 的 T>C event burden 作为正式判定依据。"
    else:
        threshold_text = f"至少 {new_threshold} 个 T>C 事件"
        threshold_rule_text = f"T>C 事件数 >= {new_threshold}"
        read_length_tc_note = "read length 与 T>C 事件数的关系；当前 fixed_threshold 模式直接使用 T>C event count 判定 new/old。"

    report_stage = args.report_stage.resolve() if args.report_stage else sample_root / "05_sample_report"
    read_flow_dir = args.read_flow_dir.resolve() if args.read_flow_dir else report_stage / "04_tool_output" / "read_flow"
    figures_dir = report_stage / "02_figures"
    metadata_dir = report_stage / "05_metadata"
    for child in ("01_tables", "02_figures", "03_logs", "04_tool_output", "05_metadata"):
        (report_stage / child).mkdir(parents=True, exist_ok=True)

    total_reads = int_or_zero(trim_qc.get("total_reads"))
    trimmed_reads = int_or_zero(trim_qc.get("trimmed_reads"))
    forward_reads = int_or_zero(trim_qc.get("forward_reads"))
    reverse_reads = int_or_zero(trim_qc.get("reverse_reads"))
    assigned = int_or_zero(slam_qc.get("assigned_gene_reads", slam_qc.get("assigned_reads")))
    new_reads = int_or_zero(slam_qc.get("new_reads"))
    old_reads = int_or_zero(slam_qc.get("old_reads"))
    ambiguous_reads = int_or_zero(slam_qc.get("ambiguous_reads", slam_qc.get("ambiguous_read_count")))
    genes_detected = int_or_zero(slam_qc.get("genes_detected"))
    gene_reads = int_or_zero(slam_qc.get("gene_reads"))
    intergenic_reads = int_or_zero(slam_qc.get("intergenic_reads"))
    region_ambiguous_reads = int_or_zero(slam_qc.get("region_ambiguous_reads"))
    bg_tc_rate = float_or_zero(slam_qc.get("bg_tc_rate"))
    tc_rate = float_or_zero(slam_qc.get("global_tc_rate"))
    t_covered = int_or_zero(slam_qc.get("t_covered"))
    bg_rate_source = str(slam_qc.get("bg_rate_source", slam_qc.get("bg_tc_rate_source", "not_set")))
    tc_background_fold = tc_rate / bg_tc_rate if bg_tc_rate else math.nan
    assignment_note = {
        "gene_span": "gene span",
        "exon_gene": "exon_gene",
        "transcript_aware_gene": "transcript_aware_gene",
    }.get(assignment_mode, assignment_mode)

    figures = {
        "composition": composition_three_svg(new_reads, old_reads, ambiguous_reads),
        "formal_rates": t_enrichment_svg(slam_qc),
        "region": assignment_composition_svg(assignment_rows),
        "genes": horizontal_ratio_chart("Top genes 新生比例", genes),
    }

    for name, svg in figures.items():
        if svg:
            save_svg(figures_dir, f"{args.sample_id}.{name}.svg", svg)

    read_flow_fragment_path = read_flow_dir / "reports" / "read_flow_sankey_section.html"
    if not read_flow_fragment_path.is_file() or read_flow_fragment_path.stat().st_size == 0:
        raise SystemExit(f"Missing validated read-flow fragment: {read_flow_fragment_path}")
    if not args.plotly_js.is_file() or args.plotly_js.stat().st_size == 0:
        raise SystemExit(f"Missing local Plotly asset: {args.plotly_js}")
    read_flow_fragment = read_flow_fragment_path.read_text(encoding="utf-8")
    plotly_js = args.plotly_js.read_text(encoding="utf-8")
    report_title = "三代 Bulk SLAM-seq 分析报告"
    fold_display = f"{tc_background_fold:.2f}×" if math.isfinite(tc_background_fold) else "NA"
    grand_status = str(grand_slam_compatibility.get("status", "not_evaluated"))
    content = f"""<!doctype html>
<html lang="zh-CN">
<head>
<meta charset="utf-8">
<meta name="viewport" content="width=device-width, initial-scale=1.0">
<title>{esc(report_title)} - {esc(args.sample_id)}</title>
<style>
*,*::before,*::after{{box-sizing:border-box}}
body{{margin:0;background:#f7f9fc;color:#172033;font-family:"Segoe UI","Helvetica Neue",Arial,"Noto Sans SC",sans-serif;line-height:1.65}}
.hero{{background:#10233f;color:#fff;padding:44px 0 36px;margin-bottom:34px}}
.container{{max-width:1180px;margin:0 auto;padding:0 24px}}
h1{{font-size:2rem;margin:0 0 8px;font-weight:800}}
h2{{font-size:1.35rem;margin:0 0 18px;color:#1e40af;border-bottom:2px solid #1e40af;padding-bottom:10px}}
h3{{font-size:1.02rem;margin:0 0 12px;color:#172033}}
section{{margin-bottom:34px}}
.meta{{display:grid;grid-template-columns:repeat(auto-fit,minmax(180px,1fr));gap:0;border-top:1px solid rgba(255,255,255,.25);border-bottom:1px solid rgba(255,255,255,.25);margin-top:24px}}
.meta div{{padding:12px 18px;border-left:1px solid rgba(255,255,255,.18)}}.meta div:first-child{{border-left:none}}
.meta span{{display:block;font-size:.76rem;opacity:.68}}.meta strong{{display:block;font-size:1rem;white-space:nowrap;overflow:hidden;text-overflow:ellipsis}}
.kpi-grid{{display:grid;grid-template-columns:repeat(auto-fit,minmax(165px,1fr));gap:14px}}
.kpi,.card,.figure{{background:#fff;border:1px solid #dfe7f1;border-radius:8px;box-shadow:0 1px 3px rgba(15,23,42,.06);padding:18px}}
.kpi-val{{font-size:1.35rem;font-weight:800;font-family:"SF Mono","Consolas",monospace;color:#1e40af}}
.kpi-label{{font-size:.82rem;color:#526173;margin-top:4px}}.kpi-note{{font-size:.76rem;color:#748094;margin-top:3px}}
.grid-2{{display:grid;grid-template-columns:1fr 1fr;gap:18px}}@media(max-width:850px){{.grid-2{{grid-template-columns:1fr}}}}
.figure{{position:relative;cursor:zoom-in}}.figure::after{{content:"点击放大";position:absolute;right:14px;top:12px;background:rgba(15,23,42,.78);color:#fff;font-size:.72rem;font-weight:700;border-radius:999px;padding:3px 9px;opacity:.0;transition:opacity .15s}}
.figure:hover::after{{opacity:1}}.figure svg{{width:100%;height:auto;display:block}}.figure p{{margin:10px 0 0;color:#526173;font-size:.9rem}}
.zoom-modal{{position:fixed;inset:0;background:rgba(15,23,42,.82);display:none;align-items:center;justify-content:center;padding:26px;z-index:20}}
.zoom-modal.active{{display:flex}}.zoom-panel{{background:#fff;border-radius:8px;max-width:min(96vw,1280px);max-height:92vh;width:100%;padding:18px;box-shadow:0 18px 60px rgba(0,0,0,.34);overflow:auto;cursor:default}}
.zoom-panel svg{{width:100%;height:auto;min-width:760px}}.zoom-close{{position:fixed;right:22px;top:18px;border:0;border-radius:999px;background:#fff;color:#172033;font-size:1.5rem;line-height:1;width:40px;height:40px;cursor:pointer;box-shadow:0 4px 16px rgba(0,0,0,.22)}}
table{{width:100%;border-collapse:collapse;background:#fff;font-size:.9rem}}th,td{{padding:9px 10px;border-bottom:1px solid #e8eef6;text-align:left}}th{{background:#eef4fb;color:#405066;font-weight:700}}.num{{text-align:right;font-family:"SF Mono","Consolas",monospace}}.mono{{font-family:"SF Mono","Consolas",monospace;font-size:.82rem}}
.note{{color:#526173;font-size:.9rem}}.method-list li{{margin-bottom:8px}}
.refs li{{margin-bottom:8px}}a{{color:#1d4ed8;text-decoration:none}}a:hover{{text-decoration:underline}}
.sidebar{{position:fixed;inset:0 auto 0 0;width:248px;background:#10233f;color:#fff;padding:28px 20px;overflow:auto;z-index:8}}
.sidebar strong{{display:block;margin-bottom:4px}}.sidebar .sidebar-sample{{font-size:.82rem;color:#b8c7da;margin-bottom:24px;word-break:break-all}}
.sidebar a{{display:block;color:#e8f0fa;text-decoration:none;padding:9px 10px;border-radius:6px;margin:3px 0;font-size:.9rem}}.sidebar a:hover{{background:#1c3a62;color:#fff;text-decoration:none}}
.report-main{{margin-left:248px;padding:30px;max-width:1500px}}.report-main .hero{{border-radius:10px;padding:32px;margin-bottom:24px}}
.report-section{{background:#fff;border:1px solid #dfe7f1;border-radius:10px;padding:24px;margin-bottom:24px;box-shadow:0 1px 4px rgba(15,23,42,.05)}}
.footer{{margin-left:248px;padding:20px 30px;background:#10233f;color:#fff}}
.block-gap{{margin-top:18px}}.empty{{border:1px dashed #cbd5e1;background:#f8fafc;color:#64748b;border-radius:8px;padding:20px;text-align:center}}
.toolbar{{display:flex;gap:12px;align-items:center;justify-content:space-between;flex-wrap:wrap;margin-bottom:12px}}
.search-input{{min-width:240px;max-width:360px;width:100%;border:1px solid #cbd5e1;border-radius:6px;padding:8px 10px;font-size:.9rem}}
.badge{{display:inline-flex;align-items:center;gap:6px;padding:4px 9px;border-radius:999px;background:#e0f2fe;color:#075985;font-size:.78rem;font-weight:700}}
.pager{{display:flex;gap:10px;align-items:center;justify-content:flex-end;margin-top:14px}}.pager button{{border:1px solid #cbd5e1;background:#fff;border-radius:6px;padding:6px 12px;cursor:pointer}}.pager button:disabled{{opacity:.45;cursor:not-allowed}}
details.card summary{{cursor:pointer;font-weight:700}}details.card[open] summary{{margin-bottom:10px}}
.svg-title{{font-size:20px;font-weight:700;fill:#172033}}.svg-label{{font-size:13px;fill:#334155}}.svg-small{{font-size:11px;fill:#475569}}.svg-inbar{{font-size:13px;fill:#fff;font-weight:700}}.axis{{stroke:#94a3b8;stroke-width:1}}
@media(max-width:900px){{.sidebar{{position:static;width:auto}}.report-main,.footer{{margin-left:0}}.report-main{{padding:18px}}.grid-2{{grid-template-columns:1fr}}}}
@media print{{.sidebar{{display:none}}.report-main,.footer{{margin-left:0}}.report-section{{break-inside:avoid;box-shadow:none}}}}
</style>
<script>{plotly_js}</script>
<script>
var genePage = 1;
var genePageSize = 20;
function filterGeneTable() {{
  var query = document.getElementById('geneSearch').value.toLowerCase();
  var rows = Array.from(document.querySelectorAll('#geneTable tbody tr'));
  var matches = rows.filter(function(row) {{ return row.textContent.toLowerCase().indexOf(query) >= 0; }});
  var pages = Math.max(1, Math.ceil(matches.length / genePageSize));
  genePage = Math.min(Math.max(1, genePage), pages);
  rows.forEach(function(row) {{ row.hidden = true; }});
  matches.slice((genePage - 1) * genePageSize, genePage * genePageSize).forEach(function(row) {{ row.hidden = false; }});
  document.getElementById('geneShown').textContent = matches.length;
  document.getElementById('genePageInfo').textContent = '第 ' + genePage + ' / ' + pages + ' 页';
  document.getElementById('genePrev').disabled = genePage <= 1;
  document.getElementById('geneNext').disabled = genePage >= pages;
}}
function resetGeneFilter() {{ genePage = 1; filterGeneTable(); }}
function turnGenePage(step) {{ genePage += step; filterGeneTable(); }}
function openFigureZoom(figure) {{
  var svg = figure.querySelector('svg');
  if (!svg) {{ return; }}
  document.getElementById('zoomContent').innerHTML = svg.outerHTML;
  document.getElementById('figureZoom').classList.add('active');
}}
function closeFigureZoom() {{
  document.getElementById('figureZoom').classList.remove('active');
  document.getElementById('zoomContent').innerHTML = '';
}}
document.addEventListener('click', function(event) {{
  var figure = event.target.closest('.figure');
  if (figure) {{
    openFigureZoom(figure);
  }}
  if (event.target.id === 'figureZoom') {{
    closeFigureZoom();
  }}
}});
document.addEventListener('keydown', function(event) {{
  if (event.key === 'Escape') {{
    closeFigureZoom();
  }}
}});
document.addEventListener('DOMContentLoaded', filterGeneTable);
</script>
</head>
<body>
<nav class="sidebar">
<strong>ISO-SLAM 单样本报告</strong>
<div class="sidebar-sample">{esc(args.sample_id)}</div>
<a href="#overview">1 样本概览</a>
<a href="#processing">2 Reads处理流与质量</a>
<a href="#slam">3 SLAM转换信号</a>
<a href="#genes">4 基因层面结果</a>
<a href="#method">5 方法与解释边界</a>
</nav>
<main class="report-main">
<header class="hero">
<h1>{esc(report_title)}</h1>
<div>样本 {esc(args.sample_id)} · 人源 bulk 三代 SLAM-seq · T>C 转换新生 RNA 分析</div>
<div class="meta">
<div><span>项目</span><strong>{esc(args.project_name)}</strong></div>
<div><span>样本</span><strong>{esc(args.sample_id)}</strong></div>
<div><span>分组</span><strong>{esc(sample_group or "未提供")}</strong></div>
<div><span>参考 / 报告</span><strong>hg38 · schema v3</strong></div>
</div>
</header>

<section id="overview" class="report-section">
<h2>1. 样本概览</h2>
<div class="kpi-grid">
{metric_card("可用 insert reads", fmt_int(trimmed_reads), fmt_pct(trimmed_reads / total_reads if total_reads else 0))}
{metric_card("MAPQ≥20 保留率", fmt_pct(primary_reads / trimmed_reads if trimmed_reads else 0), f"{fmt_int(primary_reads)} reads")}
{metric_card("基因分配率", fmt_pct(assigned / primary_reads if primary_reads else 0), f"{fmt_int(assigned)} reads")}
{metric_card("T>C / 背景富集", fold_display, f"{fmt_pct(tc_rate)} vs {fmt_pct(bg_tc_rate)}")}
{metric_card("New fraction", fmt_pct(new_reads / assigned if assigned else 0), "描述性 hard call")}
{metric_card("检测基因数", fmt_int(genes_detected), f"strand: {library_strand_mode}")}
</div>
<p class="note block-gap">方向模型：<span class="mono">{esc(library_strand_mode)}</span>；GRAND-SLAM 兼容性：<span class="mono">{esc(grand_status)}</span>。NC 为有标记生物学对照，背景值仅用于诊断。</p>
</section>

<section id="processing" class="report-section">
<h2>2. Reads处理流与长读长质量</h2>
<p class="note">流图以本样本 Raw reads 为根；所有保留和损失分支均来自明确互斥计数，带宽、节点高度和明细百分比使用同一组已验证数据。</p>
{read_flow_fragment}
<details class="card block-gap"><summary>处理流计数口径</summary><ul class="method-list"><li>Insert 分支来自 trim summary，文库结构未识别为其互斥损失。</li><li>Primary mapped、unmapped 和 MAPQ&lt;20 来自 samtools flagstat 的 primary 口径；supplementary alignment records 不重复计作 reads。</li><li>基因分配末端由 assigned、intergenic、shared antisense ambiguous、shared gene ambiguous 和 antisense only 五类闭合。</li></ul></details>
<h3 class="block-gap">长读长质控</h3>
<div class="card" style="margin-top:18px;overflow-x:auto">
<h3>核心 FASTQ 指标</h3>
{nanoplot_table(nanoplot_rows)}
<p class="note">客户报告不展示 NanoPlot 图墙。表格优先使用流程内置 FASTQ QC 的一致口径；NanoPlot 成功时，其完整 NanoStats 与 HTML 仅保留在技术输出目录。</p>
</div>
</section>

<section id="slam" class="report-section">
<h2>3. SLAM 转换信号</h2>
<div class="grid-2">
<div class="figure">{figures["composition"]}<p>Hard call 仅作描述性 QC，不能视为单条 read 的确定新旧标签。</p></div>
<div class="figure">{figures["formal_rates"]}<p>T>C 与诊断背景使用同一线性坐标；背景来自 {esc(bg_rate_source)}，不是 no-4sU 实测背景。</p></div>
</div>
<div class="card block-gap" style="overflow-x:auto"><h3>替换事件精确汇总</h3>{substitution_table(substitution_rows)}</div>
<div class="figure" style="margin-top:18px">{figures["region"]}<p>分母为 MAPQ≥20 primary reads；详细 assignment 原始分类保留在 TSV，不再用九个尺度悬殊的柱子并列展示。</p></div>
</section>

<section id="genes" class="report-section">
<h2>4. 基因层面结果</h2>
<div class="figure">{figures["genes"]}<p>灰色为旧 RNA 比例，蓝色为新生 RNA 比例；右侧显示新生比例和该基因总 reads 数。</p></div>
<div class="card" style="margin-top:18px;overflow-x:auto">
<div class="toolbar">
<h3>Top genes 明细</h3>
<span class="badge">当前显示 <span id="geneShown">{len(genes)}</span> / {len(genes)}</span>
<input id="geneSearch" class="search-input" type="search" placeholder="搜索基因名或 Gene ID" oninput="resetGeneFilter()">
</div>
<table id="geneTable" data-page-size="20">
<thead><tr><th>基因</th><th>Gene ID</th><th class="num">总 reads</th><th class="num">新生 reads</th><th class="num">旧 reads</th><th class="num">ambiguous</th><th class="num">新生比例</th><th class="num">gene span reads</th><th class="num">T>C 事件数</th><th class="num">T>C 转化率</th></tr></thead>
<tbody>
{table_rows(genes)}
</tbody>
</table>
<div class="pager"><button id="genePrev" type="button" onclick="turnGenePage(-1)">上一页</button><span id="genePageInfo"></span><button id="geneNext" type="button" onclick="turnGenePage(1)">下一页</button></div>
</div>
</section>

<section id="method" class="report-section">
<h2>5. 方法说明与解释边界</h2>
<div class="card">
<ul class="method-list">
<li>预处理：按用户提供的文库结构截取 insert，支持反向互补方向，并剔除短于配置阈值的 insert。</li>
<li>长读长 QC：客户报告仅展示 reads、median length、N50、median Q 和 Q≥20 保留率；内置 QC 的 Q≥20 为碱基口径，NanoStats 回退时为 reads 口径并在表中标注。完整 NanoPlot 报告仅作技术附件，不把 reads 与 bases 混画。</li>
<li>比对：使用 minimap2 splice 模式比对到人 hg38，保留非 secondary/supplementary 的 primary reads；正式 calling 中低于配置 MAPQ 阈值的 reads 标记为 ambiguous 或不参与 site evidence。</li>
<li>基因分配：根据 GENCODE 基因区间重叠给 read 分配唯一最佳基因；并按基因链方向统计等价 T>C 转换。</li>
<li>区域分类：`transcript_aware_gene` 模式优先做 transcript-compatible 分配；其次通过 gene-level salvage 保留跨 transcript 但仍落在同一 gene 的 reads；仅 gene overlap 未满足 transcript 兼容的 reads 降级为 gene_overlap_unresolved 诊断类但保留 gene-level assignment；反义链 reads 优先尝试 gene-level salvage；无基因重叠的 reads 标记为 intergenic。</li>
<li>Strand 处理：显式使用 <span class="mono">library_strand_mode={esc(library_strand_mode)}</span>；方向审计必须先于 transcript-aware assignment，避免以大规模 antisense salvage 掩盖方向错误。</li>
<li>T 碱基替换频率：在转录本链方向统一后，只统计 query base quality 达到配置阈值的参考 T 位点；低质量碱基不计入 T 覆盖或 mismatch 计数。</li>
<li>site-level evidence：按每个参考 T 位点统计 T>C、T>A、T>G 和 deletion；T>C 比例高且 T>A/T>G 低的位点标记为 high_conf_tc_site，杂乱错配高的位点标记为 noisy_site。</li>
<li>新旧 RNA 判定：本次主模式为 <span class="mono">{esc(main_call)}</span>，即每条 read 使用 {threshold_rule_text} 判定 new RNA；不满足 new 且不满足 ambiguous 的 reads 判定为 old RNA。</li>
<li>背景率：本次诊断背景为 <span class="mono">{fmt_pct(bg_tc_rate)}</span>，来源为 <span class="mono">{esc(bg_rate_source)}</span>；NC 是有标记生物学对照，不能充当 no-4sU 背景。</li>
<li>GRAND-SLAM：本样本兼容性状态为 <span class="mono">{esc(grand_status)}</span>；只有工具、BAM 标签和合成真值门禁通过后，NTR posterior 才能进入正式统计。</li>
<li>解释边界：本报告不估计绝对合成率、降解率或半衰期；new fraction 仅作为本样本描述性指标。</li>
</ul>
<ol class="refs">
<li>TimeLapse-seq：通过代谢标记引发的 U/T-to-C 转换推断新生 RNA。</li>
<li>SLAM-seq：以 4-thiouridine 标记和 T>C 转换解析 RNA 动态。</li>
<li>GRAND-SLAM / grandR：用于评估基因级 NTR posterior 与下游 nucleotide-conversion RNA-seq 统计。</li>
<li>SLAM-DUNK：其 SNP、mutation spectrum 和位置偏倚 QC 思路可借鉴，但 3′UTR/QuantSeq 假设不直接替代本 PacBio 主流程。</li>
</ol>
</div>
<div class="card note" style="margin-top:18px">
<h3>文件路径</h3>
<p>配置文件：<span class="mono">{esc(args.config_path)}</span></p>
<p>样本表：<span class="mono">{esc(args.samples_path)}</span></p>
<p>结果目录：<span class="mono">{esc(results_dir)}</span></p>
</div>
</section>
</main>
<footer class="footer">ISO-SLAM · example_sample · 中文离线交付报告</footer>
<div id="figureZoom" class="zoom-modal" aria-hidden="true">
<button class="zoom-close" type="button" onclick="closeFigureZoom()" aria-label="关闭放大图">×</button>
<div id="zoomContent" class="zoom-panel"></div>
</div>
</body>
</html>
"""
    report_path = report_stage / f"{args.sample_id}.report.html"
    report_path.write_text(content, encoding="utf-8")
    module_artifacts = [
        ("assignment_summary", "全量 gene assignment 汇总", slam_dir / f"{args.sample_id}.assignment_summary.tsv"),
        ("substitution_summary", "全量 T 碱基替换汇总", slam_dir / f"{args.sample_id}.substitution_summary.tsv"),
        ("ntr_sufficient_stats", "基因级 NTR 充分统计", slam_dir / f"{args.sample_id}.ntr_sufficient_stats.tsv"),
        ("grand_slam_compatibility", "GRAND-SLAM 兼容性审计", slam_dir / f"{args.sample_id}.grand_slam_compatibility.json"),
        ("read_flow_nodes", "守恒处理流节点表", read_flow_dir / "tables" / "read_flow_nodes.tsv"),
        ("read_flow_links", "守恒处理流连接表", read_flow_dir / "tables" / "read_flow_links.tsv"),
        ("read_flow_payload", "静态和交互式流图共享数据", read_flow_dir / "tables" / "read_flow_payload.json"),
        ("read_flow_sources", "处理流源文件哈希清单", read_flow_dir / "tables" / "read_flow_source_manifest.tsv"),
        ("read_flow_svg", "处理流静态 SVG", read_flow_dir / "figures" / "read_flow_sankey" / "sample_processing.svg"),
        ("read_flow_validation", "Sankey 输出验证", read_flow_dir / "reports" / "sankey_validation.json"),
        ("read_flow_adapter_validation", "ISO-SLAM 计数守恒验证", read_flow_dir / "reports" / "read_flow_adapter_validation.json"),
    ]
    missing_artifacts = [str(path) for _, _, path in module_artifacts if not path.is_file() or path.stat().st_size == 0]
    if missing_artifacts:
        raise SystemExit(f"Missing read-flow delivery artifacts: {', '.join(missing_artifacts)}")
    manifest_path = metadata_dir / "delivery_manifest.tsv"
    with manifest_path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=["artifact", "description", "target", "status"],
            delimiter="\t",
            lineterminator="\n",
        )
        writer.writeheader()
        for artifact, description, path in module_artifacts:
            writer.writerow(
                {"artifact": artifact, "description": description, "target": str(path.resolve()), "status": "up_to_date"}
            )
        writer.writerow(
            {
                "artifact": "final_sample_report",
                "description": f"{args.sample_id} final sample report",
                "target": str(report_path.resolve()),
                "status": "up_to_date",
            }
        )
    build_metadata = {
        "schema_version": "iso-slam-sample-report-v3",
        "sample_id": args.sample_id,
        "report": str(report_path.resolve()),
        "manifest": str(manifest_path.resolve()),
        "read_flow": str(read_flow_dir.resolve()),
        "sections": ["overview", "processing", "slam", "genes", "method"],
        "plotly_inline_count": 1,
    }
    (metadata_dir / "report_build.json").write_text(
        json.dumps(build_metadata, ensure_ascii=False, indent=2) + "\n", encoding="utf-8"
    )
    legacy_report = report_stage / "report.html"
    if legacy_report.is_file() and legacy_report.resolve() != report_path.resolve():
        legacy_report.unlink()
    print(report_path)


if __name__ == "__main__":
    main()
