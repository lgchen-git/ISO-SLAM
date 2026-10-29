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


COLORS = ["#2563eb", "#16a34a", "#f97316", "#7c3aed", "#0f766e", "#dc2626"]
QUALITY_THRESHOLDS = [10, 15, 20, 25, 30]


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Build bulk SLAM-seq report.")
    parser.add_argument("--results-dir", required=True)
    parser.add_argument("--sample-id", required=True)
    parser.add_argument("--project-name", required=True)
    parser.add_argument("--report-title", required=True)
    parser.add_argument("--config-path", required=True)
    parser.add_argument("--samples-path", required=True)
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
        f'<rect x="92" y="120" width="{new_w:.1f}" height="42" fill="#2563eb" rx="6"/>',
        f'<rect x="{92 + new_w:.1f}" y="120" width="{old_w:.1f}" height="42" fill="#94a3b8" rx="0"/>',
        f'<rect x="{92 + new_w + old_w:.1f}" y="120" width="{ambiguous_w:.1f}" height="42" fill="#f97316" rx="6"/>',
        f'<text x="92" y="190" class="svg-label">new {fmt_int(new_reads)} ({fmt_pct(new_reads / total)})</text>',
        f'<text x="342" y="190" class="svg-label">old {fmt_int(old_reads)} ({fmt_pct(old_reads / total)})</text>',
        f'<text x="592" y="190" class="svg-label">ambiguous {fmt_int(ambiguous_reads)} ({fmt_pct(ambiguous_reads / total)})</text>',
    ]
    parts.append("</svg>")
    return "\n".join(parts)


def region_svg(slam_qc: dict[str, object], assignment_mode: str) -> str:
    if assignment_mode == "transcript_aware_gene":
        values = [
            ("tx_unique", int_or_zero(slam_qc.get("tx_unique_reads"))),
            ("gene_multi_tx", int_or_zero(slam_qc.get("gene_unique_multi_tx_reads"))),
            ("intron_gene", int_or_zero(slam_qc.get("intron_only_unique_gene_reads"))),
            ("gene_overlap", int_or_zero(slam_qc.get("gene_overlap_unresolved_reads"))),
            ("shared_resolved", int_or_zero(slam_qc.get("shared_gene_resolved_reads"))),
            ("antisense_salvage", int_or_zero(slam_qc.get("antisense_gene_salvaged_reads"))),
            ("shared_ambiguous", int_or_zero(slam_qc.get("shared_gene_ambiguous_reads"))),
            ("antisense_only", int_or_zero(slam_qc.get("antisense_only_reads"))),
            ("intergenic", int_or_zero(slam_qc.get("intergenic_reads"))),
        ]
        return bar_chart("transcript_aware gene assignment distribution", values, "reads")
    return bar_chart(
        "基因分配 reads 分布",
        [
            ("gene span", int_or_zero(slam_qc.get("gene_reads"))),
            ("基因间", int_or_zero(slam_qc.get("intergenic_reads"))),
            ("基因分配歧义", int_or_zero(slam_qc.get("region_ambiguous_reads"))),
        ],
        "reads",
    )


def t_mutation_svg(slam_qc: dict[str, object]) -> str:
    t_covered = int_or_zero(slam_qc.get("t_covered"))
    data = []
    for label, key in (("T>A", "t_to_a"), ("T>C", "t_to_c"), ("T>G", "t_to_g")):
        value = int_or_zero(slam_qc.get(key))
        data.append((label, (value / t_covered * 100) if t_covered else 0.0))
    return bar_chart("T 碱基替换频率", data, "% of covered T")


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


def read_nanoplot_rows(qc_dir: Path, sample_id: str) -> list[dict[str, object]]:
    rows: list[dict[str, object]] = []
    for label in ("raw", "insert"):
        status_path = qc_dir / "nanoplot" / label / "nanoplot_status.json"
        status = read_json(status_path)
        if not status:
            continue
        stats_path = Path(str(status.get("nanostats", "")))
        report_path = Path(str(status.get("report", "")))
        metrics = parse_nanostats(stats_path) if status.get("status") == "complete" else {"thresholds": {}, "raw_rows": []}
        rows.append(
            {
                "label": label,
                "status": status.get("status", "missing"),
                "version": status.get("version", ""),
                "reason": status.get("reason", ""),
                "metrics": metrics,
                "stats_path": stats_path,
                "report_path": report_path,
                "report_href": f"../samples/{sample_id}/00_longread_qc/nanoplot/{label}/{report_path.name}" if report_path.name else "",
            }
        )
    return rows


def nanoplot_row_map(rows: list[dict[str, object]]) -> dict[str, dict[str, object]]:
    return {str(row.get("label", "")): row for row in rows}


def available_nanoplot_rows(rows: list[dict[str, object]]) -> list[dict[str, object]]:
    return [row for row in rows if row.get("status") == "complete"]


def metric_from_row(row: dict[str, object], key: str) -> float:
    metrics = row.get("metrics", {})
    if not isinstance(metrics, dict):
        return 0.0
    return float_or_zero(metrics.get(key))


def nanoplot_table(rows: list[dict[str, object]]) -> str:
    if not rows:
        return "<p class='note'>NanoPlot QC 未启用或尚未生成。</p>"
    out = [
        "<table>",
        "<thead><tr><th>FASTQ</th><th>Status</th><th>NanoPlot</th><th class='num'>Reads</th><th class='num'>Bases</th><th class='num'>N50</th><th class='num'>Mean len</th><th class='num'>Mean Q</th><th class='num'>>Q20</th><th class='num'>>Q30</th><th>输出</th></tr></thead><tbody>",
    ]
    for row in rows:
        metrics = row.get("metrics", {})
        metrics = metrics if isinstance(metrics, dict) else {}
        report_href = str(row.get("report_href", ""))
        report_link = f"<a href='{esc(report_href)}'>NanoPlot report</a>" if report_href else ""
        q20 = threshold_metric(metrics, 20, "fraction")
        q30 = threshold_metric(metrics, 30, "fraction")
        out.append(
            "<tr>"
            f"<td>{esc(row.get('label', ''))}</td>"
            f"<td>{esc(row.get('status', ''))}</td>"
            f"<td class='mono'>{esc(str(row.get('version', '')).splitlines()[-1] if row.get('version') else '')}</td>"
            f"<td class='num'>{fmt_int(metrics.get('number_of_reads', 0)) if row.get('status') == 'complete' else ''}</td>"
            f"<td class='num'>{fmt_int(metrics.get('number_of_bases', 0)) if row.get('status') == 'complete' else ''}</td>"
            f"<td class='num'>{fmt_int(metrics.get('n50', 0)) if row.get('status') == 'complete' else ''}</td>"
            f"<td class='num'>{fmt_float(metrics.get('mean_read_length', 0), 1) if row.get('status') == 'complete' else ''}</td>"
            f"<td class='num'>{fmt_float(metrics.get('mean_qual', 0), 1) if row.get('status') == 'complete' else ''}</td>"
            f"<td class='num'>{fmt_pct(q20) if row.get('status') == 'complete' else ''}</td>"
            f"<td class='num'>{fmt_pct(q30) if row.get('status') == 'complete' else ''}</td>"
            f"<td>{report_link or esc(row.get('reason', ''))}</td>"
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
    read_rows = read_tsv_rows(slam_dir / f"{args.sample_id}.bulk_read_calls.tsv.gz", limit=200000)
    site_rows = read_tsv_rows(slam_dir / f"{args.sample_id}.bulk_t_sites.tsv.gz", limit=200000)
    primary_reads = read_flagstat_total(align_dir / f"{args.sample_id}.primary.flagstat.txt")

    new_threshold = read_config_int(Path(args.config_path), "new_threshold", 1)
    long_tuc_formal_threshold = read_config_int(Path(args.config_path), "long_tuc_formal_threshold", int_or_zero(slam_qc.get("long_tuc_formal_threshold")) or 20)
    long_tuc_report_thresholds = read_config_list_int(Path(args.config_path), "long_tuc_report_thresholds", [6, 20, 30])
    main_call = read_config_scalar(Path(args.config_path), "main_call", str(slam_qc.get("main_call", "site_binomial")))
    assignment_mode = read_config_scalar(Path(args.config_path), "assignment_mode", "gene_span")
    tc_event_extra_levels = max(0, read_config_int(Path(args.config_path), "tc_event_extra_levels", 3))

    if main_call == "site_binomial":
        threshold_text = "read-level binomial + site evidence"
        threshold_rule_text = "binomial p-value 和 site-level T>C evidence"
        read_length_tc_note = "read length 与 T>C 事件数的关系；当前正式判定使用 site_binomial，T>C event burden 仅作为辅助观察指标。"
        tc_bin_threshold = new_threshold
        tc_bin_extra_levels = tc_event_extra_levels
        strict_tc_tail = False
    elif main_call == "long_tuc_threshold":
        threshold_text = f"Long-TUC 文献阈值：T>C 事件数 > {long_tuc_formal_threshold}"
        threshold_rule_text = f"T>C 事件数 > {long_tuc_formal_threshold}"
        read_length_tc_note = "read length 与 T>C 事件数的关系；Long-TUC 对齐模式使用每条 read 的 T>C event burden 作为正式判定依据。"
        tc_bin_threshold = min(long_tuc_report_thresholds) if long_tuc_report_thresholds else long_tuc_formal_threshold
        tc_bin_extra_levels = max(0, long_tuc_formal_threshold - tc_bin_threshold + 1)
        strict_tc_tail = True
    else:
        threshold_text = f"至少 {new_threshold} 个 T>C 事件"
        threshold_rule_text = f"T>C 事件数 >= {new_threshold}"
        read_length_tc_note = "read length 与 T>C 事件数的关系；当前 fixed_threshold 模式直接使用 T>C event count 判定 new/old。"
        tc_bin_threshold = new_threshold
        tc_bin_extra_levels = tc_event_extra_levels
        strict_tc_tail = False

    tc_bins = read_tc_event_bins(
        slam_dir / f"{args.sample_id}.bulk_read_calls.tsv.gz",
        tc_bin_threshold,
        tc_bin_extra_levels,
        strict_tc_tail,
    )

    reports_dir = sample_root / "05_sample_report"
    deliverables_dir = sample_root / "05_sample_report" / "04_tool_output" / "delivery"
    figures_dir = sample_root / "05_sample_report" / "02_figures"
    reports_dir.mkdir(parents=True, exist_ok=True)
    deliverables_dir.mkdir(parents=True, exist_ok=True)

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
    assignment_note = {
        "gene_span": "gene span",
        "exon_gene": "exon_gene",
        "transcript_aware_gene": "transcript_aware_gene",
    }.get(assignment_mode, assignment_mode)

    read_tc = [int_or_zero(row.get("tc_events")) for row in read_rows]
    read_t_covered = [int_or_zero(row.get("t_covered")) for row in read_rows]
    read_other_rate = [float_or_zero(row.get("other_mismatch_rate")) for row in read_rows]
    read_scatter = [(float_or_zero(row.get("read_length")), float_or_zero(row.get("tc_events"))) for row in read_rows[:8000]]

    tc_rate_bins = [("0-1%", 0.0, 0.01), ("1-5%", 0.01, 0.05), ("5-10%", 0.05, 0.10), ("10-20%", 0.10, 0.20), (">=20%", 0.20, 1.01)]
    site_tc_rate = [float_or_zero(row.get("tc_rate")) for row in site_rows]
    site_other_rate = [float_or_zero(row.get("other_mismatch_rate")) for row in site_rows]
    site_call_counts: dict[str, int] = {}
    for row in site_rows:
        site_call = row.get("site_call", "")
        if site_call:
            site_call_counts[site_call] = site_call_counts.get(site_call, 0) + 1

    qc_note = available_qc_figure_note(available_qc_rows)

    figures = {
        "funnel": funnel_chart(
            [
                ("Raw reads", total_reads),
                ("Insert reads", trimmed_reads),
                ("Primary Q20 reads", primary_reads),
                ("Assigned gene reads", assigned),
            ],
            "数据保留漏斗",
        ),
        "qc_counts": qc_count_chart(available_qc_rows) if available_qc_rows else "",
        "qc_lengths": qc_length_chart(available_qc_rows) if available_qc_rows else "",
        "qc_quality": qc_quality_chart(available_qc_rows) if available_qc_rows else "",
        "qc_thresholds": qc_threshold_chart(available_qc_rows) if available_qc_rows else "",
        "composition": composition_three_svg(new_reads, old_reads, ambiguous_reads),
        "formal_rates": t_mutation_svg(slam_qc),
        "region": region_svg(slam_qc, assignment_mode),
        "site_tc_rate": histogram_chart("site-level T>C rate 分布", site_tc_rate, tc_rate_bins, "sites"),
        "site_other_rate": histogram_chart("site-level other mismatch rate 分布", site_other_rate, tc_rate_bins, "sites"),
        "site_calls": bar_chart("site-level calls", [(key, value) for key, value in sorted(site_call_counts.items())], "sites"),
        "read_t_covered": histogram_chart(
            "read-level T coverage 分布",
            [float(value) for value in read_t_covered],
            [("0-9", 0, 10), ("10-19", 10, 20), ("20-49", 20, 50), ("50-99", 50, 100), ("100-199", 100, 200), (">=200", 200, 10**9)],
            "reads",
        ),
        "read_tc_count": bar_chart("read-level T>C 事件数分布", [(label, value) for label, value in tc_bins], "reads"),
        "read_other_rate": histogram_chart(
            "read-level other mismatch rate 分布",
            read_other_rate,
            [("0-1%", 0.0, 0.01), ("1-5%", 0.01, 0.05), ("5-10%", 0.05, 0.10), ("10-20%", 0.10, 0.20), (">=20%", 0.20, 1.01)],
            "reads",
        ),
        "read_length_tc": scatter_chart("read length vs tc_count", read_scatter, "read length", "tc_count"),
        "genes": horizontal_ratio_chart("Top genes 新生比例", genes),
    }

    for name, svg in figures.items():
        if svg:
            save_svg(figures_dir, f"{args.sample_id}.{name}.svg", svg)

    report_title = "三代 Bulk SLAM-seq 分析报告"
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
.tab-bar{{display:flex;gap:6px;flex-wrap:wrap;margin:0 0 18px;padding:6px;background:#e8eef6;border:1px solid #d7e1ee;border-radius:8px;position:sticky;top:0;z-index:3}}
.tab-btn{{border:0;background:transparent;color:#405066;font-weight:700;font-size:.9rem;padding:10px 14px;border-radius:6px;cursor:pointer}}
.tab-btn:hover{{background:#f8fafc;color:#1e40af}}.tab-btn.active{{background:#1e40af;color:#fff;box-shadow:0 1px 3px rgba(30,64,175,.22)}}
.tab-panel{{display:none}}.tab-panel.active{{display:block}}
.toolbar{{display:flex;gap:12px;align-items:center;justify-content:space-between;flex-wrap:wrap;margin-bottom:12px}}
.search-input{{min-width:240px;max-width:360px;width:100%;border:1px solid #cbd5e1;border-radius:6px;padding:8px 10px;font-size:.9rem}}
.badge{{display:inline-flex;align-items:center;gap:6px;padding:4px 9px;border-radius:999px;background:#e0f2fe;color:#075985;font-size:.78rem;font-weight:700}}
.svg-title{{font-size:20px;font-weight:700;fill:#172033}}.svg-label{{font-size:13px;fill:#334155}}.svg-small{{font-size:11px;fill:#475569}}.svg-inbar{{font-size:13px;fill:#fff;font-weight:700}}.axis{{stroke:#94a3b8;stroke-width:1}}
</style>
<script>
function switchReportTab(event, id) {{
  document.querySelectorAll('.tab-btn').forEach(function(btn) {{ btn.classList.remove('active'); }});
  document.querySelectorAll('.tab-panel').forEach(function(panel) {{ panel.classList.remove('active'); }});
  event.currentTarget.classList.add('active');
  document.getElementById(id).classList.add('active');
}}
function filterGeneTable() {{
  var query = document.getElementById('geneSearch').value.toLowerCase();
  var shown = 0;
  document.querySelectorAll('#geneTable tbody tr').forEach(function(row) {{
    var match = row.textContent.toLowerCase().indexOf(query) >= 0;
    row.style.display = match ? '' : 'none';
    if (match) {{ shown += 1; }}
  }});
  document.getElementById('geneShown').textContent = shown;
}}
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
</script>
</head>
<body>
<header class="hero">
<div class="container">
<h1>{esc(report_title)}</h1>
<div>样本 {esc(args.sample_id)} · 人源 bulk 三代 SLAM-seq · T>C 转换新生 RNA 分析</div>
<div class="meta">
<div><span>项目</span><strong>{esc(args.project_name)}</strong></div>
<div><span>样本</span><strong>{esc(args.sample_id)}</strong></div>
<div><span>参考基因组</span><strong>hg38 / GENCODE</strong></div>
<div><span>报告类型</span><strong>中文图文版</strong></div>
</div>
</div>
</header>
<main class="container">
<div class="sidebar"><div class="tab-bar" role="tablist" aria-label="报告结果模块">
<button class="tab-btn active" onclick="switchReportTab(event,'tab-overview')">样本概览</button>
<button class="tab-btn" onclick="switchReportTab(event,'tab-qc')">NanoPlot 质控</button>
<button class="tab-btn" onclick="switchReportTab(event,'tab-slam')">SLAM 转换</button>
<button class="tab-btn" onclick="switchReportTab(event,'tab-genes')">基因结果</button>
<button class="tab-btn" onclick="switchReportTab(event,'tab-method')">方法与文件</button>
</div></div>

<section id="tab-overview" class="tab-panel active">
<h2>1. 样本概览</h2>
<div class="kpi-grid">
{metric_card("原始 reads", fmt_int(total_reads), "trim summary")}
{metric_card("成功提取 insert", f"{fmt_int(trimmed_reads)} ({fmt_pct(trimmed_reads / total_reads if total_reads else 0)})", "文库结构识别摘要")}
{metric_card("MAPQ>=20 primary reads", fmt_int(primary_reads), "比对后高质量 reads")}
{metric_card("基因分配 reads", f"{fmt_int(assigned)} ({fmt_pct(assigned / primary_reads if primary_reads else 0)})", assignment_note)}
{metric_card("新生 RNA reads", f"{fmt_int(new_reads)} ({fmt_pct(new_reads / assigned if assigned else 0)})", threshold_text)}
{metric_card("ambiguous reads", f"{fmt_int(ambiguous_reads)} ({fmt_pct(ambiguous_reads / (new_reads + old_reads + ambiguous_reads) if (new_reads + old_reads + ambiguous_reads) else 0)})", "formal caller")}
{metric_card("Raw NanoPlot N50", fmt_int(raw_metrics.get("n50", 0)), "NanoPlot")}
{metric_card("Insert NanoPlot N50", fmt_int(insert_metrics.get("n50", 0)), "NanoPlot")}
{metric_card("Insert mean length", fmt_float(insert_metrics.get("mean_read_length", 0), 1), "NanoPlot")}
{metric_card("Insert mean Q", fmt_float(insert_metrics.get("mean_qual", 0), 1), "NanoPlot")}
{metric_card("Insert >Q20 reads", f"{fmt_int(threshold_metric(insert_metrics, 20, 'count'))} ({fmt_pct(threshold_metric(insert_metrics, 20, 'fraction'))})", "NanoPlot")}
{metric_card("Insert >Q30 reads", f"{fmt_int(threshold_metric(insert_metrics, 30, 'count'))} ({fmt_pct(threshold_metric(insert_metrics, 30, 'fraction'))})", "NanoPlot")}
{metric_card("背景 T>C rate", fmt_pct(bg_tc_rate), bg_rate_source)}
{metric_card("T>C 转化率", fmt_pct(tc_rate), f"covered T: {fmt_int(t_covered)}")}
{metric_card("检测基因数", fmt_int(genes_detected), "gene-level summary")}
</div>
<div class="grid-2" style="margin-top:18px">
<div class="figure">{figures["funnel"]}<p>从原始 reads 到基因分配 reads 的逐级保留情况，用于定位主要数据损失环节。</p></div>
<div class="card">
<h3>结果摘要</h3>
<ul class="method-list">
<li>本样本共读取 {fmt_int(total_reads)} 条三代 reads，成功提取 insert {fmt_int(trimmed_reads)} 条；其中正向 {fmt_int(forward_reads)} 条，反向 {fmt_int(reverse_reads)} 条。</li>
<li>NanoPlot raw QC 读取 {fmt_int(raw_metrics.get("number_of_reads", 0))} 条 reads，insert QC 读取 {fmt_int(insert_metrics.get("number_of_reads", 0))} 条 reads。</li>
<li>高质量 primary 比对 reads 为 {fmt_int(primary_reads)}，其中 {fmt_int(assigned)} 条可分配到唯一基因。</li>
<li>本次按{threshold_text}判定新生 RNA，new reads 为 {fmt_int(new_reads)}，old reads 为 {fmt_int(old_reads)}，ambiguous reads 为 {fmt_int(ambiguous_reads)}；所有覆盖的转录本 T 碱基中，T>C 转化率为 {fmt_pct(tc_rate)}。</li>
</ul>
</div>
</div>
</section>

<section id="tab-qc" class="tab-panel">
<h2>2. NanoPlot 长读长质控</h2>
<div class="card">
<p class="note">{esc(qc_note)}</p>
<p class="note">本模块所有图和表只使用 NanoPlot 结果。primer/adaptor 命中、insert 提取率、比对保留率和基因分配率不在此处作为 QC 绘图数据源展示。</p>
</div>
{f"<div class='grid-2'><div class='figure'>{figures['qc_counts']}<p>Raw / insert 的 reads 数和总碱基数，全部直接来自 NanoPlot NanoStats。</p></div><div class='figure'>{figures['qc_lengths']}<p>Raw / insert 的 mean length、median length 和 N50，全部直接来自 NanoPlot NanoStats。</p></div></div>" if available_qc_rows else ""}
{f"<div class='grid-2' style='margin-top:18px'><div class='figure'>{figures['qc_quality']}<p>Raw / insert 的 mean Q、median Q，全部直接来自 NanoPlot NanoStats。</p></div><div class='figure'>{figures['qc_thresholds']}<p>Reads >Q10 / >Q15 / >Q20 / >Q25 / >Q30 的比例，全部直接来自 NanoPlot NanoStats。</p></div></div>" if available_qc_rows else ""}
<div class="card" style="margin-top:18px;overflow-x:auto">
<h3>NanoPlot QC 明细</h3>
{nanoplot_table(nanoplot_rows)}
</div>
<div class="card" style="margin-top:18px">
<h3>质控解读</h3>
<ul class="method-list">
<li>raw NanoPlot 状态：{esc(raw_nanoplot.get("status", "missing"))}；insert NanoPlot 状态：{esc(insert_nanoplot.get("status", "missing"))}。</li>
<li>如果某个 FASTQ 的 NanoPlot 状态为 failed 或 skipped，本模块不会为其生成伪 0 图，只在表格中保留状态和错误信息。</li>
<li>文库结构识别率为 {fmt_pct(trimmed_reads / total_reads if total_reads else 0)}，但该指标来源于 trim summary，只作为流程摘要，不属于 NanoPlot QC 图口径。</li>
</ul>
</div>
</section>

<section id="tab-slam" class="tab-panel">
<h2>3. SLAM 转换信号</h2>
<div class="grid-2">
<div class="figure">{figures["composition"]}<p>正式版 caller 输出 new / old / ambiguous 三类；ambiguous 主要来自低 MAPQ、T 覆盖不足、T>A/T>G 杂乱错配过高或命中 noisy site。</p></div>
<div class="figure">{figures["formal_rates"]}<p>T>C 是 SLAM 标记关注的主信号，T>A 和 T>G 用作背景错配参考；未提供 control 时，背景 T>C rate 使用 T>A/T>G 近似。</p></div>
</div>
<div class="figure" style="margin-top:18px">{figures["region"]}<p>当前 assignment 模式为 {esc(assignment_note)}。tx_unique / gene_multi_tx / intron_gene / gene_overlap / antisense_salvage 均计入 gene summary（gene-level assigned）。shared_ambiguous / antisense_only / intergenic 不计入 gene summary（unresolved）。</p></div>
<div class="grid-2" style="margin-top:18px">
<div class="figure">{figures["site_tc_rate"]}<p>每个参考 T 位点的 T>C read-supported conversion rate 分布，用于评估转换位点是否集中在高可信区间。</p></div>
<div class="figure">{figures["site_other_rate"]}<p>每个参考 T 位点的 T>A/T>G 杂乱错配比例；高比例位点更可能来自测序错误、错比对、SNP 或 RNA editing。</p></div>
</div>
<div class="grid-2" style="margin-top:18px">
<div class="figure">{figures["site_calls"]}<p>site-level 分类数量：high_conf_tc_site 支持 new RNA，noisy_site 不作为 new RNA 证据。</p></div>
<div class="figure">{figures["read_t_covered"]}<p>read-level 高质量参考 T 位点覆盖数；低于 min_t_covered 的 reads 进入 ambiguous。</p></div>
</div>
<div class="grid-2" style="margin-top:18px">
<div class="figure">{figures["read_tc_count"]}<p>每条 read 的高质量 T>C 事件数分布，用于观察转换信号强度。</p></div>
<div class="figure">{figures["read_other_rate"]}<p>每条 read 的 T>A/T>G 杂乱错配比例；超过阈值的 reads 标记为 ambiguous。</p></div>
</div>
<div class="figure" style="margin-top:18px">{figures["read_length_tc"]}<p>{esc(read_length_tc_note)}</p></div>
</section>

<section id="tab-genes" class="tab-panel">
<h2>4. 基因层面结果</h2>
<div class="figure">{figures["genes"]}<p>灰色为旧 RNA 比例，蓝色为新生 RNA 比例；右侧显示新生比例和该基因总 reads 数。</p></div>
<div class="card" style="margin-top:18px;overflow-x:auto">
<div class="toolbar">
<h3>Top genes 明细</h3>
<span class="badge">当前显示 <span id="geneShown">{len(genes)}</span> / {len(genes)}</span>
<input id="geneSearch" class="search-input" type="search" placeholder="搜索基因名或 Gene ID" oninput="filterGeneTable()">
</div>
<table id="geneTable" data-page-size="20">
<thead><tr><th>基因</th><th>Gene ID</th><th class="num">总 reads</th><th class="num">新生 reads</th><th class="num">旧 reads</th><th class="num">ambiguous</th><th class="num">新生比例</th><th class="num">gene span reads</th><th class="num">T>C 事件数</th><th class="num">T>C 转化率</th></tr></thead>
<tbody>
{table_rows(genes)}
</tbody>
</table>
</div>
</section>

<section id="tab-method" class="tab-panel">
<h2>5. 方法说明与文件</h2>
<div class="card">
<ul class="method-list">
<li>预处理：按用户提供的文库结构截取 insert，支持反向互补方向，并剔除短于配置阈值的 insert。</li>
<li>长读长 QC：bulk 报告中的 QC 图和 QC KPI 只读取 NanoPlot 结果，即 <span class="mono">nanoplot_status.json</span> 与 <span class="mono">NanoStats.txt</span>。</li>
<li>如果配置 <span class="mono">enable_nanoplot_qc: true</span>，流程会对 raw 和/或 insert FASTQ 运行 NanoPlot，并将 NanoStats 与 NanoPlot HTML report 放入样本 <span class="mono">00_longread_qc/nanoplot/</span> 目录。</li>
<li>比对：使用 minimap2 splice 模式比对到人 hg38，保留非 secondary/supplementary 的 primary reads；正式 calling 中低于配置 MAPQ 阈值的 reads 标记为 ambiguous 或不参与 site evidence。</li>
<li>基因分配：根据 GENCODE 基因区间重叠给 read 分配唯一最佳基因；并按基因链方向统计等价 T>C 转换。</li>
<li>区域分类：`transcript_aware_gene` 模式优先做 transcript-compatible 分配；其次通过 gene-level salvage 保留跨 transcript 但仍落在同一 gene 的 reads；仅 gene overlap 未满足 transcript 兼容的 reads 降级为 gene_overlap_unresolved 诊断类但保留 gene-level assignment；反义链 reads 优先尝试 gene-level salvage；无基因重叠的 reads 标记为 intergenic。</li>
<li>Strand 处理：`assignment_require_same_strand` 控制是否要求 read alignment 链向与基因注释链向一致；bulk PacBio cDNA 数据建议关闭此选项。</li>
<li>T 碱基替换频率：在转录本链方向统一后，只统计 query base quality 达到配置阈值的参考 T 位点；低质量碱基不计入 T 覆盖或 mismatch 计数。</li>
<li>site-level evidence：按每个参考 T 位点统计 T>C、T>A、T>G 和 deletion；T>C 比例高且 T>A/T>G 低的位点标记为 high_conf_tc_site，杂乱错配高的位点标记为 noisy_site。</li>
<li>新旧 RNA 判定：本次主模式为 <span class="mono">{esc(main_call)}</span>，即每条 read 使用 {threshold_rule_text} 判定 new RNA；不满足 new 且不满足 ambiguous 的 reads 判定为 old RNA。</li>
<li>背景率：本次 bg_tc_rate 为 <span class="mono">{fmt_pct(bg_tc_rate)}</span>，来源为 <span class="mono">{esc(bg_rate_source)}</span>。如果来源不是 control-derived background，应视为 T>A/T>G 近似背景。</li>
</ul>
<ol class="refs">
<li>TimeLapse-seq：通过代谢标记引发的 U/T-to-C 转换推断新生 RNA。</li>
<li>SLAM-seq：以 4-thiouridine 标记和 T>C 转换解析 RNA 动态。</li>
<li>SLAM-DUNK：提供 T>C 转换 reads 和 gene/sample 层面 labeled fraction 的分析框架。</li>
<li>grandR：强调 nucleotide-conversion RNA-seq 的 QC、可视化和基因层面解释。</li>
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
<footer class="footer"><div class="container">ISO-SLAM · example_sample · 离线交付报告</div></footer>
<div id="figureZoom" class="zoom-modal" aria-hidden="true">
<button class="zoom-close" type="button" onclick="closeFigureZoom()" aria-label="关闭放大图">×</button>
<div id="zoomContent" class="zoom-panel"></div>
</div>
</body>
</html>
"""
    report_path = reports_dir / f"{args.sample_id}.bulk_slam_report.html"
    report_path.write_text(content, encoding="utf-8")
    (deliverables_dir / report_path.name).write_text(content, encoding="utf-8")


if __name__ == "__main__":
    main()
