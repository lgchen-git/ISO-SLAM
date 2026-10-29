#!/usr/bin/env python3
"""Summarize T-site mutation spectra for bulk long-read SLAM-seq runs."""

from __future__ import annotations

import argparse
import csv
import gzip
import html
import json
from collections import Counter, defaultdict
from pathlib import Path


RATE_BINS = [
    ("0", 0.0, 0.0),
    ("<0.1%", 0.0, 0.001),
    ("0.1-0.5%", 0.001, 0.005),
    ("0.5-1%", 0.005, 0.01),
    ("1-5%", 0.01, 0.05),
    ("5-10%", 0.05, 0.10),
    (">=10%", 0.10, None),
]
TC_COUNT_BINS = ["0", "1", "2", "3-5", "6-10", ">10"]
COLORS = {
    "T>C": "#2563eb",
    "T>A": "#f97316",
    "T>G": "#0f766e",
    "Other mismatch": "#7c3aed",
    "Background T>C": "#64748b",
}


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Build mutation QC tables and SVGs for lr_bulk_slam outputs.")
    parser.add_argument("--results-dir", required=True, type=Path)
    parser.add_argument("--samples", required=True, help="Comma-separated sample IDs")
    return parser.parse_args()


def esc(value: object) -> str:
    return html.escape("" if value is None else str(value))


def num(value: object) -> float:
    try:
        return float(value)
    except (TypeError, ValueError):
        return 0.0


def integer(value: object) -> int:
    return int(num(value))


def fmt_pct(value: object) -> str:
    return f"{num(value) * 100:.2f}%"


def fmt_rate(value: object) -> str:
    value_num = num(value)
    if value_num == 0:
        return "0"
    if abs(value_num) < 0.001:
        return f"{value_num:.2e}"
    return f"{value_num:.4f}"


def opener(path: Path):
    return gzip.open if path.suffix == ".gz" else open


def read_first_tsv(path: Path) -> dict[str, str]:
    if not path.is_file():
        return {}
    with opener(path)(path, "rt", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        return next(reader, {})


def read_json(path: Path) -> dict[str, object]:
    if not path.is_file():
        return {}
    return json.loads(path.read_text(encoding="utf-8"))


def rate_bin(value: object) -> str:
    val = num(value)
    if val == 0:
        return "0"
    for label, low, high in RATE_BINS[1:]:
        if high is None:
            if val >= low:
                return label
        elif low < val <= high:
            return label
    return "0"


def tc_count_bin(value: object) -> str:
    val = integer(value)
    if val <= 2:
        return str(max(val, 0))
    if val <= 5:
        return "3-5"
    if val <= 10:
        return "6-10"
    return ">10"


def iter_rows(path: Path, limit: int):
    if not path.is_file():
        return
    with opener(path)(path, "rt", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        for idx, row in enumerate(reader, start=1):
            if limit and idx > limit:
                break
            yield row


def write_tsv(path: Path, headers: list[str], rows: list[dict[str, object]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, delimiter="\t", fieldnames=headers)
        writer.writeheader()
        for row in rows:
            writer.writerow({key: row.get(key, "") for key in headers})


def bar_svg(chart_id: str, title: str, rows: list[tuple[str, float, str]], y_label: str, percent: bool) -> str:
    width, height = 920, 360
    left, right, top, bottom = 88, 38, 58, 82
    plot_w = width - left - right
    plot_h = height - top - bottom
    max_value = max([value for _, value, _ in rows] + [1.0 if not percent else 0.01])
    slot = plot_w / max(len(rows), 1)
    bar_w = slot * 0.58
    parts = [
        f"<svg class='svg-chart' id='{esc(chart_id)}' xmlns='http://www.w3.org/2000/svg' viewBox='0 0 {width} {height}' role='img' aria-label='{esc(title)}'>",
        "<rect width='100%' height='100%' fill='#fff'/>",
        f"<text x='{width/2:.1f}' y='31' text-anchor='middle' class='svg-title'>{esc(title)}</text>",
        f"<line x1='{left}' y1='{top + plot_h}' x2='{width - right}' y2='{top + plot_h}' class='axis'/>",
        f"<line x1='{left}' y1='{top}' x2='{left}' y2='{top + plot_h}' class='axis'/>",
        f"<text x='24' y='{top + plot_h/2:.1f}' transform='rotate(-90 24 {top + plot_h/2:.1f})' text-anchor='middle' class='svg-label'>{esc(y_label)}</text>",
    ]
    for idx, (label, value, color) in enumerate(rows):
        value = max(float(value), 0.0)
        x = left + idx * slot + (slot - bar_w) / 2
        h = plot_h * value / max_value if max_value else 0
        y = top + plot_h - h
        shown = fmt_pct(value) if percent else f"{value:,.0f}"
        parts.append(f"<g data-tip='{esc(label)}: {esc(shown)}'>")
        parts.append(f"<rect x='{x:.1f}' y='{y:.1f}' width='{bar_w:.1f}' height='{h:.1f}' fill='{color}' rx='4'/>")
        parts.append(f"<text x='{x + bar_w/2:.1f}' y='{max(y - 8, 48):.1f}' text-anchor='middle' class='svg-small'>{esc(shown)}</text>")
        parts.append(
            f"<text x='{x + bar_w/2:.1f}' y='{top + plot_h + 24}' text-anchor='middle' "
            f"transform='rotate(30 {x + bar_w/2:.1f} {top + plot_h + 24})' class='svg-small'>{esc(label)}</text>"
        )
        parts.append("</g>")
    parts.append("</svg>")
    return "".join(parts)


def horizontal_bar_svg(chart_id: str, title: str, rows: list[tuple[str, float, str]], x_label: str, percent: bool) -> str:
    width, height = 920, 360
    left, right, top, bottom = 170, 52, 58, 48
    plot_w = width - left - right
    row_h = (height - top - bottom) / max(len(rows), 1)
    bar_h = min(28.0, row_h * 0.52)
    max_value = max([value for _, value, _ in rows] + [0.01 if percent else 1.0])
    parts = [
        f"<svg class='svg-chart' id='{esc(chart_id)}' xmlns='http://www.w3.org/2000/svg' viewBox='0 0 {width} {height}' role='img' aria-label='{esc(title)}'>",
        "<rect width='100%' height='100%' fill='#fff'/>",
        f"<text x='{width/2:.1f}' y='31' text-anchor='middle' class='svg-title'>{esc(title)}</text>",
        f"<line x1='{left}' y1='{height - bottom}' x2='{width - right}' y2='{height - bottom}' class='axis'/>",
        f"<text x='{left + plot_w/2:.1f}' y='{height - 14}' text-anchor='middle' class='svg-label'>{esc(x_label)}</text>",
    ]
    for idx, (label, value, color) in enumerate(rows):
        value = max(float(value), 0.0)
        y = top + idx * row_h + (row_h - bar_h) / 2
        w = plot_w * value / max_value if max_value else 0.0
        shown = fmt_pct(value) if percent else f"{value:,.0f}"
        parts.append(f"<text x='{left - 12}' y='{y + bar_h/2 + 4:.1f}' text-anchor='end' class='svg-label'>{esc(label)}</text>")
        parts.append(f"<g data-tip='{esc(label)}: {esc(shown)}'>")
        parts.append(f"<rect x='{left}' y='{y:.1f}' width='{w:.1f}' height='{bar_h:.1f}' fill='{color}' rx='4'/>")
        parts.append(f"<text x='{left + w + 8:.1f}' y='{y + bar_h/2 + 4:.1f}' class='svg-small'>{esc(shown)}</text>")
        parts.append("</g>")
    parts.append("</svg>")
    return "".join(parts)


def grouped_bar_svg(chart_id: str, title: str, samples: list[dict[str, object]]) -> str:
    width, height = 980, 390
    left, right, top, bottom = 82, 36, 60, 86
    plot_w = width - left - right
    plot_h = height - top - bottom
    metrics = [("T>C", "global_tc_rate"), ("T>A", "global_ta_rate"), ("T>G", "global_tg_rate")]
    max_value = max([num(sample.get(key)) for sample in samples for _, key in metrics] + [0.01])
    group_w = plot_w / max(len(samples), 1)
    bar_w = min(28.0, group_w / 5.0)
    parts = [
        f"<svg class='svg-chart' id='{esc(chart_id)}' xmlns='http://www.w3.org/2000/svg' viewBox='0 0 {width} {height}' role='img' aria-label='{esc(title)}'>",
        "<rect width='100%' height='100%' fill='#fff'/>",
        f"<text x='{width/2:.1f}' y='31' text-anchor='middle' class='svg-title'>{esc(title)}</text>",
        f"<line x1='{left}' y1='{top + plot_h}' x2='{width - right}' y2='{top + plot_h}' class='axis'/>",
        f"<line x1='{left}' y1='{top}' x2='{left}' y2='{top + plot_h}' class='axis'/>",
        f"<text x='23' y='{top + plot_h/2:.1f}' transform='rotate(-90 23 {top + plot_h/2:.1f})' text-anchor='middle' class='svg-label'>rate</text>",
    ]
    for sidx, sample in enumerate(samples):
        gx = left + sidx * group_w + group_w / 2
        parts.append(f"<text x='{gx:.1f}' y='{top + plot_h + 26}' text-anchor='middle' class='svg-label'>{esc(sample.get('sample'))}</text>")
        for midx, (label, key) in enumerate(metrics):
            value = num(sample.get(key))
            x = gx - bar_w * 1.7 + midx * bar_w * 1.25
            h = plot_h * value / max_value if max_value else 0.0
            y = top + plot_h - h
            color = COLORS[label]
            parts.append(f"<g data-tip='{esc(sample.get('sample'))} {label}: {fmt_pct(value)}'>")
            parts.append(f"<rect x='{x:.1f}' y='{y:.1f}' width='{bar_w:.1f}' height='{h:.1f}' fill='{color}' rx='3'/>")
            parts.append("</g>")
    legend_y = height - 24
    for idx, (label, _) in enumerate(metrics):
        x = left + idx * 115
        parts.append(f"<rect x='{x}' y='{legend_y - 12}' width='15' height='15' fill='{COLORS[label]}' rx='3'/>")
        parts.append(f"<text x='{x + 21}' y='{legend_y}' class='svg-label'>{label}</text>")
    parts.append("</svg>")
    return "".join(parts)


def stacked_bins_svg(
    chart_id: str,
    title: str,
    bins: list[str],
    series: list[tuple[str, Counter[str], str]],
    y_label: str,
) -> str:
    width, height = 980, 390
    left, right, top, bottom = 88, 38, 60, 94
    plot_w = width - left - right
    plot_h = height - top - bottom
    max_value = max([counter.get(label, 0) for label in bins for _, counter, _ in series] + [1])
    group_w = plot_w / max(len(bins), 1)
    bar_w = min(28.0, group_w / (len(series) + 1.3))
    parts = [
        f"<svg class='svg-chart' id='{esc(chart_id)}' xmlns='http://www.w3.org/2000/svg' viewBox='0 0 {width} {height}' role='img' aria-label='{esc(title)}'>",
        "<rect width='100%' height='100%' fill='#fff'/>",
        f"<text x='{width/2:.1f}' y='31' text-anchor='middle' class='svg-title'>{esc(title)}</text>",
        f"<line x1='{left}' y1='{top + plot_h}' x2='{width - right}' y2='{top + plot_h}' class='axis'/>",
        f"<line x1='{left}' y1='{top}' x2='{left}' y2='{top + plot_h}' class='axis'/>",
        f"<text x='24' y='{top + plot_h/2:.1f}' transform='rotate(-90 24 {top + plot_h/2:.1f})' text-anchor='middle' class='svg-label'>{esc(y_label)}</text>",
    ]
    for bidx, label in enumerate(bins):
        gx = left + bidx * group_w + group_w / 2
        parts.append(
            f"<text x='{gx:.1f}' y='{top + plot_h + 24}' text-anchor='middle' "
            f"transform='rotate(25 {gx:.1f} {top + plot_h + 24})' class='svg-small'>{esc(label)}</text>"
        )
        for sidx, (series_label, counter, color) in enumerate(series):
            value = counter.get(label, 0)
            x = gx - (len(series) * bar_w * 0.6) + sidx * bar_w * 1.25
            h = plot_h * value / max_value if max_value else 0
            y = top + plot_h - h
            parts.append(f"<g data-tip='{esc(series_label)} {esc(label)}: {value:,}'>")
            parts.append(f"<rect x='{x:.1f}' y='{y:.1f}' width='{bar_w:.1f}' height='{h:.1f}' fill='{color}' rx='3'/>")
            parts.append("</g>")
    legend_y = height - 20
    for idx, (series_label, _, color) in enumerate(series):
        x = left + idx * 145
        parts.append(f"<rect x='{x}' y='{legend_y - 12}' width='15' height='15' fill='{color}' rx='3'/>")
        parts.append(f"<text x='{x + 21}' y='{legend_y}' class='svg-label'>{esc(series_label)}</text>")
    parts.append("</svg>")
    return "".join(parts)


def summarize_sample(results_dir: Path, sample: str) -> dict[str, object]:
    slam_dir = results_dir / "samples" / sample / "03_bulk_slam"
    summary = read_first_tsv(slam_dir / f"{sample}.bulk_slam_summary.tsv")
    qc = read_json(slam_dir / f"{sample}.bulk_slam_qc.json")
    sample_row = {
        "sample": sample,
        "total_t_covered": integer(summary.get("total_t_covered") or qc.get("t_covered")),
        "total_tc_count": integer(summary.get("total_tc_count") or qc.get("t_to_c")),
        "total_ta_count": integer(summary.get("total_ta_count") or qc.get("t_to_a")),
        "total_tg_count": integer(summary.get("total_tg_count") or qc.get("t_to_g")),
        "global_tc_rate": num(summary.get("global_tc_rate") or qc.get("global_tc_rate")),
        "global_ta_rate": num(summary.get("global_ta_rate") or qc.get("t_to_a_rate")),
        "global_tg_rate": num(summary.get("global_tg_rate") or qc.get("t_to_g_rate")),
        "global_other_mismatch_rate": num(summary.get("global_other_mismatch_rate") or qc.get("global_other_mismatch_rate")),
        "bg_tc_rate": num(summary.get("bg_tc_rate") or qc.get("bg_tc_rate")),
        "bg_rate_source": summary.get("bg_rate_source") or qc.get("bg_rate_source", ""),
    }

    return sample_row


def write_svgs(fig_dir: Path, sample_rows: list[dict[str, object]]) -> None:
    fig_dir.mkdir(parents=True, exist_ok=True)
    (fig_dir / "sample_mutation_spectrum.svg").write_text(
        grouped_bar_svg("sample-mutation-spectrum", "样本级 T 位点突变谱", sample_rows),
        encoding="utf-8",
    )
    bg_rows: list[tuple[str, float, str]] = []
    for sample in sample_rows:
        sid = str(sample["sample"])
        bg_rows.append((f"{sid} T>C", num(sample["global_tc_rate"]), COLORS["T>C"]))
        bg_rows.append((f"{sid} bg", num(sample["bg_tc_rate"]), COLORS["Background T>C"]))
    (fig_dir / "sample_tc_vs_background.svg").write_text(
        bar_svg("sample-tc-vs-background", "T>C 与背景 T>C 对比", bg_rows, "rate", True),
        encoding="utf-8",
    )
    for sample in sample_rows:
        sid = str(sample["sample"])
        (fig_dir / f"{sid}.mutation_spectrum.svg").write_text(
            horizontal_bar_svg(
                f"{sid}-mutation-spectrum",
                f"{sid} T 位点突变谱",
                [
                    ("T>C", num(sample["global_tc_rate"]), COLORS["T>C"]),
                    ("T>A", num(sample["global_ta_rate"]), COLORS["T>A"]),
                    ("T>G", num(sample["global_tg_rate"]), COLORS["T>G"]),
                    ("Other mismatch", num(sample["global_other_mismatch_rate"]), COLORS["Other mismatch"]),
                    ("Background T>C", num(sample["bg_tc_rate"]), COLORS["Background T>C"]),
                ],
                "rate",
                True,
            ),
            encoding="utf-8",
        )


def main() -> None:
    args = parse_args()
    samples = [item.strip() for item in args.samples.split(",") if item.strip()]
    if not samples:
        raise SystemExit("No samples provided")
    table_dir = args.results_dir / "tables" / "lr_bulk_slam" / "mutation_qc"
    fig_dir = args.results_dir / "figures" / "lr_bulk_slam" / "mutation_qc"
    sample_rows: list[dict[str, object]] = []
    for sample in samples:
        sample_row = summarize_sample(
            args.results_dir,
            sample,
        )
        sample_rows.append(sample_row)

    write_tsv(
        table_dir / "lr_bulk_slam_mutation_summary.tsv",
        [
            "sample",
            "total_t_covered",
            "total_tc_count",
            "total_ta_count",
            "total_tg_count",
            "global_tc_rate",
            "global_ta_rate",
            "global_tg_rate",
            "global_other_mismatch_rate",
            "bg_tc_rate",
            "bg_rate_source",
        ],
        sample_rows,
    )
    write_svgs(fig_dir, sample_rows)


if __name__ == "__main__":
    main()
