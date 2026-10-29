#!/usr/bin/env python3
"""Summarize long-read QC and long-read plus SLAM joint metrics."""

from __future__ import annotations

import argparse
import csv
import gzip
import html
import json
from collections import Counter, defaultdict
from pathlib import Path


COLORS = {
    "new": "#2563eb",
    "old": "#94a3b8",
    "ambiguous": "#f97316",
    "assignment": "#7c3aed",
    "slam": "#dc2626",
    "other": "#0f766e",
    "raw": "#64748b",
    "insert": "#2563eb",
}
LENGTH_BINS = [("lt_500", "<500"), ("500_999", "500-999"), ("1000_1499", "1000-1499"), ("1500_1999", "1500-1999"), ("ge_2000", ">=2000")]
T_COV_BINS = [("0", "0"), ("1_9", "1-9"), ("10_19", "10-19"), ("20_49", "20-49"), ("ge_50", ">=50")]
ASSIGNMENT_REASONS = {"no_gene", "shared_antisense_ambiguous", "gene_assignment_tie", "antisense_only", "shared_gene_ambiguous"}


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Build long-read and long-read+SLAM QC summaries.")
    parser.add_argument("--results-dir", required=True, type=Path)
    parser.add_argument("--samples", required=True)
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


def read_json(path: Path) -> dict[str, object]:
    if not path.is_file():
        return {}
    return json.loads(path.read_text(encoding="utf-8"))


def read_tsv(path: Path) -> list[dict[str, str]]:
    if not path.is_file():
        return []
    opener = gzip.open if path.suffix == ".gz" else open
    with opener(path, "rt", encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def read_first_tsv(path: Path) -> dict[str, str]:
    rows = read_tsv(path)
    return rows[0] if rows else {}


def write_tsv(path: Path, headers: list[str], rows: list[dict[str, object]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, delimiter="\t", fieldnames=headers)
        writer.writeheader()
        for row in rows:
            writer.writerow({key: row.get(key, "") for key in headers})


def parse_flagstat(path: Path) -> dict[str, int]:
    stats = {"total": 0, "primary": 0, "secondary": 0, "supplementary": 0, "mapped": 0}
    if not path.is_file():
        return stats
    for line in path.read_text(encoding="utf-8", errors="ignore").splitlines():
        value = integer(line.split()[0]) if line.split() else 0
        if " in total " in line:
            stats["total"] = value
        elif " primary" in line and "primary duplicates" not in line and "primary mapped" not in line:
            stats["primary"] = value
        elif " secondary" in line:
            stats["secondary"] = value
        elif " supplementary" in line:
            stats["supplementary"] = value
        elif " mapped " in line and "primary mapped" not in line:
            stats["mapped"] = value
    return stats


def length_bin(length: object) -> str:
    val = integer(length)
    if val < 500:
        return "lt_500"
    if val < 1000:
        return "500_999"
    if val < 1500:
        return "1000_1499"
    if val < 2000:
        return "1500_1999"
    return "ge_2000"


def t_cov_bin(value: object) -> str:
    val = integer(value)
    if val == 0:
        return "0"
    if val < 10:
        return "1_9"
    if val < 20:
        return "10_19"
    if val < 50:
        return "20_49"
    return "ge_50"


def svg_bar(chart_id: str, title: str, rows: list[tuple[str, float, str]], y_label: str = "reads") -> str:
    width, height = 900, 340
    left, right, top, bottom = 78, 34, 56, 78
    plot_w = width - left - right
    plot_h = height - top - bottom
    is_rate = y_label.lower() == "rate"
    max_value = max([value for _, value, _ in rows] + ([0.01] if is_rate else [1.0]))
    slot = plot_w / max(len(rows), 1)
    bar_w = slot * 0.58
    parts = [
        f"<svg class='svg-chart' id='{esc(chart_id)}' xmlns='http://www.w3.org/2000/svg' viewBox='0 0 {width} {height}' role='img' aria-label='{esc(title)}'>",
        "<rect width='100%' height='100%' fill='#fff'/>",
        f"<text x='{width/2:.1f}' y='30' text-anchor='middle' class='svg-title'>{esc(title)}</text>",
        f"<line x1='{left}' y1='{top + plot_h}' x2='{width - right}' y2='{top + plot_h}' class='axis'/>",
        f"<line x1='{left}' y1='{top}' x2='{left}' y2='{top + plot_h}' class='axis'/>",
        f"<text x='22' y='{top + plot_h/2:.1f}' transform='rotate(-90 22 {top + plot_h/2:.1f})' text-anchor='middle' class='svg-label'>{esc(y_label)}</text>",
    ]
    for idx, (label, value, color) in enumerate(rows):
        value = max(float(value), 0.0)
        x = left + idx * slot + (slot - bar_w) / 2
        h = plot_h * value / max_value if max_value else 0.0
        y = top + plot_h - h
        shown = fmt_pct(value) if is_rate or value < 1 else f"{value:,.0f}"
        parts.append(f"<g data-tip='{esc(label)}: {esc(shown)}'>")
        parts.append(f"<rect x='{x:.1f}' y='{y:.1f}' width='{bar_w:.1f}' height='{h:.1f}' fill='{color}' rx='4'/>")
        parts.append(f"<text x='{x + bar_w/2:.1f}' y='{max(y - 8, 48):.1f}' text-anchor='middle' class='svg-small'>{esc(shown)}</text>")
        parts.append(
            f"<text x='{x + bar_w/2:.1f}' y='{top + plot_h + 24}' text-anchor='middle' "
            f"transform='rotate(28 {x + bar_w/2:.1f} {top + plot_h + 24})' class='svg-small'>{esc(label)}</text>"
        )
        parts.append("</g>")
    parts.append("</svg>")
    return "".join(parts)


def svg_horizontal(chart_id: str, title: str, rows: list[tuple[str, float, str]], x_label: str = "reads") -> str:
    width, height = 920, max(320, 68 + 34 * len(rows))
    left, right, top, bottom = 230, 52, 52, 44
    plot_w = width - left - right
    row_h = (height - top - bottom) / max(len(rows), 1)
    bar_h = min(24.0, row_h * 0.58)
    max_value = max([value for _, value, _ in rows] + [1.0])
    parts = [
        f"<svg class='svg-chart' id='{esc(chart_id)}' xmlns='http://www.w3.org/2000/svg' viewBox='0 0 {width} {height}' role='img' aria-label='{esc(title)}'>",
        "<rect width='100%' height='100%' fill='#fff'/>",
        f"<text x='{width/2:.1f}' y='29' text-anchor='middle' class='svg-title'>{esc(title)}</text>",
        f"<line x1='{left}' y1='{height - bottom}' x2='{width - right}' y2='{height - bottom}' class='axis'/>",
        f"<text x='{left + plot_w/2:.1f}' y='{height - 12}' text-anchor='middle' class='svg-label'>{esc(x_label)}</text>",
    ]
    for idx, (label, value, color) in enumerate(rows):
        value = max(float(value), 0.0)
        y = top + idx * row_h + (row_h - bar_h) / 2
        w = plot_w * value / max_value if max_value else 0.0
        shown = f"{value:,.0f}" if value >= 1 else fmt_pct(value)
        parts.append(f"<text x='{left - 12}' y='{y + bar_h/2 + 4:.1f}' text-anchor='end' class='svg-label'>{esc(label)}</text>")
        parts.append(f"<g data-tip='{esc(label)}: {esc(shown)}'>")
        parts.append(f"<rect x='{left}' y='{y:.1f}' width='{w:.1f}' height='{bar_h:.1f}' fill='{color}' rx='4'/>")
        parts.append(f"<text x='{left + w + 8:.1f}' y='{y + bar_h/2 + 4:.1f}' class='svg-small'>{esc(shown)}</text>")
        parts.append("</g>")
    parts.append("</svg>")
    return "".join(parts)


def svg_stacked_percent(chart_id: str, title: str, bins: list[tuple[str, str]], counts: dict[str, Counter[str]]) -> str:
    width, height = 920, 360
    left, right, top, bottom = 118, 44, 56, 78
    plot_w = width - left - right
    plot_h = height - top - bottom
    row_w = plot_w / max(len(bins), 1)
    parts = [
        f"<svg class='svg-chart' id='{esc(chart_id)}' xmlns='http://www.w3.org/2000/svg' viewBox='0 0 {width} {height}' role='img' aria-label='{esc(title)}'>",
        "<rect width='100%' height='100%' fill='#fff'/>",
        f"<text x='{width/2:.1f}' y='30' text-anchor='middle' class='svg-title'>{esc(title)}</text>",
        f"<line x1='{left}' y1='{top + plot_h}' x2='{width - right}' y2='{top + plot_h}' class='axis'/>",
        f"<line x1='{left}' y1='{top}' x2='{left}' y2='{top + plot_h}' class='axis'/>",
    ]
    for idx, (key, label) in enumerate(bins):
        counter = counts.get(key, Counter())
        total = max(sum(counter.values()), 1)
        x = left + idx * row_w + row_w * 0.18
        bar_w = row_w * 0.52
        y = top + plot_h
        for call in ("new", "old", "ambiguous"):
            value = counter.get(call, 0)
            frac = value / total
            h = plot_h * frac
            y -= h
            parts.append(f"<rect data-tip='{esc(label)} {call}: {fmt_pct(frac)} ({value:,})' x='{x:.1f}' y='{y:.1f}' width='{bar_w:.1f}' height='{h:.1f}' fill='{COLORS[call]}' rx='3'/>")
        parts.append(f"<text x='{x + bar_w/2:.1f}' y='{top + plot_h + 24}' text-anchor='middle' class='svg-small'>{esc(label)}</text>")
    for idx, call in enumerate(("new", "old", "ambiguous")):
        lx = left + idx * 126
        parts.append(f"<rect x='{lx}' y='{height - 22}' width='15' height='15' fill='{COLORS[call]}' rx='3'/>")
        parts.append(f"<text x='{lx + 21}' y='{height - 10}' class='svg-label'>{call}</text>")
    parts.append("</svg>")
    return "".join(parts)


def summarize_sample(results_dir: Path, sample: str) -> tuple[list[dict[str, object]], list[dict[str, object]], list[dict[str, object]], dict[str, str]]:
    sample_dir = results_dir / "samples" / sample
    qc_dir = sample_dir / "00_longread_qc" / "fastq"
    trim = read_json(sample_dir / "01_trim" / f"{sample}.trim_summary.json")
    raw_qc = read_json(qc_dir / f"{sample}.raw_fastq_qc.json")
    insert_qc = read_json(qc_dir / f"{sample}.insert_fastq_qc.json")
    flag = parse_flagstat(sample_dir / "02_alignment" / f"{sample}.flagstat.txt")
    primary_flag = parse_flagstat(sample_dir / "02_alignment" / f"{sample}.primary.flagstat.txt")
    slam = read_first_tsv(sample_dir / "03_bulk_slam" / f"{sample}.bulk_slam_summary.tsv")

    total_reads = integer(trim.get("total_reads") or raw_qc.get("read_count"))
    insert_reads = integer(trim.get("trimmed_reads") or insert_qc.get("read_count"))
    primary_q20 = primary_flag["total"]
    assigned = integer(slam.get("assigned_reads"))
    longread_rows = [
        {
            "sample": sample,
            "raw_reads": total_reads,
            "insert_reads": insert_reads,
            "raw_qc_reads": integer(raw_qc.get("read_count")),
            "insert_qc_reads": integer(insert_qc.get("read_count")),
            "fastq_qc_max_reads": integer(raw_qc.get("max_reads")),
            "insert_rate": insert_reads / total_reads if total_reads else 0.0,
            "raw_n50": integer(raw_qc.get("n50")),
            "insert_n50": integer(insert_qc.get("n50")),
            "raw_mean_length": num(raw_qc.get("mean_length")),
            "insert_mean_length": num(insert_qc.get("mean_length")),
            "raw_q20_base_fraction": num(raw_qc.get("q20_base_fraction")),
            "insert_q20_base_fraction": num(insert_qc.get("q20_base_fraction")),
            "raw_q30_base_fraction": num(raw_qc.get("q30_base_fraction")),
            "insert_q30_base_fraction": num(insert_qc.get("q30_base_fraction")),
            "raw_gc_fraction": num(raw_qc.get("gc_fraction")),
            "insert_gc_fraction": num(insert_qc.get("gc_fraction")),
            "raw_n_fraction": num(raw_qc.get("n_fraction")),
            "insert_n_fraction": num(insert_qc.get("n_fraction")),
            "align_total": flag["total"],
            "align_primary": flag["primary"],
            "align_supplementary": flag["supplementary"],
            "primary_q20_reads": primary_q20,
            "primary_q20_rate": primary_q20 / insert_reads if insert_reads else 0.0,
            "assigned_reads": assigned,
            "assigned_rate": assigned / primary_q20 if primary_q20 else 0.0,
        }
    ]

    reason_counts: Counter[str] = Counter()
    assignment_unresolved = 0
    slam_ambiguous = 0
    length_call_counts: dict[str, Counter[str]] = defaultdict(Counter)
    tcov_call_counts: dict[str, Counter[str]] = defaultdict(Counter)
    mutation_by_len: dict[str, Counter[str]] = defaultdict(Counter)

    with gzip.open(sample_dir / "03_bulk_slam" / f"{sample}.bulk_read_calls.tsv.gz", "rt", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        for row in reader:
            call = row.get("final_call", "unknown") or "unknown"
            lb = length_bin(row.get("read_length"))
            tb = t_cov_bin(row.get("t_covered"))
            length_call_counts[lb][call] += 1
            tcov_call_counts[tb][call] += 1
            covered = integer(row.get("t_covered"))
            mutation_by_len[lb]["t_covered"] += covered
            mutation_by_len[lb]["tc"] += integer(row.get("tc_count") or row.get("t_to_c"))
            mutation_by_len[lb]["ta"] += integer(row.get("ta_count") or row.get("t_to_a"))
            mutation_by_len[lb]["tg"] += integer(row.get("tg_count") or row.get("t_to_g"))
            if call == "ambiguous":
                reasons = [item for item in (row.get("ambiguity_reason", "") or "unknown").split(",") if item]
                reason_counts.update(reasons)
                if any(reason in ASSIGNMENT_REASONS for reason in reasons):
                    assignment_unresolved += 1
                else:
                    slam_ambiguous += 1

    reason_rows = [
        {
            "sample": sample,
            "reason": reason,
            "reads": count,
            "reason_class": "assignment_unresolved" if reason in ASSIGNMENT_REASONS else "slam_ambiguous",
        }
        for reason, count in sorted(reason_counts.items(), key=lambda item: (-item[1], item[0]))
    ]
    joint_rows: list[dict[str, object]] = []
    for key, label in LENGTH_BINS:
        counter = length_call_counts.get(key, Counter())
        total = sum(counter.values())
        mut = mutation_by_len.get(key, Counter())
        covered = mut.get("t_covered", 0)
        joint_rows.append(
            {
                "sample": sample,
                "metric": "read_length_bin",
                "bin": label,
                "total_reads": total,
                "new_reads": counter.get("new", 0),
                "old_reads": counter.get("old", 0),
                "ambiguous_reads": counter.get("ambiguous", 0),
                "tc_rate": mut.get("tc", 0) / covered if covered else 0.0,
                "ta_rate": mut.get("ta", 0) / covered if covered else 0.0,
                "tg_rate": mut.get("tg", 0) / covered if covered else 0.0,
            }
        )
    for key, label in T_COV_BINS:
        counter = tcov_call_counts.get(key, Counter())
        total = sum(counter.values())
        joint_rows.append(
            {
                "sample": sample,
                "metric": "t_covered_bin",
                "bin": label,
                "total_reads": total,
                "new_reads": counter.get("new", 0),
                "old_reads": counter.get("old", 0),
                "ambiguous_reads": counter.get("ambiguous", 0),
                "tc_rate": "",
                "ta_rate": "",
                "tg_rate": "",
            }
        )

    callable_total = integer(slam.get("new_read_count")) + integer(slam.get("old_read_count")) + slam_ambiguous
    summary = {
        "assignment_unresolved": str(assignment_unresolved),
        "slam_ambiguous": str(slam_ambiguous),
        "callable_total": str(callable_total),
        "callable_ambiguous_fraction": f"{slam_ambiguous / callable_total if callable_total else 0.0:.8f}",
    }

    fig_qc = results_dir / "figures" / "lr_bulk_slam" / "longread_qc"
    fig_joint = results_dir / "figures" / "lr_bulk_slam" / "longread_slam"
    fig_qc.mkdir(parents=True, exist_ok=True)
    fig_joint.mkdir(parents=True, exist_ok=True)

    length_rows = []
    for label, color, qc in (("raw", COLORS["raw"], raw_qc), ("insert", COLORS["insert"], insert_qc)):
        bins = qc.get("length_bins", {}) if isinstance(qc.get("length_bins"), dict) else {}
        for key, shown in [("lt_200", "<200"), ("200_499", "200-499"), ("500_999", "500-999"), ("1000_1499", "1000-1499"), ("1500_1999", "1500-1999"), ("ge_2000", ">=2000")]:
            length_rows.append((f"{label} {shown}", integer(bins.get(key)), color))
    (fig_qc / f"{sample}.read_length_bins.svg").write_text(svg_bar(f"{sample}-read-length-bins", f"{sample} 三代 read length 分布", length_rows, "reads"), encoding="utf-8")

    q_rows = []
    for label, color, qc in (("raw", COLORS["raw"], raw_qc), ("insert", COLORS["insert"], insert_qc)):
        bins = qc.get("mean_q_bins", {}) if isinstance(qc.get("mean_q_bins"), dict) else {}
        for key, shown in [("lt_7", "<7"), ("7_9", "7-9"), ("10_14", "10-14"), ("15_19", "15-19"), ("20_29", "20-29"), ("ge_30", ">=30")]:
            q_rows.append((f"{label} Q{shown}", integer(bins.get(key)), color))
    (fig_qc / f"{sample}.quality_bins.svg").write_text(svg_bar(f"{sample}-quality-bins", f"{sample} read quality 分布", q_rows, "reads"), encoding="utf-8")

    (fig_qc / f"{sample}.alignment_assignment.svg").write_text(
        svg_bar(
            f"{sample}-alignment-assignment",
            f"{sample} 三代数据保留与分配",
            [
                ("Raw", total_reads, COLORS["raw"]),
                ("Insert", insert_reads, COLORS["insert"]),
                ("Primary Q20", primary_q20, COLORS["other"]),
                ("Gene assigned", assigned, COLORS["new"]),
            ],
            "reads",
        ),
        encoding="utf-8",
    )
    (fig_joint / f"{sample}.ambiguous_reason.svg").write_text(
        svg_horizontal(
            f"{sample}-ambiguous-reason",
            f"{sample} ambiguous reason",
            [(row["reason"], float(row["reads"]), COLORS["assignment"] if row["reason_class"] == "assignment_unresolved" else COLORS["slam"]) for row in reason_rows[:12]],
            "reads",
        ),
        encoding="utf-8",
    )
    (fig_joint / f"{sample}.read_length_call.svg").write_text(
        svg_stacked_percent(f"{sample}-read-length-call", f"{sample} read length 与 SLAM call", LENGTH_BINS, length_call_counts),
        encoding="utf-8",
    )
    (fig_joint / f"{sample}.t_covered_call.svg").write_text(
        svg_stacked_percent(f"{sample}-t-covered-call", f"{sample} T coverage 与 SLAM call", T_COV_BINS, tcov_call_counts),
        encoding="utf-8",
    )
    length_mut_rows = []
    for key, label in LENGTH_BINS:
        mut = mutation_by_len.get(key, Counter())
        covered = max(mut.get("t_covered", 0), 1)
        length_mut_rows.append((f"{label} T>C", mut.get("tc", 0) / covered, COLORS["new"]))
        length_mut_rows.append((f"{label} T>A", mut.get("ta", 0) / covered, COLORS["old"]))
        length_mut_rows.append((f"{label} T>G", mut.get("tg", 0) / covered, COLORS["other"]))
    (fig_joint / f"{sample}.length_mutation_rates.svg").write_text(
        svg_bar(f"{sample}-length-mutation-rates", f"{sample} read length bin 突变率", length_mut_rows, "rate"),
        encoding="utf-8",
    )

    return longread_rows, reason_rows, joint_rows, summary


def main() -> None:
    args = parse_args()
    samples = [item.strip() for item in args.samples.split(",") if item.strip()]
    longread_rows: list[dict[str, object]] = []
    reason_rows: list[dict[str, object]] = []
    joint_rows: list[dict[str, object]] = []
    summary_rows: list[dict[str, object]] = []
    for sample in samples:
        lr, reasons, joint, summary = summarize_sample(args.results_dir, sample)
        longread_rows.extend(lr)
        reason_rows.extend(reasons)
        joint_rows.extend(joint)
        summary_rows.append({"sample": sample, **summary})

    write_tsv(
        args.results_dir / "tables" / "lr_bulk_slam" / "longread_qc" / "lr_bulk_slam_longread_qc.tsv",
        [
            "sample",
            "raw_reads",
            "insert_reads",
            "raw_qc_reads",
            "insert_qc_reads",
            "fastq_qc_max_reads",
            "insert_rate",
            "raw_n50",
            "insert_n50",
            "raw_mean_length",
            "insert_mean_length",
            "raw_q20_base_fraction",
            "insert_q20_base_fraction",
            "raw_q30_base_fraction",
            "insert_q30_base_fraction",
            "raw_gc_fraction",
            "insert_gc_fraction",
            "raw_n_fraction",
            "insert_n_fraction",
            "align_total",
            "align_primary",
            "align_supplementary",
            "primary_q20_reads",
            "primary_q20_rate",
            "assigned_reads",
            "assigned_rate",
        ],
        longread_rows,
    )
    write_tsv(
        args.results_dir / "tables" / "lr_bulk_slam" / "longread_slam" / "lr_bulk_slam_ambiguous_reasons.tsv",
        ["sample", "reason", "reads", "reason_class"],
        reason_rows,
    )
    write_tsv(
        args.results_dir / "tables" / "lr_bulk_slam" / "longread_slam" / "lr_bulk_slam_longread_slam_bins.tsv",
        ["sample", "metric", "bin", "total_reads", "new_reads", "old_reads", "ambiguous_reads", "tc_rate", "ta_rate", "tg_rate"],
        joint_rows,
    )
    write_tsv(
        args.results_dir / "tables" / "lr_bulk_slam" / "longread_slam" / "lr_bulk_slam_callable_summary.tsv",
        ["sample", "assignment_unresolved", "slam_ambiguous", "callable_total", "callable_ambiguous_fraction"],
        summary_rows,
    )


if __name__ == "__main__":
    main()
