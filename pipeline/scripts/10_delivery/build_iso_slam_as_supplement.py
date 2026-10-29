#!/usr/bin/env python3
"""Build a standalone, gene-annotated alternative-splicing supplement."""

from __future__ import annotations

import argparse
import csv
import getpass
import json
import os
import re
import shutil
import subprocess
import sys
import tempfile
from collections import Counter
from pathlib import Path

ROOT = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(ROOT / "pipeline/scripts"))
from gene_annotation import fields_with_gene_name, load_gtf_gene_map  # noqa: E402
from iso_slam_v4_common import read_config  # noqa: E402


TABLES = (
    "differential_new_rna_splicing.tsv",
    "sample_event_psi.tsv",
    "novel_structure_inventory.tsv",
)
ENSEMBL_ID = re.compile(r"ENSG\d+(?:\.\d+)?")
BASE_FIGURES = (
    ("AS_contrast_overview.png", "比较检验数和显著事件数概览（PNG）"),
    ("AS_contrast_overview.pdf", "比较检验数和显著事件数概览（PDF）"),
    ("AS_significant_delta_psi.png", "全部显著比较事件的 ΔPSI（PNG）"),
    ("AS_significant_delta_psi.pdf", "全部显著比较事件的 ΔPSI（PDF）"),
)


def resolve_gene_name(gene_id: str, gene_map: dict[str, str]) -> tuple[str, bool]:
    value = gene_id.strip()
    if not value:
        return "NA", False
    if value in gene_map:
        return gene_map[value], True
    components = ENSEMBL_ID.findall(value)
    if len(components) > 1:
        return ";".join(gene_map.get(gene, gene) for gene in components), all(gene in gene_map for gene in components)
    return gene_map.get(value, value), False


def annotate_table(source: Path, target: Path, gene_map: dict[str, str]) -> dict[str, int]:
    counts = Counter()
    with source.open(encoding="utf-8", newline="") as source_handle, target.open("w", encoding="utf-8", newline="") as target_handle:
        reader = csv.DictReader(source_handle, delimiter="\t")
        if not reader.fieldnames or "gene_id" not in reader.fieldnames:
            raise RuntimeError(f"Expected a gene_id field in {source}")
        fields = fields_with_gene_name(reader.fieldnames)
        writer = csv.DictWriter(target_handle, fieldnames=fields, delimiter="\t", lineterminator="\n", extrasaction="ignore")
        writer.writeheader()
        for row in reader:
            gene_id = row.get("gene_id", "")
            gene_name, annotated = resolve_gene_name(gene_id, gene_map)
            row["gene_name"] = gene_name
            counts["rows"] += 1
            counts["unmapped_gene_ids"] += int(not annotated)
            counts["composite_gene_ids"] += int(len(ENSEMBL_ID.findall(gene_id)) > 1)
            counts["empty_gene_names"] += int(not gene_name)
            writer.writerow(row)
    if counts["rows"] == 0 or counts["empty_gene_names"]:
        raise RuntimeError(f"Gene-name annotation failed for {source}: {dict(counts)}")
    return dict(counts)


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open(encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def write_summary(path: Path, validation: dict[str, object], differential: list[dict[str, str]], psi_page_count: int) -> None:
    significant = [row for row in differential if row.get("significant", "").lower() == "true"]
    candidate_event_count = len({row["event_id"] for row in significant})
    significant.sort(key=lambda row: (float(row.get("FDR") or 1), -abs(float(row.get("delta_PSI") or 0))))
    lines = [
        "# 新生 RNA 可变剪切分析补充结果",
        "",
        "本目录汇总新生 RNA 可变剪切分析结果、逐样本 PSI、差异统计和图形。",
        "",
        "## 分析概况",
        "",
        f"- 参考注释事件数：{int(validation.get('reference_event_count', 0)):,}。",
        f"- 事件类型：{', '.join(validation.get('event_types', []))}。",
        f"- 新生 RNA 支持门槛：每样本至少 {validation.get('support_threshold_expected_new_read_equivalents_per_sample', 'NA')} 个期望 reads 等价数。",
        f"- 显著阈值：{validation.get('significance_threshold', 'NA')}。",
        f"- 差异方法：{validation.get('method', 'NA')}。",
        f"- 每组生物学重复数：{validation.get('biological_replicates_per_group', 'NA')}；结果用于探索性筛查。",
        "",
        "## 各比较结果",
        "",
        "| 比较 | 进入检验事件数 | 达到双阈值事件数 |",
        "|---|---:|---:|",
    ]
    tested = validation.get("tested_event_count_by_contrast", {})
    sig_counts = validation.get("significant_event_count_by_contrast", {})
    for contrast in tested:
        lines.append(f"| {contrast.replace('_vs_', ' vs ')} | {tested[contrast]} | {sig_counts.get(contrast, 0)} |")
    lines.extend(["", "## 显著候选事件", "", "| 比较 | 事件 | 基因名 | Gene ID | ΔPSI | FDR |", "|---|---|---|---|---:|---:|"])
    for row in significant[:20]:
        lines.append("| {contrast} | {event} | {name} | {gene} | {delta:.3f} | {fdr:.3g} |".format(
            contrast=row.get("contrast", ""), event=row.get("event_id", ""), name=row.get("gene_name", ""),
            gene=row.get("gene_id", ""), delta=float(row.get("delta_PSI") or 0), fdr=float(row.get("FDR") or 1)))
    if not significant:
        lines.append("| 无 | — | — | — | — | — |")
    lines.extend([
        "",
        "## 结果图",
        "",
        "### 比较概览",
        "",
        "![各比较进入检验事件数和显著事件数](figures/AS_contrast_overview.png)",
        "",
        "[矢量 PDF](figures/AS_contrast_overview.pdf)。上方面板为进入检验事件数，下方面板为达到双阈值的事件数。",
        "",
        "### 显著事件 ΔPSI",
        "",
        "![全部显著比较事件的 ΔPSI](figures/AS_significant_delta_psi.png)",
        "",
        "[矢量 PDF](figures/AS_significant_delta_psi.pdf)。ΔPSI = group 1 PSI − group 2 PSI；颜色表示事件类型，候选详情见差异结果表。",
        "",
        "### 候选事件的样本级 PSI",
        "",
        f"共展示 {candidate_event_count} 个独立候选事件。每个点代表一个样本，形状区分重复；NA 未绘制且未按 0 处理。",
        "",
    ])
    for page in range(1, psi_page_count + 1):
        lines.extend([f"![候选事件样本级 PSI，第 {page} 页](figures/AS_candidate_sample_psi_page{page}.png)", ""])
    lines.extend([
        f"[{psi_page_count} 页矢量 PDF](figures/AS_candidate_sample_psi.pdf)。",
        "",
        "## 解读边界",
        "",
        "每组 n=2，显著候选需结合重复一致性及后续实验验证。新发现转录本结构见 `novel_structure_inventory.tsv`；现有结果无法可靠估计这些结构逐读段的新生 RNA 后验，因此未将其纳入新生 RNA PSI 差异检验。",
        "基因名按参考 GTF 的 gene_name 字段展示；未提供区别于 Ensembl ID 的符号时保留该 ID，不推断名称。无法匹配的 novel/复合 ID 也保留原 ID。",
        "",
        "`sample_event_psi.tsv` 提供逐样本事件 PSI 与支持量；`differential_new_rna_splicing.tsv` 提供差异统计；`DELIVERY_VALIDATION.json` 记录本补充包检查结果。",
    ])
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def render_figures(differential_path: Path, sample_psi_path: Path, output_dir: Path, rscript: str) -> None:
    script = Path(__file__).with_name("plot_iso_slam_as_supplement.R")
    subprocess.run(
        [rscript, str(script), str(differential_path), str(sample_psi_path), str(output_dir)],
        check=True,
    )


def ensure_customer_package_desensitized(root: Path) -> None:
    forbidden = [b"/internal/", b"/work/", b"/internal/", b"/media/", b"/tmp/", b"file://", getpass.getuser().encode()]
    forbidden.extend([b"wangjiahe", b"z-wangjiahe_20260930"])
    path_pattern = re.compile(rb"(?i)(?:/(?:home|work|mnt|media|tmp|scratch|Users)/[^\s\"'<>]+)")
    for path in root.rglob("*"):
        if not path.is_file():
            continue
        relative = str(path.relative_to(root))
        name_bytes = relative.encode("utf-8").lower()
        if any(token.lower() in name_bytes for token in forbidden if token):
            raise RuntimeError(f"Delivery filename contains internal identity: {relative}")
        tail = b""
        with path.open("rb") as handle:
            while chunk := handle.read(1024 * 1024):
                scan = tail + chunk
                lowered = scan.lower()
                found = [token.decode("utf-8", errors="replace") for token in forbidden if token and token.lower() in lowered]
                if found or path_pattern.search(scan):
                    reason = ", ".join(found) if found else "absolute path"
                    raise RuntimeError(f"Delivery file contains internal identity or path: {reason} ({relative})")
                tail = scan[-512:]


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--config", required=True, type=Path)
    parser.add_argument("--outdir", required=True, type=Path)
    args = parser.parse_args()
    config = read_config(args.config)
    results_root = Path(config["results_root"])
    source = results_root / "_project/08_project_qc/02_group_analysis/tables/differential/new_rna_splicing"
    source_validation = source / "new_rna_splicing_validation.json"
    if not source_validation.is_file():
        raise RuntimeError(f"Completed new-RNA splicing results are missing: {source_validation}")
    validation = json.loads(source_validation.read_text(encoding="utf-8"))
    if validation.get("status") != "PASS":
        raise RuntimeError("New-RNA splicing validation status is not PASS")
    missing = [name for name in TABLES if not (source / name).is_file() or (source / name).stat().st_size == 0]
    if missing:
        raise RuntimeError(f"Missing or empty AS result tables: {', '.join(missing)}")
    if args.outdir.exists():
        raise RuntimeError(f"Refusing to overwrite existing supplement: {args.outdir}")
    args.outdir.parent.mkdir(parents=True, exist_ok=True)
    tmp_root = Path(config["tmp_root"]).resolve()
    tmp_root.mkdir(parents=True, exist_ok=True)
    stage = Path(tempfile.mkdtemp(prefix="as_supplement_", dir=tmp_root))
    try:
        gene_map = load_gtf_gene_map(Path(config["gtf"]))
        annotation_counts = {}
        for filename in TABLES:
            annotation_counts[filename] = annotate_table(source / filename, stage / filename, gene_map)
        shutil.copy2(source_validation, stage / "new_rna_splicing_validation.json")
        differential_rows = read_tsv(stage / TABLES[0])
        figure_dir = stage / "figures"
        figure_dir.mkdir()
        render_figures(stage / TABLES[0], stage / TABLES[1], figure_dir, str(config.get("rscript_bin", "Rscript")))
        page_files = sorted(figure_dir.glob("AS_candidate_sample_psi_page*.png"))
        figures = [*BASE_FIGURES]
        figures.extend((path.name, f"候选事件逐样本 PSI，第 {index} 页（PNG）") for index, path in enumerate(page_files, start=1))
        if page_files:
            figures.append(("AS_candidate_sample_psi.pdf", f"候选事件逐样本 PSI（{len(page_files)} 页 PDF）"))
        missing_figures = [name for name, _ in figures if not (figure_dir / name).is_file() or (figure_dir / name).stat().st_size == 0]
        if missing_figures:
            raise RuntimeError(f"Missing generated figures: {', '.join(missing_figures)}")
        write_summary(stage / "AS_结果摘要.md", validation, differential_rows, len(page_files))
        summary = {
            "status": "PASS",
            "package_type": "alternative_splicing_supplement",
            "source_result_count": len(TABLES),
            "gene_annotation_source": "GRCh38-2024-A gene annotation",
            "gene_name_column": "gene_name",
            "gene_annotation_counts": annotation_counts,
            "previous_full_delivery_rebuilt": False,
            "figure_count": len(figures),
            "significant_comparison_row_count": sum(row.get("significant", "").lower() == "true" for row in differential_rows),
            "distinct_significant_event_count": len({row["event_id"] for row in differential_rows if row.get("significant", "").lower() == "true"}),
            "privacy_scan_status": "PASS",
        }
        (stage / "DELIVERY_VALIDATION.json").write_text(json.dumps(summary, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")
        manifest_rows = [
            {"file": filename, "description": description, "status": "included"}
            for filename, description in [
                ("AS_结果摘要.md", "可变剪切结果、阈值和解读限制"),
                (TABLES[0], "组间新生 RNA 可变剪切事件差异结果，含基因名"),
                (TABLES[1], "逐样本新生 RNA PSI 与事件支持量，含基因名"),
                (TABLES[2], "新发现转录本结构清单，含基因名"),
                ("new_rna_splicing_validation.json", "分析统计与筛选阈值"),
                ("DELIVERY_VALIDATION.json", "补充交付包验证信息"),
            ]
        ]
        manifest_rows.extend(
            {"file": f"figures/{filename}", "description": description, "status": "included"}
            for filename, description in figures
        )
        with (stage / "delivery_manifest.tsv").open("w", encoding="utf-8", newline="") as handle:
            writer = csv.DictWriter(handle, fieldnames=["file", "description", "status"], delimiter="\t", lineterminator="\n")
            writer.writeheader()
            writer.writerows(manifest_rows)
        ensure_customer_package_desensitized(stage)
        os.replace(stage, args.outdir)
    except Exception:
        shutil.rmtree(stage, ignore_errors=True)
        raise
    print(json.dumps({"status": "PASS", "delivery_root": str(args.outdir.resolve()), "files": len(manifest_rows), "gene_annotation_counts": annotation_counts}, ensure_ascii=False))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
