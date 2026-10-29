#!/usr/bin/env python3
"""Build an atomic, copy-only customer delivery package from completed results."""

from __future__ import annotations

import argparse
import csv
import getpass
import gzip
import hashlib
import json
import os
import re
import shutil
import subprocess
import sys
from collections import Counter
from datetime import datetime
from html.parser import HTMLParser
from pathlib import Path


ROOT = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(ROOT / "pipeline/scripts"))
sys.path.insert(0, str(ROOT / "pipeline/scripts/10_delivery"))
from iso_slam_v4_common import as_bool, read_config, write_json  # noqa: E402
from gene_annotation import add_gene_names, fields_with_gene_name, load_gtf_gene_map  # noqa: E402
from build_iso_slam_as_supplement import (  # noqa: E402
    BASE_FIGURES as AS_BASE_FIGURES,
    TABLES as AS_TABLES,
    annotate_table as annotate_as_table,
    read_tsv as read_as_tsv,
    render_figures as render_as_figures,
    write_summary as write_as_summary,
)


DETAILED_TABLES = (
    ("isoform_catalog.tsv", "04_isoform_catalog/01_tables/{sample}.isoform_catalog.tsv", "02_Transcript_Catalog_and_Structure_QC/{sample}/isoform_catalog.tsv"),
    ("gene_rna_labeling.tsv", "05_slam_model/01_tables/{sample}.grand3_gene_ntr.tsv", "03_Expression_Quantification_and_Matrices/{sample}/gene_rna_labeling.tsv"),
    ("global_rna_labeling_parameters.tsv", "05_slam_model/01_tables/{sample}.grand3_global_parameters.tsv", "03_Expression_Quantification_and_Matrices/{sample}/global_rna_labeling_parameters.tsv"),
    ("gene_dynamics.tsv", "06_isoform_dynamics/01_tables/{sample}.gene_dynamics.tsv", "03_Expression_Quantification_and_Matrices/{sample}/gene_dynamics.tsv"),
    ("isoform_dynamics.tsv", "06_isoform_dynamics/01_tables/{sample}.isoform_dynamics.tsv", "03_Expression_Quantification_and_Matrices/{sample}/isoform_dynamics.tsv"),
    ("isoform_usage_candidates.tsv", "06_isoform_dynamics/01_tables/{sample}.isoform_switch_candidates.tsv", "05_Isoform_Usage/{sample}/isoform_usage_candidates.tsv"),
)

LAYOUT = {
    "00": "00_Report_and_Figures",
    "01": "01_QC_and_Preprocessing",
    "02": "02_Transcript_Catalog_and_Structure_QC",
    "03": "03_Expression_Quantification_and_Matrices",
    "04": "04_Differential_Analysis",
    "05": "05_Isoform_Usage",
    "06": "06_Methods_and_Run_Metadata",
}


def group_delivery_target(relative: Path) -> Path | None:
    if relative.parts[0] == "figures":
        return Path(LAYOUT["00"]) / "figures" / relative.name
    if relative.name == "group_analysis_validation.json":
        return None
    if relative.parts[0] != "tables":
        return None
    if relative.name == "contrast_summary.tsv":
        return Path(LAYOUT["04"]) / relative.name
    if "differential" in relative.parts:
        return Path(LAYOUT["04"]) / Path(*relative.parts[1:])
    return Path(LAYOUT["01"]) / relative.name

PROJECT_QC_FIELDS = (
    "sample_id",
    "raw_reads",
    "insert_reads",
    "insert_rate",
    "isoform_assigned_reads",
    "dynamics_eligible_reads",
    "no_t_sites_reads",
    "assignment_conflict_reads",
    "assignment_conflict_rate",
    "catalog_isoforms",
    "filtered_isoforms",
    "grand3_global_binom_ntr",
    "posterior_mean",
    "posterior_grand3_absolute_difference",
    "switch_candidates",
)


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open(encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def write_tsv(path: Path, fields: list[str], rows: list[dict[str, object]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", lineterminator="\n", extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)


def digest(path: Path) -> str:
    value = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            value.update(block)
    return value.hexdigest()


def copy_file(source: Path, target: Path) -> None:
    if source.is_symlink() or not source.is_file() or source.stat().st_size == 0:
        raise RuntimeError(f"Delivery source must be a nonempty regular file: {source}")
    target.parent.mkdir(parents=True, exist_ok=True)
    shutil.copy2(source, target, follow_symlinks=False)
    if target.is_symlink() or not target.is_file() or target.stat().st_size != source.stat().st_size:
        raise RuntimeError(f"Delivery copy failed: {source} -> {target}")


def copy_result_table(source: Path, target: Path, gene_map: dict[str, str]) -> None:
    rows = read_tsv(source)
    if not rows or "gene_id" not in rows[0]:
        copy_file(source, target)
        return
    fields = fields_with_gene_name(list(rows[0]))
    add_gene_names(rows, gene_map)
    write_tsv(target, fields, rows)


def copy_tree(source: Path, target: Path) -> int:
    copied = 0
    for path in sorted(source.rglob("*")):
        if path.is_symlink():
            raise RuntimeError(f"Delivery source tree contains a symbolic link: {path}")
        if path.is_file():
            copy_file(path, target / path.relative_to(source))
            copied += 1
    return copied


class VisibleText(HTMLParser):
    def __init__(self) -> None:
        super().__init__(convert_charrefs=True)
        self.hidden = 0
        self.parts: list[str] = []

    def handle_starttag(self, tag: str, attrs: list[tuple[str, str | None]]) -> None:
        if tag in {"script", "style"}:
            self.hidden += 1

    def handle_endtag(self, tag: str) -> None:
        if tag in {"script", "style"} and self.hidden:
            self.hidden -= 1

    def handle_data(self, data: str) -> None:
        if not self.hidden:
            self.parts.append(data)


def ensure_customer_text(path: Path) -> None:
    parser = VisibleText()
    parser.feed(path.read_text(encoding="utf-8"))
    text = " ".join(parser.parts)
    forbidden = ("V4", "v4", "No_T_sites", "Assignment conflict", "No T sites", "修复分析")
    found = [item for item in forbidden if item in text]
    if found:
        raise RuntimeError(f"Customer report contains internal wording: {', '.join(found)} ({path})")


REVIEW_REQUIRED_FIELDS = {
    "sample_id", "raw_reads", "raw_reads_status", "fastq_q20_base_fraction", "fastq_q20_status",
    "primary_mapq20_fraction", "primary_mapq20_status", "genes_detected", "genes_detected_status",
    "insert_retention_rate", "insert_mapping_rate", "overall_review_status",
}


def validate_review_table(path: Path, sample_ids: list[str]) -> dict[str, int]:
    rows = read_tsv(path)
    if not rows:
        raise RuntimeError(f"QC review table is empty: {path}")
    missing_fields = REVIEW_REQUIRED_FIELDS - set(rows[0])
    if missing_fields:
        raise RuntimeError(f"QC review table is missing fields: {sorted(missing_fields)}")
    observed = [row["sample_id"] for row in rows]
    if len(observed) != len(set(observed)) or set(observed) != set(sample_ids):
        raise RuntimeError("QC review table sample IDs do not match the sample sheet")
    value_fields = REVIEW_REQUIRED_FIELDS - {
        "sample_id", "raw_reads_status", "fastq_q20_status", "primary_mapq20_status",
        "genes_detected_status", "overall_review_status",
    }
    for row in rows:
        for key in value_fields:
            if not row.get(key, "").strip():
                raise RuntimeError(f"QC review table has an empty {key} value for {row['sample_id']}")
    return {
        "sample_count": len(rows),
        "pass_sample_count": sum(row["overall_review_status"] == "PASS" for row in rows),
        "below_threshold_sample_count": sum(row["overall_review_status"] == "BELOW_THRESHOLD" for row in rows),
    }


def number_range(values: list[float], digits: int = 0, percent: bool = False) -> str:
    low, high = min(values), max(values)
    if percent:
        return f"{low * 100:.{digits}f}%–{high * 100:.{digits}f}%"
    return f"{low:,.{digits}f}–{high:,.{digits}f}"


def build_analysis_summary(
    overview: list[dict[str, object]], calls_rows: list[dict[str, object]], candidates: list[dict[str, object]],
    group_rows: list[dict[str, str]], group_validation: dict[str, object], include_sample_reports: bool,
    splicing_validation: dict[str, object] | None = None,
) -> str:
    splicing_validation = splicing_validation or {}
    raw_reads = [float(row["raw_reads"]) for row in overview]
    insert_rates = [float(row["insert_rate"]) for row in overview]
    assigned_reads = [float(row["isoform_assigned_reads"]) for row in overview]
    isoforms = [float(row["filtered_isoforms"]) for row in overview]
    new_rna = [float(row["global_new_rna_fraction"]) for row in overview]
    screen_rates = [float(row["transcript_assignment_screen_rate"]) for row in overview]
    calls = {key: sum(int(float(row.get(key, 0) or 0)) for row in calls_rows) for key in (
        "new_rna_reads", "old_rna_reads", "ambiguous_reads", "not_included_after_assignment_screen", "reads_without_usable_t_sites",
    )}
    no_t_note = ""
    if calls["reads_without_usable_t_sites"]:
        no_t_note = f"\n- {calls['reads_without_usable_t_sites']:,} 条 read 未具备可用于标记判定的 T 位点，未纳入相应的判定。"
    group_text = ""
    sample_reports_line = "- `sample_reports/`：各样本的独立离线报告。\n" if include_sample_reports else ""
    if group_rows:
        by_contrast: dict[str, int] = {}
        for row in group_rows:
            by_contrast[row["contrast"]] = by_contrast.get(row["contrast"], 0) + int(row["significant_features"])
        contrasts = group_validation.get("contrasts", list(by_contrast))
        lines = "\n".join(f"- {contrast.replace('_vs_', ' vs ')}：各结果层级合计 {by_contrast.get(contrast, 0):,} 条显著记录。" for contrast in contrasts)
        groups = group_validation.get("groups", [])
        replicate_counts = group_validation.get("biological_replicates_per_group", {})
        replicate_summary = "、".join(f"{group} n={replicate_counts.get(group, '未知')}" for group in groups)
        design_summary = f"{len(groups)} 个独立组（{replicate_summary}），完成 {len(contrasts)} 项预设比较。"
        splice_counts = splicing_validation.get("significant_event_count_by_contrast", {})
        splice_lines = "\n".join(f"- {contrast.replace('_vs_', ' vs ')}：{count} 个新生 RNA 可变剪切事件达到统计与效应阈值。" for contrast, count in splice_counts.items())
        group_text = f"""

## 组内与组间分析

{design_summary}组内一致性使用 RNA 标记信号、基因丰度相关性和 PCA 评估；组间结果覆盖基因丰度、基因标记比例、转录本丰度、总 RNA 转录本使用和新生 RNA 转录本使用。

{lines}

显著记录数按比较和分析层级累加，同一特征可能在多个比较中重复出现；完整效应量、p 值和 FDR 见 `04_Differential_Analysis/differential/`。

新生 RNA 可变剪切分析纳入 {splicing_validation.get('reference_event_count', 0):,} 个参考注释事件；每组 n=2，结果为探索性筛查。每组比较的显著事件数：

{splice_lines}

新发现转录本结构列于 `04_Differential_Analysis/differential/new_rna_splicing/novel_structure_inventory.tsv`；由于现有流程结果不能可靠分配逐读段新生 RNA 后验，这些结构未纳入新生 RNA PSI 差异检验。

可变剪切结果摘要、带基因名的结果表及图形见 `04_Differential_Analysis/differential/new_rna_splicing/AS_结果摘要.md`。
"""
    return f"""# 长读长 RNA 标记分析结果摘要

主报告：[ISO-SLAM 离线项目报告](00_Report_and_Figures/ISO_SLAM_project_report.html)。文件清单见 [delivery_manifest.tsv](delivery_manifest.tsv)。交付内容按 00–06 编号分区，组别和样本名保留原分析标签。

## 数据范围

本交付包汇总 {len(overview)} 个独立样本的长读长 RNA 标记分析结果。各样本独立处理；结果用于技术质量和描述性信号汇总。

## 主要结果

- 原始 reads 数范围为 {number_range(raw_reads)}；插入序列保留率为 {number_range(insert_rates, 2, True)}。
- 通过筛选并分配至转录本的 reads 数范围为 {number_range(assigned_reads)}；每个样本保留的转录本数量为 {number_range(isoforms)}。
- 全局新生 RNA 比例范围为 {number_range(new_rna, 2, True)}，反映在既定分析筛选条件下的 RNA 标记信号。
- 转录本分配一致性筛查比例为 {number_range(screen_rates, 2, True)}。未通过该筛查的 read 不纳入转录本层面的新生/已有 RNA 汇总。
- 汇总后，新生 RNA read 为 {calls['new_rna_reads']:,}，已有 RNA read 为 {calls['old_rna_reads']:,}，判定不确定 read 为 {calls['ambiguous_reads']:,}；另有 {calls['not_included_after_assignment_screen']:,} 条 read 未纳入转录本层面汇总。{no_t_note}
- `05_Isoform_Usage/candidate_isoform_usage_by_sample.tsv` 汇集了 {len(candidates):,} 条按样本筛选的重点候选转录本使用变化记录；每个样本的全量候选结果见 `05_Isoform_Usage/<sample_id>/`，供后续复核与实验验证。

## 结果解读范围

- 单样本候选结果不单独用于组间显著性推断；正式组间统计见本摘要后续章节及 `04_Differential_Analysis/`。
- 本结果不估计 RNA 半衰期、绝对合成速率或降解速率。
- 候选转录本使用变化为样本内描述性筛查结果，建议结合生物学重复或正交实验进一步验证。
- 各项比例均基于对应分析步骤的筛选 read 集合；不同指标的分母可能不同，不宜直接互换比较。

## 文件说明

- `RNA_Labeling_Long_Read_Report.html`：离线项目总报告。
{sample_reports_line}- `tables/`：项目级汇总和重点结果摘录。
- `detailed_tables/`：各样本的全量结果表。
- `group_analysis/`：组内一致性、组间差异表和图形。
- `TABLE_GUIDE.md`：结果表说明。
{group_text}"""


def build_table_guide() -> str:
    return """# 结果表说明

## 01_QC_and_Preprocessing

- `01_QC_and_Preprocessing/sample_qc_review.tsv`：逐样本审核指标及阈值判断；审核标准见同目录 `REVIEW_CRITERIA.md`。
- `01_QC_and_Preprocessing/project_qc_detailed.tsv`：各样本的测序、转录本分配、RNA 标记与动态分析质量指标；不包含本地来源路径。
- `01_QC_and_Preprocessing/sample_overview.tsv`、`01_QC_and_Preprocessing/read_classification_summary.tsv`：用于快速核对样本级概览和 read 分类。

## 02–05 全量结果表

各样本结果分别位于 `02_Transcript_Catalog_and_Structure_QC/<sample_id>/`、`03_Expression_Quantification_and_Matrices/<sample_id>/` 和 `05_Isoform_Usage/<sample_id>/`：

- `isoform_catalog.tsv`：转录本结构、支持 reads 和筛选状态。
- `gene_rna_labeling.tsv`、`global_rna_labeling_parameters.tsv`：基因层面的 RNA 标记比例估计及全局模型参数。
- `gene_dynamics.tsv`、`isoform_dynamics.tsv`：基因/转录本层面的新生与已有 RNA 估计及使用量指标。
- `isoform_usage_candidates.tsv`：全量转录本使用变化候选；该结果为样本内描述性筛查，不代表组间统计结论。

## 解读提示

各表的比例和 reads 计数基于对应步骤的筛选集合，分母可能不同。单样本表不单独用于组间显著性推断；正式组间结果见 `04_Differential_Analysis/`。所有结果均不用于推断 RNA 半衰期或绝对动力学参数。

## 组内与组间结果

- `01_QC_and_Preprocessing/`：样本相关性、PCA 坐标、组内信号汇总等质量评估结果。
- `04_Differential_Analysis/`：比较汇总及五类完整差异结果。
- `04_Differential_Analysis/differential/new_rna_splicing/`：新生 RNA 事件级 PSI、差异结果和新结构清单；差异模型采用加权 limma/eBayes。
- 同目录 `AS_结果摘要.md` 和 `figures/`：可变剪切结果摘要及图形；三张 AS 结果表包含 `gene_name`。
- 同目录 `AS_结果摘要.md`：可变剪切结果摘要和图表；三张 AS 结果表包含基因名。
- `00_Report_and_Figures/figures/`：组内分布、相关性、PCA 和差异汇总图；PNG 用于浏览，PDF 用于矢量输出。
- 丰度差异使用 edgeR；比例差异使用 read 支持加权的 limma。完整阈值和模型限制见项目总报告。
"""


def ensure_delivery_desensitized(root: Path) -> None:
    forbidden = (b"/internal/", b"/work/", b"/internal/", b"/media/", b"/tmp/", b"file://", getpass.getuser().encode())
    path_pattern = re.compile(rb"(?i)(?:/(?:home|work|mnt|media|tmp|scratch|Users)/[^\s\"'<>]+)")
    data_uri_pattern = re.compile(rb"(?i)data:[^;,\s\"']+;base64,[a-z0-9+/=\r\n]+")
    text_suffixes = {".html", ".md", ".tsv", ".json", ".csv", ".txt", ".yaml", ".yml"}
    for path in root.rglob("*"):
        if not path.is_file():
            continue
        if path.suffix == ".gz":
            with gzip.open(path, "rb") as handle:
                content = handle.read()
        else:
            content = path.read_bytes()
        if path.suffix.lower() == ".html":
            content = data_uri_pattern.sub(b"", content)
        if path.suffix.lower() in text_suffixes or path.suffix.lower() == ".gz":
            scan_chunks = [content]
        else:
            strings = shutil.which("strings")
            if strings:
                process = subprocess.Popen([strings, "-a", "-n", "6", str(path)], stdout=subprocess.PIPE, stderr=subprocess.DEVNULL)
                assert process.stdout is not None
                scan_chunks = process.stdout
            else:
                scan_chunks = [content[:1024 * 1024]]
        found: list[str] = []
        for chunk in scan_chunks:
            lowered = chunk.lower()
            found.extend(item.decode("utf-8", errors="replace") for item in forbidden if item.lower() in lowered)
            if path_pattern.search(chunk):
                found.append("absolute path")
        if path.suffix.lower() not in text_suffixes and path.suffix.lower() != ".gz" and strings:
            process.stdout.close()
            process.wait()
        if found:
            raise RuntimeError(f"Delivery file contains local identity or paths: {', '.join(found)} ({path})")


AS_VALIDATION_FIELDS = (
    "status", "method", "event_types", "reference_event_count", "tested_event_count_by_contrast",
    "significant_event_count_by_contrast", "novel_catalog_records", "new_RNA_novel_event_differential",
    "support_threshold_expected_new_read_equivalents_per_sample", "significance_threshold",
    "biological_replicates_per_group", "interpretation",
)


def package_as_results(
    group_source: Path,
    staging: Path,
    gene_map: dict[str, str],
    rscript: str,
) -> dict[str, object]:
    as_source = group_source / "tables/differential/new_rna_splicing"
    validation_path = as_source / "new_rna_splicing_validation.json"
    if not validation_path.is_file():
        raise RuntimeError(f"Completed new-RNA splicing results are missing: {validation_path}")
    validation = json.loads(validation_path.read_text(encoding="utf-8"))
    if validation.get("status") != "PASS":
        raise RuntimeError(f"New-RNA splicing validation did not pass: {validation_path}")
    missing = [name for name in AS_TABLES if not (as_source / name).is_file() or not (as_source / name).stat().st_size]
    if missing:
        raise RuntimeError(f"Missing or empty AS result tables: {', '.join(missing)}")

    as_target = staging / LAYOUT["04"] / "differential/new_rna_splicing"
    annotation_counts = {}
    for filename in AS_TABLES:
        annotated_path = as_target / f".{filename}.annotating"
        annotation_counts[filename] = annotate_as_table(as_target / filename, annotated_path, gene_map)
        os.replace(annotated_path, as_target / filename)
    safe_validation = {key: validation[key] for key in AS_VALIDATION_FIELDS if key in validation}
    write_json(as_target / "new_rna_splicing_validation.json", safe_validation)

    differential_rows = read_as_tsv(as_target / AS_TABLES[0])
    figure_dir = as_target / "figures"
    figure_dir.mkdir(parents=True, exist_ok=True)
    render_as_figures(as_target / AS_TABLES[0], as_target / AS_TABLES[1], figure_dir, rscript)
    page_files = sorted(figure_dir.glob("AS_candidate_sample_psi_page*.png"))
    figures = [*AS_BASE_FIGURES]
    figures.extend((path.name, f"候选事件逐样本 PSI，第 {index} 页（PNG）") for index, path in enumerate(page_files, start=1))
    if page_files:
        figures.append(("AS_candidate_sample_psi.pdf", f"候选事件逐样本 PSI（{len(page_files)} 页 PDF）"))
    missing_figures = [name for name, _ in figures if not (figure_dir / name).is_file() or not (figure_dir / name).stat().st_size]
    if missing_figures:
        raise RuntimeError(f"Missing generated AS figures: {', '.join(missing_figures)}")
    write_as_summary(as_target / "AS_结果摘要.md", safe_validation, differential_rows, len(page_files))
    significant = [row for row in differential_rows if row.get("significant", "").lower() == "true"]
    return {
        "status": "PASS",
        "gene_annotation_source": "configured gene annotation GTF",
        "gene_annotation_counts": annotation_counts,
        "figure_count": len(figures),
        "significant_comparison_row_count": len(significant),
        "distinct_significant_event_count": len({row["event_id"] for row in significant}),
        "validation_status": safe_validation["status"],
    }


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--config", required=True, type=Path)
    parser.add_argument("--replace", action="store_true")
    args = parser.parse_args()
    config = read_config(args.config)
    gene_map = load_gtf_gene_map(Path(config["gene_annotation_gtf"])) if config.get("gene_annotation_gtf") else {}
    results_root, tmp_root, destination = Path(config.get("report_results_root", config["results_root"])), Path(config["tmp_root"]), Path(config["delivery_root"])
    complete = results_root / "_project/RUN_COMPLETE.tsv"
    if not complete.is_file() or "completed" not in complete.read_text(encoding="utf-8"):
        raise RuntimeError(f"Completed run marker is required: {complete}")
    if destination.exists() and not args.replace:
        raise RuntimeError(f"Delivery destination already exists; rerun with --replace: {destination}")
    sample_rows = read_tsv(Path(config["samples_tsv"]))
    sample_ids = [row["sample_id"] for row in sample_rows]
    sample_group = {row["sample_id"]: row.get("group", "") for row in sample_rows}
    include_sample_reports = as_bool(config, "include_sample_reports", True)
    token = datetime.now().strftime("%Y%m%d%H%M%S%f")
    staging = tmp_root / f"delivery_package_stage_{token}"
    if staging.exists():
        raise RuntimeError(f"Refusing to reuse staging directory: {staging}")
    staging.mkdir(parents=True)
    project_report = results_root / "_project/09_project_report/report.html"
    copy_file(project_report, staging / "RNA_Labeling_Long_Read_Report.html")
    ensure_customer_text(staging / "RNA_Labeling_Long_Read_Report.html")
    review_table_source = results_root / "_project/09_project_report/01_tables/sample_qc_review.tsv"
    review_table_target = staging / f"{LAYOUT['01']}/sample_qc_review.tsv"
    copy_file(review_table_source, review_table_target)
    review_summary = validate_review_table(review_table_target, sample_ids)
    copy_file(ROOT / "docs/result_review_criteria.md", staging / f"{LAYOUT['01']}/REVIEW_CRITERIA.md")
    report_parser = VisibleText()
    report_parser.feed((staging / "RNA_Labeling_Long_Read_Report.html").read_text(encoding="utf-8"))
    report_text = " ".join(report_parser.parts)
    required_review_labels = ("逐样本审核指标", "FASTQ 碱基 Q20", "MAPQ≥20", "检测基因数")
    if any(label not in report_text for label in required_review_labels):
        raise RuntimeError("Project HTML report does not expose all required QC review metrics")
    qc_rows = read_tsv(results_root / "_project/08_project_qc/01_tables/sample_qc_summary.tsv")
    qc_by_sample = {row["sample_id"]: row for row in qc_rows}
    overview: list[dict[str, object]] = []
    calls_rows: list[dict[str, object]] = []
    top_isoforms: list[dict[str, object]] = []
    candidates: list[dict[str, object]] = []
    detailed_table_count = 0
    for sample in sample_ids:
        source_root = Path(next(row for row in sample_rows if row["sample_id"] == sample).get("source_results_root") or results_root) / sample
        root = source_root
        if include_sample_reports:
            report = results_root / sample / f"07_sample_report/{sample}_report.html"
            copied = staging / f"sample_reports/{sample}.html"
            copy_file(report, copied)
            ensure_customer_text(copied)
        for delivered_name, source_pattern, target_pattern in DETAILED_TABLES:
            copy_result_table(root / source_pattern.format(sample=sample), staging / target_pattern.format(sample=sample), gene_map)
            detailed_table_count += 1
        row = qc_by_sample[sample]
        overview.append({"sample_id": sample, "group": sample_group[sample], "raw_reads": row["raw_reads"], "insert_rate": row["insert_rate"], "isoform_assigned_reads": row["isoform_assigned_reads"], "filtered_isoforms": row["filtered_isoforms"], "global_new_rna_fraction": row["grand3_global_binom_ntr"], "transcript_assignment_screen_rate": row["assignment_conflict_rate"], "candidate_isoform_usage_changes": row["switch_candidates"]})
        validation = json.loads((root / "06_isoform_dynamics/05_metadata/dynamics_validation.json").read_text(encoding="utf-8"))
        calls = validation["posterior_call_counts"]
        calls_row: dict[str, object] = {"sample_id": sample, "new_rna_reads": calls.get("New", 0), "old_rna_reads": calls.get("Old", 0), "ambiguous_reads": calls.get("Ambiguous", 0), "not_included_after_assignment_screen": calls.get("Assignment_conflict", 0)}
        if int(calls.get("No_T_sites", 0)):
            calls_row["reads_without_usable_t_sites"] = calls["No_T_sites"]
        calls_rows.append(calls_row)
        dynamics = read_tsv(root / f"06_isoform_dynamics/01_tables/{sample}.isoform_dynamics.tsv")
        for item in sorted(dynamics, key=lambda value: float(value["expected_new_reads"]), reverse=True)[:20]:
            add_gene_names([item], gene_map)
            top_isoforms.append({"sample_id": sample, **item})
        switches = read_tsv(root / f"06_isoform_dynamics/01_tables/{sample}.isoform_switch_candidates.tsv")
        for item in sorted(switches, key=lambda value: abs(float(value["delta_isoform_usage"])), reverse=True)[:20]:
            add_gene_names([item], gene_map)
            candidates.append({"sample_id": sample, **item})
    write_tsv(staging / f"{LAYOUT['01']}/sample_overview.tsv", list(overview[0]), overview)
    call_fields = ["sample_id", "new_rna_reads", "old_rna_reads", "ambiguous_reads", "not_included_after_assignment_screen"]
    if any("reads_without_usable_t_sites" in row for row in calls_rows):
        call_fields.append("reads_without_usable_t_sites")
    write_tsv(staging / f"{LAYOUT['01']}/read_classification_summary.tsv", call_fields, calls_rows)
    isoform_fields = ["sample_id", "gene_id", "gene_name", "isoform_id", "total_reads", "expected_new_reads", "expected_old_reads", "posterior_new_fraction", "technical_bootstrap_ci_low", "technical_bootstrap_ci_high", "total_isoform_usage", "new_isoform_usage", "delta_isoform_usage", "switch_status"]
    write_tsv(staging / f"{LAYOUT['05']}/top_isoforms_by_sample.tsv", isoform_fields, top_isoforms)
    write_tsv(staging / f"{LAYOUT['05']}/candidate_isoform_usage_by_sample.tsv", isoform_fields, candidates)
    project_qc_rows = [{field: row.get(field, "") for field in PROJECT_QC_FIELDS} for row in qc_rows]
    write_tsv(staging / f"{LAYOUT['01']}/project_qc_detailed.tsv", list(PROJECT_QC_FIELDS), project_qc_rows)
    group_source = results_root / "_project/08_project_qc/02_group_analysis"
    group_rows: list[dict[str, str]] = []
    group_validation: dict[str, object] = {}
    group_file_count = 0
    as_delivery_validation: dict[str, object] = {}
    if group_source.is_dir():
        validation_path = group_source / "metadata/group_analysis_validation.json"
        if not validation_path.is_file():
            raise RuntimeError(f"Group analysis validation is required: {validation_path}")
        group_validation = json.loads(validation_path.read_text(encoding="utf-8"))
        if group_validation.get("status") != "PASS":
            raise RuntimeError(f"Group analysis validation did not pass: {validation_path}")
        observed_groups = dict(Counter(sample_group.values()))
        if observed_groups != group_validation.get("biological_replicates_per_group", {}):
            raise RuntimeError("Sample sheet group counts do not match validated group analysis")
        group_rows = read_tsv(group_source / "tables/contrast_summary.tsv")
        for source in sorted(group_source.rglob("*")):
            if not source.is_file():
                continue
            relative = source.relative_to(group_source)
            target_relpath = group_delivery_target(relative)
            if target_relpath is None:
                continue
            target = staging / target_relpath
            if source.suffix.lower() == ".tsv":
                copy_result_table(source, target, gene_map)
            else:
                copy_file(source, target)
            group_file_count += 1
        as_delivery_validation = package_as_results(
            group_source, staging, gene_map, str(config.get("rscript_bin", "Rscript"))
        )
        figure_dir = staging / LAYOUT["00"] / "figures"
        if figure_dir.exists():
            # The report figure collection is the single browseable location for project figures.
            pass
        write_json(staging / LAYOUT["06"] / "group_analysis_validation.json", {
            key: group_validation[key] for key in (
                "status", "sample_count", "groups", "biological_replicates_per_group",
                "label_control_group", "transcript_comparison_scope", "usage_denominator",
                "batch_covariate", "contrast_details", "contrasts", "layers",
                "gene_features_union", "isoform_features_union", "gene_features_present_in_all_samples",
                "isoform_features_present_in_all_samples", "gene_annotation_coverage", "thresholds",
            ) if key in group_validation
        })
    qc_validation_path = results_root / "_project/08_project_qc/05_metadata/project_qc_validation.json"
    if qc_validation_path.is_file():
        raw_qc_validation = json.loads(qc_validation_path.read_text(encoding="utf-8"))
        safe_qc_validation = {key: raw_qc_validation[key] for key in ("status", "sample_count", "checks", "summary") if key in raw_qc_validation}
        write_json(staging / LAYOUT["06"] / "project_qc_validation.json", safe_qc_validation)
    groups = group_validation.get("groups", [])
    contrasts = group_validation.get("contrasts", [])
    replicate_counts = group_validation.get("biological_replicates_per_group", {})
    replicate_summary = "、".join(f"{group} n={replicate_counts.get(group, '未知')}" for group in groups)
    group_summary = (
        f"组间分析覆盖 {len(groups)} 个组（{replicate_summary}），完成 {len(contrasts)} 项预设比较；"
        "不估计 RNA 半衰期或绝对动力学参数。\n"
        if group_validation else ""
    )
    sample_report_line = "- `sample_reports/`：各样本的独立离线报告。\n" if include_sample_reports else ""
    report_path = f"{LAYOUT['00']}/ISO_SLAM_project_report.html"
    readme = f"# {config.get('delivery_name', '长读长 RNA 标记分析结果交付包')}\n\n本文件夹按编号分区组织，可离线打开。所有文件均为实体文件，不包含原始数据、BAM、read 级中间结果或软链接。\n\n- `{report_path}`：项目总报告，含逐样本审核指标。\n- `RESULTS_ANALYSIS_CN.md`：中文结果摘要。\n- `{LAYOUT['00']}/figures/`：项目图形。\n- `{LAYOUT['01']}/sample_qc_review.tsv`：逐样本审核指标数据。\n- `{LAYOUT['01']}/REVIEW_CRITERIA.md`：审核阈值和统计口径。\n- `{LAYOUT['01']}/`：样本质控与组内一致性汇总。\n- `{LAYOUT['02']}/`：逐样本转录本目录表。\n- `{LAYOUT['03']}/`：逐样本 RNA 标记与动态定量表。\n- `{LAYOUT['04']}/`：组间差异分析结果。\n- `{LAYOUT['05']}/`：逐样本转录本使用候选。\n- `{LAYOUT['06']}/`：脱敏后的分析设计及验证元数据。\n- `TABLE_GUIDE.md`：结果表说明。\n- `delivery_manifest.tsv`：文件清单。\n\n{group_summary}"
    readme = readme.replace("- `TABLE_GUIDE.md`", f"- `{LAYOUT['04']}/differential/new_rna_splicing/AS_结果摘要.md`：新生 RNA 可变剪切摘要和结果图。\n- `TABLE_GUIDE.md`")
    (staging / "README.md").write_text(readme, encoding="utf-8")
    splicing_validation_path = results_root / "_project/08_project_qc/02_group_analysis/tables/differential/new_rna_splicing/new_rna_splicing_validation.json"
    splicing_validation = json.loads(splicing_validation_path.read_text(encoding="utf-8")) if splicing_validation_path.is_file() else {}
    (staging / "RESULTS_ANALYSIS_CN.md").write_text(build_analysis_summary(overview, calls_rows, candidates, group_rows, group_validation, include_sample_reports, splicing_validation), encoding="utf-8")
    (staging / "TABLE_GUIDE.md").write_text(build_table_guide(), encoding="utf-8")
    manifest_rows: list[dict[str, object]] = []
    report_target = staging / report_path
    report_target.parent.mkdir(parents=True, exist_ok=True)
    os.replace(staging / "RNA_Labeling_Long_Read_Report.html", report_target)
    for path in sorted(staging.rglob("*")):
        if path.is_symlink():
            raise RuntimeError(f"Delivery must not contain symbolic links: {path}")
        if path.is_file():
            manifest_rows.append({"relative_path": str(path.relative_to(staging)), "bytes": path.stat().st_size, "sha256": digest(path), "kind": path.suffix.lstrip(".") or "text"})
    final_report_relpath = report_path
    report_manifest_rows = [row for row in manifest_rows if row["relative_path"] == final_report_relpath]
    manifest_rows = [row for row in manifest_rows if row["relative_path"] != final_report_relpath] + report_manifest_rows
    write_tsv(staging / "delivery_manifest.tsv", ["category", "source", "target", "status", "description"], [
        {"category": path.split("/", 1)[0], "source": "ISO-SLAM analysis outputs", "target": path, "status": "copied", "description": Path(path).name}
        for path in (row["relative_path"] for row in manifest_rows)
    ])
    write_json(staging / "DELIVERY_VALIDATION.json", {"status": "PASS", "package_name": config.get("delivery_name", "ISO-SLAM 结果交付包"), "sample_count": len(sample_ids), "group_count": len(groups), "contrasts": contrasts, "sample_reports_included": include_sample_reports, "qc_review_metrics_included": True, "qc_review_sample_count": review_summary["sample_count"], "qc_review_pass_sample_count": review_summary["pass_sample_count"], "qc_review_below_threshold_sample_count": review_summary["below_threshold_sample_count"], "detailed_table_count": detailed_table_count, "group_analysis_file_count": group_file_count, "alternative_splicing": as_delivery_validation, "file_count": len(manifest_rows) + 2, "symlink_count": 0, "internal_customer_terms": False})
    try:
        ensure_delivery_desensitized(staging)
    except Exception:
        shutil.rmtree(staging)
        raise
    destination.parent.mkdir(parents=True, exist_ok=True)
    archived = ""
    if destination.exists():
        archive = tmp_root / "delivery_archive" / f"{destination.name}_{token}"
        archive.parent.mkdir(parents=True, exist_ok=True)
        shutil.move(str(destination), str(archive))
        archived = str(archive)
    os.replace(staging, destination)
    print(json.dumps({"status": "PASS", "delivery_root": str(destination), "sample_count": len(sample_ids), "group_count": len(groups), "contrast_count": len(contrasts), "sample_reports_included": include_sample_reports, "archived_previous_delivery": archived}, ensure_ascii=False))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
