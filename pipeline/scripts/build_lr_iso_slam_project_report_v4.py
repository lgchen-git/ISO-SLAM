#!/usr/bin/env python3
"""Build the customer-facing offline project report with group analysis."""

from __future__ import annotations

import argparse
import base64
import csv
import html
import json
from pathlib import Path

from iso_slam_v4_common import ensure_stage, read_config, write_json, write_tsv

GROUP_PALETTE = ["#3C5488", "#E64B35", "#00A087", "#4DBBD5", "#F39B7F", "#8491B4", "#91D1C2", "#DC0000"]
GROUP_COLORS = {"A": "#E64B35", "AB": "#4DBBD5", "B": "#00A087", "NC": "#3C5488"}
LAYER_LABELS = {"gene_abundance": "基因丰度", "gene_labeling": "基因标记比例", "isoform_abundance": "转录本丰度", "total_isoform_usage": "总 RNA 转录本使用", "new_isoform_usage": "新生 RNA 转录本使用"}


def esc(value: object) -> str:
    return html.escape(str(value))


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open(encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def pct(value: object, digits: int = 2) -> str:
    return f"{float(value) * 100:.{digits}f}%"


def data_uri(path: Path) -> str:
    return "data:image/png;base64," + base64.b64encode(path.read_bytes()).decode("ascii")


def finite_float(value: object, default: float) -> float:
    try:
        number = float(value)
        return number if number == number and abs(number) != float("inf") else default
    except (TypeError, ValueError):
        return default


def figure(path: Path, title: str, caption: str) -> str:
    uri = data_uri(path)
    return f'<figure class="figure-card"><img src="{uri}" alt="{esc(title)}"><figcaption><strong>{esc(title)}</strong><span>{esc(caption)}</span></figcaption></figure>'


def static_table(rows: list[dict[str, str]], columns: list[tuple[str, str]]) -> str:
    head = "".join(f"<th>{esc(label)}</th>" for _, label in columns)
    body = "".join("<tr>" + "".join(f"<td>{esc(row.get(key, ''))}</td>" for key, _ in columns) + "</tr>" for row in rows)
    return f'<div class="table-wrap"><table><thead><tr>{head}</tr></thead><tbody>{body}</tbody></table></div>'


REVIEW_FIELDS = [
    "sample_id", "group", "raw_reads", "raw_reads_threshold", "raw_reads_status",
    "fastq_q20_base_fraction", "fastq_q20_threshold", "fastq_q20_status", "fastq_qc_reads",
    "primary_mapped_reads", "primary_mapq20_reads", "primary_mapq20_fraction",
    "primary_mapq20_threshold", "primary_mapq20_status", "genes_detected",
    "genes_detected_threshold", "genes_detected_status", "insert_reads",
    "insert_retention_rate", "insert_mapping_rate", "overall_review_status",
]


def build_qc_review_rows(results_root: Path, qc_rows: list[dict[str, str]], config: dict[str, str]) -> list[dict[str, object]]:
    thresholds = {
        "raw_reads": int(config.get("qc_review_min_raw_reads", "4000000")),
        "fastq_q20": float(config.get("qc_review_min_fastq_q20_fraction", "0.8")),
        "primary_mapq20": float(config.get("qc_review_min_primary_mapq20_fraction", "0.9")),
        "genes_detected": int(config.get("qc_review_min_detected_genes", "10000")),
    }
    if thresholds["raw_reads"] <= 0 or thresholds["genes_detected"] <= 0 or not 0 <= thresholds["fastq_q20"] <= 1 or not 0 <= thresholds["primary_mapq20"] <= 1:
        raise ValueError(f"Invalid QC review thresholds: {thresholds}")

    rows: list[dict[str, object]] = []
    for qc in qc_rows:
        sample = qc["sample_id"]
        sample_root = results_root / sample
        raw_qc_path = sample_root / f"01_raw_qc/04_tool_output/longread_qc/fastq/{sample}.raw_fastq_qc.json"
        flow_path = sample_root / "07_sample_report/04_tool_output/read_flow/tables/read_flow_nodes.tsv"
        caller_path = sample_root / f"05_slam_model/04_tool_output/caller/{sample}.bulk_slam_qc.json"
        for path in (raw_qc_path, flow_path, caller_path):
            if not path.is_file() or path.stat().st_size == 0:
                raise RuntimeError(f"QC review source is missing for {sample}: {path}")

        raw_qc = json.loads(raw_qc_path.read_text(encoding="utf-8"))
        flow_rows = read_tsv(flow_path)
        nodes = {row["node_id"]: int(row["value"]) for row in flow_rows if row.get("view_id") == "sample_reads"}
        caller_qc = json.loads(caller_path.read_text(encoding="utf-8"))
        needed_nodes = {"raw", "insert", "mapped", "mapq20"}
        if not needed_nodes.issubset(nodes):
            raise RuntimeError(f"Read-flow QC nodes are incomplete for {sample}: missing {sorted(needed_nodes - nodes.keys())}")

        q20_reads = int(raw_qc.get("read_count", 0))
        q20_fraction = float(raw_qc.get("q20_base_fraction", float("nan")))
        mapped = nodes["mapped"]
        if q20_reads <= 0 or not 0 <= q20_fraction <= 1 or mapped <= 0:
            raise RuntimeError(f"Invalid FASTQ or mapping QC values for {sample}")
        genes_detected = int(caller_qc.get("genes_detected", -1))
        if genes_detected < 0:
            raise RuntimeError(f"Caller QC does not contain genes_detected for {sample}")

        statuses = {
            "raw_reads": int(qc["raw_reads"]) >= thresholds["raw_reads"],
            "fastq_q20": q20_fraction >= thresholds["fastq_q20"],
            "primary_mapq20": nodes["mapq20"] / mapped >= thresholds["primary_mapq20"],
            "genes_detected": genes_detected >= thresholds["genes_detected"],
        }
        insert = nodes["insert"]
        values: dict[str, object] = {
            "sample_id": sample,
            "group": qc.get("group_label_metadata_only", ""),
            "raw_reads": int(qc["raw_reads"]),
            "raw_reads_threshold": thresholds["raw_reads"],
            "raw_reads_status": "PASS" if statuses["raw_reads"] else "BELOW_THRESHOLD",
            "fastq_q20_base_fraction": q20_fraction,
            "fastq_q20_threshold": thresholds["fastq_q20"],
            "fastq_q20_status": "PASS" if statuses["fastq_q20"] else "BELOW_THRESHOLD",
            "fastq_qc_reads": q20_reads,
            "primary_mapped_reads": mapped,
            "primary_mapq20_reads": nodes["mapq20"],
            "primary_mapq20_fraction": nodes["mapq20"] / mapped,
            "primary_mapq20_threshold": thresholds["primary_mapq20"],
            "primary_mapq20_status": "PASS" if statuses["primary_mapq20"] else "BELOW_THRESHOLD",
            "genes_detected": genes_detected,
            "genes_detected_threshold": thresholds["genes_detected"],
            "genes_detected_status": "PASS" if statuses["genes_detected"] else "BELOW_THRESHOLD",
            "insert_reads": insert,
            "insert_retention_rate": insert / nodes["raw"] if nodes["raw"] else 0.0,
            "insert_mapping_rate": mapped / insert if insert else 0.0,
            "overall_review_status": "PASS" if all(statuses.values()) else "BELOW_THRESHOLD",
        }
        if int(values["raw_reads"]) != nodes["raw"]:
            raise RuntimeError(f"Raw read count differs between project QC and read-flow for {sample}")
        rows.append(values)
    return rows


def review_html_rows(rows: list[dict[str, object]]) -> list[dict[str, str]]:
    result = []
    for row in rows:
        result.append({
            "sample_id": str(row["sample_id"]),
            "group": str(row["group"]),
            "raw_reads_display": f'{int(row["raw_reads"]):,} / {int(row["raw_reads_threshold"]):,} ({row["raw_reads_status"]})',
            "q20_display": f'{float(row["fastq_q20_base_fraction"]):.1%} / {float(row["fastq_q20_threshold"]):.0%} ({row["fastq_q20_status"]}; n={int(row["fastq_qc_reads"]):,})',
            "mapq20_display": f'{float(row["primary_mapq20_fraction"]):.1%} / {float(row["primary_mapq20_threshold"]):.0%} ({row["primary_mapq20_status"]})',
            "genes_display": f'{int(row["genes_detected"]):,} / {int(row["genes_detected_threshold"]):,} ({row["genes_detected_status"]})',
            "insert_retention_display": f'{float(row["insert_retention_rate"]):.2%}',
            "insert_mapping_display": f'{float(row["insert_mapping_rate"]):.2%}',
            "overall_review_status": str(row["overall_review_status"]),
        })
    return result


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--results-root", required=True, type=Path)
    parser.add_argument("--config", type=Path)
    args = parser.parse_args()
    qc_path = args.results_root / "_project/08_project_qc/01_tables/sample_qc_summary.tsv"
    group_root = args.results_root / "_project/08_project_qc/02_group_analysis"
    validation_path = group_root / "metadata/group_analysis_validation.json"
    if not validation_path.is_file():
        raise RuntimeError(f"Completed group analysis is required: {validation_path}")
    validation = json.loads(validation_path.read_text(encoding="utf-8"))
    if validation.get("status") != "PASS":
        raise RuntimeError("Group analysis validation did not pass")
    global GROUP_COLORS
    GROUP_COLORS = dict(zip(validation["groups"], GROUP_PALETTE))
    replicate_counts = validation.get("biological_replicates_per_group", {})
    replicate_summary = "、".join(f"{group} n={replicate_counts.get(group, '未知')}" for group in validation["groups"])
    contrast_directions = "；".join(value.replace("_vs_", " 减 ") for value in validation["contrasts"])
    label_control_group = validation.get("label_control_group", validation.get("control_group", ""))
    qc_rows = read_tsv(qc_path)
    run_config = read_config(args.config) if args.config else {}
    review_rows = build_qc_review_rows(args.results_root, qc_rows, run_config)
    review_path = args.results_root / "_project/09_project_report/01_tables/sample_qc_review.tsv"
    write_tsv(review_path, review_rows, REVIEW_FIELDS)
    review_display = review_html_rows(review_rows)
    signal_rows = read_tsv(group_root / "tables/sample_labeling_signal.tsv")
    within_rows = read_tsv(group_root / "tables/within_group_signal_summary.tsv")
    contrast_rows = read_tsv(group_root / "tables/contrast_summary.tsv")
    splicing_root = group_root / "tables/differential/new_rna_splicing"
    splicing_validation = json.loads((splicing_root / "new_rna_splicing_validation.json").read_text(encoding="utf-8"))
    splicing_rows = read_tsv(splicing_root / "differential_new_rna_splicing.tsv")
    splicing_rows.sort(key=lambda row: (finite_float(row.get("FDR"), 1), -abs(finite_float(row.get("delta_PSI"), 0))))
    splicing_display = [{**row, "delta_PSI": f'{finite_float(row.get("delta_PSI"), 0):.3f}',
                         "mean_PSI_group_1": f'{finite_float(row.get("mean_PSI_group_1"), 0):.3f}',
                         "mean_PSI_group_2": f'{finite_float(row.get("mean_PSI_group_2"), 0):.3f}',
                         "FDR": f'{finite_float(row.get("FDR"), 1):.3g}',
                         "significant": "是" if row.get("significant", "").lower() == "true" else "否"}
                        for row in splicing_rows[:20]]
    splicing_significant = sum(1 for row in splicing_rows if row.get("significant", "").lower() == "true")
    for row in qc_rows:
        row["group"] = row.get("group_label_metadata_only", "")
        row["insert_rate_display"] = pct(row["insert_rate"])
        row["ntr_display"] = pct(row["grand3_global_binom_ntr"])
        row["screen_display"] = pct(row["assignment_conflict_rate"])
    for row in signal_rows:
        row["tc_display"] = pct(row["tc_rate"], 3)
        row["background_display"] = f"{float(row['diagnostic_background_rate']) * 1000:.3f}‰"
    for row in within_rows:
        row["tc_mean_display"] = pct(row["tc_rate_mean"], 3)
        row["tc_cv_display"] = f"{float(row['tc_rate_cv']):.3f}"
        row["background_mean_display"] = f"{float(row['diagnostic_background_mean']) * 1000:.3f}‰"
        row["background_cv_display"] = f"{float(row['diagnostic_background_cv']):.3f}"
    empty_uri = "data:image/svg+xml;base64," + base64.b64encode('<svg xmlns="http://www.w3.org/2000/svg" width="760" height="120"><text x="20" y="60" font-size="20">没有满足展示条件的特征。</text></svg>'.encode()).decode()
    result_payload: dict[str, dict[str, object]] = {}
    for contrast in validation["contrasts"]:
        result_payload[contrast] = {}
        for layer in validation["layers"]:
            rows = read_tsv(group_root / f"tables/differential/{contrast}/{layer}.tsv")
            rows.sort(key=lambda row: (float(row.get("FDR", 1) or 1), -abs(float(row.get("effect", 0) or 0))))
            compact = [{"feature_id": row.get("feature_id", ""), "gene_id": row.get("gene_id", ""), "gene_name": row.get("gene_name", "") or row.get("gene_id", ""), "effect": finite_float(row.get("effect"), 0), "effect_type": row.get("effect_type", ""), "fdr": finite_float(row.get("FDR"), 1), "mean_group_1": finite_float(row.get("mean_group_1"), 0), "mean_group_2": finite_float(row.get("mean_group_2"), 0), "significant": str(row.get("significant", "")).upper() == "TRUE"} for row in rows]
            heat_rows = read_tsv(group_root / f"figures/differential/{contrast}/{layer}_heatmap_data.tsv")
            heatmap = [{"feature_id": row.get("feature_id", ""), "gene_id": row.get("gene_id", ""), "gene_name": row.get("gene_name", "") or row.get("gene_id", ""), "sample_id": row.get("sample_id", ""), "group": row.get("group", ""), "z_score": finite_float(row.get("z_score"), 0)} for row in heat_rows]
            volcano_path = group_root / f"figures/differential/{contrast}/{layer}_volcano.png"
            heatmap_path = group_root / f"figures/differential/{contrast}/{layer}_heatmap.png"
            result_payload[contrast][layer] = {"label": LAYER_LABELS[layer], "rows": compact, "heatmap_rows": heatmap,
                "volcano_uri": data_uri(volcano_path) if rows and volcano_path.is_file() else empty_uri,
                "heatmap_uri": data_uri(heatmap_path) if heat_rows and heatmap_path.is_file() else empty_uri}
    total_significant = sum(int(row["significant_features"]) for row in contrast_rows)
    sample_count = len(qc_rows)
    group_count = len(validation["groups"])
    contrast_count = len(validation["contrasts"])
    contrast_options = "".join(f'<option value="{esc(value)}">{esc(value.replace("_vs_", " vs "))}</option>' for value in validation["contrasts"])
    layer_options = "".join(f'<option value="{esc(value)}">{esc(LAYER_LABELS[value])}</option>' for value in validation["layers"])
    group_legend = "".join(f'<span><i style="background:{color}"></i>{esc(group)}</span>' for group, color in GROUP_COLORS.items())
    plotly_b64 = base64.b64encode((Path(__file__).resolve().parents[1] / "assets/plotly-3.5.1.min.js").read_bytes()).decode("ascii")
    logo_uri = data_uri(Path(__file__).resolve().parents[1] / "assets/epi-logo.png")
    stage = args.results_root / "_project/09_project_report"
    ensure_stage(stage)
    content = f'''<!doctype html><html lang="zh-CN"><head><meta charset="utf-8"><meta name="viewport" content="width=device-width,initial-scale=1"><title>长读长 RNA 标记分析结果报告</title><style>
*{{box-sizing:border-box}}html{{scroll-behavior:smooth}}body{{margin:0;background:#f3f6fa;color:#182230;font:14px/1.65 Arial,"Noto Sans SC",sans-serif}}.sidebar{{position:fixed;inset:0 auto 0 0;width:250px;background:#17324d;color:#fff;padding:24px 18px;overflow:auto}}.sidebar strong{{font-size:17px}}.sidebar a{{display:block;color:#e9f0f8;text-decoration:none;padding:7px 0}}main{{margin-left:250px;padding:28px;max-width:1540px}}section{{background:#fff;border:1px solid #dce5ef;border-radius:10px;padding:24px;margin-bottom:20px}}h1{{font-size:30px;margin:0 0 8px}}h2{{color:#284d73;border-bottom:2px solid #dce6f0;padding-bottom:8px}}h3{{color:#385b7b}}.subtitle,.notice{{background:#f4f8fb;color:#506477;padding:12px 14px;border-left:4px solid #4DBBD5;border-radius:5px}}.warning{{border-left-color:#E64B35}}.kpis{{display:grid;grid-template-columns:repeat(4,1fr);gap:12px;margin:18px 0}}.kpi{{background:#f7faff;border-top:4px solid #3C5488;border-radius:7px;padding:13px}}.kpi span,.kpi small{{display:block;color:#667085}}.kpi strong{{display:block;font-size:22px;color:#17324d;margin:4px 0}}.grid{{display:grid;grid-template-columns:1fr 1fr;gap:18px}}.concepts{{display:grid;grid-template-columns:repeat(3,1fr);gap:12px}}.concept{{border:1px solid #dde6ef;border-radius:8px;padding:13px}}.concept strong{{display:block;color:#284d73;margin-bottom:4px}}.figure-card{{margin:0;border:1px solid #dde6ef;border-radius:9px;padding:10px;background:#fff}}.figure-card img{{display:block;width:100%;height:auto}}.plot{{width:100%;height:520px}}.static-fallback{{display:none!important}}figcaption{{padding:9px 4px 2px}}figcaption strong,figcaption span{{display:block}}figcaption span{{color:#607184}}.legend{{display:flex;gap:18px;flex-wrap:wrap;margin:10px 0}}.legend i{{display:inline-block;width:12px;height:12px;margin-right:6px;border-radius:2px}}.table-wrap{{overflow:auto}}table{{width:100%;border-collapse:collapse;font-size:12.5px}}th,td{{padding:7px 8px;border-bottom:1px solid #e5eaf0;text-align:left;white-space:nowrap}}th{{background:#eef4f9}}.selectors{{display:flex;gap:12px;flex-wrap:wrap;margin-bottom:14px}}select,input,button{{font:inherit}}select,input{{padding:7px 9px;border:1px solid #bdc8d5;border-radius:5px;background:#fff}}.pager{{display:flex;align-items:center;gap:10px;margin-top:10px}}.pager button{{padding:6px 10px;border:1px solid #bdc8d5;background:#fff;border-radius:5px}}.references li{{margin-bottom:8px}}.footer{{margin-left:250px;background:#fff;color:#526477;padding:22px 30px;text-align:center;border-top:1px solid #dce5ef}}.footer img{{display:block;width:min(320px,80vw);height:auto;margin:0 auto 10px}}.modal{{display:none;position:fixed;inset:0;background:#000b;z-index:20;padding:3vh 3vw}}.modal.open{{display:flex;align-items:center;justify-content:center}}.modal-box{{background:#fff;max-width:94vw;max-height:94vh;padding:12px;border-radius:8px;overflow:auto}}.modal img{{max-width:90vw;max-height:84vh}}.modal-close{{float:right;border:0;background:#17324d;color:#fff;padding:6px 10px;border-radius:4px}}@media(max-width:1000px){{.kpis{{grid-template-columns:repeat(2,1fr)}}.grid,.concepts{{grid-template-columns:1fr}}}}@media(max-width:760px){{.sidebar{{position:static;width:auto}}main,.footer{{margin-left:0}}}}@media print{{.sidebar,.selectors,.pager,.modal,.plot{{display:none!important}}.static-fallback{{display:block!important}}main,.footer{{margin-left:0}}section,.figure-card{{break-inside:avoid}}#result-body tr{{display:table-row!important}}}}
</style><script>eval(atob("{plotly_b64}"))</script></head><body><nav class="sidebar"><strong>项目总报告</strong><a href="#overview">1. 项目概览</a><a href="#concepts">2. 概念说明</a><a href="#quality">3. 数据质量与标记信号</a><a href="#within">4. 组内一致性</a><a href="#between">5. 组间整体结构</a><a href="#differential">6. 差异分析</a><a href="#new-rna-splicing">7. 新生 RNA 可变剪切</a><a href="#interpretation">8. 结果解读</a><a href="#methods">9. 方法与限制</a><a href="#references">10. 参考文献</a></nav><main>
<section id="overview"><h1>长读长 RNA 标记分析结果报告</h1><p class="subtitle">{esc("、".join(validation["groups"]))} 共 {group_count} 个实验组；完成 {contrast_count} 项指定比较。</p><div class="legend">{group_legend}</div><div class="kpis"><div class="kpi"><span>样本数</span><strong>{sample_count}</strong><small>{group_count} 个实验组</small></div><div class="kpi"><span>组间比较</span><strong>{contrast_count}</strong><small>配置指定比较</small></div><div class="kpi"><span>差异结果层级</span><strong>{len(validation["layers"])}</strong><small>基因、标记与转录本</small></div><div class="kpi"><span>显著结果记录</span><strong>{total_significant:,}</strong><small>各比较与层级合计</small></div></div><p class="notice warning">显著记录数是不同比较和分析层级的结果条目总和，不能理解为互不重复的基因数。</p></section>
<section id="concepts"><h2>2. 概念说明</h2><div class="concepts"><div class="concept"><strong>T&gt;C 转化率</strong>在通过筛选并覆盖参考 T 位点的碱基中，观察到 T&gt;C 的比例；分母是纳入统计的参考 T 碱基。</div><div class="concept"><strong>样本内诊断背景</strong>背景/低深度 T 位点中 T&gt;A 与 T&gt;G 比率的较大值，用于模型诊断，不等同于实测未标记对照。</div><div class="concept"><strong>一致性筛查</strong>检查 read 的转录本分配是否与后续转录本汇总一致；未通过者不进入转录本层面的新生/已有 RNA 和使用变化结果。</div><div class="concept"><strong>NTR</strong>模型估计的新生 RNA 占总 RNA 的比例；它是校正后的比例估计，不是绝对合成速率。</div><div class="concept"><strong>log2FC 与 Δ比例</strong>计数型丰度使用 log2 倍数变化；标记和使用比例使用两组原始比例均值之差。</div><div class="concept"><strong>p 值与 FDR</strong>p 值描述零假设下的证据；FDR 为同一比较和结果层级内的 BH 多重检验校正值。</div></div></section>
<section id="quality"><h2>3. 数据质量与标记信号</h2><h3>逐样本审核指标</h3>{static_table(review_display, [("sample_id","样本"),("group","组别"),("raw_reads_display","原始 reads / 阈值"),("q20_display","FASTQ 碱基 Q20 / 阈值"),("mapq20_display","MAPQ≥20 / Primary mapped"),("genes_display","检测基因数 / 阈值"),("insert_retention_display","Insert 保留率"),("insert_mapping_display","Insert 比对率"),("overall_review_status","审核状态")])}<p class="notice">FASTQ 碱基 Q20 统计每样本原始 FASTQ 前 200,000 条 reads；MAPQ≥20 比例以 Primary mapped reads 为分母。Insert 保留率和比对率为报告指标，不设本表阈值。阈值来自本次运行配置；PASS 表示达到表列审核线，不替代流程验证或生物学判断。</p><h3>原有项目 QC 汇总</h3>{static_table(qc_rows, [("sample_id","样本"),("group","组别"),("raw_reads","原始 reads"),("insert_rate_display","插入序列保留"),("isoform_assigned_reads","转录本分配 reads"),("ntr_display","全局新生 RNA"),("screen_display","一致性筛查")])}<h3>T&gt;C 与诊断背景</h3>{static_table(signal_rows, [("sample_id","样本"),("group","组别"),("tc_display","T>C 转化率"),("background_display","诊断背景")])}<p class="notice">不同指标来自不同筛选链和分母；T&gt;C 观测率、诊断背景和模型估计 NTR 不应互相替代。</p></section>
<section id="within"><h2>4. 组内一致性</h2><div class="grid">{figure(group_root / "figures/within_group_labeling_signal.png", "组内 RNA 标记信号", "散点为独立生物学重复，横线为组均值；两个面板采用各自量程。")}{figure(group_root / "figures/gene_abundance_correlation.png", "样本相关性", "基于过滤后基因 logCPM 的 Spearman 相关性。")}</div>{static_table(within_rows, [("group","组别"),("tc_mean_display","平均 T>C"),("tc_cv_display","T>C CV"),("background_mean_display","平均背景"),("background_cv_display","背景 CV")])}</section>
<section id="between"><h2>5. 组间整体结构</h2><div class="grid">{figure(group_root / "figures/gene_abundance_pca.png", "基因丰度 PCA", "基于 TMM 标准化后过滤基因的 logCPM；点间距离反映全局丰度结构。")}{figure(group_root / "figures/differential_feature_counts.png", "显著差异结果数量", "按比较和分析层级汇总；丰度和比例结果使用各自预设效应阈值。")}</div></section>
<section id="differential"><h2>6. 差异分析</h2><div class="selectors"><label>比较 <select id="contrast-select">{contrast_options}</select></label><label>结果层级 <select id="layer-select">{layer_options}</select></label><label>表格筛选 <input id="result-filter" type="search" placeholder="基因名、基因或转录本 ID"></label></div><div class="grid"><article class="figure-card"><div id="volcano-plot" class="plot" aria-label="交互式差异火山图"></div><img id="volcano-static" class="static-fallback" alt="差异火山图静态版本"><figcaption><strong id="volcano-title">差异火山图</strong><span>横轴为效应量，纵轴为 -log10(FDR)；鼠标悬停可查看基因信息。</span></figcaption></article><article class="figure-card"><div id="heatmap-plot" class="plot" aria-label="交互式 Top 特征热图"></div><img id="heatmap-static" class="static-fallback" alt="Top 特征热图静态版本"><figcaption><strong id="heatmap-title">Top 特征热图</strong><span>每行按样本进行 z-score 标准化；鼠标悬停可查看基因、样本和组别。</span></figcaption></article></div><div class="table-wrap"><table><thead><tr><th>特征</th><th>基因名</th><th>基因 ID</th><th>效应量</th><th>效应类型</th><th>FDR</th><th>显著</th></tr></thead><tbody id="result-body"></tbody></table></div><div class="pager"><button id="prev-page" type="button">上一页</button><span id="page-status"></span><button id="next-page" type="button">下一页</button></div><p class="notice">交互图展示全部通过分析过滤的特征；表格分页浏览相同结果，完整字段见交付包中的 TSV 文件。</p></section>
<section id="new-rna-splicing"><h2>7. 新生 RNA 可变剪切</h2><p class="notice">参考注释事件 {int(splicing_validation.get('reference_event_count', 0)):,} 个；各比较中满足新生 RNA 支持要求并进入检验的事件 {esc(json.dumps(splicing_validation.get('tested_event_count_by_contrast', {}), ensure_ascii=False))}；达到 FDR&lt;0.05 且 |ΔPSI|≥0.10 的事件 {splicing_significant} 个。</p>{static_table(splicing_display, [("contrast","比较"),("event_id","事件"),("gene_name","基因名"),("gene_id","基因 ID"),("mean_PSI_group_1","实验组 PSI"),("mean_PSI_group_2","对照组 PSI"),("delta_PSI","ΔPSI"),("FDR","FDR"),("significant","显著")])}<p class="notice warning">PSI 由转录本层面的期望新生 reads 等价数计算；事件限定于参考注释的 SE、A3、A5、MX、RI。新发现转录本结构另列清单，但现有结果无法可靠估计其逐读段新生 RNA 后验，因此不将其纳入新生 RNA 事件差异检验。每组 n=2，统计结果仅作探索性筛查。完整样本 PSI、差异结果及新结构清单见交付包 TSV。</p></section>
<section id="interpretation"><h2>8. 结果解读</h2><ul><li>组内分析用于判断独立生物学重复的技术和全局信号一致性；PCA 分离本身不等同于统计显著。</li><li>基因和转录本丰度差异使用计数模型；标记比例与转录本使用差异使用 read 支持加权的比例模型。</li><li>差异结果同时要求 FDR 和预设效应量阈值，完整表仍保留所有通过分析过滤的特征。</li><li>{esc(label_control_group)} 是已标记的生物学对照，不作为未标记背景对照。</li></ul></section>
<section id="methods"><h2>9. 方法与限制</h2><details open><summary>统计方法</summary><ul><li>{esc("、".join(validation["groups"]))} 多组独立设计（{esc(replicate_summary)}）；未提供批次协变量。比较效应方向依次为 {esc(contrast_directions)}。</li><li>基因和转录本丰度：edgeR TMM、filterByExpr、负二项准似然模型；显著阈值 FDR&lt;0.05 且 |log2FC|≥1。</li><li>基因标记比例：支持 reads≥20 且比较两组各至少 2 个样本满足条件后，经验 logit 转换、read 支持加权 limma；阈值 FDR&lt;0.05 且 |ΔNTR|≥0.05。</li><li>转录本使用：同样使用加权 limma；阈值 FDR&lt;0.05 且 |Δusage|≥0.10。</li><li>转录本仅比较能唯一对应参考注释的完整剪接匹配结构；使用比例的分母为同基因内可比较参考转录本的 reads。独立发现的新转录本保留单样本结果。</li><li>新生 RNA 可变剪切：以每样本期望新生 reads 等价数计算参考注释事件 PSI；每样本支持≥20后进行经验 logit、read 支持加权 limma 和 eBayes，并按比较进行 BH 校正，阈值 FDR&lt;0.05 且 |ΔPSI|≥0.10。</li></ul></details><details><summary>分析限制</summary><ul><li>未提供批次变量，因此模型不校正未知批次；若后续确认批次，需要重新拟合。</li><li>单次标记且无时间序列，不估计绝对合成速率、降解速率或半衰期。</li><li>SQANTI3 未运行 ORF 预测；coding 状态未评估，不能据此把转录本解释为编码或非编码。</li></ul></details></section>
<section id="references"><h2>10. 参考文献</h2><ol class="references"><li>Jürges C, et al. GRAND-SLAM enables genome-wide profiling of RNA synthesis and degradation kinetics for single-cell and bulk RNA sequencing. <em>Bioinformatics</em>. 2018;34:i218–i226. DOI: 10.1093/bioinformatics/bty256.</li><li>Robinson MD, McCarthy DJ, Smyth GK. edgeR: a Bioconductor package for differential expression analysis of digital gene expression data. <em>Bioinformatics</em>. 2010;26:139–140. DOI: 10.1093/bioinformatics/btp616.</li><li>Ritchie ME, et al. limma powers differential expression analyses for RNA-sequencing and microarray studies. <em>Nucleic Acids Research</em>. 2015;43:e47. DOI: 10.1093/nar/gkv007.</li><li>Benjamini Y, Hochberg Y. Controlling the false discovery rate: a practical and powerful approach to multiple testing. <em>Journal of the Royal Statistical Society Series B</em>. 1995;57:289–300. DOI: 10.1111/j.2517-6161.1995.tb02031.x.</li><li>Prjibelski AD, et al. IsoQuant: a tool for accurate novel isoform discovery with long reads. <em>Nature Biotechnology</em>. 2023;41:915–918. DOI: 10.1038/s41587-022-01565-y.</li><li>De Coster W, et al. NanoPack: visualizing and processing long-read sequencing data. <em>Bioinformatics</em>. 2018;34:2666–2669. DOI: 10.1093/bioinformatics/bty149.</li></ol></section></main><footer class="footer"><img src="{logo_uri}" alt="EPIBIOTEK 麦观生物"><span>长读长 RNA 标记分析 · 离线项目报告</span></footer><script>
const resultData={json.dumps(result_payload, ensure_ascii=False, separators=(',', ':'))};
const pageSize=10;let page=1,printing=false;
const contrast=document.getElementById('contrast-select'),layer=document.getElementById('layer-select'),filter=document.getElementById('result-filter'),body=document.getElementById('result-body'),status=document.getElementById('page-status');
function current(){{return resultData[contrast.value][layer.value]}}
function volcanoPlot(d){{
  if(!d.rows.length){{document.getElementById("volcano-plot").innerHTML="<p>没有满足检验条件的特征。</p>";return}}
  const traces=[false,true].map(sig=>{{const rows=d.rows.filter(r=>r.significant===sig);return {{type:'scattergl',mode:'markers',name:sig?'Significant':'Not significant',x:rows.map(r=>r.effect),y:rows.map(r=>-Math.log10(Math.max(r.fdr,1e-300))),customdata:rows.map(r=>[r.gene_name,r.gene_id,r.feature_id,r.effect_type,r.fdr,r.mean_group_1,r.mean_group_2,sig?'Yes':'No']),marker:{{size:sig?7:5,color:sig?'#E64B35':'#9AA4B2',opacity:sig?.82:.55}},hovertemplate:'Gene name: %{{customdata[0]}}<br>Gene ID: %{{customdata[1]}}<br>Feature ID: %{{customdata[2]}}<br>%{{customdata[3]}}: %{{x:.5f}}<br>FDR: %{{customdata[4]:.3e}}<br>Group 1 mean: %{{customdata[5]:.5f}}<br>Group 2 mean: %{{customdata[6]:.5f}}<br>Significant: %{{customdata[7]}}<extra></extra>'}}}});
  Plotly.react('volcano-plot',traces,{{margin:{{l:65,r:20,t:30,b:60}},xaxis:{{title:d.rows[0]?.effect_type||'Effect'}},yaxis:{{title:'-log10(FDR)'}},legend:{{orientation:'h'}},hovermode:'closest',paper_bgcolor:'#fff',plot_bgcolor:'#fff'}},{{responsive:true,displaylogo:false,modeBarButtonsToRemove:['lasso2d','select2d']}});
}}
function heatmapPlot(d){{
  if(!d.heatmap_rows.length){{document.getElementById("heatmap-plot").innerHTML="<p>没有可展示的特征热图。</p>";return}}
  const rows=d.heatmap_rows,features=[...new Set(rows.map(r=>r.feature_id))],samples=[...new Set(rows.map(r=>r.sample_id))];
  const lookup=new Map(rows.map(r=>[r.feature_id+'\t'+r.sample_id,r]));
  const z=features.map(f=>samples.map(s=>lookup.get(f+'\t'+s)?.z_score??0));
  const custom=features.map(f=>samples.map(s=>{{const r=lookup.get(f+'\t'+s)||{{}};return [r.gene_name||'',r.gene_id||'',f,s,r.group||'']}}));
  const labels=features.map(f=>{{const r=rows.find(x=>x.feature_id===f);return r?.gene_name&&r.gene_name!==r.gene_id?r.gene_name:f}});
  Plotly.react('heatmap-plot',[{{type:'heatmap',z,x:samples,y:labels,customdata:custom,colorscale:[[0,'#3C5488'],[.5,'#F7F7F7'],[1,'#E64B35']],zmid:0,colorbar:{{title:'Row z-score'}},hovertemplate:'Gene name: %{{customdata[0]}}<br>Gene ID: %{{customdata[1]}}<br>Feature ID: %{{customdata[2]}}<br>Sample: %{{customdata[3]}}<br>Group: %{{customdata[4]}}<br>Row z-score: %{{z:.3f}}<extra></extra>'}}],{{margin:{{l:120,r:25,t:30,b:85}},xaxis:{{tickangle:-45}},yaxis:{{automargin:true}},paper_bgcolor:'#fff'}},{{responsive:true,displaylogo:false}});
}}
function render(){{
  const d=current();document.getElementById('volcano-static').src=d.volcano_uri;document.getElementById('heatmap-static').src=d.heatmap_uri;const q=filter.value.trim().toLowerCase();let rows=d.rows.filter(r=>!q||r.feature_id.toLowerCase().includes(q)||r.gene_id.toLowerCase().includes(q)||r.gene_name.toLowerCase().includes(q));const pages=Math.max(1,Math.ceil(rows.length/pageSize));page=Math.min(page,pages);const shown=printing?d.rows:rows.slice((page-1)*pageSize,page*pageSize);
  body.innerHTML=shown.map(r=>`<tr><td>${{r.feature_id}}</td><td>${{r.gene_name}}</td><td>${{r.gene_id}}</td><td>${{r.effect.toFixed(5)}}</td><td>${{r.effect_type}}</td><td>${{r.fdr.toExponential(3)}}</td><td>${{r.significant?'是':'否'}}</td></tr>`).join('');status.textContent=`第 ${{page}} / ${{pages}} 页，共 ${{rows.length}} 条`;
  document.getElementById('volcano-title').textContent=`${{contrast.value.replace('_vs_',' vs ')}} · Volcano plot`;document.getElementById('heatmap-title').textContent=`${{contrast.value.replace('_vs_',' vs ')}} · Top-feature heatmap`;volcanoPlot(d);heatmapPlot(d);
}}
contrast.onchange=layer.onchange=()=>{{page=1;render()}};filter.oninput=()=>{{page=1;render()}};document.getElementById('prev-page').onclick=()=>{{page=Math.max(1,page-1);render()}};document.getElementById('next-page').onclick=()=>{{page++;render()}};
render();
window.addEventListener('beforeprint',()=>{{printing=true;render()}});window.addEventListener('afterprint',()=>{{printing=false;render()}});
</script></body></html>'''
    for source, target in {
        "长读长 RNA 标记分析结果报告": "ISO-SLAM分析报告",
        "长读长 RNA 标记分析 · 离线项目报告": "ISO-SLAM分析报告",
        ">差异火山图<": ">Volcano plot<",
        ">Top 特征热图<": ">Top-feature heatmap<",
    }.items():
        content = content.replace(source, target)
    content = content.replace('<section id=', '<section class="report-section" id=')
    content = content.replace('<figure class="figure-card">', '<article class="figure-card">').replace('</figure>', '</article>')
    first = next(iter(next(iter(result_payload.values())).values()))
    content = content.replace('<img id="volcano-static"', '<img src="' + first["volcano_uri"] + '" id="volcano-static"')
    content = content.replace('<img id="heatmap-static"', '<img src="' + first["heatmap_uri"] + '" id="heatmap-static"')
    launch_path = args.results_root / "run_launch_metadata.json"
    if launch_path.is_file():
        launch = json.loads(launch_path.read_text())
        runtime = '<details><summary>运行与参考信息</summary><p>' + esc(f"参考：{launch['reference']}；CPU上限：{launch['cores']}；内存上限：{launch['memory_mb']} MB；启动时间：{launch['started_utc']}") + '</p><pre>' + esc(json.dumps(launch.get('software_versions', {}), ensure_ascii=False, indent=2)) + '</pre></details>'
        content = content.replace('<h2>8. 方法与限制</h2>', '<h2>8. 方法与限制</h2>' + runtime)
    from report_components import enhance_report
    content = content.replace("<details open>", "<details>")
    content = enhance_report(content, expandable_figures=False)
    report = stage / "report.html"
    report.write_text(content, encoding="utf-8")
    manifest = stage / "05_metadata/delivery_manifest.tsv"
    manifest.parent.mkdir(parents=True, exist_ok=True)
    with manifest.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(["artifact", "description", "target", "status"])
        writer.writerow(["sample_qc_summary", "Technical QC summary", str(qc_path.resolve()), "up_to_date"])
        writer.writerow(["sample_qc_review", "Per-sample QC review metrics", str(review_path.resolve()), "up_to_date"])
        writer.writerow(["group_analysis", "Within- and between-group analysis", str(validation_path.resolve()), "up_to_date"])
        writer.writerow(["new_rna_splicing", "New-RNA alternative-splicing event analysis", str((splicing_root / "differential_new_rna_splicing.tsv").resolve()), "up_to_date"])
        writer.writerow(["new_rna_splicing_validation", "New-RNA alternative-splicing validation", str((splicing_root / "new_rna_splicing_validation.json").resolve()), "up_to_date"])
        writer.writerow(["final_project_report", "Final project report", str(report.resolve()), "up_to_date"])
    write_json(stage / "05_metadata/report_build.json", {"schema_version": "group-analysis-report-1", "status": "PASS", "sample_count": len(qc_rows), "group_inference": True, "contrast_count": len(validation["contrasts"]), "report": str(report)})
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
