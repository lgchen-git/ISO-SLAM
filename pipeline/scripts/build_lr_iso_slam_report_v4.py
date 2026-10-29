#!/usr/bin/env python3
"""Build one self-contained customer-facing long-read RNA labeling report."""

from __future__ import annotations

import argparse
import base64
import csv
import html
import json
import re
from pathlib import Path

from iso_slam_v4_common import ensure_stage, write_json
from gene_annotation import add_gene_names, load_gtf_gene_map


SCHEMA = "iso-slam-sample-report-v7"
NPG = ["#3C5488", "#00A087", "#E64B35", "#4DBBD5", "#8491B4", "#F39B7F"]


def read_json(path: Path) -> dict[str, object]:
    return json.loads(path.read_text(encoding="utf-8"))


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open(encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def esc(value: object) -> str:
    return html.escape(str(value))


def pct(value: float) -> str:
    return f"{value:.1%}"


def bars(title: str, data: list[tuple[str, float]], percent: bool = False) -> str:
    data = sorted(data, key=lambda item: item[1], reverse=True)
    width, height = 760, 86 + 38 * len(data)
    maximum = max((value for _, value in data), default=1) or 1
    parts = [f'<svg viewBox="0 0 {width} {height}" role="img" aria-label="{esc(title)}"><text x="20" y="27" font-size="17" font-weight="700">{esc(title)}</text>']
    for index, (label, value) in enumerate(data):
        y = 50 + index * 38
        bar = 430 * value / maximum
        shown = pct(value) if percent else (f"{value:,.3f}" if value < 10 else f"{value:,.0f}")
        parts.append(f'<text x="20" y="{y+15}" font-size="13">{esc(label)}</text><rect x="200" y="{y}" width="{bar:.2f}" height="20" rx="4" fill="{NPG[index % len(NPG)]}"/><text x="740" y="{y+15}" text-anchor="end" font-size="12">{esc(shown)}</text>')
    return "".join(parts) + "</svg>"


def mismatch_spectrum(rows: list[dict[str, str]]) -> str:
    values = {(row["reference_base"], row["observed_base"]): float(row["probability"]) for row in rows}
    counts = {(row["reference_base"], row["observed_base"]): int(row["count"]) for row in rows}
    maximum = max(values.values(), default=0.0) or 1.0
    cells = ['<svg viewBox="0 0 760 300" role="img" aria-label="全碱基错配概率谱"><text x="20" y="28" font-size="17" font-weight="700">全碱基错配概率谱</text><text x="20" y="48" font-size="11" fill="#596b7d">分母为 high-quality primary alignment 中对应参考碱基的覆盖位点；单位为 ‰</text>']
    for row, ref in enumerate("ACGT"):
        y = 72 + row * 52
        cells.append(f'<text x="22" y="{y+20}" font-size="13">{ref} →</text>')
        col = 0
        for query in "ACGT":
            if query == ref:
                continue
            value = values.get((ref, query), 0.0)
            count = counts.get((ref, query), 0)
            x = 92 + col * 205
            opacity = 0.16 + 0.84 * value / maximum
            cells.append(f'<rect x="{x}" y="{y}" width="185" height="38" rx="5" fill="#3C5488" fill-opacity="{opacity:.3f}"/><text x="{x+9}" y="{y+16}" font-size="12" fill="#172033">{ref}&gt;{query}</text><text x="{x+9}" y="{y+31}" font-size="11" fill="#172033">{value * 1000:.3f}‰ · {count:,}</text>')
            col += 1
    return "".join(cells) + "</svg>"


def compact_table(rows: list[dict[str, object]], columns: list[tuple[str, str]], limit: int | None = None) -> str:
    if not rows:
        return '<div class="empty">该模块没有满足展示阈值的结果。</div>'
    head = "".join(f"<th>{esc(label)}</th>" for _, label in columns)
    body = "".join("<tr>" + "".join(f"<td>{esc(row.get(key, ''))}</td>" for key, _ in columns) + "</tr>" for row in (rows if limit is None else rows[:limit]))
    return f'<table data-page-size="10"><thead><tr>{head}</tr></thead><tbody>{body}</tbody></table>'


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--sample-id", required=True)
    parser.add_argument("--sample-root", required=True, type=Path)
    parser.add_argument("--v3-caller", required=True, type=Path)
    parser.add_argument("--report-title", required=True)
    parser.add_argument("--gtf", type=Path)
    parser.add_argument("--report-root", type=Path)
    parser.add_argument("--config", type=Path)
    args = parser.parse_args()
    sample = args.sample_id
    stage = (args.report_root or args.sample_root) / sample / "07_sample_report" if args.report_root else args.sample_root / "07_sample_report"
    ensure_stage(stage)
    trim = read_json(args.sample_root / f"02_trim/04_tool_output/trim/{sample}.trim_summary.json")
    catalog_val = read_json(args.sample_root / "04_isoform_catalog/05_metadata/catalog_validation.json")
    grand_val = read_json(args.sample_root / "05_slam_model/05_metadata/grand3_validation.json")
    dyn_val = read_json(args.sample_root / "06_isoform_dynamics/05_metadata/dynamics_validation.json")
    slam_qc = read_json(args.v3_caller / f"{sample}.bulk_slam_qc.json")
    raw, insert = int(trim["total_reads"]), int(trim["trimmed_reads"])
    joined = int(dyn_val["joined_reads"])
    q20 = int(slam_qc.get("total_reads", slam_qc.get("read_pass_reads", 0)))
    background_table = args.v3_caller / f"{sample}.substitution_summary.tsv"
    substitutions = read_tsv(background_table)
    tc_rows = [row for row in substitutions if row["substitution"] == "T>C"]
    background_rows = [row for row in substitutions if row["substitution"] == "diagnostic_background"]
    if len(tc_rows) != 1 or len(background_rows) != 1 or background_rows[0].get("source") != "max_ta_tg_rate_background_site_filtered":
        raise ValueError(f"{sample}: invalid T>C or diagnostic background summary")
    tc_rate = float(tc_rows[0]["rate"])
    background_rate = float(background_rows[0]["rate"])
    ntr = float(grand_val["global_parameters"]["binom_ntr"])
    catalog = read_tsv(args.sample_root / f"04_isoform_catalog/01_tables/{sample}.isoform_catalog.tsv")
    dynamics = read_tsv(args.sample_root / f"06_isoform_dynamics/01_tables/{sample}.isoform_dynamics.tsv")
    switches = read_tsv(args.sample_root / f"06_isoform_dynamics/01_tables/{sample}.isoform_switch_candidates.tsv")
    genes = read_tsv(args.sample_root / f"05_slam_model/01_tables/{sample}.grand3_gene_ntr.tsv")
    gene_map = load_gtf_gene_map(args.gtf) if args.gtf else {}
    add_gene_names(catalog, gene_map)
    add_gene_names(dynamics, gene_map)
    add_gene_names(switches, gene_map)
    add_gene_names(genes, gene_map)
    calls = dyn_val["posterior_call_counts"]
    screened = int(calls.get("Assignment_conflict", 0))
    no_t = int(calls.get("No_T_sites", 0))
    class_counts: dict[str, int] = {}
    for row in catalog:
        class_counts[row["structural_category"]] = class_counts.get(row["structural_category"], 0) + 1
    call_bars = [("New RNA", float(calls.get("New", 0))), ("Existing RNA", float(calls.get("Old", 0))), ("Ambiguous", float(calls.get("Ambiguous", 0)))]
    if no_t:
        call_bars.append(("No usable T sites", float(no_t)))
    if screened:
        call_bars.append(("Excluded after assignment screen", float(screened)))
    kpis = [("插入序列保留", pct(insert / raw if raw else 0), f"{insert:,} / {raw:,} reads"), ("高质量比对", pct(q20 / insert if insert else 0), f"{q20:,} reads"), ("转录本分配", pct(joined / q20 if q20 else 0), f"{joined:,} reads"), ("通过筛选的转录本", f"{int(catalog_val['filtered_isoforms']):,}", f"目录共 {int(catalog_val['catalog_isoforms']):,} 条"), ("全局新生 RNA 比例", pct(ntr), "基因层面估计"), ("一致性筛查", f"{screened:,}", pct(screened / joined if joined else 0))]
    kpi_html = "".join(f"<div class='kpi'><span>{esc(label)}</span><strong>{esc(value)}</strong><small>{esc(note)}</small></div>" for label, value, note in kpis)
    top_isoforms = sorted(dynamics, key=lambda row: float(row["expected_new_reads"]), reverse=True)[:10]
    top_genes = sorted(genes, key=lambda row: float(row["read_count"]), reverse=True)[:10]
    fragment = (args.sample_root / "07_sample_report/04_tool_output/read_flow/reports/read_flow_sankey_section.html").read_text(encoding="utf-8")
    fragment = fragment.replace("ISO-SLAM V4", "ISO-SLAM").replace("Assignment conflict", "Excluded after assignment screen")
    fragment = fragment.replace("No T sites", "No usable T sites" if no_t else "")
    fragment = fragment.replace("No_T_sites", "No usable T sites" if no_t else "")
    for source, target in {"原始 reads": "Raw reads", "插入序列 reads": "Insert reads", "高质量比对": "High-quality alignments", "转录本分配": "Transcript assignment", "新合成 RNA": "New RNA", "已有 RNA": "Existing RNA", "不确定": "Ambiguous", "未纳入转录本汇总": "Excluded after assignment screen"}.items():
        fragment = fragment.replace(source, target)
    fragment = fragment.replace("V4", "")
    plotly = "eval(atob(" + json.dumps(base64.b64encode(Path("pipeline/assets/plotly-3.5.1.min.js").read_bytes()).decode()) + "));"
    content = f'''<!doctype html><html lang="zh-CN"><head><meta charset="utf-8"><meta name="viewport" content="width=device-width"><title>{esc(sample)} · 长读长 RNA 标记分析</title><style>*{{box-sizing:border-box}}body{{margin:0;background:#f4f7fb;color:#182230;font:14px/1.65 Arial,"Noto Sans SC",sans-serif}}.sidebar{{position:fixed;inset:0 auto 0 0;width:230px;background:#17324d;color:#fff;padding:25px 18px;overflow:auto}}.sidebar a{{display:block;color:#e9f0f8;text-decoration:none;padding:8px 0}}main{{margin-left:230px;padding:30px;max-width:1540px}}section{{background:#fff;border:1px solid #dce5ef;border-radius:12px;padding:24px;margin-bottom:20px;box-shadow:0 2px 7px #17324d0d}}h1{{margin:0 0 4px;font-size:29px}}h2{{color:#284d73;border-bottom:2px solid #dce6f0;padding-bottom:8px}}.subtitle,.note{{color:#596b7d;background:#f5f8fb;padding:11px 13px;border-radius:7px}}.kpis{{display:grid;grid-template-columns:repeat(4,1fr);gap:12px;margin-top:18px}}.kpi{{border-top:4px solid #3C5488;background:#f7faff;padding:13px;border-radius:7px}}.kpi span,.kpi small{{display:block;color:#667085}}.kpi strong{{display:block;font-size:22px;margin:5px 0;color:#17324d}}.grid{{display:grid;grid-template-columns:1fr 1fr;gap:18px}}.figure{{border:1px solid #dde6ef;border-radius:9px;padding:12px;overflow:auto;background:#fff}}svg{{width:100%;height:auto}}details{{border:1px solid #dde6ef;border-radius:8px;padding:10px 13px;margin-top:14px}}summary{{cursor:pointer;font-weight:700;color:#385b7b}}table{{width:100%;border-collapse:collapse;font-size:12.5px;margin-top:9px}}th,td{{padding:7px 8px;border-bottom:1px solid #e5eaf0;text-align:left}}th{{background:#eef4f9}}.footer{{margin-left:230px;background:#17324d;color:#fff;padding:18px 30px}}@media(max-width:1000px){{.kpis{{grid-template-columns:repeat(2,1fr)}}.grid{{grid-template-columns:1fr}}}}@media(max-width:760px){{.sidebar{{position:static;width:auto}}main,.footer{{margin-left:0}}}}</style><script>__ISO_SLAM_PLOTLY__</script></head><body><nav class="sidebar"><strong>{esc(sample)}</strong><a href="#overview">1. 样本概览</a><a href="#flow">2. 数据处理与长读长质量</a><a href="#catalog">3. 转录本结构概览</a><a href="#signal">4. RNA 标记信号</a><a href="#dynamics">5. 转录本层面结果</a></nav><main><section id="overview"><h1>{esc(args.report_title)}</h1><p class="subtitle">样本 <strong>{esc(sample)}</strong> · 独立样本分析 · 齐碳全长 cDNA（按 ONT cDNA 类数据处理）</p><div class="kpis">{kpi_html}</div></section><section id="flow"><h2>2. 数据处理与长读长质量</h2>{fragment}<p class="note">流程图以 reads 为单位，展示从原始数据到转录本层面汇总的保留与未纳入分支。</p></section><section id="catalog"><h2>3. 转录本结构概览</h2><div class="grid"><div class="figure">{bars("转录本结构类别", class_counts.items())}</div><div class="figure">{bars("支持 reads 最高的基因", [(row['gene_id'], float(row['read_count'])) for row in top_genes])}</div></div></section><section id="signal"><h2>4. RNA 标记信号</h2><div class="grid"><div class="figure">{bars("T>C 与诊断背景转化率", [("T>C", tc_rate), ("诊断背景", background_rate)], True)}</div><div class="figure">{bars("read 分类", call_bars)}</div></div><p class="note">诊断背景由通过筛选的低深度背景位点 T&gt;A 与 T&gt;G 转化率较大值估计，仅用于样本内诊断，不代表未标记生物学对照。</p><p class="note">read 分类中的 assignment screen 用于核对转录本分配与后续汇总的一致性；未通过者不进入转录本层面的新生/已有 RNA 和使用变化结果。</p></section><section id="dynamics"><h2>5. 转录本层面结果</h2><div class="grid"><div class="figure">{bars("预期新生 reads 最高的转录本", [(row['isoform_id'], float(row['expected_new_reads'])) for row in top_isoforms])}</div><div class="figure">{bars("候选转录本使用变化", [(row['isoform_id'], abs(float(row['delta_isoform_usage']))) for row in switches[:10]], True)}</div></div><details><summary>查看 Top 转录本精简表</summary>{compact_table(dynamics, [('isoform_id','转录本'),('gene_name','基因名'),('gene_id','基因 ID'),('total_reads','Reads'),('posterior_new_fraction','新生比例')])}</details></section></main><footer class="footer">长读长 RNA 标记分析 · 离线结果报告</footer></body></html>'''
    visible_replacements = {
        f"{sample} · 长读长 RNA 标记分析": "ISO-SLAM分析报告",
        "长读长 RNA 标记分析 · 离线结果报告": "ISO-SLAM分析报告",
        "转录本结构类别": "Transcript structure categories",
        "支持 reads 最高的基因": "Top genes by supporting reads",
        "T&gt;C 与诊断背景转化率": "T&gt;C conversion and diagnostic background",
        "诊断背景": "Diagnostic background",
        "read 分类": "Read classification",
        "预期新生 reads 最高的转录本": "Top transcripts by expected new RNA reads",
        "候选转录本使用变化": "Candidate transcript usage changes",
    }
    for source, target in visible_replacements.items():
        content = content.replace(source, target)
    labels = {esc(row["gene_id"]): esc(f'{row["gene_name"]} ({row["gene_id"]})') for row in top_genes}
    labels.update({esc(row["isoform_id"]): esc(f'{row["gene_name"]} ({row["isoform_id"]})') for row in top_isoforms + switches[:10]})
    if labels:
        pattern = r"(?<![A-Za-z0-9_.-])(?:" + "|".join(re.escape(key) for key in labels) + r")(?![A-Za-z0-9_.-])"
        content = re.sub(pattern, lambda match: labels[match.group()], content)
    content = content.replace(
        "本报告展示独立样本的技术质量与描述性结果，不进行组间差异或显著性检验。",
        "本报告展示独立样本的技术质量与描述性结果；组间差异和显著性结果见项目总报告。",
    )
    from report_components import enhance_report
    if args.config:
        from iso_slam_v4_common import read_config
        run_config = read_config(args.config)
        content = content.replace("NC 是已标记的生物学对照，不作为未标记背景对照。", esc(run_config["control_group"]) + " 是已标记的生物学对照，不作为未标记背景对照。")
    content = enhance_report(content)
    content = content.replace("__ISO_SLAM_PLOTLY__", plotly)
    report = stage / f"{sample}_report.html"
    report.write_text(content, encoding="utf-8")
    artifacts = [args.sample_root / f"04_isoform_catalog/01_tables/{sample}.isoform_catalog.tsv", background_table, args.sample_root / f"06_isoform_dynamics/01_tables/{sample}.isoform_dynamics.tsv", args.sample_root / f"06_isoform_dynamics/01_tables/{sample}.isoform_switch_candidates.tsv", report]
    manifest = stage / "05_metadata/delivery_manifest.tsv"
    with manifest.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(["artifact", "description", "target", "status"])
        for path in artifacts[:-1]:
            writer.writerow([path.stem, path.name, str(path.resolve()), "up_to_date"])
        writer.writerow(["final_sample_report", f"{sample} final report", str(report.resolve()), "up_to_date"])
    write_json(stage / "05_metadata/report_build.json", {"schema_version": SCHEMA, "status": "PASS", "sample_id": sample, "report": str(report), "single_final_report": True, "offline": True, "no_t_sites_displayed": bool(no_t)})
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
