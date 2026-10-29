#!/usr/bin/env python3
"""Build the self-contained Chinese ISO-SLAM project summary report."""
from __future__ import annotations
import argparse, base64, csv, html
from pathlib import Path

p=argparse.ArgumentParser(); p.add_argument("--results-root",required=True); p.add_argument("--samples",required=True); a=p.parse_args()
root=Path(a.results_root); cohort=root/"_project"/"06_cohort_qc"; de=root/"_project"/"07_differential_expression"; stage=root/"_project"/"08_project_report"
for d in ("01_tables","02_figures","03_logs","04_tool_output","05_metadata"): (stage/d).mkdir(parents=True,exist_ok=True)

def tsv(path):
    with path.open(encoding="utf-8",newline="") as h:return list(csv.DictReader(h,delimiter="\t"))
def fmt(x,n=2):
    try:return f"{float(x):,.{n}f}"
    except:return html.escape(str(x))
def data_image(path):
    if not path.is_file(): return '<div class="empty">图尚未生成</div>'
    return 'data:image/png;base64,'+base64.b64encode(path.read_bytes()).decode()
metrics=tsv(cohort/"01_tables"/"sample_metrics.tsv")
comparisons=[]
for kind in ("total","new"):
    idx=de/"01_tables"/f"{kind}_contrast_index.tsv"
    if idx.is_file(): comparisons.extend(tsv(idx))

metric_rows=''.join(f'<tr><td>{html.escape(r["sample_id"])}</td><td>{html.escape(r["group"])}</td><td>{int(r["raw_reads"]):,}</td><td>{float(r["trimmed_rate"]):.1%}</td><td>{float(r["primary_rate"]):.1%}</td><td>{float(r["assigned_rate"]):.1%}</td><td>{int(r["new_reads"]):,}</td><td>{float(r["tc_rate"]):.3%}</td><td>{int(r["genes_detected"]):,}</td></tr>' for r in metrics)
contrast_rows=''.join(f'<tr><td>{html.escape(r["rna_type"])}</td><td>{html.escape(r["comparison"])}</td><td>{int(r["tested"]):,}</td><td>{int(r["significant"]):,}</td></tr>' for r in comparisons)
sample_reports=''.join(f'<tr><td>{html.escape(r["sample_id"])}</td><td>{html.escape(r["group"])}</td><td class="mono">{html.escape(r["sample_report"])}</td></tr>' for r in metrics)
total_pca=data_image(cohort/"02_figures"/"total_pca.png"); new_pca=data_image(cohort/"02_figures"/"new_pca.png")
report=stage/"report.html"; manifest=stage/"05_metadata"/"delivery_manifest.tsv"
content=f'''<!doctype html><html lang="zh-CN"><head><meta charset="utf-8"><meta name="viewport" content="width=device-width"><title>ISO-SLAM 项目汇总报告</title>
<style>*{{box-sizing:border-box}}body{{margin:0;background:#f7f9fc;color:#172033;font:15px/1.65 Arial,"Noto Sans SC",sans-serif}}.sidebar{{position:fixed;left:0;top:0;bottom:0;width:220px;background:#10233f;color:white;padding:24px 18px}}.sidebar a{{display:block;color:white;text-decoration:none;padding:8px 0}}main{{margin-left:220px;padding:30px;max-width:1500px}}section{{background:white;border:1px solid #dde5ef;border-radius:8px;padding:22px;margin:0 0 22px}}h1{{font-size:27px}}h2{{color:#3C5488;border-bottom:2px solid #3C5488;padding-bottom:8px}}table{{width:100%;border-collapse:collapse;font-size:13px}}th,td{{padding:8px;border-bottom:1px solid #e6ebf2;text-align:left}}th{{background:#eef3f8}}.grid{{display:grid;grid-template-columns:1fr 1fr;gap:18px}}.figure{{cursor:zoom-in}}.figure img{{width:100%;background:white}}.mono{{font-family:monospace;font-size:12px}}.note,.empty{{color:#667085}}.footer{{margin-left:220px;padding:20px 30px;background:#10233f;color:white}}@media(max-width:850px){{.sidebar{{position:static;width:auto}}main,.footer{{margin-left:0}}.grid{{grid-template-columns:1fr}}}}</style></head><body>
<nav class="sidebar"><strong>ISO-SLAM 汇总</strong><a href="#overview">1 项目概览</a><a href="#qc">2 跨样本质控</a><a href="#de">3 组间差异</a><a href="#reports">4 单样本报告</a><a href="#limits">5 解释边界</a></nav><main>
<section id="overview"><h1>ISO-SLAM 项目汇总报告</h1><h2>1 项目概览</h2><p>人源 bulk 三代 SLAM-seq，共 4 组、每组 3 个生物学重复。正式 caller 为 site_binomial；long-TUC 6/20/30 仅作敏感性展示。</p></section>
<section id="qc"><h2>2 跨样本质控</h2><table data-page-size="12"><thead><tr><th>样本</th><th>组别</th><th>Raw reads</th><th>Insert率</th><th>MAPQ≥20率</th><th>基因分配率</th><th>New reads</th><th>T>C率</th><th>基因数</th></tr></thead><tbody>{metric_rows}</tbody></table><div class="grid"><div class="figure"><h3>Total RNA PCA</h3><img src="{total_pca}"></div><div class="figure"><h3>New RNA PCA</h3><img src="{new_pca}"></div></div></section>
<section id="de"><h2>3 组间差异</h2><p>total RNA 与 new RNA 分别使用 edgeR TMM、robust quasi-likelihood 模型；显著标准为 BH FDR&lt;0.05。</p><table data-page-size="12"><thead><tr><th>RNA</th><th>比较</th><th>检验基因</th><th>显著基因</th></tr></thead><tbody>{contrast_rows}</tbody></table></section>
<section id="reports"><h2>4 单样本报告</h2><table data-page-size="12"><thead><tr><th>样本</th><th>组别</th><th>报告路径</th></tr></thead><tbody>{sample_reports}</tbody></table></section>
<section id="limits"><h2>5 解释边界</h2><ul><li>NC 是有标记生物学对照，不用于估计无标记背景。</li><li>AB 不解释为已验证的协同或交互效应。</li><li>new fraction 仅报告描述性效应量，不进行正式显著性检验。</li><li>不报告绝对合成率、降解率或半衰期；transcript_id 仅作为基因分配证据。</li></ul></section></main><footer class="footer">ISO-SLAM · example_sample · 中文离线报告</footer></body></html>'''
report.write_text(content,encoding="utf-8")
with manifest.open("w",encoding="utf-8",newline="") as h:
    w=csv.writer(h,delimiter="\t",lineterminator="\n"); w.writerow(["description","target","status"])
    for r in metrics:w.writerow([f'{r["sample_id"]} sample report',r["sample_report"],"up_to_date"])
    w.writerow(["final project report",str(report.resolve()),"up_to_date"])
print(report)
