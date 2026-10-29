#!/usr/bin/env python3
from __future__ import annotations
import csv, subprocess, sys, tempfile
from pathlib import Path

ROOT=Path(__file__).resolve().parents[1]
GROUPS=["A","B","AB","NC"]

def main():
  with tempfile.TemporaryDirectory(prefix="iso_slam_project_") as raw:
    tmp=Path(raw); results=tmp/"results"; cohort=results/"_project"/"06_cohort_qc"; tables=cohort/"01_tables"; tables.mkdir(parents=True)
    samples=[]
    for group in GROUPS:
      for rep in range(1,4): samples.append((f"{group}-{rep}",group))
    sheet=tmp/"samples.tsv"
    with sheet.open("w",encoding="utf-8",newline="") as h:
      w=csv.writer(h,delimiter="\t",lineterminator="\n"); w.writerow(["sample_id","group"]); w.writerows(samples)
    for kind in ("total","new"):
      with (tables/f"{kind}_rna_counts.tsv").open("w",encoding="utf-8",newline="") as h:
        w=csv.writer(h,delimiter="\t",lineterminator="\n"); w.writerow(["gene_id","gene_name",*[s for s,_ in samples]])
        for g in range(1,101):
          vals=[]
          for i,(_,grp) in enumerate(samples):
            effect=60 if grp=="A" and g<=15 else 0
            base=30+g+(i%3)*2+effect
            vals.append(max(1,base if kind=="total" else base//3))
          w.writerow([f"ENSG{g:06d}",f"GENE{g}",*vals])
    sample_metrics=[]
    for sample,group in samples:
      report=results/sample/"05_sample_report"/f"{sample}.report.html"; report.parent.mkdir(parents=True); report.write_text("sample report",encoding="utf-8")
      sample_metrics.append([sample,group,200000,0.92,180000,0.90,170000,0.94,80000,85000,5000,0.012,15000,str(report)])
    with (tables/"sample_metrics.tsv").open("w",encoding="utf-8",newline="") as h:
      w=csv.writer(h,delimiter="\t",lineterminator="\n"); w.writerow(["sample_id","group","raw_reads","trimmed_rate","trimmed_reads","primary_rate","primary_reads","assigned_rate","new_reads","old_reads","ambiguous_reads","tc_rate","genes_detected","sample_report"]); w.writerows(sample_metrics)
    subprocess.run(["Rscript",str(ROOT/"pipeline/scripts/run_project_edger.R"),str(results),str(sheet),"full"],check=True,cwd=ROOT)
    subprocess.run([sys.executable,str(ROOT/"pipeline/scripts/build_project_report.py"),"--results-root",str(results),"--samples",str(sheet)],check=True,cwd=ROOT)
    report=results/"_project"/"08_project_report"/"report.html"; manifest=results/"_project"/"08_project_report"/"05_metadata"/"delivery_manifest.tsv"
    subprocess.run([sys.executable,str(ROOT/"pipeline/scripts/00_shared/validate_delivery_report.py"),"--report",str(report),"--manifest",str(manifest),"--require-section","1 项目概览"],check=True,cwd=ROOT)

if __name__=="__main__": main()
