# 05_sample_report

Stage wrapper directory. The reusable implementation remains in the top-level scripts directory and is invoked by `workflow/run.sh`.

Report schema V3 is generated from compact exact tables. It contains six overview KPIs, the conserved reads Sankey, table-only NanoPlot QC, a T>C-versus-diagnostic-background chart, and a five-category assignment composition. The single final HTML is named `<sample_id>.report.html`, is self-contained, and must pass the delivery validator; do not hand-edit it.
