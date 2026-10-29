#!/usr/bin/env python3
"""Calculate independent-sample SUPPA2 AS effect sizes and record APA availability."""

from __future__ import annotations

import argparse
import csv
from pathlib import Path

from iso_slam_v4_common import ensure_stage, executable, nonempty, run_logged, write_json, write_tsv


SCHEMA = "iso-slam-isoform-events-v4"


def arguments() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--sample-id", required=True)
    parser.add_argument("--filtered-gtf", required=True, type=Path)
    parser.add_argument("--isoform-dynamics", required=True, type=Path)
    parser.add_argument("--outdir", required=True, type=Path)
    parser.add_argument("--suppa-bin", default="suppa.py")
    parser.add_argument("--polya-table", type=Path)
    parser.add_argument("--polya-evaluable-fraction", type=float, default=0.0)
    parser.add_argument("--polya-min-fraction", type=float, default=0.70)
    parser.add_argument("--skip-suppa", action="store_true")
    return parser.parse_args()


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open(encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def main() -> int:
    args = arguments()
    ensure_stage(args.outdir)
    dynamics = read_tsv(args.isoform_dynamics)
    expression = args.outdir / f"04_tool_output/suppa2/{args.sample_id}.isoform_expression.tsv"
    expression.parent.mkdir(parents=True, exist_ok=True)
    with expression.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(["total", "new"])
        for row in dynamics:
            writer.writerow([row["isoform_id"].split("|", 1)[-1], row["total_reads"], row["expected_new_reads"]])

    prefix = args.outdir / f"04_tool_output/suppa2/{args.sample_id}"
    as_rows: list[dict[str, object]] = []
    suppa_status = "skipped_by_request" if args.skip_suppa else "complete"
    if not args.skip_suppa:
        suppa = executable(args.suppa_bin)
        event_marker = Path(str(prefix) + "_SE_strict.ioe")
        if not nonempty(event_marker):
            run_logged(
                [suppa, "generateEvents", "-i", str(args.filtered_gtf), "-o", str(prefix), "-f", "ioe", "-e", "SE", "SS", "MX", "RI", "FL"],
                args.outdir / f"03_logs/{args.sample_id}.suppa_generate.log",
            )
        for ioe in sorted(prefix.parent.glob(prefix.name + "_*.ioe")):
            psi_prefix = Path(str(ioe).removesuffix(".ioe") + ".psi")
            psi_file = Path(str(psi_prefix) + ".psi")
            if not nonempty(psi_file):
                run_logged(
                    [suppa, "psiPerEvent", "-i", str(ioe), "-e", str(expression), "-o", str(psi_prefix)],
                    args.outdir / f"03_logs/{args.sample_id}.suppa_psi_{ioe.stem}.log",
                )
            if not nonempty(psi_file):
                continue
            with psi_file.open(encoding="utf-8", newline="") as handle:
                for row in csv.DictReader(handle, delimiter="\t"):
                    event_id = row.get("event_id") or row.get("EventID") or next(iter(row.values()), "")
                    try:
                        total_psi = float(row.get("total", "nan"))
                        new_psi = float(row.get("new", "nan"))
                    except ValueError:
                        continue
                    as_rows.append(
                        {
                            "sample_id": args.sample_id,
                            "event_id": event_id,
                            "event_type": ioe.stem.split("_")[-2] if "_" in ioe.stem else "unknown",
                            "total_psi": total_psi,
                            "new_psi": new_psi,
                            "delta_psi": new_psi - total_psi,
                            "inference": "effect_size_only",
                        }
                    )
    as_fields = ["sample_id", "event_id", "event_type", "total_psi", "new_psi", "delta_psi", "inference"]
    write_tsv(args.outdir / f"01_tables/{args.sample_id}.alternative_splicing.tsv", as_rows, as_fields)

    apa_rows: list[dict[str, object]] = []
    apa_status = "not_evaluable"
    apa_closure_error = 0.0
    if args.polya_evaluable_fraction >= args.polya_min_fraction:
        apa_status = "complete"
        gene_total: dict[str, float] = {}
        gene_new: dict[str, float] = {}
        for row in dynamics:
            gene = row["gene_id"]
            gene_total[gene] = gene_total.get(gene, 0.0) + float(row["total_reads"])
            gene_new[gene] = gene_new.get(gene, 0.0) + float(row["expected_new_reads"])
        for row in dynamics:
            gene = row["gene_id"]
            total_usage = float(row["total_reads"]) / gene_total[gene] if gene_total[gene] else 0.0
            new_usage = float(row["expected_new_reads"]) / gene_new[gene] if gene_new[gene] else 0.0
            apa_rows.append(
                {
                    "sample_id": args.sample_id,
                    "gene_id": gene,
                    "site_id": row["isoform_id"],
                    "total_usage": total_usage,
                    "new_usage": new_usage,
                    "delta_usage": new_usage - total_usage,
                    "site_definition": "IsoQuant_transcript_end",
                }
            )
        for gene in gene_total:
            gene_rows = [row for row in apa_rows if row["gene_id"] == gene]
            apa_closure_error = max(apa_closure_error, abs(sum(float(row["total_usage"]) for row in gene_rows) - 1.0))
            if gene_new[gene] > 0:
                apa_closure_error = max(apa_closure_error, abs(sum(float(row["new_usage"]) for row in gene_rows) - 1.0))
    apa_fields = ["sample_id", "gene_id", "site_id", "total_usage", "new_usage", "delta_usage", "site_definition"]
    write_tsv(args.outdir / f"01_tables/{args.sample_id}.apa_usage.tsv", apa_rows, apa_fields)
    write_json(
        args.outdir / "05_metadata/events_validation.json",
        {
            "schema_version": SCHEMA,
            "status": "PASS",
            "sample_id": args.sample_id,
            "suppa_status": suppa_status,
            "as_events": len(as_rows),
            "as_inference": "effect_size_only_no_p_values",
            "apa_status": apa_status,
            "polya_evaluable_fraction": args.polya_evaluable_fraction,
            "polya_min_fraction": args.polya_min_fraction,
            "apa_usage_closure_error": apa_closure_error,
            "apa_usage_closure_status": "PASS" if apa_closure_error <= 1e-6 else "FAIL",
        },
    )
    if apa_closure_error > 1e-6:
        raise RuntimeError("APA usage denominators do not close")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
