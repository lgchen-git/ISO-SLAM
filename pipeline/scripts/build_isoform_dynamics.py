#!/usr/bin/env python3
"""Join IsoQuant assignments with T>C counts and calculate read/isoform posteriors."""

from __future__ import annotations

import argparse
import csv
import gzip
import hashlib
import math
import random
from collections import Counter, defaultdict
from pathlib import Path

from iso_slam_v4_common import ensure_stage, write_json, write_tsv


SCHEMA = "iso-slam-isoform-dynamics-v6"
CALLER_REQUIRED_COLUMNS = {"read_id", "gene_id", "region_class", "assignment_level", "t_covered"}
ASSIGNMENT_REQUIRED_COLUMNS = {"read_id", "isoform_id", "gene_id"}


def arguments() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--sample-id", required=True)
    parser.add_argument("--read-calls", required=True, type=Path)
    parser.add_argument("--read-assignments", required=True, type=Path)
    parser.add_argument("--catalog", required=True, type=Path)
    parser.add_argument("--grand3-global", required=True, type=Path)
    parser.add_argument("--grand3-genes", required=True, type=Path)
    parser.add_argument("--outdir", required=True, type=Path)
    parser.add_argument("--new-threshold", type=float, default=0.50)
    parser.add_argument("--old-threshold", type=float, default=0.30)
    parser.add_argument("--bootstrap-iterations", type=int, default=500)
    parser.add_argument("--bootstrap-seed", type=int, default=20260818)
    parser.add_argument("--gene-min-reads", type=int, default=20)
    parser.add_argument("--isoform-min-reads", type=int, default=5)
    parser.add_argument("--switch-min-delta", type=float, default=0.10)
    parser.add_argument("--concordance-threshold", type=float, default=0.05)
    parser.add_argument("--baseline-read-posterior", type=Path)
    return parser.parse_args()


def open_text(path: Path):
    return gzip.open(path, "rt", encoding="utf-8", newline="") if path.suffix == ".gz" else path.open(encoding="utf-8", newline="")


def tsv(path: Path) -> list[dict[str, str]]:
    with open_text(path) as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def required_columns(rows: list[dict[str, str]], required: set[str], label: str) -> None:
    if not rows:
        raise RuntimeError(f"{label} is empty")
    missing = sorted(required - set(rows[0]))
    if missing:
        raise RuntimeError(f"{label} is missing required columns: {', '.join(missing)}")


def normalized_gene(value: str | None) -> str:
    cleaned = (value or "").strip()
    return "" if cleaned.lower() in {"", ".", "na", "none", "null"} else cleaned


def conflict_reasons(caller: dict[str, str], assignment: dict[str, str]) -> list[str]:
    reasons: list[str] = []
    region = (caller.get("region_class") or "").strip().lower()
    level = (caller.get("assignment_level") or "").strip().lower()
    if region == "intergenic":
        reasons.append("Caller_intergenic")
    if region == "ambiguous":
        reasons.append("Caller_ambiguous")
    if level == "unresolved":
        reasons.append("Unresolved_caller_assignment")
    if normalized_gene(caller.get("gene_id")) != normalized_gene(assignment.get("gene_id")):
        reasons.append("Gene_id_conflict")
    return reasons


def binom_logpmf(k: int, n: int, p: float) -> float:
    if n < 0 or k < 0 or k > n:
        return -math.inf
    if p <= 0:
        return 0.0 if k == 0 else -math.inf
    if p >= 1:
        return 0.0 if k == n else -math.inf
    return math.lgamma(n + 1) - math.lgamma(k + 1) - math.lgamma(n - k + 1) + k * math.log(p) + (n - k) * math.log1p(-p)


def read_posterior(k: int, n: int, prior: float, p_err: float, p_conv: float) -> float:
    prior = min(max(prior, 1e-12), 1 - 1e-12)
    p_old = min(max(p_err, 1e-12), 1 - 1e-12)
    p_new = 1 - (1 - p_err) * (1 - p_conv)
    p_new = min(max(p_new, 1e-12), 1 - 1e-12)
    if n == 0:
        return prior
    old_log = math.log1p(-prior) + binom_logpmf(k, n, p_old)
    new_log = math.log(prior) + binom_logpmf(k, n, p_new)
    pivot = max(old_log, new_log)
    return math.exp(new_log - pivot) / (math.exp(new_log - pivot) + math.exp(old_log - pivot))


def percentile(values: list[float], q: float) -> float:
    if not values:
        return math.nan
    ordered = sorted(values)
    position = (len(ordered) - 1) * q
    lower = int(math.floor(position))
    upper = int(math.ceil(position))
    if lower == upper:
        return ordered[lower]
    return ordered[lower] * (upper - position) + ordered[upper] * (position - lower)


def bootstrap_ci(values: list[float], iterations: int, seed: int) -> tuple[float, float]:
    if not values:
        return math.nan, math.nan
    if len(values) == 1 or iterations <= 0:
        return values[0], values[0]
    rng = random.Random(seed)
    means = [sum(rng.choice(values) for _ in values) / len(values) for _ in range(iterations)]
    return percentile(means, 0.025), percentile(means, 0.975)


def main() -> int:
    args = arguments()
    if not 0 <= args.old_threshold < args.new_threshold <= 1:
        raise SystemExit("posterior thresholds must satisfy 0 <= old < new <= 1")
    ensure_stage(args.outdir)
    global_rows = tsv(args.grand3_global)
    if len(global_rows) != 1:
        raise RuntimeError("Grand3 global parameter table must contain exactly one row")
    params = global_rows[0]
    prior = float(params["binom_ntr"])
    p_err = float(params["binom_p_err"])
    p_conv = float(params["binom_p_conv"])
    assignment_rows = tsv(args.read_assignments)
    required_columns(assignment_rows, ASSIGNMENT_REQUIRED_COLUMNS, "IsoQuant read assignments")
    assignments = {row["read_id"]: row for row in assignment_rows}
    catalog_rows = tsv(args.catalog)
    passed = {row["source_isoform_id"]: row for row in catalog_rows if row["filter_status"] == "PASS"}
    grand_genes = {row["gene_id"]: row for row in tsv(args.grand3_genes)}
    expected_join_reads = sum(1 for row in assignments.values() if row.get("isoform_id") in passed)

    read_rows: list[dict[str, object]] = []
    iso_values: dict[str, list[float]] = defaultdict(list)
    iso_gene: dict[str, str] = {}
    call_counts: Counter = Counter({"New": 0, "Old": 0, "Ambiguous": 0, "No_T_sites": 0, "Assignment_conflict": 0})
    conflict_subtype_counts: Counter = Counter()
    unmatched = 0
    with open_text(args.read_calls) as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        missing = sorted(CALLER_REQUIRED_COLUMNS - set(reader.fieldnames or []))
        if missing:
            raise RuntimeError(f"SLAM read calls are missing required columns: {', '.join(missing)}")
        for row in reader:
            assignment = assignments.get(row["read_id"])
            if assignment is None:
                unmatched += 1
                continue
            isoform = assignment["isoform_id"]
            if isoform not in passed:
                continue
            n = int(float(row.get("t_covered", "0") or 0))
            k = int(float(row.get("tc_count", row.get("tc_events", "0")) or 0))
            posterior = read_posterior(k, n, prior, p_err, p_conv)
            reasons = conflict_reasons(row, assignment)
            if reasons:
                call = "Assignment_conflict"
                conflict_subtype_counts[reasons[0]] += 1
            elif n == 0:
                call = "No_T_sites"
            elif posterior >= args.new_threshold:
                call = "New"
            elif posterior <= args.old_threshold:
                call = "Old"
            else:
                call = "Ambiguous"
            gene = assignment.get("gene_id") or passed[isoform].get("gene_id", ".")
            prefixed = f"{args.sample_id}|{isoform}"
            if not reasons:
                iso_values[prefixed].append(posterior)
                iso_gene[prefixed] = gene
            call_counts[call] += 1
            read_rows.append(
                {
                    "sample_id": args.sample_id,
                    "read_id": row["read_id"],
                    "gene_id": gene,
                    "caller_gene_id": row.get("gene_id", "."),
                    "caller_region_class": row.get("region_class", "."),
                    "caller_assignment_level": row.get("assignment_level", "."),
                    "assignment_conflict_subtype": reasons[0] if reasons else "",
                    "assignment_conflict_reasons": ";".join(reasons),
                    "isoform_id": prefixed,
                    "structural_category": passed[isoform].get("structural_category", "unknown"),
                    "t_covered": n,
                    "tc_count": k,
                    "grand3_binom_ntr_prior": prior,
                    "grand3_binom_p_err": p_err,
                    "grand3_binom_p_conv": p_conv,
                    "p_new": posterior,
                    "posterior_call": call,
                }
            )
    if not read_rows:
        raise RuntimeError("No filtered isoform assignments joined to SLAM read calls")
    if not iso_values:
        raise RuntimeError("No consistency-qualified reads remain for isoform dynamics")
    read_path = args.outdir / f"01_tables/{args.sample_id}.read_new_posterior.tsv.gz"
    read_tmp = read_path.with_suffix(read_path.suffix + ".tmp")
    with gzip.open(read_tmp, "wt", encoding="utf-8", newline="") as handle:
        fields = list(read_rows[0])
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(read_rows)
    read_tmp.replace(read_path)

    transition_rows: list[dict[str, object]] = []
    if args.baseline_read_posterior:
        baseline_rows = tsv(args.baseline_read_posterior)
        required_columns(baseline_rows, {"read_id", "posterior_call"}, "baseline read posterior")
        baseline_calls = {row["read_id"]: row["posterior_call"] for row in baseline_rows}
        transitions: Counter = Counter((baseline_calls.get(str(row["read_id"]), "Not_joined"), str(row["posterior_call"])) for row in read_rows)
        transition_rows = [
            {"sample_id": args.sample_id, "baseline_posterior_call": baseline, "revised_posterior_call": revised, "read_count": count}
            for (baseline, revised), count in sorted(transitions.items())
        ]
        write_tsv(
            args.outdir / f"01_tables/{args.sample_id}.posterior_call_reclassification.tsv",
            transition_rows,
            ["sample_id", "baseline_posterior_call", "revised_posterior_call", "read_count"],
        )

    gene_total: Counter = Counter()
    gene_new: Counter = Counter()
    for isoform, values in iso_values.items():
        gene_total[iso_gene[isoform]] += len(values)
        gene_new[iso_gene[isoform]] += sum(values)
    isoform_rows: list[dict[str, object]] = []
    switch_rows: list[dict[str, object]] = []
    for isoform, values in sorted(iso_values.items(), key=lambda item: len(item[1]), reverse=True):
        gene = iso_gene[isoform]
        expected_new = sum(values)
        expected_old = len(values) - expected_new
        fraction = expected_new / len(values)
        seed_offset = int(hashlib.sha256(isoform.encode()).hexdigest()[:8], 16)
        ci_low, ci_high = bootstrap_ci(values, args.bootstrap_iterations, args.bootstrap_seed + seed_offset)
        total_usage = len(values) / gene_total[gene] if gene_total[gene] else 0.0
        new_usage = expected_new / gene_new[gene] if gene_new[gene] else 0.0
        delta = new_usage - total_usage
        status = (
            "candidate"
            if gene_total[gene] >= args.gene_min_reads and len(values) >= args.isoform_min_reads and abs(delta) >= args.switch_min_delta
            else "below_threshold"
        )
        item = {
            "sample_id": args.sample_id,
            "gene_id": gene,
            "isoform_id": isoform,
            "total_reads": len(values),
            "expected_new_reads": expected_new,
            "expected_old_reads": expected_old,
            "posterior_new_fraction": fraction,
            "technical_bootstrap_ci_low": ci_low,
            "technical_bootstrap_ci_high": ci_high,
            "total_isoform_usage": total_usage,
            "new_isoform_usage": new_usage,
            "delta_isoform_usage": delta,
            "switch_status": status,
        }
        isoform_rows.append(item)
        if status == "candidate":
            switch_rows.append(item)
    write_tsv(args.outdir / f"01_tables/{args.sample_id}.isoform_dynamics.tsv", isoform_rows)
    write_tsv(
        args.outdir / f"01_tables/{args.sample_id}.isoform_switch_candidates.tsv",
        switch_rows,
        list(isoform_rows[0]),
    )
    gene_rows = []
    for gene, total in gene_total.most_common():
        expected_new = gene_new[gene]
        gene_rows.append(
            {
                "sample_id": args.sample_id,
                "gene_id": gene,
                "assigned_reads": total,
                "posterior_expected_new_reads": expected_new,
                "posterior_expected_old_reads": total - expected_new,
                "posterior_new_fraction": expected_new / total if total else 0.0,
                "grand3_binom_ntr_map": grand_genes.get(gene, {}).get("grand3_binom_ntr_map", "NA"),
            }
        )
    write_tsv(args.outdir / f"01_tables/{args.sample_id}.gene_dynamics.tsv", gene_rows)
    eligible_rows = [row for row in read_rows if row["posterior_call"] != "Assignment_conflict"]
    posterior_mean = sum(float(row["p_new"]) for row in eligible_rows) / len(eligible_rows)
    no_t_sites_consistency = all(
        row["posterior_call"] != "No_T_sites"
        or (
            int(row["t_covered"]) == 0
            and not str(row["assignment_conflict_subtype"])
            and normalized_gene(str(row["caller_gene_id"])) == normalized_gene(str(row["gene_id"]))
        )
        for row in read_rows
    )
    validation = {
        "schema_version": SCHEMA,
        "status": "PASS",
        "sample_id": args.sample_id,
        "joined_reads": len(read_rows),
        "dynamics_eligible_reads": len(eligible_rows),
        "assignment_conflict_reads": int(call_counts["Assignment_conflict"]),
        "assignment_conflict_subtype_counts": dict(conflict_subtype_counts),
        "unmatched_slam_reads": unmatched,
        "posterior_call_counts": dict(call_counts),
        "posterior_mean": posterior_mean,
        "grand3_global_binom_ntr": prior,
        "posterior_grand3_absolute_difference": abs(posterior_mean - prior),
        "posterior_concordance_threshold": args.concordance_threshold,
        "posterior_concordance_status": "PASS" if abs(posterior_mean - prior) <= args.concordance_threshold else "FAIL",
        "no_t_sites_consistency_status": "PASS" if no_t_sites_consistency else "FAIL",
        "expected_filtered_assignment_reads": expected_join_reads,
        "read_id_join_fraction": len(read_rows) / expected_join_reads if expected_join_reads else 0.0,
        "weighted_closure_error": max(abs(float(row["expected_new_reads"]) + float(row["expected_old_reads"]) - int(row["total_reads"])) for row in isoform_rows),
        "technical_ci_only": True,
        "biological_p_values": False,
        "switch_candidates": len(switch_rows),
        "baseline_reclassification_rows": len(transition_rows),
    }
    write_json(args.outdir / "05_metadata/dynamics_validation.json", validation)
    if validation["weighted_closure_error"] > 1e-6 or validation["posterior_concordance_status"] != "PASS" or validation["no_t_sites_consistency_status"] != "PASS" or validation["read_id_join_fraction"] < 0.99:
        raise RuntimeError("Isoform dynamics validation failed")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
