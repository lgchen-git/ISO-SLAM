#!/usr/bin/env python3
"""Run the pinned Grand3 compatibility pilot and compare sample-level NTR estimates."""

from __future__ import annotations

import argparse
import csv
import gzip
import json
import os
import re
import subprocess
from pathlib import Path


def read_config(path: str) -> dict[str, str]:
    config: dict[str, str] = {}
    with Path(path).open(encoding="utf-8") as handle:
        for raw in handle:
            if not raw.strip() or raw.lstrip().startswith("#") or raw[:1].isspace() or ":" not in raw:
                continue
            key, value = raw.split(":", 1)
            config[key.strip()] = value.split("#", 1)[0].strip().strip("'\"")
    return config


def file_nonempty(path: Path) -> bool:
    return path.is_file() and path.stat().st_size > 0


def run_logged(command: list[str], log_path: Path, env: dict[str, str]) -> None:
    log_path.parent.mkdir(parents=True, exist_ok=True)
    with log_path.open("w", encoding="utf-8") as log:
        log.write("$ " + " ".join(command) + "\n")
        log.flush()
        subprocess.run(command, check=True, stdout=log, stderr=subprocess.STDOUT, env=env)


def parse_single_row(path: Path) -> dict[str, str]:
    with path.open(encoding="utf-8", newline="") as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    if len(rows) != 1:
        raise RuntimeError(f"Expected one row in {path}, found {len(rows)}")
    return rows[0]


def select_column(header: list[str], sample_id: str, suffix: str) -> str:
    matches = [name for name in header if sample_id in name and name.endswith(suffix)]
    if len(matches) != 1:
        raise RuntimeError(f"Could not uniquely select {suffix!r} for {sample_id}: {matches}")
    return matches[0]


def read_global_binomial_ntr(
    model_parameters: Path,
    sample_ids: list[str],
    strandness: str,
) -> dict[str, float]:
    """Read sample-level global NTRs from Grand3's Separate binomial models."""
    values: dict[str, float] = {}
    with model_parameters.open(encoding="utf-8", newline="") as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            if row.get("Estimator") != "Separate" or row.get("Label") != "4sU":
                continue
            if row.get("Subread") != strandness:
                continue
            condition_name = Path(row.get("Condition", "")).name
            for sample_id in sample_ids:
                if condition_name == f"{sample_id}.4sU":
                    if sample_id in values:
                        raise RuntimeError(f"Duplicate Grand3 global model for {sample_id}")
                    values[sample_id] = float(row["Binom ntr"])
    missing = [sample_id for sample_id in sample_ids if sample_id not in values]
    if missing:
        raise RuntimeError(f"Missing Grand3 global binomial NTR for: {', '.join(missing)}")
    return values


def compare_sample_ntr(
    grand3_table: Path,
    model_parameters: Path,
    results_root: Path,
    sample_ids: list[str],
    strandness: str,
    threshold: float,
    output: Path,
) -> tuple[bool, list[dict[str, object]]]:
    global_ntr = read_global_binomial_ntr(model_parameters, sample_ids, strandness)
    opener = gzip.open if grand3_table.suffix == ".gz" else open
    with opener(grand3_table, "rt", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        header = reader.fieldnames or []
        grand_rows = list(reader)
    comparisons: list[dict[str, object]] = []
    for sample_id in sample_ids:
        count_col = select_column(header, sample_id, "Read count")
        ntr_col = select_column(header, sample_id, "4sU Binom NTR MAP")
        numerator = 0.0
        denominator = 0.0
        genes = 0
        for row in grand_rows:
            try:
                count = float(row[count_col])
                ntr = float(row[ntr_col])
            except (TypeError, ValueError):
                continue
            if count <= 0:
                continue
            numerator += count * ntr
            denominator += count
            genes += 1
        gene_weighted_ntr = numerator / denominator if denominator else float("nan")
        ours = parse_single_row(
            results_root
            / sample_id
            / "04_slam_calling/04_tool_output/caller"
            / f"{sample_id}.bulk_slam_summary.tsv"
        )
        ours_ntr = float(ours["new_fraction_excluding_ambiguous"])
        grand3_ntr = global_ntr[sample_id]
        absolute_error = abs(grand3_ntr - ours_ntr)
        comparisons.append(
            {
                "sample_id": sample_id,
                "our_new_fraction_excluding_ambiguous": ours_ntr,
                "grand3_global_binomial_ntr": grand3_ntr,
                "grand3_gene_weighted_binomial_ntr_diagnostic": gene_weighted_ntr,
                "absolute_error": absolute_error,
                "grand3_genes_used": genes,
                "status": "PASS" if absolute_error <= threshold else "FAIL",
            }
        )
    output.parent.mkdir(parents=True, exist_ok=True)
    with output.open("w", encoding="utf-8", newline="") as handle:
        fields = list(comparisons[0]) if comparisons else ["sample_id", "status"]
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(comparisons)
    return bool(comparisons) and all(row["status"] == "PASS" for row in comparisons), comparisons


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--config", required=True)
    parser.add_argument("--samples", required=True)
    parser.add_argument("--results-root", required=True)
    parser.add_argument("--threads", type=int, default=8)
    args = parser.parse_args()

    config = read_config(args.config)
    gedi = Path(config["grand_slam_bin"])
    java_bin = config.get("grand_slam_java_bin", "")
    genomic = Path(config["grand_slam_genomic"])
    genome_fa = config["genome_fa"]
    gtf = config["gtf"]
    tmp_root = Path(config["tmp_root"])
    results_root = Path(args.results_root)
    threshold = float(config.get("grand_slam_ntr_mae_max", "0.05"))

    with Path(args.samples).open(encoding="utf-8", newline="") as handle:
        sample_ids = [row["sample_id"] for row in csv.DictReader(handle, delimiter="\t")]
    if not sample_ids:
        raise SystemExit("No pilot samples")

    project_stage = results_root / "_project/06_cohort_qc"
    output_root = project_stage / "04_tool_output/grand3"
    metadata_root = project_stage / "05_metadata"
    table_root = project_stage / "01_tables"
    work_root = tmp_root / "grand3_pilot/work"
    log_root = tmp_root / "grand3_pilot/logs"
    gedi_home = tmp_root / "grand3_pilot/gedi_home"
    for path in (output_root, metadata_root, table_root, work_root, log_root, gedi_home / ".gedi", genomic.parent):
        path.mkdir(parents=True, exist_ok=True)

    env = os.environ.copy()
    env["PATH"] = os.pathsep.join(
        part for part in (java_bin, str(gedi.parent), env.get("PATH", "")) if part
    )
    common = [str(gedi), "-Duser.home=" + str(gedi_home)]

    if not file_nonempty(genomic):
        index_cmd = common + [
            "-mem",
            "16G",
            "-e",
            "IndexGenome",
            "-s",
            genome_fa,
            "-a",
            gtf,
            "-f",
            str(genomic.parent),
            "-n",
            genomic.stem,
            "-o",
            str(genomic),
            "-nomapping",
            "-p",
            "-D",
        ]
        run_logged(index_cmd, log_root / "01_index_genome.log", env)
    if not file_nonempty(genomic):
        raise RuntimeError(f"GEDI genome index was not created: {genomic}")

    cit_paths: list[Path] = []
    for sample_id in sample_ids:
        bam = (
            results_root
            / sample_id
            / "03_alignment/04_tool_output/minimap2_samtools"
            / f"{sample_id}.primary.q20.bam"
        )
        if not file_nonempty(bam):
            raise RuntimeError(f"Pilot BAM is missing: {bam}")
        bam_link = work_root / f"{sample_id}.4sU.bam"
        if bam_link.is_symlink() and bam_link.resolve() != bam.resolve():
            raise RuntimeError(f"Pilot BAM link points elsewhere: {bam_link}")
        if not bam_link.exists():
            bam_link.symlink_to(bam)
        cit = Path(str(bam_link) + ".cit")
        if not file_nonempty(cit):
            run_logged(
                common + ["-mem", "12G", "-e", "Bam2CIT", str(cit), str(bam_link)],
                log_root / f"02_bam2cit_{sample_id}.log",
                env,
            )
            run_logged(
                common + ["-mem", "12G", "-e", "CorrectCIT", str(cit)],
                log_root / f"03_correct_cit_{sample_id}.log",
                env,
            )
        if not file_nonempty(cit):
            raise RuntimeError(f"CIT conversion failed for {sample_id}: {cit}")
        cit_paths.append(cit)

    combined_cit = work_root / "pilot.4sU.cit"
    if not file_nonempty(combined_cit):
        run_logged(
            # Keep per-sample CITs so a resumed pilot does not have to reconvert BAMs.
            common + ["-mem", "12G", "-e", "MergeCIT", str(combined_cit), *map(str, cit_paths)],
            log_root / "04_merge_cit.log",
            env,
        )
        run_logged(
            common + ["-mem", "12G", "-e", "ReadCount", str(combined_cit)],
            log_root / "05_read_count.log",
            env,
        )
    if not file_nonempty(combined_cit):
        raise RuntimeError(f"Merged CIT was not created: {combined_cit}")

    selected = json.loads((tmp_root / "strand_orientation/selected_mode.json").read_text(encoding="utf-8"))
    mode_map = {"forward": "Sense", "reverse": "Antisense", "unstranded": "Unspecific"}
    strandness = mode_map.get(str(selected.get("selected_mode", "")))
    if not strandness:
        raise RuntimeError("Grand3 pilot requires a passed project strand selection")

    prefix = output_root / "pilot"
    grand3_table = Path(str(prefix) + ".targets/data.tsv.gz")
    model_parameters = Path(str(prefix) + ".model.parameters.tsv")
    if not file_nonempty(grand3_table):
        grand3_cmd = common + [
            "-mem",
            "24G",
            "-t",
            str(tmp_root / "grand3_pilot/runtime_tmp"),
            "-e",
            "Grand3",
            "-auto",
            "-reads",
            str(combined_cit),
            "-genomic",
            str(genomic),
            "-prefix",
            str(prefix),
            "-nthreads",
            str(args.threads),
            "-strandness",
            strandness,
            "-targets",
            "genes",
            "-noplot",
            "-nosnps",
            "-output-table",
            "-D",
        ]
        run_logged(grand3_cmd, log_root / "06_grand3.log", env)
    if not file_nonempty(grand3_table) or not file_nonempty(model_parameters):
        raise RuntimeError(
            f"Grand3 output tables were not created: {grand3_table}, {model_parameters}"
        )

    comparison_path = table_root / "grand3_ntr_concordance.tsv"
    passed, comparisons = compare_sample_ntr(
        grand3_table, model_parameters, results_root, sample_ids, strandness, threshold, comparison_path
    )
    mean_absolute_error = sum(float(row["absolute_error"]) for row in comparisons) / len(comparisons)
    version_text = ""
    version_log = log_root / "00_grand3_version.log"
    version_cmd = common + ["-mem", "2G", "-e", "Grand3", "-h"]
    run_logged(version_cmd, version_log, env)
    version_text = version_log.read_text(encoding="utf-8", errors="replace")
    match = re.search(r"GRAND3 version\s+([^\s]+)", version_text, re.IGNORECASE)
    validation = {
        "schema_version": "iso-slam-grand3-pilot-v1",
        "status": "PASS" if passed else "FAIL",
        "grand3_version": match.group(1) if match else "",
        "gedi_executable": str(gedi),
        "genomic_index": str(genomic),
        "reads": str(combined_cit),
        "strandness": strandness,
        "sample_ids": sample_ids,
        "ntr_absolute_error_max": threshold,
        "ntr_mean_absolute_error": mean_absolute_error,
        "comparison_table": str(comparison_path),
        "grand3_table": str(grand3_table),
        "grand3_model_parameters": str(model_parameters),
        "comparison_basis": "Grand3 sample-level global Separate Binom ntr versus read-level hard-call fraction; gene-weighted target MAP is diagnostic only.",
        "note": "NC is treated as a labeled biological control, not as a no-4sU error-control library.",
    }
    (metadata_root / "grand3_pilot_validation.json").write_text(
        json.dumps(validation, indent=2) + "\n", encoding="utf-8"
    )
    if not passed:
        raise SystemExit(1)


if __name__ == "__main__":
    main()
