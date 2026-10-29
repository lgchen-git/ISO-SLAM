#!/usr/bin/env python3
"""Run Grand3 for one labelled sample and normalize gene/global Binom NTR outputs."""

from __future__ import annotations

import argparse
import csv
import gzip
import os
import re
from pathlib import Path

from iso_slam_v4_common import ensure_stage, executable, nonempty, read_config, run_logged, write_json, write_tsv


SCHEMA = "iso-slam-grand3-sample-v4"


def arguments() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--config", required=True)
    parser.add_argument("--sample-id", required=True)
    parser.add_argument("--bam", required=True, type=Path)
    parser.add_argument("--outdir", required=True, type=Path)
    parser.add_argument("--tmp-dir", required=True, type=Path)
    parser.add_argument("--threads", type=int, default=8)
    parser.add_argument("--skip-tools", action="store_true")
    parser.add_argument("--grand3-table", type=Path)
    parser.add_argument("--model-parameters", type=Path)
    return parser.parse_args()


def strandness(config: dict[str, str]) -> str:
    mapping = {"forward": "Sense", "reverse": "Antisense", "unstranded": "Unspecific"}
    value = mapping.get(config.get("library_strand_mode", ""))
    if not value:
        raise RuntimeError("V4 requires a fixed forward, reverse, or unstranded library_strand_mode")
    return value


def read_global(path: Path, sample_id: str, strand: str) -> dict[str, object]:
    candidates: list[dict[str, str]] = []
    with path.open(encoding="utf-8", newline="") as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            if row.get("Estimator") != "Separate" or row.get("Label") != "4sU" or row.get("Subread") != strand:
                continue
            if Path(row.get("Condition", "")).name != f"{sample_id}.4sU":
                continue
            candidates.append(row)
    if len(candidates) != 1:
        raise RuntimeError(f"Expected one Grand3 Separate Binom model for {sample_id}, found {len(candidates)}")
    row = candidates[0]
    required = ("Binom ntr", "Binom p.err", "Binom p.conv")
    missing = [key for key in required if row.get(key, "") in {"", "NA", "NaN"}]
    if missing:
        raise RuntimeError(f"Grand3 model lacks required parameters: {missing}")
    return {
        "sample_id": sample_id,
        "estimator": "Separate",
        "label": "4sU",
        "strandness": strand,
        "binom_ntr": float(row["Binom ntr"]),
        "binom_p_err": float(row["Binom p.err"]),
        "binom_p_conv": float(row["Binom p.conv"]),
    }


def select(header: list[str], sample_id: str, suffix: str) -> str:
    matches = [name for name in header if sample_id in name and name.endswith(suffix)]
    if len(matches) != 1:
        raise RuntimeError(f"Could not select Grand3 column {suffix!r}: {matches}")
    return matches[0]


def normalize_genes(path: Path, sample_id: str) -> list[dict[str, object]]:
    opener = gzip.open if path.suffix == ".gz" else open
    with opener(path, "rt", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        header = reader.fieldnames or []
        count_col = select(header, sample_id, "Read count")
        ntr_col = select(header, sample_id, "4sU Binom NTR MAP")
        rows = []
        for row in reader:
            try:
                count = float(row[count_col])
                ntr = float(row[ntr_col])
            except (TypeError, ValueError):
                continue
            gene = row.get("Gene") or row.get("gene") or row.get("Name") or ""
            if gene and count > 0:
                rows.append({"gene_id": gene, "read_count": count, "grand3_binom_ntr_map": ntr})
    if not rows:
        raise RuntimeError(f"No usable Grand3 gene NTR rows in {path}")
    return rows


def main() -> int:
    args = arguments()
    config = read_config(args.config)
    ensure_stage(args.outdir)
    strand = strandness(config)
    output_root = args.outdir / "04_tool_output/grand3"
    output_root.mkdir(parents=True, exist_ok=True)
    args.tmp_dir.mkdir(parents=True, exist_ok=True)
    prefix = output_root / args.sample_id
    table = args.grand3_table or Path(str(prefix) + ".targets/data.tsv.gz")
    models = args.model_parameters or Path(str(prefix) + ".model.parameters.tsv")

    if not args.skip_tools and (not nonempty(table) or not nonempty(models)):
        gedi = executable(config["grand_slam_bin"])
        genomic = Path(config["grand_slam_genomic"])
        if not nonempty(genomic):
            raise RuntimeError(f"Grand3 genomic index is missing: {genomic}")
        bam_link = args.tmp_dir / f"{args.sample_id}.4sU.bam"
        if bam_link.is_symlink() and bam_link.resolve() != args.bam.resolve():
            raise RuntimeError(f"Existing Grand3 BAM link points elsewhere: {bam_link}")
        if not bam_link.exists():
            bam_link.symlink_to(args.bam.resolve())
        bai_link = Path(str(bam_link) + ".bai")
        if not bai_link.exists():
            bai_link.symlink_to(Path(str(args.bam.resolve()) + ".bai"))
        java_bin = config.get("grand_slam_java_bin", "")
        env_path = os.pathsep.join(part for part in (java_bin, str(Path(gedi).parent), os.environ.get("PATH", "")) if part)
        for directory in (args.tmp_dir / "gedi_home/.gedi", args.tmp_dir / "runtime_tmp"):
            directory.mkdir(parents=True, exist_ok=True)
        command_prefix = [gedi, "-Duser.home=" + str(args.tmp_dir / "gedi_home")]
        cit = Path(str(bam_link) + ".cit")
        if not nonempty(cit):
            env = os.environ.copy()
            env["PATH"] = env_path
            os.environ.update({"PATH": env_path})
            run_logged(command_prefix + ["-mem", "12G", "-e", "Bam2CIT", str(cit), str(bam_link)], args.outdir / f"03_logs/{args.sample_id}.bam2cit.log")
            run_logged(command_prefix + ["-mem", "12G", "-e", "CorrectCIT", str(cit)], args.outdir / f"03_logs/{args.sample_id}.correctcit.log")
        command = command_prefix + [
            "-mem", "24G", "-t", str(args.tmp_dir / "runtime_tmp"), "-e", "Grand3", "-auto",
            "-reads", str(cit), "-genomic", str(genomic), "-prefix", str(prefix),
            "-nthreads", str(args.threads), "-strandness", strand, "-targets", "genes",
            "-noplot", "-output-table", "-D",
        ]
        # Deliberately omit -nosnps: V4 retains Grand3's SNP/error modelling.
        run_logged(command, args.outdir / f"03_logs/{args.sample_id}.grand3.log")
    if not nonempty(table) or not nonempty(models):
        raise RuntimeError(f"Grand3 outputs are missing: {table}, {models}")

    global_row = read_global(models, args.sample_id, strand)
    genes = normalize_genes(table, args.sample_id)
    write_tsv(args.outdir / f"01_tables/{args.sample_id}.grand3_global_parameters.tsv", [global_row])
    write_tsv(args.outdir / f"01_tables/{args.sample_id}.grand3_gene_ntr.tsv", genes)
    valid = (
        0 <= float(global_row["binom_ntr"]) <= 1
        and 0 <= float(global_row["binom_p_err"]) < 1
        and 0 < float(global_row["binom_p_conv"]) <= 1
    )
    metadata = {
        "schema_version": SCHEMA,
        "status": "PASS" if valid else "FAIL",
        "sample_id": args.sample_id,
        "grand3_version": config.get("grand_slam_version", ""),
        "targets": "genes",
        "snp_modeling": "enabled",
        "nosnps_flag_used": False,
        "gene_rows": len(genes),
        "global_parameters": global_row,
        "source_table": str(table),
        "source_model_parameters": str(models),
    }
    write_json(args.outdir / "05_metadata/grand3_validation.json", metadata)
    if not valid:
        raise RuntimeError("Grand3 parameters are outside valid probability ranges")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
