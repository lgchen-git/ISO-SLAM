#!/usr/bin/env bash
set -euo pipefail
export NUMBA_CACHE_DIR="${NUMBA_CACHE_DIR:-/runtime/numba}"
mkdir -p "$NUMBA_CACHE_DIR"
/opt/conda/bin/snakemake --version
python -c 'import pysam,numpy,scipy,yaml; print("Python imports PASS",pysam.__version__,numpy.__version__,scipy.__version__,yaml.__version__)'
minimap2 --version
samtools --version
isoquant --help >/runtime/isoquant_help.txt
sqanti3_qc.py --help >/runtime/sqanti3_help.txt
suppa.py --help >/runtime/suppa_help.txt
/opt/project/software/micromamba-root/envs/iso_slam/bin/python -c 'import nanoplot.NanoPlot; print("NanoPlot import PASS")'
/opt/project/software/micromamba-root/envs/r/bin/Rscript -e 'suppressPackageStartupMessages({library(edgeR);library(limma);library(data.table);library(ggplot2);library(jsonlite)}); cat("R imports PASS\n")'
java -version
/opt/project/software/gedi/ddcfa6444a9d61195e371c2ca9e49c193a816d0e/bin/gedi -h
