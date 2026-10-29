SHELL := /usr/bin/env bash
PYTHON_BIN ?= $(CURDIR)/.envs/iso_slam_v4/bin/python
CONFIG ?= config/run.yaml
export TMPDIR := $(CURDIR)/tmp/test_runs
export MPLCONFIGDIR := $(CURDIR)/tmp/test_runs/mplconfig
export XDG_CACHE_HOME := $(CURDIR)/tmp/test_runs/cache
export PYTHONPYCACHEPREFIX := $(CURDIR)/tmp/pycache
export PATH := /opt/project/software/micromamba-root/envs/r/bin:$(PATH)
.PHONY: lint validate smoke dry-run full docker fixture audit
lint:
	bash -n workflow/run.sh pipeline/docker/*.sh
	"$(PYTHON_BIN)" -m py_compile pipeline/scripts/*.py pipeline/scripts/00_shared/*.py pipeline/scripts/10_delivery/*.py pipeline/docker/*.py tests/*.py
validate:
	samples_tsv="$$("$(PYTHON_BIN)" -c 'import sys,yaml; print(yaml.safe_load(open(sys.argv[1]))["samples_tsv"])' "$(CONFIG)")"; \
	"$(PYTHON_BIN)" pipeline/scripts/check_samples.py --samples "$$samples_tsv" --analysis-type lr_bulk_slam
	"$(PYTHON_BIN)" pipeline/scripts/validate_run_config.py --config "$(CONFIG)"
smoke:
	"$(PYTHON_BIN)" tests/test_lr_bulk_slam_formal.py
	"$(PYTHON_BIN)" tests/test_lr_bulk_slam_assignment_mode.py
	"$(PYTHON_BIN)" tests/test_lr_bulk_slam_nanoplot_report.py
	"$(PYTHON_BIN)" tests/test_lr_bulk_slam_runner_helpers.py
	"$(PYTHON_BIN)" tests/test_lr_strand_audit.py
	"$(PYTHON_BIN)" tests/test_grand3_pilot.py
	"$(PYTHON_BIN)" pipeline/scripts/test_build_read_flow_sankey.py
	"$(PYTHON_BIN)" tests/test_lr_bulk_read_flow.py
	"$(PYTHON_BIN)" tests/test_sample_report_v2.py
	"$(PYTHON_BIN)" tests/test_project_delivery.py
	"$(PYTHON_BIN)" tests/test_iso_slam_v4.py
	"$(PYTHON_BIN)" tests/test_iso_slam_v4_report.py
	"$(PYTHON_BIN)" tests/test_iso_slam_delivery_package.py
	"$(PYTHON_BIN)" tests/test_iso_slam_as_delivery.py

	"$(PYTHON_BIN)" tests/test_production_groups.py
	"$(PYTHON_BIN)" tests/test_isoform_catalog_unique_reads.py

dry-run:
	bash workflow/run.sh "$(CONFIG)" --dry-run
full:
	bash workflow/run.sh "$(CONFIG)"
docker:
	bash pipeline/docker/build_export.sh
fixture:
	bash pipeline/docker/test_fixture.sh
audit:
	"$(PYTHON_BIN)" pipeline/scripts/00_shared/audit_snake_layout.py --config config/analysis_layout.yaml --project-root . --outputs-root outputs --json
