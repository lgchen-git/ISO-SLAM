#!/usr/bin/env python3
"""Validate the run's inputs and reference without starting analysis."""
import argparse
import csv
from pathlib import Path
import yaml
parser=argparse.ArgumentParser(description=__doc__)
parser.add_argument('--config',required=True)
args=parser.parse_args()
root=Path.cwd().resolve()
cfg=yaml.safe_load(Path(args.config).read_text())
for key in ['results_root','tmp_root','log_root']:
    path=Path(cfg[key]).resolve()
    expected={'results_root':'outputs','tmp_root':'tmp','log_root':'log'}[key]
    if not path.is_relative_to(root/expected):raise SystemExit(f'{key} must be under {expected}/')
if cfg.get('delivery_root'):
    delivery_path=Path(cfg['delivery_root']).resolve()
    if not delivery_path.is_relative_to(root/'outputs'):raise SystemExit('delivery_root must be under outputs/')
for key, default in {
    'qc_review_min_raw_reads': 4_000_000,
    'qc_review_min_fastq_q20_fraction': 0.80,
    'qc_review_min_primary_mapq20_fraction': 0.90,
    'qc_review_min_detected_genes': 10_000,
}.items():
    value=float(cfg.get(key, default))
    if value <= 0 or (key.endswith('_fraction') and value > 1):
        raise SystemExit(f'Invalid QC review threshold {key}={value}')
with open(cfg['samples_tsv']) as handle:rows=list(csv.DictReader(handle,delimiter='\t'))
if not rows or len({row['sample_id'] for row in rows})!=len(rows):raise SystemExit('Empty or duplicate sample identifiers')
groups=[row.get('group','').strip() for row in rows]
if any(not group for group in groups):raise SystemExit('Every sample must have a group')
group_counts={group:groups.count(group) for group in dict.fromkeys(groups)}
if len(group_counts)<2 or any(count<2 for count in group_counts.values()):raise SystemExit('At least two samples per group are required')
if cfg.get('control_group') and cfg['control_group'] not in group_counts:raise SystemExit('control_group is absent from the sample sheet')
expected_replicates=cfg.get('biological_replicates_per_group')
if expected_replicates is not None and any(count!=int(expected_replicates) for count in group_counts.values()):
    raise SystemExit(f'Group replicate counts do not match biological_replicates_per_group={expected_replicates}: {group_counts}')
configured_contrasts=cfg.get('differential_contrasts')
if configured_contrasts is None:
    configured_contrasts=[{'treatment':cfg.get('treatment_group'),'control':cfg.get('control_group')}]
if not isinstance(configured_contrasts,list) or not configured_contrasts:
    raise SystemExit('differential_contrasts must be a non-empty list')
contrast_names=[]
for item in configured_contrasts:
    if not isinstance(item,dict):raise SystemExit('Each differential contrast must define treatment and control')
    treatment=item.get('treatment');control=item.get('control')
    if not treatment or not control or treatment==control:raise SystemExit('Each differential contrast must name two distinct groups')
    if treatment not in group_counts or control not in group_counts:raise SystemExit(f'Contrast refers to a group absent from the sample sheet: {treatment} vs {control}')
    contrast_names.append(f'{treatment}_vs_{control}')
if len(set(contrast_names))!=len(contrast_names):raise SystemExit('Differential contrast names must be unique')
for row in rows:
    path=Path(row['fastq'])
    if not path.is_file() or not path.stat().st_size:raise SystemExit(f'Missing input: {path}')
    if not path.resolve().is_relative_to(root/'inputs'):raise SystemExit('Input is outside inputs/')
for key in ['genome_fa','gtf','mmi','junc_bed']:
    if not cfg.get(key):continue
    path=Path(cfg[key])
    if not path.is_file() or not path.stat().st_size:raise SystemExit(f'Missing reference: {path}')
print(f'Run configuration OK: {len(rows)} samples')
