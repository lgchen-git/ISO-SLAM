#!/usr/bin/env python3
"""Record host configuration and Git state before the Docker runner starts."""
import argparse
import datetime
import json
import subprocess
from pathlib import Path
import yaml
parser=argparse.ArgumentParser(description=__doc__)
parser.add_argument('--config',required=True)
args=parser.parse_args()
cfg=yaml.safe_load(Path(args.config).read_text())
out=Path(cfg['results_root']);out.mkdir(parents=True,exist_ok=True)
state=subprocess.run(['git','status','--short'],text=True,capture_output=True,check=True)
(out/'git_state.txt').write_text(state.stdout)
metadata={'started_utc':datetime.datetime.now(datetime.timezone.utc).isoformat(),'image':'iso-slam:20260930-v1.0.0',
          'config':args.config,'cores':cfg.get('cores',32),'memory_mb':cfg.get('memory_mb',196608),
          'reference':cfg.get('reference_label', str(Path(cfg['genome_fa']).parent)),
          'label_control_group':cfg.get('control_group'),
          'differential_contrasts':cfg.get('differential_contrasts') or [{'treatment':cfg['treatment_group'],'control':cfg['control_group']}],
          'biological_replicates_per_group':cfg['biological_replicates_per_group'],'labeling_design':cfg['labeling_design'],
          'analysis_scope':'planned_group_contrasts','git_dirty':bool(state.stdout)}
software={}
for key in ['python_bin','rscript_bin','nanoplot_bin','grand_slam_java_bin']:
    path=Path(cfg[key]);prefix=path.parent.parent if key!='grand_slam_java_bin' else path.parent
    for package in (prefix/'conda-meta').glob('*.json'):
        entry=json.loads(package.read_text())
        name=entry.get('name')
        if name in {'python','isoquant','sqanti3','suppa','samtools','minimap2','pysam','scipy','nanoplot','openjdk','bioconductor-edger','bioconductor-limma','r-base','r-ggplot2','r-data.table','r-jsonlite'}:
            software[key+':'+name]=entry.get('version')
metadata['software_versions']=software
(out/'run_launch_metadata.json').write_text(json.dumps(metadata,ensure_ascii=False,indent=2)+'\n')
