#!/usr/bin/env python3
"""Stage installed runtime environments only, preserving their software prefixes."""
import shutil
import sys
from pathlib import Path
target=Path(sys.argv[1]).resolve()
repo=Path(__file__).resolve().parents[2]
if not target.is_relative_to(repo/'tmp'):
    raise SystemExit('Build staging must be beneath project tmp/')
sources=[repo/'.envs/iso_slam_v4',Path('/opt/project/software/micromamba-root/envs/iso_slam'),
         Path('/opt/project/software/micromamba-root/envs/r'),Path('/opt/project/software/iso-slam-gedi-env'),
         Path('/opt/project/software/gedi/ddcfa6444a9d61195e371c2ca9e49c193a816d0e')]
for source in sources:
    if not source.is_dir():raise SystemExit(f'Missing runtime: {source}')
    dest=target/str(source).lstrip('/')
    print(f'Staging {source}',flush=True)
    shutil.copytree(source,dest,symlinks=True,dirs_exist_ok=True,ignore=shutil.ignore_patterns('__pycache__','.cache'))
