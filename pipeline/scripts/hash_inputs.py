#!/usr/bin/env python3
"""Write a deterministic SHA256 and size manifest for sample FASTQ inputs."""
import argparse, csv, hashlib
from pathlib import Path
p=argparse.ArgumentParser(); p.add_argument("--samples",required=True); p.add_argument("--output",required=True); a=p.parse_args()
with Path(a.samples).open(encoding="utf-8",newline="") as h: rows=list(csv.DictReader(h,delimiter="\t"))
target=Path(a.output)
inputs=[]; seen=set()
for row in rows:
    path=Path(row.get("fastq") or row["input_1"]).resolve()
    if path in seen: continue
    seen.add(path); inputs.append((row["sample_id"],row["group"],path))

# A completed manifest remains valid while every path/size matches and no input
# is newer. This makes a failed downstream resume avoid rereading ~170 GB.
if target.is_file() and target.stat().st_size > 0:
    with target.open(encoding="utf-8",newline="") as h:
        cached=list(csv.DictReader(h,delimiter="\t"))
    expected=[(sample,group,str(path),str(path.stat().st_size)) for sample,group,path in inputs]
    observed=[(row.get("sample_id",""),row.get("group",""),row.get("path",""),row.get("bytes","")) for row in cached]
    hashes_ok=all(len(row.get("sha256",""))==64 for row in cached)
    mtimes_ok=all(path.stat().st_mtime_ns <= target.stat().st_mtime_ns for _sample,_group,path in inputs)
    if observed == expected and hashes_ok and mtimes_ok:
        print(f"Input SHA256 manifest is current: {target}")
        raise SystemExit(0)

out=[]
for sample_id,group,path in inputs:
    digest=hashlib.sha256()
    with path.open("rb") as h:
        for block in iter(lambda:h.read(8*1024*1024),b""): digest.update(block)
    out.append([sample_id,group,str(path),path.stat().st_size,digest.hexdigest()])
target.parent.mkdir(parents=True,exist_ok=True)
with target.open("w",encoding="utf-8",newline="") as h:
    w=csv.writer(h,delimiter="\t",lineterminator="\n"); w.writerow(["sample_id","group","path","bytes","sha256"]); w.writerows(out)
