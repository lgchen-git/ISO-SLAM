#!/usr/bin/env python3
"""Create deterministic synthetic cDNA/4sU libraries in an isolated project."""
import csv
import os
import random
import shutil
import sys
from pathlib import Path
import yaml

def main():
    repo=Path(__file__).resolve().parents[2]
    fixture_name=os.environ.get('ISO_SLAM_FIXTURE_NAME','docker_fixture_20260930')
    if not fixture_name.startswith('docker_fixture_') or Path(fixture_name).name!=fixture_name:
        raise SystemExit('ISO_SLAM_FIXTURE_NAME must be a child name beginning with docker_fixture_')
    root=(repo/'tmp'/fixture_name).resolve()
    if root.exists():
        raise SystemExit(f'Refusing to overwrite existing fixture: {root}')
    root.mkdir(parents=True,exist_ok=True)
    for name in ['config','inputs','log','outputs','reference','tmp']:(root/name).mkdir(exist_ok=True)
    shutil.copytree(repo/'pipeline',root/'pipeline',ignore=shutil.ignore_patterns('__pycache__'))
    rng=random.Random(20260930)
    genome=[];gtf=[];transcripts=[];position=1
    for gene in range(1,41):
        chunk=''.join(rng.choices('ACGT',k=3200))
        # Two canonical splice structures share terminal exons.
        bounds=[(position,position+299),(position+900,position+1099),(position+2100,position+2399)]
        chars=list(chunk)
        for _left,right in bounds[:-1]:chars[right-position+1:right-position+3]='GT'
        for left,_right in bounds[1:]:chars[left-position-2:left-position]='AG'
        chunk=''.join(chars);genome.append(chunk)
        gid=f'G{gene}'; attrs=f'gene_id "{gid}"; gene_name "Gene{gene}"; gene_type "protein_coding";'
        gtf.append(f'chrFixture\tfixture\tgene\t{position}\t{position+2399}\t.\t+\t.\t{attrs}')
        for iso,exons in [(1,bounds),(2,[bounds[0],bounds[2]])]:
            tid=f'TX{gene}.{iso}';ta=attrs+f' transcript_id "{tid}"; transcript_type "protein_coding";'
            gtf.append(f'chrFixture\tfixture\ttranscript\t{exons[0][0]}\t{exons[-1][1]}\t.\t+\t.\t{ta}')
            for number,(left,right) in enumerate(exons,1):gtf.append(f'chrFixture\tfixture\texon\t{left}\t{right}\t.\t+\t.\t{ta} exon_number "{number}";')
            seq=''.join(chunk[left-position:right-position+1] for left,right in exons)
            transcripts.append((gene,iso,seq))
        position+=len(chunk)
    ref=root/'reference/fixture';ref.mkdir(exist_ok=True)
    sequence=''.join(genome)
    (ref/'genome.fa').write_text('>chrFixture\n'+'\n'.join(sequence[i:i+80] for i in range(0,len(sequence),80))+'\n')
    import pysam
    pysam.faidx(str(ref/'genome.fa'))
    import gzip
    with gzip.open(ref/'genes.gtf.gz','wt') as handle:handle.write('\n'.join(gtf)+'\n')
    cfg=yaml.safe_load((repo/'config/example_sample.yaml').read_text())
    cfg.update(project_id='synthetic_fixture',results_root='outputs/synthetic_fixture',tmp_root='tmp/synthetic_fixture',log_root='log/synthetic_fixture',
               samples_tsv='config/samples.tsv',genome_fa='reference/fixture/genome.fa',gtf='reference/fixture/genes.gtf.gz',mmi='',junc_bed='',
               grand_slam_genomic='tmp/synthetic_fixture/reference/fixture.oml',threads=4,cores=16,memory_mb=131072,
               enable_nanoplot_qc=False,longread_qc_max_reads=10000,posterior_bootstrap_iterations=20,
               report_title='合成长读长 RNA 标记流程验证报告')
    (root/'config/fixture.yaml').write_text(yaml.safe_dump(cfg,sort_keys=False,allow_unicode=True))
    rows=[];complement=str.maketrans('ACGT','TGCA')
    for group in ['DMSO','Dtag-short']:
        for replicate in [1,2]:
            sid=f'{group}-{replicate}';dest=root/'inputs'/f'{sid}.fastq'
            ntr=(0.35 if group=='DMSO' else 0.65)+(replicate-1)*0.015
            with dest.open('w') as handle:
                for gene,iso,seq in transcripts:
                    n=400+(gene%7)*60+replicate*7
                    if group=='Dtag-short' and gene%4==0:n=int(n*(1.8 if iso==1 else 0.6))
                    for read in range(n):
                        labeled=rng.random()<ntr;bases=list(seq)
                        for pos,base in enumerate(bases):
                            if base=='T' and labeled and rng.random()<0.055:bases[pos]='C'
                            elif rng.random()<0.001:bases[pos]=rng.choice([x for x in 'ACGT' if x!=base])
                        cdna='T'*20+''.join(bases).translate(complement)[::-1]
                        raw=cfg['primer_5p']+cdna+cfg['adapter_3p']
                        handle.write(f'@{sid}_g{gene}_i{iso}_r{read}\n{raw}\n+\n'+('I'*len(raw))+'\n')
            rows.append([sid,group,f'inputs/{sid}.fastq',f'inputs/{sid}.fastq',f'inputs/{sid}.fastq','human',cfg['genome_fa'],cfg['gtf'],4,cfg['results_root'],cfg['results_root']])
    with (root/'config/samples.tsv').open('w') as handle:
        writer=csv.writer(handle,delimiter='\t',lineterminator='\n')
        writer.writerow(['sample_id','group','input_1','input_2','fastq','species','genome_fa','gtf','threads','source_results_root','caller_results_root']);writer.writerows(rows)
    print(root)

if __name__=='__main__':main()
