#!/usr/bin/env python3
"""Verify planned multi-group contrasts, direction and reference-only IDs."""
import csv
import json
import os
import sys
import subprocess
import tempfile
from pathlib import Path

ROOT=Path(__file__).resolve().parents[1]
def write(path,fields,rows):
    path.parent.mkdir(parents=True,exist_ok=True)
    with path.open('w') as handle:
        writer=csv.DictWriter(handle,fields,delimiter='\t',lineterminator='\n');writer.writeheader();writer.writerows(rows)

def main():
    with tempfile.TemporaryDirectory(dir=ROOT/'tmp/test_runs',prefix='two-groups-') as directory:
        root=Path(directory);results=root/'outputs';samples=[]
        genes=range(1,41)
        gtf=[]
        for gene in genes:
            attrs=f'gene_id "G{gene}"; gene_name "Gene{gene}"; gene_type "protein_coding";'
            gtf.append(f'chr1\ttest\tgene\t{gene*1000}\t{gene*1000+900}\t.\t+\t.\t{attrs}')
            for iso in [1,2]:gtf.append(f'chr1\ttest\ttranscript\t{gene*1000}\t{gene*1000+900}\t.\t+\t.\t{attrs} transcript_id "TX{gene}.{iso}";')
        import gzip
        reference=root/'genes.gtf.gz'
        with gzip.open(reference,'wt') as handle:handle.write('\n'.join(gtf)+'\n')
        groups=['DMSO','Dtag-short','Dtag-long','Dtag-FL','Dtag-PR']
        count_factor={'DMSO':1,'Dtag-short':2,'Dtag-long':3,'Dtag-FL':4,'Dtag-PR':5}
        ntr_value={'DMSO':.20,'Dtag-short':.35,'Dtag-long':.50,'Dtag-FL':.65,'Dtag-PR':.80}
        usage_value={'DMSO':.20,'Dtag-short':.32,'Dtag-long':.45,'Dtag-FL':.60,'Dtag-PR':.80}
        for group in groups:
            for rep in [1,2]:
                sid=f'{group}-{rep}';sample_root=results/sid;samples.append(dict(sample_id=sid,group=group,source_results_root=str(results),caller_results_root=str(results)))
                gene_dyn=[];gene_ntr=[];iso_dyn=[];catalog=[]
                for gene in genes:
                    count=(300+gene*13+rep*7)*count_factor[group]
                    ntr=ntr_value[group]+rep*.005
                    gene_dyn.append(dict(gene_id=f'G{gene}',assigned_reads=count))
                    gene_ntr.append(dict(gene_id=f'G{gene}',read_count=count,grand3_binom_ntr_map=ntr))
                    for iso in [1,2]:
                        iid=f'{sid}|local{gene}_{iso}'
                        usage=usage_value[group] if iso==1 else 1-usage_value[group]
                        iso_dyn.append(dict(isoform_id=iid,gene_id=f'G{gene}',total_reads=round(count*usage),expected_new_reads=count*usage*ntr,
                                            posterior_new_fraction=ntr,total_isoform_usage=usage,new_isoform_usage=usage))
                        catalog.append(dict(isoform_id=iid,associated_transcript=f'TX{gene}.{iso}',filter_status='PASS',structural_category='full-splice_match'))
                    # The same local novel ID in each sample is deliberately incomparable.
                    iso_dyn.append(dict(isoform_id=f'{sid}|novel_collision_{gene}',gene_id=f'G{gene}',total_reads=40,expected_new_reads=15,
                                        posterior_new_fraction=.375,total_isoform_usage=.05,new_isoform_usage=.05))
                    catalog.append(dict(isoform_id=f'{sid}|novel_collision_{gene}',associated_transcript=f'TX{gene}.1',filter_status='PASS',structural_category='novel_in_catalog'))
                for name,rows,folder in [('gene_dynamics',gene_dyn,'06_isoform_dynamics'),('isoform_dynamics',iso_dyn,'06_isoform_dynamics'),('grand3_gene_ntr',gene_ntr,'05_slam_model'),('isoform_catalog',catalog,'04_isoform_catalog')]:
                    write(sample_root/f'{folder}/01_tables/{sid}.{name}.tsv',list(rows[0]),rows)
                write(sample_root/f'05_slam_model/04_tool_output/caller/{sid}.substitution_summary.tsv', ['substitution','rate','source'],
                      [dict(substitution='T>C',rate=.02,source='observed'),dict(substitution='diagnostic_background',rate=.001,source='max_ta_tg_rate_background_site_filtered')])
                fastq_qc=sample_root/f'01_raw_qc/04_tool_output/longread_qc/fastq/{sid}.raw_fastq_qc.json'
                fastq_qc.parent.mkdir(parents=True,exist_ok=True)
                fastq_qc.write_text(json.dumps({'sample_id':sid,'read_count':1000,'q20_base_fraction':.85}))
                flow=sample_root/'07_sample_report/04_tool_output/read_flow/tables/read_flow_nodes.tsv'
                write(flow,['view_id','node_id','stage','label','value','unit','scope','node_kind'],[
                    dict(view_id='sample_reads',node_id=node,label=node,value=value,unit='reads',scope='sample',node_kind='count')
                    for node,value in [('raw',1000),('insert',900),('mapped',850),('mapq20',800)]])
                caller_qc=sample_root/f'05_slam_model/04_tool_output/caller/{sid}.bulk_slam_qc.json'
                caller_qc.parent.mkdir(parents=True,exist_ok=True)
                caller_qc.write_text(json.dumps({'sample_id':sid,'genes_detected':12000}))
        sheet=root/'samples.tsv';write(sheet,list(samples[0]),samples)
        dest=root/'group_analysis'
        subprocess.run(['/opt/project/software/micromamba-root/envs/r/bin/Rscript',str(ROOT/'pipeline/scripts/run_iso_slam_group_analysis.R'),
                        '--samples',str(sheet),'--gtf',str(reference),'--outdir',str(dest),'--label-control-group','DMSO',
                        '--contrasts-json',json.dumps([
                            {'treatment':'Dtag-short','control':'DMSO'},
                            {'treatment':'Dtag-long','control':'DMSO'},
                            {'treatment':'Dtag-FL','control':'DMSO'},
                            {'treatment':'Dtag-PR','control':'Dtag-FL'},
                            {'treatment':'Dtag-long','control':'Dtag-short'},
                        ])],check=True)
        metadata=json.loads((dest/'metadata/group_analysis_validation.json').read_text())
        expected_contrasts=['Dtag-short_vs_DMSO','Dtag-long_vs_DMSO','Dtag-FL_vs_DMSO','Dtag-PR_vs_Dtag-FL','Dtag-long_vs_Dtag-short']
        assert metadata['sample_count']==10 and metadata['groups']==groups and metadata['contrasts']==expected_contrasts,metadata
        assert metadata['biological_replicates_per_group']=={group:2 for group in groups},metadata
        assert metadata['label_control_group']=='DMSO'
        assert len(list(csv.DictReader((dest/'tables/contrast_summary.tsv').open(),delimiter='\t')))==25
        for contrast_name in expected_contrasts:
            contrast=dest/'tables/differential'/contrast_name
            with (contrast/'gene_labeling.tsv').open() as handle:labeling=list(csv.DictReader(handle,delimiter='\t'))
            assert labeling and all(float(row['effect'])>0 for row in labeling),contrast_name
        with (contrast/'isoform_abundance.tsv').open() as handle:isoforms=list(csv.DictReader(handle,delimiter='\t'))
        assert isoforms and all(row['feature_id'].startswith('TX') for row in isoforms)
        assert len({row['feature_id'] for row in isoforms})==80
        # The production HTML must be offline and validate with the same manifest contract.
        import shutil
        group_root=results/'_project/08_project_qc/02_group_analysis'
        shutil.copytree(dest,group_root)
        splicing_root=group_root/'tables/differential/new_rna_splicing'
        splicing_root.mkdir(parents=True,exist_ok=True)
        (splicing_root/'new_rna_splicing_validation.json').write_text(json.dumps({
            'status':'PASS','reference_event_count':1,
            'tested_event_count_by_contrast':{'Dtag-short_vs_DMSO':1},
        }))
        write(splicing_root/'differential_new_rna_splicing.tsv',
              ['contrast','event_id','gene_name','gene_id','mean_PSI_group_1','mean_PSI_group_2','delta_PSI','FDR','significant'],
              [dict(contrast='Dtag-short_vs_DMSO',event_id='SE:test',gene_name='Gene1',gene_id='G1',mean_PSI_group_1=.6,mean_PSI_group_2=.3,delta_PSI=.3,FDR=.02,significant='TRUE')])
        qc=[]
        for row in samples:
            qc.append(dict(sample_id=row['sample_id'],group_label_metadata_only=row['group'],raw_reads=1000,insert_reads=900,insert_rate=.9,isoform_assigned_reads=800,dynamics_eligible_reads=800,no_t_sites_reads=0,assignment_conflict_reads=0,assignment_conflict_rate=0,catalog_isoforms=80,filtered_isoforms=80,grand3_global_binom_ntr=.35,posterior_mean=.35,posterior_grand3_absolute_difference=0,switch_candidates=10))
        write(results/'_project/08_project_qc/01_tables/sample_qc_summary.tsv',list(qc[0]),qc)
        subprocess.run([sys.executable,str(ROOT/'pipeline/scripts/build_lr_iso_slam_project_report_v4.py'),'--results-root',str(results)],check=True,cwd=ROOT)
        review_path=results/'_project/09_project_report/01_tables/sample_qc_review.tsv'
        review_rows=list(csv.DictReader(review_path.open(),delimiter='\t'))
        assert len(review_rows)==10 and review_rows[0]['overall_review_status']=='BELOW_THRESHOLD'
        assert review_rows[0]['fastq_q20_base_fraction']=='0.85' and review_rows[0]['genes_detected']=='12000'
        report=results/'_project/09_project_report/report.html'
        subprocess.run([sys.executable,str(ROOT/'pipeline/scripts/00_shared/validate_bio_report.py'),'--report',str(report),'--manifest',str(report.parent/'05_metadata/delivery_manifest.tsv'),'--require-print-css','--json'],check=True)
        html=report.read_text()
        assert 'Gene1' in html and '基因名' in html
        assert all(name.replace('_vs_',' 减 ') in html for name in expected_contrasts)
        assert '五个实验组' not in html and '配置指定比较' in html and 'DMSO 是已标记的生物学对照' in html
        assert all(f'{group} n=2' in html for group in groups)
        assert 'data-pager-for="project-result-table"' in html and 'printing?d.rows' in html
        assert 'SQANTI3 未运行 ORF 预测' in html
        assert '逐样本审核指标' in html and 'FASTQ 碱基 Q20' in html and 'MAPQ≥20' in html and '检测基因数' in html
        # A valid gene analysis can have no comparable reference isoforms.
        for sample in samples:
            path=results/sample['sample_id']/f"04_isoform_catalog/01_tables/{sample['sample_id']}.isoform_catalog.tsv"
            with path.open() as handle:rows=list(csv.DictReader(handle,delimiter='\t'))
            for row in rows:row['structural_category']='novel_in_catalog'
            write(path,list(rows[0]),rows)
        empty=root/'group_empty_reference'
        subprocess.run(['/opt/project/software/micromamba-root/envs/r/bin/Rscript',str(ROOT/'pipeline/scripts/run_iso_slam_group_analysis.R'),
                        '--samples',str(sheet),'--gtf',str(reference),'--outdir',str(empty),'--label-control-group','DMSO',
                        '--contrasts-json',json.dumps([
                            {'treatment':'Dtag-short','control':'DMSO'},
                            {'treatment':'Dtag-long','control':'DMSO'},
                            {'treatment':'Dtag-FL','control':'DMSO'},
                            {'treatment':'Dtag-PR','control':'Dtag-FL'},
                            {'treatment':'Dtag-long','control':'Dtag-short'},
                        ])],check=True)
        with (empty/'tables/contrast_summary.tsv').open() as handle:summary=list(csv.DictReader(handle,delimiter='\t'))
        assert len(summary)==25 and all(int(row['tested_features'])==0 for row in summary if 'isoform' in row['layer'])
        legacy_sheet=root/'legacy_samples.tsv';write(legacy_sheet,list(samples[0]),samples[:4])
        legacy=root/'legacy_group_analysis'
        subprocess.run(['/opt/project/software/micromamba-root/envs/r/bin/Rscript',str(ROOT/'pipeline/scripts/run_iso_slam_group_analysis.R'),
                        '--samples',str(legacy_sheet),'--gtf',str(reference),'--outdir',str(legacy),
                        '--control-group','DMSO','--treatment-group','Dtag-short'],check=True)
        legacy_meta=json.loads((legacy/'metadata/group_analysis_validation.json').read_text())
        assert legacy_meta['contrasts']==['Dtag-short_vs_DMSO']
        print('Five planned contrasts, legacy two-group compatibility and empty transcript layers PASS')

if __name__=='__main__':main()
