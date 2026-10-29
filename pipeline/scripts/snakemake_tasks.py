#!/usr/bin/env python3
"""Execute exactly one stage selected by the Snakemake DAG."""
from __future__ import annotations
import argparse
import csv
import json
import os
import shutil
import subprocess
import sys
from pathlib import Path
import yaml

HERE = Path(__file__).resolve().parent

def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--config', required=True)
    parser.add_argument('--task', required=True)
    parser.add_argument('--sample', default='')
    parser.add_argument('--threads', type=int, default=1)
    args = parser.parse_args()
    config = yaml.safe_load(Path(args.config).read_text())
    out, tmp, log = (Path(config[key]) for key in ['results_root', 'tmp_root', 'log_root'])
    with open(config['samples_tsv']) as handle:
        samples = list(csv.DictReader(handle, delimiter='\t'))
    sample = args.sample
    sample_root = out / sample
    work = tmp / (sample or '_project') / args.task
    work.mkdir(parents=True, exist_ok=True)
    for name in ['runtime_tmp', 'cache', 'mplconfig', 'home']:
        (work / name).mkdir(exist_ok=True)
    os.environ.update(TMPDIR=str((work/'runtime_tmp').resolve()), TMP=str((work/'runtime_tmp').resolve()),
                      TEMP=str((work/'runtime_tmp').resolve()), XDG_CACHE_HOME=str((work/'cache').resolve()),
                      MPLCONFIGDIR=str((work/'mplconfig').resolve()), HOME=str((work/'home').resolve()),
                      ISO_LOG_ROOT=str(log.resolve()))
    os.environ['PATH'] = os.pathsep.join([str(Path(config['python_bin']).parent), config['grand_slam_java_bin'], os.environ['PATH']])
    def run(name, *params):
        subprocess.run([sys.executable, str(HERE/name), *map(str, params)], check=True)
    def stage(name):
        p = sample_root/name
        p.mkdir(parents=True, exist_ok=True)
        return p
    bam = sample_root/f'03_alignment/04_tool_output/minimap2_samtools/{sample}.primary.q20.bam'
    insert = sample_root/f'02_trim/04_tool_output/trim/{sample}.insert.fastq.gz'
    caller = sample_root/'05_slam_model/04_tool_output/caller'
    def selected_config():
        selection = json.loads((tmp/'strand/selected_mode.json').read_text())
        if selection['status'] != 'PASS':
            raise RuntimeError('Strand selection failed')
        changed = dict(config, library_strand_mode=selection['selected_mode'])
        path = work/'selected.yaml'
        path.write_text(yaml.safe_dump(changed, allow_unicode=True))
        return path
    if args.task == 'raw_qc':
        row = next(row for row in samples if row['sample_id']==sample)
        dest = stage('01_raw_qc')/'04_tool_output/longread_qc/fastq'
        run('lr_longread_fastq_qc.py', '--fastq', row['fastq'], '--sample-id', sample, '--outdir', dest,
            '--label', 'raw', '--max-reads', config.get('longread_qc_max_reads', 200000))
        if config.get('enable_nanoplot_qc', True):
            nano = work/'nanoplot'
            subprocess.run([config['nanoplot_bin'], '--fastq', row['fastq'], '--outdir', str(nano),
                            '--threads', '1', '--downsample', '200000', '--no_static', '--tsv_stats'], check=True)
            for p in nano.glob('*Stats*'):
                shutil.copy2(p, dest/p.name)
    elif args.task == 'trim':
        row = next(row for row in samples if row['sample_id']==sample)
        run('lr_bulk_trim_reads.py', '--fastq', row['fastq'], '--sample-id', sample,
            '--outdir', stage('02_trim')/'04_tool_output/trim', '--primer-5p', config['primer_5p'],
            '--adapter-3p', config['adapter_3p'], '--seed-len', config.get('adapter_seed_len', 6),
            '--search-window', config.get('adapter_search_window', 400), '--min-insert-len', config.get('min_insert_len', 50),
            '--adapter-min-overlap', config.get('adapter_min_overlap', 12))
    elif args.task == 'alignment':
        params = ['--fastq', insert, '--reference', config.get('mmi') or config['genome_fa'],
                  '--outdir', stage('03_alignment')/'04_tool_output/minimap2_samtools', '--sample-id', sample,
                  '--threads', args.threads, '--minimap2-bin', config['minimap2_bin'], '--samtools-bin', config['samtools_bin']]
        if config.get('junc_bed'): params += ['--junc-bed', config['junc_bed']]
        run('lr_align_reads.py', *params)
    elif args.task == 'strand_audit':
        run('audit_lr_library_strand.py', '--bam', bam, '--gtf', config['gtf'], '--sample-id', sample,
            '--outdir', tmp/sample/'strand')
    elif args.task == 'strand_selection':
        audits = [json.loads((tmp/row['sample_id']/'strand'/f"{row['sample_id']}.strand_orientation_audit.json").read_text()) for row in samples]
        modes = {audit['selected_mode'] for audit in audits}
        passed = all(audit['status']=='PASS' for audit in audits) and len(modes)==1
        if not passed: raise RuntimeError('Samples disagree on library orientation; inspect strand audits')
        dest = tmp/'strand/selected_mode.json';dest.parent.mkdir(parents=True, exist_ok=True)
        dest.write_text(json.dumps({'status':'PASS','selected_mode':next(iter(modes)), 'sample_count':len(audits)}, indent=2)+'\n')
    elif args.task == 'grand3_reference':
        index = Path(config['grand_slam_genomic']);index.parent.mkdir(parents=True, exist_ok=True)
        subprocess.run([config['grand_slam_bin'], '-Duser.home='+str((work/'home').resolve()), '-mem', '16G', '-t', str(work/'runtime_tmp'),
                        '-e', 'IndexGenome', '-s', str(Path(config['genome_fa']).resolve()),
                        '-a', str(Path(config['gtf']).resolve()), '-f', str(index.parent.resolve()),
                        '-n', index.stem, '-o', str(index.resolve()), '-nomapping', '-p', '-D'], check=True)
    elif args.task == 'read_calling':
        selected = json.loads((tmp/'strand/selected_mode.json').read_text())['selected_mode']
        params = ['--bam', bam, '--fasta', config['genome_fa'], '--gtf', config['gtf'], '--sample-id', sample,
                  '--outdir', caller, '--library-strand-mode', selected, '--assignment-mode', config['assignment_mode'],
                  '--main-call', config['main_call'], '--site-evidence-mode', 'supportive',
                  '--assignment-end-distance-tiebreak', '--assignment-allow-gene-level-salvage',
                  '--assignment-count-shared-resolved-reads', '--low-t-zero-tc-as-old', '--low-t-high-conf-tc-as-new',
                  '--bg-from-background-sites']
        for key in ['min_mapq','min_tc_baseq','min_t_covered','min_site_depth','min_site_tc_count','min_site_tc_rate',
                    'max_site_other_mismatch_rate','max_read_other_mismatch_rate','alpha']:
            params += ['--'+key.replace('_','-'), config[key]]
        run('lr_call_bulk_slam.py', *params)
    elif args.task == 'isoform_catalog':
        run('run_isoform_catalog.py', '--config', args.config, '--sample-id', sample, '--bam', bam,
            '--fastq', insert, '--fasta', config['genome_fa'], '--gtf', config['gtf'],
            '--outdir', stage('04_isoform_catalog'), '--tmp-dir', work, '--threads', args.threads)
    elif args.task == 'slam_model':
        run('run_grand3_sample.py', '--config', selected_config(), '--sample-id', sample, '--bam', bam,
            '--outdir', stage('05_slam_model'), '--tmp-dir', work, '--threads', args.threads)
    elif args.task == 'isoform_dynamics':
        catalog = sample_root/'04_isoform_catalog/01_tables'; model = sample_root/'05_slam_model/01_tables'
        run('build_isoform_dynamics.py', '--sample-id', sample, '--read-calls', caller/f'{sample}.bulk_read_calls.tsv.gz',
            '--read-assignments', catalog/f'{sample}.read_isoform_assignments.tsv.gz', '--catalog', catalog/f'{sample}.isoform_catalog.tsv',
            '--grand3-global', model/f'{sample}.grand3_global_parameters.tsv', '--grand3-genes', model/f'{sample}.grand3_gene_ntr.tsv',
            '--outdir', stage('06_isoform_dynamics'), '--new-threshold', config['posterior_new_threshold'],
            '--old-threshold', config['posterior_old_threshold'], '--bootstrap-iterations', config['posterior_bootstrap_iterations'],
            '--bootstrap-seed', config['posterior_bootstrap_seed'], '--concordance-threshold', config['grand_slam_ntr_mae_max'])
    elif args.task == 'sample_report':
        dest = stage('07_sample_report')
        run('build_lr_iso_slam_read_flow.py', '--sample-id', sample, '--sample-root', sample_root,
            '--outdir', dest/'04_tool_output/read_flow', '--plotly-js', 'pipeline/assets/plotly-3.5.1.min.js')
        run('build_lr_iso_slam_report_v4.py', '--sample-id', sample, '--sample-root', sample_root,
            '--v3-caller', caller, '--report-title', config['report_title'], '--gtf', config['gtf'], '--config', args.config)
        validate(dest/f'{sample}_report.html', dest/'05_metadata/delivery_manifest.tsv', dest/'05_metadata/report_validation.json')
    elif args.task == 'project_qc':
        run('build_lr_iso_slam_project_qc_v4.py', '--results-root', out, '--samples', config['samples_tsv'])
    elif args.task == 'group_analysis':
        command = [config['rscript_bin'], str(HERE/'run_iso_slam_group_analysis.R'), '--samples', config['samples_tsv'],
                   '--gtf', config['gtf'], '--outdir', str(out/'_project/08_project_qc/02_group_analysis')]
        if config.get('differential_contrasts'):
            command.extend(['--contrasts-json', json.dumps(config['differential_contrasts'], ensure_ascii=False)])
        else:
            command.extend(['--control-group', config['control_group'], '--treatment-group', config['treatment_group']])
        if config.get('control_group'):
            command.extend(['--label-control-group', config['control_group']])
        subprocess.run(command, check=True)
    elif args.task == 'new_rna_splicing':
        contrasts = config.get('differential_contrasts') or [
            {'treatment': config['treatment_group'], 'control': config['control_group']}
        ]
        run('run_new_rna_splicing.py', '--samples', config['samples_tsv'], '--gtf', config['gtf'],
            '--results-root', out, '--outdir', out/'_project/08_project_qc/02_group_analysis/tables/differential/new_rna_splicing',
            '--workdir', tmp/'_project/new_rna_splicing',
            '--suppa', config['suppa_bin'], '--rscript', config['rscript_bin'],
            '--contrasts-json', json.dumps(contrasts, ensure_ascii=False))
    elif args.task == 'project_report':
        run('build_lr_iso_slam_project_report_v4.py', '--results-root', out, '--config', args.config)
        dest = out/'_project/09_project_report'
        validate(dest/'report.html', dest/'05_metadata/delivery_manifest.tsv', dest/'05_metadata/report_validation.json', require_expandable=False)
    elif args.task == 'customer_delivery':
        package_config = dict(config)
        package_config.setdefault('delivery_root', str(out/'_project/10_delivery/ISO_SLAM_Delivery'))
        package_config.setdefault('delivery_name', 'ISO-SLAM 长读长 RNA 标记分析结果交付包')
        package_config.setdefault('include_sample_reports', False)
        package_config.setdefault('gene_annotation_gtf', config.get('gtf', ''))
        package_config['samples_tsv'] = str(Path(config['samples_tsv']).resolve())
        package_config_path = work/'customer_delivery_config.yaml'
        package_config_path.write_text(yaml.safe_dump(package_config, allow_unicode=True, sort_keys=False), encoding='utf-8')
        run('10_delivery/build_iso_slam_delivery_package.py', '--config', package_config_path, '--replace')
    elif args.task == 'delivery':
        complete = out/'_project/RUN_COMPLETE.tsv'
        if complete.exists(): complete.unlink()
        artifacts = []
        for row in samples:
            sid=row['sample_id'];root=out/sid
            for folder in ['04_isoform_catalog/01_tables','05_slam_model/01_tables','06_isoform_dynamics/01_tables']:
                artifacts.extend(sorted((root/folder).glob('*.tsv')))
            artifacts.append(root/f'07_sample_report/{sid}_report.html')
        artifacts.extend(sorted((out/'_project/08_project_qc/02_group_analysis/tables').rglob('*.tsv')))
        report=out/'_project/09_project_report/report.html'
        artifacts.append(report)
        manifest=out/'delivery_manifest.tsv'
        with manifest.open('w') as handle:
            writer=csv.writer(handle,delimiter='\t',lineterminator='\n');writer.writerow(['artifact','description','target','status'])
            for p in artifacts:
                if not p.is_file() or not p.stat().st_size: raise RuntimeError(f'Missing deliverable: {p}')
                writer.writerow(['final_report' if p==report else ('sample_report' if p.suffix=='.html' else p.stem),p.name,str(p.resolve()),'up_to_date'])
        validate(report, manifest, out/'delivery_validation.json', require_expandable=False)
        shutil.copy2(args.config,out/'run_config.yaml')
        shutil.copy2(config['samples_tsv'],out/'samples.tsv')
        if not (out/'run_launch_metadata.json').is_file():
            raise RuntimeError('Host launch metadata is missing')
        group_validation = json.loads((out/'_project/08_project_qc/02_group_analysis/metadata/group_analysis_validation.json').read_text())
        contrasts = group_validation.get('contrasts', [])
        if group_validation.get('status') != 'PASS' or not contrasts:
            raise RuntimeError('Validated group contrasts are required before marking delivery complete')
        complete.write_text('status\tsample_count\treport_count\tcontrast\treport_schema\n'+
                            f"completed\t{len(samples)}\t{len(samples)+1}\t{';'.join(contrasts)}\t{config['report_schema']}\n")
    else:
        raise SystemExit('Unknown Snakemake stage: '+args.task)

def validate(report, manifest, output, require_expandable=True):
    command = [sys.executable,str(HERE/'00_shared/validate_bio_report.py'), '--report',str(report),
                             '--manifest',str(manifest),'--require-print-css',
                             '--forbid-text','single_sample_report','--json']
    if require_expandable:
        command.insert(command.index('--require-print-css'), '--require-expandable-figures')
    result = subprocess.run(command,text=True,capture_output=True)
    output.parent.mkdir(parents=True,exist_ok=True);output.write_text(result.stdout)
    if result.returncode: raise RuntimeError(result.stdout+result.stderr)

if __name__=='__main__':
    main()
