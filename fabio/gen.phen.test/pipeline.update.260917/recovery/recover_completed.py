#!/usr/bin/env python3
"""Offline assembly of one stopped run. Never invokes Nextflow or CAAStools."""
import argparse
from collections import Counter
import re
import shlex
import shutil
import subprocess
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT/'scripts'))
from common import checked, identifier, read_json, read_tsv, require, sha, write_json, write_tsv
from summarize_hypothesis import summarize
from common_background import build
from functional import run as functional

TAG = re.compile(r'CAAS_POOLED \(([^/()]+)/([^:()]+):([^()]+)\)')
HASH = re.compile(r'\[([0-9a-f]{2}/[0-9a-f]+)\]')


def parse_log(path):
    """Only explicit terminal records count; scheduler polling ERROR is not one."""
    with Path(path).open() as handle:
        return parse_lines(handle)


def parse_lines(lines):
    outcomes = {}
    params_path = None
    for line in lines:
        if '$> nextflow ' in line:
            require(params_path is None, 'Multiple launches in source log')
            args = shlex.split(line.split('$> ', 1)[1])
            require('-params-file' in args, 'Log has no explicit launch parameters')
            params_path = Path(args[args.index('-params-file')+1])
        tag = TAG.search(line)
        if not tag: continue
        strategy, replicate, gene = tag.groups()
        key = (strategy+'_'+replicate, gene)
        row = outcomes.setdefault(key, dict(hypothesis_id=key[0], strategy_id=strategy,
                                  replicate_id=replicate, gene=gene, status='pending', work_dir='', task_hash=''))
        if 'Cached process >' in line or 'Submitted process >' in line or 'Error is ignored' in line:
            match = HASH.search(line)
            if match: row['task_hash'] = match[1]
        if 'Cached process >' in line:
            row['status'] = 'cached'
        elif 'Task completed > TaskHandler' in line:
            exit_status = re.search(r'; exit: ([^;]+);', line)
            require(exit_status is not None, 'Unrecognized task completion record')
            row['status'] = 'completed' if exit_status[1] == '0' else 'failed'
            work = re.search(r'workDir: (.*?) started:', line)
            if work: row['work_dir'] = work[1]
        elif 'Error is ignored' in line:
            row['status'] = 'submission_failed' if 'Error submitting process' in line else 'failed'
    require(params_path is not None, 'No launch header found in log')
    return params_path, outcomes


def validate_receipt(path, run_id, row, alignment, lock):
    receipt = read_json(path)
    expected = dict(run_id=run_id, hypothesis_id=row['hypothesis_id'], gene_id=alignment['gene_id'],
                    strategy_id=row['strategy_id'], replicate_id=row['replicate_id'], status='success',
                    alignment_sha256=alignment['sha256'], config_sha256=row['config_sha256'],
                    pool_sha256=row['pool_sha256'], settings_sha256=lock['settings_sha256'],
                    design_sha256=lock['design_sha256'], selection_sha256=lock['selection_sha256'])
    for field, value in expected.items():
        require(receipt.get(field) == value, f'Receipt provenance mismatch {field}: {path}')
    gene = alignment['gene_id']
    checked(path.parent/(gene+'.pooled.caas.events.tsv'),receipt['events_sha256'])
    checked(path.parent/(gene+'.pooled.caas.tsv'),receipt['legacy_sha256'])
    return receipt


def resolve_receipt(root, lock, run_id, row, gene, outcome, cache):
    published = root/'results'/run_id/row['strategy_id']/row['replicate_id']/'raw'/gene/(gene+'.receipt.json')
    if published.is_file(): return published
    work_root = Path(lock['work_dir']).resolve()
    if outcome['work_dir']:
        task = Path(outcome['work_dir']).resolve()
        require(task.is_relative_to(work_root), 'Task work directory outside original work root')
    else:
        require(re.fullmatch(r'[0-9a-f]{2}/[0-9a-f]+',outcome['task_hash']) is not None,'Missing safe cache hash')
        prefix, suffix = outcome['task_hash'].split('/')
        if prefix not in cache:
            cache[prefix] = list((work_root/prefix).iterdir())
        matches = [p for p in cache[prefix] if p.name.startswith(suffix)]
        require(len(matches)==1, f'Ambiguous/missing cached task directory: {outcome["task_hash"]}')
        task = matches[0]
    result = task/gene/(gene+'.receipt.json')
    require(result.is_file(), f'Log reports success but receipt missing: {result}')
    return result


def recover(root, run_id, recovery_id, log, rscript='Rscript', tables_only=False):
    identifier(run_id); identifier(recovery_id)
    root=root.resolve(); log=log.resolve()
    original_log_hash = sha(log)
    params_path, outcomes = parse_log(log)
    params = read_json(params_path)
    require(params['run_id']==run_id, 'Source log belongs to a different run')
    require(params.get('direct_discovery') is True, 'Recovery is limited to direct benchmark runs')
    lock = read_json(root/'results'/run_id/'launch.lock.json')
    require(Path(lock['results_root']).resolve()==(root/'results').resolve(), 'Nondefault results root: explicit review required')
    selection = checked(Path(params['selection_manifest']),lock['selection_sha256'])
    inventory = checked(Path(params['alignment_manifest']),lock['alignment_manifest_sha256'])
    settings = checked(Path(params['settings']),lock['settings_sha256'])
    design = Path(params['design'])
    checked(design/'design.lock.json', lock['design_sha256'])
    design_lock = read_json(design/'design.lock.json')
    rows, genes = read_tsv(selection), read_tsv(inventory)
    expected={(r['hypothesis_id'],g['gene_id']) for r in rows for g in genes}
    require(set(outcomes)==expected, 'Log does not account for exactly the selected gene/hypothesis tasks')
    require(all(r['status']!='pending' for r in outcomes.values()), 'Source log still has pending/unresolved tasks')
    output=root/'summaries'/run_id/'recovery'/recovery_id
    require(not output.exists(), 'Recovery ID already exists; choose a new ID (originals are never overwritten)')
    output.mkdir(parents=True)
    write_json(output/'recovery.started.json',dict(run_id=run_id,source_log=str(log),log_sha256=original_log_hash,
               launch_lock_sha256=sha(root/'results'/run_id/'launch.lock.json'),source_params_sha256=sha(params_path),
               expected_tasks=len(expected),outcomes=dict(Counter(x['status'] for x in outcomes.values()))))
    cache, cohort, audit = {}, [], []
    for row in rows:
        hid=row['hypothesis_id']; print(f'Verifying {hid}...',flush=True)
        receipts=[]
        anchor=dict(row,settings=read_json(settings),settings_sha256=sha(settings),
                    design_sha256=lock['design_sha256'],
                    tool_sha256=design_lock['caastools_sha256'],smoke=params.get('synthetic_inputs',False))
        for gene in genes:
            state=outcomes[(hid,gene['gene_id'])]
            record=dict(state,receipt_path='',receipt_sha256='')
            if state['status'] in ('completed','cached'):
                receipt_path=resolve_receipt(root,lock,run_id,row,gene['gene_id'],state,cache)
                receipt=validate_receipt(receipt_path,run_id,row,gene,lock)
                require(receipt['tool_sha256']==anchor['tool_sha256'],'Receipt uses a different CAAStools version')
                require(receipt['settings']==anchor['settings'],'Receipt settings content mismatch')
                receipts.append(receipt_path)
                record.update(receipt_path=str(receipt_path),receipt_sha256=sha(receipt_path))
            audit.append(record)
        folder=output/'hypotheses'/hid
        summary=summarize(hid,inventory,None,folder,True,receipt_paths=receipts,metadata=anchor)
        missing=read_tsv(folder/'non_completed_genes.tsv')
        for record in missing:
            record['reason']=outcomes[(hid,record['gene'])]['status']
        write_tsv(folder/'non_completed_genes.tsv',missing,['hypothesis_id','strategy_id','replicate_id','gene','status','reason'])
        write_json(folder/'summary.lock.json',{p.name:sha(p) for p in sorted(folder.iterdir()) if p.is_file() and p.name!='summary.lock.json'})
        cohort.append(dict(hypothesis_id=hid,summary_dir=str(folder),summary_lock_sha256=sha(folder/'summary.lock.json')))
        print(f'{hid}: {summary["completed_genes"]} completed; {summary["non_completed_genes"]} missing',flush=True)
    write_tsv(output/'source_task_audit.tsv',audit,list(audit[0]))
    write_tsv(output/'cohort.tsv',cohort,['hypothesis_id','summary_dir','summary_lock_sha256'])
    build(output/'cohort.tsv',output/'matched',allow_incomplete=True)
    functional(output/'matched',settings,output/'functional',discovery_only=True)
    if not tables_only:
        subprocess.run([rscript,str(ROOT/'scripts/plot_validation.R'),str(output/'matched/strategy_comparison.tsv'),
                        str(output/'figures'),str(output/'functional')],check=True)
        shutil.copy2(ROOT/'scripts/plot_validation.R',output/'figures/plot_validation.R')
        shutil.copy2(ROOT/'scripts/figures.md',output/'figures/README.md')
    require(sha(log)==original_log_hash,'Source log changed during recovery; stop the original driver first')
    write_json(output/'recovery.complete.json',dict(status='recovery_complete',
               biological_completeness=read_json(output/'matched/cohort.json')['execution_status'],
               expected_tasks=len(expected),non_completed_tasks=sum(x['status'] not in ('cached','completed') for x in outcomes.values()),
               figures_generated=not tables_only, recovery_script_sha256=sha(__file__),
               summary_script_sha256=sha(ROOT/'scripts/summarize_hypothesis.py')))
    print(f'Recovery complete: {output}',flush=True)
    return output


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--root',type=Path,default=ROOT)
    p.add_argument('--run-id',required=True);p.add_argument('--recovery-id',required=True)
    p.add_argument('--log',type=Path,required=True);p.add_argument('--rscript',default='Rscript')
    p.add_argument('--tables-only',action='store_true')
    a=p.parse_args();recover(a.root,a.run_id,a.recovery_id,a.log,a.rscript,a.tables_only)
