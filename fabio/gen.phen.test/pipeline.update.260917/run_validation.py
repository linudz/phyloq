#!/usr/bin/env python3
"""Explicit, collision-safe launcher. Never defaults to production or all strategies."""
import argparse
import datetime
import os
import re
import shlex
import subprocess
import sys
from pathlib import Path
ROOT = Path(__file__).resolve().parent
sys.path.insert(0, str(ROOT / 'scripts'))
from common import digest, identifier, read_json, read_tsv, require, sha, tree_hash, verify_design, write_json
from check_environment import inspect_runtime
sys.path.insert(0, str(ROOT / 'launch_support'))
from resource_resume import resume_identity

STAGES = ['prepare', 'smoke', 'benchmark', 'deterministic', 'paired', 'summarize', 'reference-rerun', 'auxiliary-p2']


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('stage', choices=STAGES)
    p.add_argument('--design', type=Path, default=ROOT / 'inputs/frozen/brain-260917-v1')
    p.add_argument('--run-id')
    p.add_argument('--selection', type=Path, help='Explicit TSV; default is the frozen selection for this stage')
    p.add_argument('--alignments', type=Path, help='Explicit alignment_manifest.tsv, not an alignment directory')
    p.add_argument('--cohort', type=Path, help='Explicit complete-summary manifest for summarize')
    p.add_argument('--settings', type=Path, default=ROOT / 'inputs/settings.json')
    p.add_argument('--approvals', type=Path)
    p.add_argument('--annotation', type=Path, help='Frozen snapshot or version-checked response-cache directory')
    p.add_argument('--profile', choices=['local', 'cluster'], default='local')
    p.add_argument('--python')
    p.add_argument('--nextflow')
    p.add_argument('--rscript')
    p.add_argument('--conda-prefix', type=Path, help='Shared installed environment; normally supplied by run_pipeline.sh')
    p.add_argument('--conda-init', type=Path, help='Conda initialization for worker jobs; normally supplied by run_pipeline.sh')
    p.add_argument('--work-dir', type=Path, default=ROOT / 'work')
    p.add_argument('--results-root', type=Path, default=ROOT / 'results')
    p.add_argument('--summaries-root', type=Path, default=ROOT / 'summaries')
    p.add_argument('--cluster-config', type=Path, help='User-reviewed resource/environment overrides')
    p.add_argument('--task-time', default='30m', help='Walltime per CAAS task; safe to change on resume')
    p.add_argument('--task-memory', default='2 GB', help='Memory per CAAS task; safe to change on resume')
    p.add_argument('--task-cpus', type=int, default=1, help='CPUs per CAAS task; safe to change on resume')
    p.add_argument('--resume', action='store_true', help='Only resume this exact run and frozen selection')
    p.add_argument('--execute', action='store_true', help='Without this flag, print/save a plan only')
    a = p.parse_args()
    require(a.task_cpus > 0, 'Task CPUs must be positive')
    require(re.fullmatch(r'[1-9][0-9]*(?:\.[0-9]+)?\s*(?:ms|s|m|h|d)', a.task_time) is not None,
            'Use a positive task time with a unit, e.g. 30m or 1h')
    require(re.fullmatch(r'[1-9][0-9]*(?:\.[0-9]+)?\s*(?:KB|MB|GB|TB)', a.task_memory, re.I) is not None,
            'Use positive task memory with a unit, e.g. 2 GB')
    if a.conda_prefix:
        a.conda_prefix = a.conda_prefix.resolve()
        require(a.conda_init is not None and a.conda_init.is_file(), 'Conda mode requires --conda-init pointing to conda.sh')
        a.conda_init = a.conda_init.resolve()
    a.python = a.python or (str(a.conda_prefix/'bin/python') if a.conda_prefix else sys.executable)
    a.nextflow = a.nextflow or (str(a.conda_prefix/'bin/nextflow') if a.conda_prefix else 'nextflow')
    a.rscript = a.rscript or (str(a.conda_prefix/'bin/Rscript') if a.conda_prefix else 'Rscript')
    design = a.design.resolve(); lock = verify_design(design)
    if a.stage == 'prepare':
        from check_design import check
        print(__import__('json').dumps(check(design), indent=2))
        return
    require(a.run_id is not None, 'An explicit --run-id is required')
    identifier(a.run_id)
    cfg = read_json(a.settings)
    # Fabio's direct cluster workflow: collect CAAS results now, without a
    # mandatory pilot or pretending that downstream statistics were approved.
    direct_discovery = a.stage == 'benchmark'
    if direct_discovery:
        require(a.selection is not None, 'benchmark requires an explicit --selection; use launch_all_caas.sh')
        require(a.annotation is None, 'benchmark is CAAS discovery only; run enrichment separately with summarize')
    selection = (a.selection or design / 'selections' / (a.stage + '.tsv')).resolve() if a.stage != 'summarize' else None
    if a.stage == 'summarize':
        require(a.cohort is not None and a.cohort.is_file(), 'summarize requires --cohort')
        selected, inventory = [], []
    else:
        require(a.alignments is not None and a.alignments.is_file(), '--alignments must name an inventory TSV')
        selected, inventory = read_tsv(selection), read_tsv(a.alignments)
        require(selected and inventory, 'Empty selection or inventory')
        require(not any(r['strategy_id'] in ('R0', 'R1') for r in selected),
                'R0/R1 randomized null series have been disabled, including custom selections.')
        allowed = {'deterministic': {'N3', 'N4', 'N5'}, 'paired': {'P1', 'N6'}, 'reference-rerun': {'P0', 'N0', 'N1', 'N2'}, 'auxiliary-p2': {'P2'}}
        if a.stage in allowed:
            require({r['strategy_id'] for r in selected} <= allowed[a.stage], 'Selection contains strategies outside the requested stage')
    # Bind code, actual selection, annotations, settings and outputs to the run.
    fingerprint = digest(dict(scripts=tree_hash(ROOT / 'scripts'), workflow=sha(ROOT / 'main.nf'),
                              config=sha(ROOT / 'nextflow.config'), cluster=tree_hash(ROOT / 'conf'), launcher=sha(__file__),
                              environment=sha(ROOT/'environment.yml'),
                              resume_helper=sha(ROOT/'launch_support/resource_resume.py')))
    implementation_fingerprint = fingerprint
    # Conda package records detect updates in-place even at the same prefix.
    package_lock = {p.name:sha(p) for p in sorted((a.conda_prefix/'conda-meta').glob('*.json'))} if a.conda_prefix else None
    identity = dict(stage=a.stage, direct_discovery=direct_discovery, design_sha256=sha(design / 'design.lock.json'),
        selection_sha256=sha(selection) if selection else None,
        alignment_manifest_sha256=sha(a.alignments) if a.alignments else None,
        cohort_sha256=sha(a.cohort) if a.cohort else None, settings_sha256=sha(a.settings),
        code_fingerprint=fingerprint, annotation_sha256=tree_hash(a.annotation) if a.annotation else None,
        cluster_config_sha256=sha(a.cluster_config) if a.cluster_config else None,
        python_command=a.python, r_command=a.rscript, nextflow_command=a.nextflow,
        conda_prefix=str(a.conda_prefix) if a.conda_prefix else None,
        conda_init_sha256=sha(a.conda_init) if a.conda_init else None,
        conda_packages_sha256=digest(package_lock) if package_lock is not None else None,
        profile=a.profile, work_dir=str(a.work_dir.resolve()),
        results_root=str(a.results_root.resolve()), summaries_root=str(a.summaries_root.resolve()))
    out = a.results_root.resolve() / a.run_id
    launch_lock = out / 'launch.lock.json'
    if launch_lock.exists():
        require(a.resume, 'Run ID already exists; use --resume for the exact same inputs, or a new run ID')
        identity = resume_identity(read_json(launch_lock), identity,
                                   read_json(ROOT/'launch_support/resource_resume_bridge.json'))
        # Keep the reviewed old task-script fingerprint to reuse valid cache.
        # Actual implementation and resources are recorded in each launch plan.
        fingerprint = identity['code_fingerprint']
    else:
        require(not a.resume, 'Cannot resume a run with no launch lock')
        require(not out.exists(), 'Output directory exists without a run lock; choose a new run ID')
    if a.execute and a.stage not in ('smoke', 'summarize', 'benchmark'):
        require(a.approvals is not None, 'No production launch without explicit reviewed approvals')
        ap = read_json(a.approvals)
        for key in ('design_sha256', 'selection_sha256', 'alignment_manifest_sha256', 'settings_sha256'):
            require(ap.get(key) == identity[key], 'Approval missing/mismatched: ' + key)
        require(all(ap.get(k) is True for k in ('pool_review', 'technical_review', 'production_review')), 'Production approval gates not satisfied')
        require(bool(ap.get('approved_by')), 'Reviewer name required')
        if any(r['strategy_id'] in ('R0', 'R1') for r in selected):
            require(ap.get('statistical_review') is True and cfg['statistics']['approved'], 'Random null requires approved statistical settings')
        if any(int(r['replicate_id']) > 19 for r in selected):
            cost = Path(ap.get('pilot_cost_report', ''))
            require(cost.is_file() and sha(cost) == ap.get('pilot_cost_report_sha256'), 'Null extension needs a checksummed reviewed pilot cost report')
    runtime = None
    if a.execute:
        runtime = inspect_runtime(a.python, a.rscript, a.nextflow, a.conda_prefix, a.profile == 'cluster')
        if a.resume:
            prior_runtime = read_json(out/'environment.json')
            require(runtime['runtime_sha256'] == prior_runtime['runtime_sha256'], 'Installed runtime changed; refusing incompatible resume')
    runtime_setup = ''
    if a.conda_prefix:
        runtime_setup = f'set +u; source {shlex.quote(str(a.conda_init))}; conda activate {shlex.quote(str(a.conda_prefix))}; set -u; '
    launch_id = datetime.datetime.now(datetime.timezone.utc).strftime('%Y%m%dT%H%M%S%fZ')
    synthetic_inputs = bool(inventory) and all(
        (Path(r['alignment_path']) if Path(r['alignment_path']).is_absolute() else a.alignments.parent/ r['alignment_path']).resolve().is_relative_to(ROOT/'tests/fixtures')
        for r in inventory)
    params = dict(stage=a.stage, direct_discovery=direct_discovery, synthetic_inputs=synthetic_inputs,
        task_time=a.task_time, task_memory=a.task_memory, task_cpus=a.task_cpus,
        run_id=a.run_id, design=str(design), settings=str(a.settings.resolve()),
        selection_manifest=str(selection) if selection else None, alignment_manifest=str(a.alignments.resolve()) if a.alignments else None,
        cohort_manifest=str(a.cohort.resolve()) if a.cohort else None, approvals=str(a.approvals.resolve()) if a.approvals else None,
        code_fingerprint=fingerprint, python_command=a.python, r_command=a.rscript, runtime_setup=runtime_setup,
        results_root=str(a.results_root.resolve()),
        summaries_root=str(a.summaries_root.resolve()), annotation_input=str(a.annotation.resolve()) if a.annotation else str(ROOT / 'inputs/annotations.disabled.txt'),
        annotation_backend=cfg['enrichment']['backend'], launch_id=launch_id)
    plan_dir = ROOT / 'launch-plans' / a.run_id
    log_dir = ROOT / 'logs' / a.run_id
    log_dir.mkdir(parents=True, exist_ok=True)
    plan_dir.mkdir(parents=True, exist_ok=True)
    param_file = plan_dir / (launch_id + '.params.json')
    write_json(param_file, params)
    cmd = [a.nextflow, '-log', str(log_dir / (launch_id + '.nextflow.log')), 'run', str(ROOT / 'main.nf'),
           '-profile', a.profile, '-params-file', str(param_file), '-work-dir', str(a.work_dir.resolve())]
    if a.cluster_config: cmd += ['-c', str(a.cluster_config.resolve())]
    if a.resume: cmd += ['-resume', 'caas_' + a.run_id.replace('-', '_')]
    else: cmd += ['-name', 'caas_' + a.run_id.replace('-', '_')]
    plan = dict(identity, run_id=a.run_id, hypotheses=len(selected), alignments=len(inventory),
                actual_implementation_fingerprint=implementation_fingerprint,
                task_resources=dict(time=a.task_time, memory=a.task_memory, cpus=a.task_cpus),
                expected_tasks=len(selected)*len(inventory), cycles=sum(int(r['selected_cycles']) for r in selected),
                command=cmd, execute=a.execute, production=a.stage not in ('smoke', 'summarize'),
                approvals_required_for_execution=a.stage not in ('smoke', 'summarize', 'benchmark'))
    write_json(plan_dir / (launch_id + '.plan.json'), plan)
    print(f"{a.stage}: {len(selected)} hypotheses × {len(inventory)} alignments = {plan['expected_tasks']} discovery tasks")
    print(shlex.join(cmd), flush=True)
    if not a.execute:
        print('PLAN ONLY: no discovery tasks submitted. Add --execute to run the selected analysis.'); return
    if not launch_lock.exists(): write_json(launch_lock, identity)
    write_json(out/'environment.json',runtime)
    env = dict(os.environ, PYTHONDONTWRITEBYTECODE='1', PYTHONHASHSEED='0', NXF_ANSI_LOG='false', NXF_DISABLE_CHECK_LATEST='true')
    env['JAVA_HOME'] = env['NXF_JAVA_HOME'] = runtime['identity']['java_home']
    if a.conda_prefix:
        env['PATH'] = str(a.conda_prefix/'bin') + os.pathsep + env.get('PATH','')
        env['PYTHONNOUSERSITE'] = '1'
    if a.stage == 'smoke': env['NXF_OFFLINE'] = 'true'
    completed = subprocess.run(cmd, cwd=ROOT, env=env)
    raise SystemExit(completed.returncode)


if __name__ == '__main__':
    main()
