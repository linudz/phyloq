"""Bounded checks of the direct cluster interface; never submits a real job."""
import os
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path
ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT/'scripts'))
from common import read_json, read_tsv, sha, write_json, write_tsv
from prepare_cluster_launch import NEW_MODES, prepare, selected_modes, select_rows
from validate_inputs import validate
from functional import run as run_functional
DESIGN = ROOT/'inputs/frozen/brain-260917-v1'
INVENTORY = ROOT/'tests/fixtures/alignment_manifest.tsv'


class UnifiedSelection(unittest.TestCase):
    def test_all_new_modes_are_full_and_unique(self):
        rows = select_rows(DESIGN, selected_modes('all'))
        self.assertEqual(len(rows), 203)
        self.assertEqual(sum(int(r['selected_cycles']) for r in rows), 20130)
        self.assertEqual({r['strategy_id'] for r in rows}, {'N3','N4','N5','R0','R1','P1','N6'})
        self.assertTrue({'R0_099','R1_099'} <= {r['hypothesis_id'] for r in rows})
        self.assertEqual(len(select_rows(DESIGN, selected_modes('all', True, True))), 208)
        self.assertEqual(len(select_rows(DESIGN, selected_modes('deterministic,r0'))), 102)

    def test_launch_inputs_reused_not_changed_on_resume(self):
        with tempfile.TemporaryDirectory() as tmp:
            out = Path(tmp)/'launch'
            report = prepare(DESIGN, out, 'unit-direct', list(NEW_MODES), inventory_path=INVENTORY)
            initial = sha(out/'selection.tsv')
            self.assertEqual(report['expected_caas_tasks'], 406)
            self.assertEqual(prepare(DESIGN, out, 'unit-direct', list(NEW_MODES), inventory_path=INVENTORY, resume=True), report)
            self.assertEqual(initial, sha(out/'selection.tsv'))
            with self.assertRaisesRegex(ValueError, 'changed'):
                prepare(DESIGN, out, 'unit-direct', ['deterministic'], inventory_path=INVENTORY, resume=True)

    def test_direct_full_null_needs_no_pilot_or_approval(self):
        rows = [r for r in read_tsv(DESIGN/'execution_manifest.tsv') if r['hypothesis_id'] in ('N3_000','R0_099','R1_099')]
        with tempfile.TemporaryDirectory() as tmp:
            tmp = Path(tmp)
            write_tsv(tmp/'selection.tsv', rows, list(rows[0]))
            lock = validate(DESIGN, tmp/'selection.tsv', INVENTORY, ROOT/'inputs/settings.json',
                            'unit-direct', tmp/'validated', direct_discovery=True)
            self.assertTrue(lock['direct_discovery'])
            self.assertFalse(read_json(ROOT/'inputs/settings.json')['statistics']['approved'])
            self.assertEqual(len(read_tsv(tmp/'validated/jobs.tsv')), 6)
            # Direct execution does not disable frozen selection validation.
            rows[0]['selected_cycles'] = '99'
            write_tsv(tmp/'selection.tsv', rows, list(rows[0]))
            with self.assertRaisesRegex(ValueError, 'differs from frozen'):
                validate(DESIGN, tmp/'selection.tsv', INVENTORY, ROOT/'inputs/settings.json',
                         'unit-direct', tmp/'invalid', direct_discovery=True)

    def test_discovery_only_does_not_run_enrichment_or_gene_list_nulls(self):
        with tempfile.TemporaryDirectory() as tmp:
            tmp = Path(tmp); matched = tmp/'matched'; matched.mkdir()
            write_json(matched/'cohort.json', {'hypotheses':['P0_000']})
            write_tsv(matched/'common_background.tsv', [{'gene':'GENE'}], ['gene'])
            write_tsv(matched/'common_queries.tsv', [{'hypothesis_id':'P0_000','gene':'GENE'}], ['hypothesis_id','gene'])
            run_functional(matched, ROOT/'inputs/settings.json', tmp/'out', discovery_only=True)
            self.assertEqual(read_json(tmp/'out/enrichment_status.json')['status'], 'not_requested_discovery_only')
            self.assertEqual(len(read_tsv(tmp/'out/gene_list_manifest.tsv')), 1)


class UnifiedShell(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        self.base = Path(self.tmp.name); self.bin = self.base/'bin'; self.bin.mkdir()
        self.log = self.base/'log'
        self.env = dict(os.environ, PATH=str(self.bin)+os.pathsep+os.environ['PATH'], CAAS_TEST_LOG=str(self.log))
        for name in ('SLURM_JOB_ID','SLURM_SUBMIT_DIR','CAAS_VALIDATION_ROOT','CAAS_LAUNCH_CALLER_DIR','CAAS_CONDA_PREFIX','CAAS_CONDA_ENV'):
            self.env.pop(name, None)

    def tearDown(self): self.tmp.cleanup()

    def fake(self, name, contents):
        p = self.bin/name; p.write_text('#!/bin/bash\n'+contents); p.chmod(0o700)

    def test_one_submission_outside_slurm(self):
        self.fake('sbatch', 'printf "%s\\n" "$@" > "$CAAS_TEST_LOG"\n')
        p = subprocess.run(['bash',str(ROOT/'launch_all_caas.sh'),'--run-id','unit-all','--inventory',str(INVENTORY)],
                           cwd=ROOT,env=self.env,capture_output=True,text=True)
        self.assertEqual(p.returncode, 0, p.stderr)
        args = self.log.read_text().splitlines()
        self.assertEqual(args.count(str(ROOT/'launch_all_caas.sh')), 1)
        self.assertIn('--export=ALL', args)

    def driver(self, plan_only=False):
        self.fake('python', 'printf "%s\\n" "$@" >> "$CAAS_TEST_LOG"\n')
        hook = self.base/'conda.sh'
        hook.write_text('conda() { export CONDA_PREFIX="'+str(self.base)+'"; return 0; }\n')
        env = dict(self.env, CAAS_CONDA_SH=str(hook))
        if not plan_only: env['SLURM_JOB_ID'] = 'synthetic-unit-test'
        args = ['bash',str(ROOT/'launch_all_caas.sh'),'--run-id','unit-all','--inventory',str(INVENTORY)]
        if plan_only: args += ['--plan-only']
        p = subprocess.run(args,cwd=ROOT,env=env,capture_output=True,text=True)
        self.assertEqual(p.returncode, 0, p.stderr)
        return self.log.read_text()

    def test_driver_activates_and_launches_one_benchmark(self):
        args = self.driver()
        self.assertIn('prepare_cluster_launch.py', args)
        self.assertIn('\nbenchmark\n', args)
        self.assertIn('\n--execute\n', args)
        self.assertNotIn('--approvals', args)

    def test_plan_only_never_submits_or_executes(self):
        args = self.driver(plan_only=True)
        self.assertIn('run_validation.py', args)
        self.assertNotIn('--execute', args)


if __name__ == '__main__': unittest.main()
