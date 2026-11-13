"""P3 input contract and real Nextflow cache extension with synthetic tasks."""
import os
import shlex
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path
ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/'scripts'))
from common import read_json, read_tsv, sha, write_json
from check_additional_design import check_additional, additional_identity
from validate_inputs import validate
DESIGN = ROOT/'inputs/frozen/brain-260917-v1'
ADDITIONAL = ROOT/'inputs/additional/venditti-upper20-261113-v1'
INVENTORY = ROOT/'tests/fixtures/alignment_manifest.tsv'


class AdditionalP3(unittest.TestCase):
    def test_exact_bundle_and_cycles(self):
        rows=check_additional(DESIGN,ADDITIONAL)
        self.assertEqual([r['hypothesis_id'] for r in rows],['P3_000'])
        self.assertEqual(int(rows[0]['possible_cycles']),10725)
        self.assertEqual(int(rows[0]['selected_cycles']),100)

    def test_reject_changed_bundle(self):
        import shutil
        with tempfile.TemporaryDirectory() as folder:
            changed=Path(folder)/'p3';shutil.copytree(ADDITIONAL,changed)
            p=changed/'configs/P3_000/pool.cfg';p.write_text(p.read_text().replace('Macaca_nigra','Homo_sapiens'))
            with self.assertRaisesRegex(ValueError,'SHA-256 mismatch'):
                check_additional(DESIGN,changed)

    def test_validation_uses_additional_pool_and_original_phenotype(self):
        with tempfile.TemporaryDirectory() as folder:
            out=Path(folder)
            validate(DESIGN,ADDITIONAL/'selection.tsv',INVENTORY,ROOT/'inputs/settings.json',
                     'unit-p3',out,direct_discovery=True,additional_design=ADDITIONAL)
            rows=read_tsv(out/'jobs.tsv')
            self.assertEqual(len(rows),2)
            job=read_json(Path(rows[0]['job_file']))
            self.assertEqual(job['hypothesis_id'],'P3_000')
            self.assertTrue(job['pool_path'].startswith(str(ADDITIONAL)))
            self.assertEqual(job['design_sha256'],sha(DESIGN/'design.lock.json'))
            self.assertEqual(job['additional_design_sha256'],sha(ADDITIONAL/'additional.lock.json'))

    def test_old_process_definitions_and_job_scientific_code_unchanged(self):
        old=(ROOT/'tests/fixtures/pre_p3_main.nf').read_text()
        new=(ROOT/'main.nf').read_text()
        for name,end in [('PREPARE_VALIDATION','VALIDATE_INPUTS'),('VALIDATE_INPUTS','CAAS_POOLED'),('CAAS_POOLED','ASSEMBLE_FILTER_AND_QUERY')]:
            original=old.split('process '+name+' {',1)[1].split('process '+end+' {',1)[0]
            current=new.split('process '+name+' {',1)[1].split('process '+end+' {',1)[0]
            if name=='VALIDATE_INPUTS':current=current.split('process VALIDATE_ADDITIONAL_INPUTS {')[0]
            self.assertEqual(original,current)
        self.assertEqual(sha(ROOT/'scripts/run_caas.py'),'5299ff0c2f457351d7018e8489ba5731f6eb0d714df9cce31dc78ac45a22ac90')

    def test_launcher_extends_same_run_without_rewriting_original_lock(self):
        from prepare_cluster_launch import select_rows,selected_modes
        from common import write_tsv
        with tempfile.TemporaryDirectory() as folder:
            base=Path(folder);rows=select_rows(DESIGN,selected_modes('all'))
            selection=base/'selection.tsv';write_tsv(selection,rows,list(rows[0]))
            run_id='unit-p3-resume'
            cmd=[sys.executable,str(ROOT/'run_validation.py'),'benchmark','--run-id',run_id,
                 '--selection',str(selection),'--alignments',str(INVENTORY),
                 '--results-root',str(base/'results'),'--summaries-root',str(base/'summaries'),
                 '--work-dir',str(base/'work')]
            def run(extra):return subprocess.run(cmd+extra,capture_output=True,text=True)
            result=run([]);self.assertEqual(result.returncode,0,result.stderr)
            plan_dir=ROOT/'launch-plans'/run_id
            plan=read_json(sorted(plan_dir.glob('*.plan.json'))[-1])
            fields=('stage direct_discovery design_sha256 selection_sha256 alignment_manifest_sha256 cohort_sha256 settings_sha256 code_fingerprint annotation_sha256 cluster_config_sha256 python_command r_command nextflow_command conda_prefix conda_init_sha256 conda_packages_sha256 profile work_dir results_root summaries_root').split()
            identity={key:plan[key] for key in fields}
            bridge=read_json(ROOT/'launch_support/resource_resume_bridge.json')
            identity['code_fingerprint']=bridge['source_code_fingerprints'][-1]
            lock=base/'results'/run_id/'launch.lock.json';write_json(lock,identity)
            before=lock.read_bytes();selection_before=selection.read_bytes()
            result=run(['--resume','--additional-design',str(ADDITIONAL)])
            self.assertEqual(result.returncode,0,result.stderr)
            extended=read_json(sorted(plan_dir.glob('*.plan.json'))[-1])
            self.assertEqual(extended['hypotheses'],6)
            self.assertEqual(extended['expected_tasks'],12)
            self.assertEqual(extended['cycles'],430)
            self.assertEqual(extended['code_fingerprint'],identity['code_fingerprint'])
            self.assertEqual(lock.read_bytes(),before)
            self.assertEqual(selection.read_bytes(),selection_before)
            write_json(lock.parent/'additional_hypotheses.lock.json',additional_identity(DESIGN,ADDITIONAL))
            result=run(['--resume'])
            self.assertNotEqual(result.returncode,0)
            self.assertIn('Additional hypothesis changed or omitted',result.stderr)

    @unittest.skipUnless(os.environ.get('NF_TEST_COMMAND'),'Set NF_TEST_COMMAND for real cache test')
    def test_five_cached_tasks_and_one_new_p3_task(self):
        with tempfile.TemporaryDirectory(prefix='p3-resume-') as folder:
            root=Path(folder);(root/'scripts').mkdir();(root/'shared-tool').mkdir()
            (root/'scripts/run_caas.py').write_text('''import json,sys
from pathlib import Path
j=json.loads(Path(sys.argv[sys.argv.index('--job')+1]).read_text())
with (Path(__file__).resolve().parents[1]/'calls.txt').open('a') as f:f.write(j['hypothesis_id']+'\\n')
out=Path(sys.argv[sys.argv.index('--output')+1]);out.mkdir();(out/'done.txt').write_text('synthetic')
''')
            ids=['N3','N4','N5','P1','N6','P3']
            for sid in ids:write_json(root/(sid+'.json'),dict(hypothesis_id=sid+'_000'))
            for name in ['alignment.phy','pool.cfg','cycles.cfg']:(root/name).write_text('immutable fixture')
            wf=(ROOT/'main.nf').read_text()
            process='process CAAS_POOLED {'+wf.split('process CAAS_POOLED {',1)[1].split('process ASSEMBLE_FILTER_AND_QUERY',1)[0]
            workflow='''
workflow {
    base = Channel.of('N3','N4','N5','P1','N6')
    selected = params.add_p3 ? base.mix(Channel.of('P3')) : base
    jobs = selected.map { sid -> tuple(sid+'_000',sid,'000','GENE',file(sid+'.json'),file('alignment.phy'),file('pool.cfg'),file('cycles.cfg')) }
    CAAS_POOLED(jobs,file('shared-tool').toAbsolutePath().toString())
}
'''
            (root/'main.nf').write_text('nextflow.enable.dsl=2\n'+process+workflow)
            (root/'nextflow.config').write_text('''params.direct_discovery=true
params.results_root="${projectDir}/results"
params.run_id='test'
params.code_fingerprint='preserved-reviewed-fingerprint'
params.python_command='python3'
params.add_p3=false
process.executor='local'
process.maxForks=1
''')
            env=dict(os.environ,NXF_ANSI_LOG='false',NXF_OFFLINE='true')
            command=shlex.split(os.environ['NF_TEST_COMMAND'])+['run','main.nf']
            def run(args):
                p=subprocess.run(command+args,cwd=root,env=env,capture_output=True,text=True,timeout=90)
                self.assertEqual(p.returncode,0,p.stdout+p.stderr)
                return p.stdout+p.stderr
            run(['-name','p3_cache_test'])
            self.assertEqual(len((root/'calls.txt').read_text().splitlines()),5)
            # A real source change adds the P3 routing; CAAS task stays byte-identical.
            with (root/'main.nf').open('a') as f:f.write('\n// Additional P3 routing version\n')
            text=run(['-resume','p3_cache_test','--add_p3','true'])
            calls=(root/'calls.txt').read_text().splitlines()
            self.assertEqual(calls.count('P3_000'),1)
            self.assertEqual(len(calls),6,text)
            self.assertEqual(text.count('Cached process > CAAS_POOLED'),5,text)
            run(['-resume','p3_cache_test','--add_p3','true'])
            self.assertEqual(len((root/'calls.txt').read_text().splitlines()),6)


if __name__=='__main__':unittest.main()
