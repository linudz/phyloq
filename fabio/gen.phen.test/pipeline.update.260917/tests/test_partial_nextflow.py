"""Real Nextflow routing/assembly with deliberately failed synthetic CAAS tasks."""
import json
import os
from pathlib import Path
import shlex
import subprocess
import sys
import tempfile
import unittest
ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/'scripts'))
from common import read_json, write_json, write_tsv


@unittest.skipUnless(os.environ.get('NF_TEST_COMMAND'), 'Set NF_TEST_COMMAND for Nextflow integration')
class PartialNextflow(unittest.TestCase):
    def test_one_failed_gene_and_one_entire_failed_hypothesis(self):
        with tempfile.TemporaryDirectory(prefix='caas-partial-') as folder:
            root=Path(folder); (root/'scripts').mkdir(); (root/'shared-tool').mkdir()
            # Symlink unmodified scientific summary code; mock only CAAS execution.
            for source in (ROOT/'scripts').iterdir():
                if source.is_file() and source.name != 'run_caas.py':
                    (root/'scripts'/source.name).symlink_to(source)
            (root/'scripts/run_caas.py').write_text('''import sys,json,hashlib
from pathlib import Path
j=json.loads(Path(sys.argv[sys.argv.index('--job')+1]).read_text())
g=j['gene_id']
if j['hypothesis_id']=='BAD_000' or g=='FAIL': sys.exit(1)
out=Path(g); out.mkdir(exist_ok=True)
legacy=out/(g+'.pooled.caas.tsv'); legacy.write_text('')
events=out/(g+'.pooled.caas.events.tsv'); events.write_text('gene\\tposition\\tevent_id\\tprimary_event\\tpositional_pvalue\\tfg_support_count\\tbg_support_count\\n')
j.update(status='success',coverage=dict(coverage_testable=True,coverage_testable_positions=1,present_fg=4,present_bg=4,missing_fg=[],missing_bg=[]),elapsed_seconds=1,cpu_seconds=1,max_rss_bytes=1,output_bytes=0,events_sha256=hashlib.sha256(events.read_bytes()).hexdigest(),legacy_sha256=hashlib.sha256(legacy.read_bytes()).hexdigest())
(out/(g+'.receipt.json')).write_text(json.dumps(j))
''')
            settings=read_json(ROOT/'inputs/settings.json')
            for strategy in ('GOOD','BAD'):
                for gene in ('PASS','FAIL'):
                    job=dict(hypothesis_id=strategy+'_000',strategy_id=strategy,replicate_id='000',
                             gene_id=gene,settings=settings,null_family='observed',smoke=True,
                             alignment_sha256='x',settings_sha256='x',config_sha256='x',pool_sha256='x',tool_sha256='x',design_sha256='x')
                    write_json(root/(strategy+'_000__'+gene+'.json'),job)
            for name in ('alignment.phy','pool.cfg','cycles.cfg'): (root/name).write_text('fixture')
            write_tsv(root/'inventory.tsv',[dict(gene_id=g,sha256='x') for g in ('PASS','FAIL')],['gene_id','sha256'])
            original=(ROOT/'main.nf').read_text()
            processes='process CAAS_POOLED {'+original.split('process CAAS_POOLED {',1)[1].split('process COMMON_BACKGROUND {',1)[0]
            grouping=original.split('        // Seed every hypothesis',1)[1].split('        // The cohort',1)[0]
            grouping='        // Seed every hypothesis'+grouping
            (root/'main.nf').write_text('nextflow.enable.dsl=2\n'+processes+'''
workflow {
    jobs = Channel.of('GOOD','BAD').flatMap { strategy -> ['PASS','FAIL'].collect { gene ->
        tuple(strategy+'_000',strategy,'000',gene,file(strategy+'_000__'+gene+'.json'),file('alignment.phy'),file('pool.cfg'),file('cycles.cfg'))
    } }
    alignments=file('inventory.tsv')
    CAAS_POOLED(jobs,file('shared-tool').toAbsolutePath().toString())
'''+grouping+'\n}\n')
            (root/'nextflow.config').write_text('''
params.direct_discovery = true
params.results_root = "${projectDir}/results"
params.run_id = 'test'
params.code_fingerprint = 'test'
params.python_command = 'python3'
process.executor = 'local'
process.maxRetries = 0
''')
            command=shlex.split(os.environ['NF_TEST_COMMAND'])+['run','main.nf']
            r=subprocess.run(command,cwd=root,env=dict(os.environ,NXF_ANSI_LOG='false',NXF_OFFLINE='true'),capture_output=True,text=True,timeout=90)
            self.assertEqual(r.returncode,0,r.stdout+r.stderr)
            for strategy,completed,missing in [('GOOD',1,1),('BAD',0,2)]:
                summary=read_json(root/f'results/test/{strategy}/000/{strategy}_000/summary.json')
                self.assertEqual((summary['status'],summary['completed_genes'],summary['non_completed_genes']),('partial',completed,missing))


if __name__ == '__main__': unittest.main()
