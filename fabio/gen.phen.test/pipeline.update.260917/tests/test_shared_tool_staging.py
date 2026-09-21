"""Run the actual CAAS process declaration with a fake runner, never real CAAS.

Set NF_TEST_COMMAND to the Nextflow launcher command to enable the local
integration test. No scheduler or cluster is used.
"""
import os
from pathlib import Path
import shlex
import subprocess
import tempfile
import unittest

ROOT = Path(__file__).resolve().parents[1]


class SharedToolStaging(unittest.TestCase):
    def test_shared_directory_is_not_a_staged_input(self):
        workflow = (ROOT/'main.nf').read_text()
        process = workflow.split('process CAAS_POOLED {', 1)[1].split('process ASSEMBLE_FILTER_AND_QUERY', 1)[0]
        self.assertIn('val caastools', process)
        self.assertNotIn('path caastools', process)
        self.assertIn('toAbsolutePath().toString()', workflow)
        self.assertIn('tree_hash(tool_dir) == job["tool_sha256"]', (ROOT/'scripts/run_caas.py').read_text())

    @unittest.skipUnless(os.environ.get('NF_TEST_COMMAND'), 'Set NF_TEST_COMMAND for local Nextflow integration')
    def test_incomplete_task_reexecution_without_tool_copy(self):
        with tempfile.TemporaryDirectory(prefix='shared-caas-') as folder:
            root = Path(folder)
            shared = root/'shared-tool'; shared.mkdir(); (shared/'ct').write_text('immutable tool fixture')
            (root/'scripts').mkdir()
            # Mock only the computation: verify the absolute shared path and emit an output.
            (root/'scripts/run_caas.py').write_text(
                'import sys\nfrom pathlib import Path\n'
                'tool=Path(sys.argv[sys.argv.index("--tool-dir")+1])\n'
                'assert tool.is_absolute() and (tool/"ct").is_file()\n'
                'assert not Path("caastools").exists()\n'
                'out=Path(sys.argv[sys.argv.index("--output")+1]); out.mkdir(exist_ok=True)\n'
                '(out/"done.txt").write_text(str(tool))\n')
            for name in ('job.json', 'alignment.phy', 'pool.cfg', 'cycles.cfg'):
                (root/name).write_text('fixture')
            wf = (ROOT/'main.nf').read_text()
            process = 'process CAAS_POOLED {' + wf.split('process CAAS_POOLED {', 1)[1].split('process ASSEMBLE_FILTER_AND_QUERY', 1)[0]
            (root/'main.nf').write_text('nextflow.enable.dsl=2\n' + process + '''
workflow {
    jobs = Channel.of(tuple('N3_000','N3','000','GENE',file('job.json'),file('alignment.phy'),file('pool.cfg'),file('cycles.cfg')))
    CAAS_POOLED(jobs, file(params.caastools_dir, checkIfExists:true).toAbsolutePath().toString())
}
''')
            (root/'nextflow.config').write_text('''
params.results_root = "${projectDir}/results"
params.run_id = 'test'
params.code_fingerprint = 'test'
params.python_command = 'python3'
params.caastools_dir = "${projectDir}/shared-tool"
process.executor = 'local'
''')
            command = shlex.split(os.environ['NF_TEST_COMMAND']) + ['run', 'main.nf']
            env = dict(os.environ, NXF_ANSI_LOG='false', NXF_OFFLINE='true')
            def run(extra):
                result = subprocess.run(command+extra, cwd=root, env=env, capture_output=True, text=True, timeout=90)
                self.assertEqual(result.returncode, 0, result.stdout+result.stderr)
                return result.stdout+result.stderr
            run([])
            task = next((root/'work').glob('*/*/.command.run')).parent
            self.assertFalse((task/'caastools').exists())
            self.assertNotIn('rm -f caastools', (task/'.command.run').read_text())
            (task/'GENE').rename(root/'saved-output')
            rerun = run(['-resume'])
            self.assertNotIn('Cached process', rerun)
            self.assertEqual((shared/'ct').read_text(), 'immutable tool fixture')


if __name__ == '__main__': unittest.main()
