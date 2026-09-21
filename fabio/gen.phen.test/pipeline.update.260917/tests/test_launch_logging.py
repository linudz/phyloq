"""Logging regression checks; no Conda, Nextflow or scheduler required."""
import os
from pathlib import Path
import shutil
import subprocess
import tempfile
import unittest

ROOT = Path(__file__).resolve().parents[1]


class LaunchLoggingTests(unittest.TestCase):
    def test_submission_and_early_conda_failure(self):
        with tempfile.TemporaryDirectory(prefix='caas logging ') as temp:
            root = Path(temp)
            (root / 'conf').mkdir()
            (root / 'main.nf').touch()
            shutil.copy(ROOT / 'launch_all_caas.sh', root)
            shutil.copy(ROOT / 'conf/logging_helpers.sh', root / 'conf')
            (root / 'conf/conda_helpers.sh').write_text(
                'caas_init_conda() { echo "TEST: conda missing" >&2; return 17; }\n')
            binary = root / 'bin'
            binary.mkdir()
            scheduler = binary / 'sbatch'
            scheduler.write_text('#!/bin/bash\nprintf "%s\\n" "$@"\n')
            scheduler.chmod(0o755)
            env = dict(os.environ, CAAS_VALIDATION_ROOT=str(root),
                       PATH=str(binary) + os.pathsep + os.environ['PATH'])
            env.pop('SLURM_JOB_ID', None)
            cmd = ['bash', str(root / 'launch_all_caas.sh'), '--run-id', 'logging-test',
                   '--alignments-pattern', '/unused/*.phy']
            submitted = subprocess.run(cmd, env=env, capture_output=True, text=True)
            self.assertEqual(submitted.returncode, 0, submitted.stderr)
            logs = root / 'logs/logging-test'
            self.assertIn('--output=' + str(logs / 'slurm-%j.out'), submitted.stdout)
            self.assertIn('--error=' + str(logs / 'slurm-%j.err'), submitted.stdout)
            self.assertEqual(len(list(logs.glob('submission-*.log'))), 1)
            env['SLURM_JOB_ID'] = '12345'
            failed = subprocess.run(cmd, env=env, capture_output=True, text=True)
            self.assertEqual(failed.returncode, 17)
            self.assertIn('TEST: conda missing', (logs / 'driver-12345.log').read_text())
            self.assertFalse(list(root.glob('*.log')))


if __name__ == '__main__':
    unittest.main()
