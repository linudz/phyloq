"""Portable shell-entry-point tests with a fake Conda, never installing packages."""
import os
import subprocess
import tempfile
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]


class EnvironmentEntryPoints(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory(prefix='caas-shell-test-')
        self.base = Path(self.tmp.name)
        self.prefix = self.base/'env'; (self.prefix/'bin').mkdir(parents=True)
        self.log = self.base/'commands.log'; self.marker = self.base/'exists'
        python = self.prefix/'bin/python'
        python.write_text('#!/bin/bash\nprintf "%s\\n" "$@"\n')
        python.chmod(0o700)
        self.hook = self.base/'conda.sh'
        self.hook.write_text('''conda() {
  printf '%s\\n' "$*" >> "$CAAS_TEST_LOG"
  if [[ "$1" == info ]]; then printf '%s\\n' "$CAAS_TEST_BASE"; return 0; fi
  if [[ "$1" == activate ]]; then
    [[ -f "$CAAS_TEST_MARKER" ]] || return 1
    export CONDA_PREFIX="$CAAS_TEST_PREFIX"
    return 0
  fi
  if [[ "$1" == env && ( "$2" == create || "$2" == update ) ]]; then
    [[ "$*" == *--dry-run* ]] || touch "$CAAS_TEST_MARKER"
    return 0
  fi
  return 1
}
''')
        self.env = dict(os.environ, CAAS_CONDA_SH=str(self.hook), CAAS_TEST_LOG=str(self.log),
                        CAAS_TEST_PREFIX=str(self.prefix),CAAS_TEST_MARKER=str(self.marker),CAAS_TEST_BASE=str(self.base))
        for var in ('CAAS_CONDA_ENV','CAAS_CONDA_PREFIX','CAAS_VALIDATION_ROOT','SLURM_SUBMIT_DIR'):
            self.env.pop(var,None)

    def tearDown(self): self.tmp.cleanup()

    def run_script(self,name,*args,env=None):
        return subprocess.run(['bash',str(ROOT/name),*args],cwd=self.base,env=env or self.env,
                              capture_output=True,text=True)

    def test_create_uses_bundled_recipe_dedicated_name(self):
        p=self.run_script('create_conda_environment.sh')
        self.assertEqual(p.returncode,0,p.stderr)
        commands=self.log.read_text()
        self.assertIn('env create --name caas-validation-260917 --file '+str(ROOT/'environment.yml'),commands)
        self.assertNotIn('--force',commands)
        self.assertIn(str(ROOT/'scripts/check_environment.py'),p.stdout)

    def test_existing_environment_is_not_updated(self):
        self.marker.touch(); p=self.run_script('create_conda_environment.sh')
        self.assertEqual(p.returncode,0,p.stderr)
        self.assertNotIn('env create',self.log.read_text()); self.assertNotIn('env update',self.log.read_text())

    def test_update_requires_flag(self):
        self.marker.touch(); p=self.run_script('create_conda_environment.sh','--update')
        self.assertEqual(p.returncode,0,p.stderr)
        self.assertIn('env update --name caas-validation-260917',self.log.read_text())
        self.assertNotIn('--prune',self.log.read_text())

    def test_dry_run_does_not_create(self):
        p=self.run_script('create_conda_environment.sh','--dry-run')
        self.assertEqual(p.returncode,0,p.stderr); self.assertFalse(self.marker.exists())
        self.assertIn('--dry-run',self.log.read_text())

    def test_shared_phyloq_is_protected(self):
        p=self.run_script('create_conda_environment.sh','--name','phyloq','--update')
        self.assertNotEqual(p.returncode,0)
        self.assertFalse(self.log.exists())

    def test_base_prefix_is_protected(self):
        p=self.run_script('create_conda_environment.sh','--prefix',str(self.base),'--update')
        self.assertNotEqual(p.returncode,0)
        self.assertNotIn('env update',self.log.read_text())

    def test_missing_environment_blocks_launch(self):
        p=self.run_script('run_pipeline.sh','prepare')
        self.assertNotEqual(p.returncode,0)
        self.assertIn('Environment unavailable',p.stderr)
        self.assertNotIn('run_validation.py',p.stdout)

    def test_driver_passes_prefix_and_preserves_arguments(self):
        self.marker.touch(); p=self.run_script('run_pipeline.sh','smoke','--run-id','test-with-spaces','--selection','a path.tsv')
        self.assertEqual(p.returncode,0,p.stderr)
        self.assertIn(str(ROOT/'run_validation.py'),p.stdout)
        self.assertIn('\na path.tsv\n',p.stdout)
        self.assertIn('--conda-prefix\n'+str(self.prefix),p.stdout)
        self.assertIn('--conda-init\n'+str(self.hook),p.stdout)

    def test_slurm_spool_uses_submit_directory(self):
        self.marker.touch()
        spool=self.base/'slurm_script'; spool.write_bytes((ROOT/'submit_validation_slurm.sh').read_bytes())
        env=dict(self.env,SLURM_SUBMIT_DIR=str(ROOT))
        p=subprocess.run(['bash',str(spool),'deterministic','--run-id','never-submitted'],cwd=self.base,env=env,capture_output=True,text=True)
        self.assertEqual(p.returncode,0,p.stderr)
        self.assertIn(str(ROOT/'run_validation.py'),p.stdout)
        self.assertIn('--profile\ncluster',p.stdout)


if __name__=='__main__': unittest.main()
