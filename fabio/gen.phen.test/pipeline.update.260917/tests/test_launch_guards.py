"""Launch failures must be visible and must occur before Nextflow submission."""
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/'scripts'))
from common import read_json, read_tsv, write_json, write_tsv
from validate_inputs import validate
from functional import cached_enrich, payload


class LaunchGuards(unittest.TestCase):
    def test_missing_production_approval(self):
        p=subprocess.run([sys.executable,str(ROOT/'run_validation.py'),'deterministic','--run-id','unit-must-not-launch',
                          '--alignments',str(ROOT/'tests/fixtures/alignment_manifest.tsv'),'--execute'],capture_output=True,text=True)
        self.assertNotEqual(p.returncode,0)
        self.assertIn('No production launch without explicit reviewed approvals',p.stderr)
        self.assertFalse((ROOT/'results/unit-must-not-launch').exists())

    def test_nonexistent_resume(self):
        p=subprocess.run([sys.executable,str(ROOT/'run_validation.py'),'smoke','--run-id','unit-never-launched',
                          '--alignments',str(ROOT/'tests/fixtures/alignment_manifest.tsv'),'--resume'],capture_output=True,text=True)
        self.assertNotEqual(p.returncode,0)
        self.assertIn('Cannot resume a run with no launch lock',p.stderr)

    def test_duplicate_inventory_rejected_before_jobs(self):
        design=ROOT/'inputs/frozen/brain-260917-v1'
        rows=read_tsv(ROOT/'tests/fixtures/alignment_manifest.tsv')
        with tempfile.TemporaryDirectory() as tmp:
            tmp=Path(tmp); write_tsv(tmp/'inventory.tsv',[rows[0],rows[0]],list(rows[0]))
            with self.assertRaisesRegex(ValueError,'Duplicate gene IDs'):
                validate(design,design/'selections/smoke.tsv',tmp/'inventory.tsv',ROOT/'inputs/settings.json',
                         'unit-invalid',tmp/'jobs',smoke=True)
            self.assertFalse((tmp/'jobs').exists())

    def test_missing_cache_is_not_zero_enrichment(self):
        cfg=read_json(ROOT/'inputs/settings.json')['enrichment']; cfg['annotation_version']='test-version'
        with tempfile.TemporaryDirectory() as tmp:
            with self.assertRaisesRegex(ValueError,'Missing frozen g:Profiler cache'):
                cached_enrich('P0',payload({'A'},{'A','B'},cfg),Path(tmp),cfg)


if __name__=='__main__': unittest.main()
