"""Preparation must work with no original project tree available."""
import shutil
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path

ROOT=Path(__file__).resolve().parents[1]


class PortablePreparation(unittest.TestCase):
    def test_bundled_sources_are_sufficient(self):
        with tempfile.TemporaryDirectory(prefix='caas-portable-') as tmp:
            tmp=Path(tmp); bundle=tmp/'standalone'
            for name in ('scripts','inputs','bin'):
                shutil.copytree(ROOT/name,bundle/name,ignore=shutil.ignore_patterns('__pycache__','*.pyc','alignments','annotation-cache'))
            shutil.copyfile(ROOT/'environment.yml',bundle/'environment.yml')
            missing=tmp/'original-project-does-not-exist'
            command=[sys.executable,'-B',str(bundle/'scripts/prepare_validation.py'),
                     '--pipeline',str(missing),'--phenotype',str(missing/'trait.tsv'),
                     '--taxonomy',str(missing/'taxonomy.tsv'),'--tree',str(missing/'tree.nwk'),
                     '--output',str(tmp/'pilot'),'--design-id','portable-pilot','--replicates','19']
            result=subprocess.run(command,cwd=tmp,capture_output=True,text=True)
            self.assertEqual(result.returncode,0,result.stderr)
            for family in ('R0','R1'):
                for i in range(1,20):
                    for filename in ('pool.cfg','cycles.cfg'):
                        rel=Path('configs')/f'{family}_{i:03d}'/filename
                        self.assertEqual((tmp/'pilot'/rel).read_bytes(),(bundle/'inputs/frozen/brain-260917-v1'/rel).read_bytes())


if __name__=='__main__': unittest.main()
