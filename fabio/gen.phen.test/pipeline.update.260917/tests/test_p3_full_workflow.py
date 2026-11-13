"""Bounded full-workflow cache migration using real synthetic pooled CAAS."""
import os
from pathlib import Path
import shlex
import shutil
import subprocess
import sys
import tempfile
import unittest
ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/'scripts'))
from common import read_json,read_tsv,write_json


@unittest.skipUnless(os.environ.get('NF_TEST_COMMAND'),'Set NF_TEST_COMMAND for full Nextflow integration')
class P3FullWorkflow(unittest.TestCase):
    def test_old_five_then_add_p3_retains_all_existing_caas_cache(self):
        with tempfile.TemporaryDirectory(prefix='p3-full-') as folder:
            root=Path(folder)
            for name in ['scripts','inputs','bin']:
                shutil.copytree(ROOT/name,root/name,ignore=shutil.ignore_patterns('__pycache__','*.pyc','launches'))
            shutil.copy2(ROOT/'nextflow.config',root/'nextflow.config')
            # Fixed synthetic alignment paths are stable across all three runs.
            inventory=read_tsv(ROOT/'tests/fixtures/alignment_manifest.tsv')[:1]
            for row in inventory:
                target=root/'tests/fixtures/alignments'/Path(row['alignment_path']).name
                target.parent.mkdir(parents=True,exist_ok=True);shutil.copy2(row['alignment_path'],target)
                row['alignment_path']=str(target)
            from common import write_tsv
            write_tsv(root/'inventory.tsv',inventory,list(inventory[0]))
            design=root/'inputs/frozen/brain-260917-v1'
            from prepare_cluster_launch import select_rows,selected_modes
            rows=select_rows(design,selected_modes('all'))
            write_tsv(root/'selection.tsv',rows,list(rows[0]))
            old=(ROOT/'tests/fixtures/pre_p3_main.nf').read_bytes()
            (root/'main.nf').write_bytes(old)
            params=dict(stage='benchmark',direct_discovery=True,synthetic_inputs=True,
                        run_id='test',design=str(design),selection_manifest=str(root/'selection.tsv'),
                        alignment_manifest=str(root/'inventory.tsv'),settings=str(root/'inputs/settings.json'),
                        code_fingerprint='retained-reviewed-fingerprint',python_command=sys.executable,r_command='Rscript',
                        caastools_dir=str(root/'bin/caastools'),results_root=str(root/'results'),
                        summaries_root=str(root/'summaries'),launch_id='old',annotation_backend='disabled',
                        annotation_input=str(root/'inputs/annotations.disabled.txt'))
            write_json(root/'params.json',params)
            cmd=shlex.split(os.environ['NF_TEST_COMMAND'])+['run','main.nf','-params-file','params.json']
            env=dict(os.environ,NXF_ANSI_LOG='false',NXF_OFFLINE='true')
            def run(args):
                p=subprocess.run(cmd+args,cwd=root,env=env,capture_output=True,text=True,timeout=150)
                self.assertEqual(p.returncode,0,p.stdout+p.stderr)
                return p.stdout+p.stderr
            first=run(['-name','p3_full_test'])
            errors='\n'.join(p.read_text() for p in (root/'work').glob('*/*/.command.err'))
            errors+='\n'+'\n'.join(p.read_text() for p in (root/'work').glob('*/*/*/caastools.log'))
            self.assertNotIn('Error is ignored',first,errors)
            raw=list((root/'results/test').glob('*/000/raw/*/*.receipt.json'))
            self.assertEqual(len(raw),5,first)
            old_bytes={str(p.relative_to(root)):p.read_bytes() for p in raw}
            # Only workflow routing changes. Historical selection and every task
            # input stay byte-identical and at the same absolute location.
            (root/'main.nf').write_bytes((ROOT/'main.nf').read_bytes())
            params.update(additional_design=str(root/'inputs/additional/venditti-upper20-261113-v1'),launch_id='extended')
            write_json(root/'params.json',params)
            second=run(['-resume','p3_full_test'])
            self.assertEqual(second.count('Cached process > CAAS_POOLED'),5,second)
            self.assertIn('Cached process > VALIDATE_INPUTS',second)
            self.assertIn('VALIDATE_ADDITIONAL_INPUTS',second)
            self.assertNotIn('Error is ignored',second)
            self.assertEqual(len(list((root/'results/test').glob('*/000/raw/*/*.receipt.json'))),6,second)
            for path,data in old_bytes.items():self.assertEqual((root/path).read_bytes(),data)
            comparisons=read_tsv(root/'summaries/test/matched/strategy_comparison.tsv')
            self.assertEqual({x['strategy_id'] for x in comparisons},{'N3','N4','N5','P1','N6','P3'})
            params['launch_id']='repeat';write_json(root/'params.json',params)
            third=run(['-resume','p3_full_test'])
            self.assertEqual(third.count('Cached process > CAAS_POOLED'),6,third)


if __name__=='__main__':unittest.main()
