import json
from pathlib import Path
import shlex
import sys
import tempfile
import unittest
ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/'scripts'));sys.path.insert(0,str(ROOT/'recovery'))
from common import read_json,read_tsv,sha,write_json,write_tsv
from recover_completed import parse_log,recover
import test_validation


class RecoveryTests(unittest.TestCase):
    def fixture(self,root):
        raw=test_validation.FilteringTests().make_results(root,[])
        row=dict(hypothesis_id='TEST_000',strategy_id='TEST',replicate_id='000',null_family='observed',
                 config_sha256='frozen',pool_sha256='frozen')
        write_tsv(root/'selection.tsv',[row],list(row))
        write_tsv(root/'inventory.tsv',[dict(gene_id='TEST',sha256='synthetic-sha'),dict(gene_id='MISSING',sha256='x')],['gene_id','sha256'])
        write_json(root/'settings.json',read_json(ROOT/'inputs/settings.json'))
        write_json(root/'design/design.lock.json',dict(caastools_sha256='frozen'))
        lock=dict(selection_sha256=sha(root/'selection.tsv'),alignment_manifest_sha256=sha(root/'inventory.tsv'),
                  settings_sha256=sha(root/'settings.json'),design_sha256=sha(root/'design/design.lock.json'),
                  work_dir=str(root/'work'),results_root=str(root/'results'))
        write_json(root/'results/run/launch.lock.json',lock)
        params=dict(run_id='run',direct_discovery=True,synthetic_inputs=True,selection_manifest=str(root/'selection.tsv'),
                    alignment_manifest=str(root/'inventory.tsv'),settings=str(root/'settings.json'),design=str(root/'design'))
        write_json(root/'params.json',params)
        receipt=read_json(raw/'TEST.receipt.json')
        receipt.update(run_id='run',settings_sha256=lock['settings_sha256'],design_sha256=lock['design_sha256'],selection_sha256=lock['selection_sha256'])
        write_json(raw/'TEST.receipt.json',receipt)
        published=root/'results/run/TEST/000/raw/TEST';published.parent.mkdir(parents=True)
        raw.rename(published)
        log=root/'run.log'
        log.write_text('$> nextflow run main.nf -params-file '+shlex.quote(str(root/'params.json'))+'\n'
                       '[aa/bbbbbb] Cached process > CAAS_POOLED (TEST/000:TEST)\n'
                       '[aa/cccccc] NOTE: Error submitting process \'CAAS_POOLED (TEST/000:MISSING)\' for execution -- Error is ignored\n')
        return log,published

    def test_recovery_separate_outputs_and_missing_report(self):
        with tempfile.TemporaryDirectory() as t:
            root=Path(t);log,pub=self.fixture(root)
            before=sha(pub/'TEST.receipt.json')
            out=recover(root,'run','rec',log,tables_only=True)
            self.assertTrue((out/'recovery.complete.json').is_file())
            self.assertEqual(sha(pub/'TEST.receipt.json'),before)
            self.assertEqual(read_tsv(out/'matched/non_completed_genes.tsv')[0]['reason'],'submission_failed')
            self.assertEqual(read_json(out/'hypotheses/TEST_000/summary.json')['completed_genes'],1)
            with self.assertRaisesRegex(ValueError,'already exists'):
                recover(root,'run','rec',log,tables_only=True)

    def test_recovery_includes_additive_p3_with_its_own_selection_digest(self):
        from check_additional_design import additional_identity
        additional=ROOT/'inputs/additional/venditti-upper20-261113-v1'
        design=ROOT/'inputs/frozen/brain-260917-v1'
        with tempfile.TemporaryDirectory() as t:
            root=Path(t);log,pub=self.fixture(root)
            params=read_json(root/'params.json');params.update(design=str(design),additional_design=str(additional))
            write_json(root/'params.json',params)
            lock=read_json(root/'results/run/launch.lock.json');lock['design_sha256']=sha(design/'design.lock.json')
            write_json(root/'results/run/launch.lock.json',lock)
            receipt=read_json(pub/'TEST.receipt.json');receipt['design_sha256']=lock['design_sha256']
            receipt['tool_sha256']=read_json(design/'design.lock.json')['caastools_sha256']
            write_json(pub/'TEST.receipt.json',receipt)
            import shutil
            p3=root/'results/run/P3/000/raw/TEST';p3.parent.mkdir(parents=True);shutil.copytree(pub,p3)
            row=read_tsv(additional/'selection.tsv')[0]
            p3_receipt=dict(receipt,**row,selection_sha256=sha(additional/'selection.tsv'))
            write_json(p3/'TEST.receipt.json',p3_receipt)
            write_json(root/'results/run/additional_hypotheses.lock.json',additional_identity(design,additional))
            with log.open('a') as f:
                f.write('[ab/123456] Cached process > CAAS_POOLED (P3/000:TEST)\n')
                f.write("[ab/654321] NOTE: Process `CAAS_POOLED (P3/000:MISSING)` terminated with an error exit status (1) -- Error is ignored\n")
            out=recover(root,'run','p3-recovery',log,tables_only=True)
            self.assertEqual({r['strategy_id'] for r in read_tsv(out/'matched/strategy_comparison.tsv')},{'TEST','P3'})
            self.assertEqual(read_json(out/'hypotheses/P3_000/summary.json')['completed_genes'],1)

    def test_wrong_receipt_run_fails(self):
        with tempfile.TemporaryDirectory() as t:
            root=Path(t);log,pub=self.fixture(root)
            r=read_json(pub/'TEST.receipt.json');r['run_id']='another-run';write_json(pub/'TEST.receipt.json',r)
            with self.assertRaisesRegex(ValueError,'provenance mismatch run_id'):
                recover(root,'run','rec',log,tables_only=True)

    def test_truncated_log_fails_before_creating_output(self):
        with tempfile.TemporaryDirectory() as t:
            root=Path(t);log,pub=self.fixture(root)
            log.write_text('\n'.join(log.read_text().splitlines()[:2])+'\n')
            with self.assertRaisesRegex(ValueError,'exactly the selected'):
                recover(root,'run','rec',log,tables_only=True)
            self.assertFalse((root/'summaries').exists())

    def test_corrupt_output_fails(self):
        with tempfile.TemporaryDirectory() as t:
            root=Path(t);log,pub=self.fixture(root)
            (pub/'TEST.pooled.caas.tsv').write_text('corrupt')
            with self.assertRaisesRegex(ValueError,'SHA-256 mismatch'):
                recover(root,'run','rec',log,tables_only=True)


if __name__=='__main__':unittest.main()
