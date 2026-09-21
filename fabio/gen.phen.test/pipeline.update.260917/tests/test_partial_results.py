"""Missing genes remain explicit missing observations, never negative results."""
import tempfile
from pathlib import Path
import sys
import unittest
ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT/'scripts'))
from common import read_json, read_tsv, write_json, write_tsv, sha
from summarize_hypothesis import summarize
from common_background import build
from test_validation import FilteringTests


class PartialResults(unittest.TestCase):
    def fixture(self, root):
        raw = FilteringTests().make_results(root, [])
        anchor = read_json(raw/'TEST.receipt.json')
        write_json(raw/'TEST_000__TEST.json', anchor)
        write_tsv(root/'inventory.tsv', [dict(gene_id='TEST', sha256='synthetic-sha'),
                                       dict(gene_id='MISSING', sha256='x')], ['gene_id','sha256'])
        return raw

    def test_missing_is_not_a_negative_gene(self):
        with tempfile.TemporaryDirectory() as folder:
            root=Path(folder); raw=self.fixture(root)
            result=summarize('TEST_000',root/'inventory.tsv',raw,root/'summary',True)
            self.assertEqual((result['status'],result['completed_genes'],result['non_completed_genes']),('partial',1,1))
            self.assertEqual(read_tsv(root/'summary/background.tsv'),[{'gene':'TEST'}])
            self.assertEqual(read_tsv(root/'summary/non_completed_genes.tsv')[0]['gene'],'MISSING')
            self.assertEqual(len(read_tsv(root/'summary/gene_inventory.tsv')),1)

    def test_all_missing_still_summarized_and_cohort_reported(self):
        with tempfile.TemporaryDirectory() as folder:
            root=Path(folder); raw=self.fixture(root)
            (raw/'TEST.receipt.json').rename(root/'saved.receipt.json')
            result=summarize('TEST_000',root/'inventory.tsv',raw,root/'summary',True)
            self.assertEqual((result['status'],result['completed_genes'],result['non_completed_genes']),('partial',0,2))
            self.assertEqual(read_tsv(root/'summary/background.tsv'),[])
            write_tsv(root/'cohort.tsv',[dict(hypothesis_id='TEST_000',summary_dir=str(root/'summary'),
                                            summary_lock_sha256=sha(root/'summary/summary.lock.json'))],
                      ['hypothesis_id','summary_dir','summary_lock_sha256'])
            build(root/'cohort.tsv',root/'matched',allow_incomplete=True)
            self.assertEqual(len(read_tsv(root/'matched/non_completed_genes.tsv')),2)
            self.assertEqual(read_json(root/'matched/cohort.json')['execution_status'],'partial')
            self.assertEqual(read_tsv(root/'matched/common_background.tsv'),[])

    def test_default_still_rejects_incomplete(self):
        with tempfile.TemporaryDirectory() as folder:
            root=Path(folder); raw=self.fixture(root)
            with self.assertRaisesRegex(ValueError,'Incomplete hypothesis'):
                summarize('TEST_000',root/'inventory.tsv',raw,root/'summary')

    def test_corrupt_success_still_rejected(self):
        with tempfile.TemporaryDirectory() as folder:
            root=Path(folder); raw=self.fixture(root)
            (raw/'TEST.pooled.caas.events.tsv').write_text('corrupted')
            with self.assertRaisesRegex(ValueError,'SHA-256 mismatch'):
                summarize('TEST_000',root/'inventory.tsv',raw,root/'summary',True)


if __name__ == '__main__': unittest.main()
