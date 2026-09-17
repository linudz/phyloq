#!/usr/bin/env python3
"""Independent acceptance checks; run with python3 -m unittest discover -s tests."""
import copy
import itertools
import math
import subprocess
import sys
import tempfile
import unittest
from unittest.mock import patch
from pathlib import Path
ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / 'scripts'))
from common import checked, cycle_read, digest, pool_read, read_json, read_tsv, sha, write_json, write_tsv
from check_design import check
from prepare_validation import random_assignments, generate_cycles
from inventory_alignments import inventory
from run_caas import coverage
from summarize_hypothesis import summarize
from functional import gene_lists, local_enrich, metrics, conditional_tail
from legacy_support import dense_cluster_positions, preferred_position_record

DESIGN = ROOT / 'inputs/frozen/brain-260917-v1'


class DesignTests(unittest.TestCase):
    def test_full_frozen_design(self):
        report = check(DESIGN)
        self.assertEqual(report['hypotheses'], 208)
        self.assertEqual(report['r0_assignments'], 4480)
        self.assertEqual(report['possible_9x10_cycles'], 26460)

    def test_stable_extension(self):
        expected = read_json(ROOT / 'inputs/expected_pools.json')
        p0 = {k: v.split() for k, v in expected['P0'].items()}
        eligible = read_tsv(DESIGN / 'eligible_species.tsv')
        for family in ('R0','R1'):
            universe = [r['species'] for r in eligible if r['randomization'] == family]
            first = random_assignments(universe, expected['quotas'], p0, 19, 260917, family)
            extended = random_assignments(universe, expected['quotas'], p0, 99, 260917, family)
            self.assertEqual(first, extended[:19])
            for i, (a,b) in enumerate(zip(first, extended), 1):
                self.assertEqual(generate_cycles(a, [], 100, 260917, f'{family}_{i:03d}'),
                                 generate_cycles(b, [], 100, 260917, f'{family}_{i:03d}'))

    def test_repeated_preparation(self):
        # An existing frozen design is verified, never silently regenerated.
        before = sha(DESIGN / 'design.lock.json')
        subprocess.run([sys.executable, str(ROOT / 'scripts/prepare_validation.py')], check=True, capture_output=True)
        self.assertEqual(before, sha(DESIGN / 'design.lock.json'))

    def test_actual_pilot_configs_byte_identical_to_extension(self):
        with tempfile.TemporaryDirectory() as tmp:
            pilot=Path(tmp)/'pilot'
            subprocess.run([sys.executable,str(ROOT/'scripts/prepare_validation.py'),'--output',str(pilot),
                            '--design-id','unit-pilot','--replicates','19'],check=True,capture_output=True)
            for family in ('R0','R1'):
                for i in range(1,20):
                    for filename in ('pool.cfg','cycles.cfg'):
                        rel=Path('configs')/f'{family}_{i:03d}'/filename
                        self.assertEqual((pilot/rel).read_bytes(),(DESIGN/rel).read_bytes())

    def test_unknown_phenotype_identifier_rejected(self):
        def edited(path):
            rows=read_tsv(path)
            return [r for r in rows if r['species']!='Macaca_nigra'] if Path(path).name=='phenotype.tsv' else rows
        with patch('check_design.read_tsv',side_effect=edited):
            with self.assertRaisesRegex(ValueError,'Unknown species'): check(DESIGN)

    def test_historical_configs_are_byte_exact(self):
        for r in read_tsv(DESIGN / 'execution_manifest.tsv'):
            if not r['historical_source']: continue
            original = next(x for x in read_tsv(DESIGN / 'sources/benchmark.configs.tsv') if x['approach'] == r['historical_source'])
            self.assertEqual(sha(DESIGN / r['hypotheses_config']), sha(DESIGN/'sources'/original['hypotheses_config']))

    def test_source_change_detected(self):
        with tempfile.TemporaryDirectory() as tmp:
            p = Path(tmp) / 'source'; p.write_text('old'); old = sha(p); p.write_text('new')
            with self.assertRaisesRegex(ValueError, 'SHA-256 mismatch'): checked(p, old)

    def test_malformed_adapters_and_pair_counts(self):
        with tempfile.TemporaryDirectory() as tmp:
            p = Path(tmp) / 'pool'; p.write_text('Macaca_nigra\t1\nMacaca_nigra\t0\n')
            with self.assertRaisesRegex(ValueError, 'conflicting'): pool_read(p)
            p.write_text('FG\tSpecies\n')
            with self.assertRaisesRegex(ValueError, 'Invalid pool'): pool_read(p)
        with self.assertRaisesRegex(ValueError, 'Insufficient linked pairs'):
            generate_cycles({}, [('g','a','b')]*3, 15, 1, 'too-few')

    def test_duplicate_cycles(self):
        with tempfile.TemporaryDirectory() as tmp:
            p = Path(tmp) / 'cycles'
            p.write_text('b1\ta,b,c,d\te,f,g,h\nb2\td,b,c,a\th,g,f,e\n')
            with self.assertRaisesRegex(ValueError, 'Duplicate cycle'):
                cycle_read(p, {'FG':list('abcd'),'BG':list('efgh')})


class AlignmentTests(unittest.TestCase):
    def test_duplicate_gene(self):
        fixture = ROOT / 'tests/fixtures/alignments/ZERO.phy'
        with tempfile.TemporaryDirectory() as tmp:
            tmp = Path(tmp)
            for name in ('SAME.a.phy','SAME.b.phy'): (tmp / name).write_bytes(fixture.read_bytes())
            with self.assertRaisesRegex(ValueError, 'Duplicate/ambiguous gene'): inventory(str(tmp/'*.phy'), tmp/'out.tsv')

    def test_malformed_alignment(self):
        with tempfile.TemporaryDirectory() as tmp:
            with self.assertRaises(ValueError): inventory(str(ROOT/'tests/fixtures/MALFORMED.phy'), Path(tmp)/'out.tsv')

    def test_zero_gene_is_coverage_testable(self):
        d = DESIGN / 'configs/P0_000'
        pools = pool_read(d/'pool.cfg'); cycles = cycle_read(d/'cycles.cfg', pools)
        cov = coverage(ROOT/'tests/fixtures/alignments/ZERO.phy', pools, cycles, read_json(ROOT/'inputs/settings.json')['discovery'])
        self.assertTrue(cov['coverage_testable'])
        self.assertEqual(cov['coverage_testable_positions'], 20)
        # A known missing alignment species is a coverage issue, not an unknown pool ID.
        pools['FG'].append('Known_but_absent')
        cov = coverage(ROOT/'tests/fixtures/alignments/ZERO.phy', pools, cycles, read_json(ROOT/'inputs/settings.json')['discovery'])
        self.assertIn('Known_but_absent', cov['missing_fg'])


class FilteringTests(unittest.TestCase):
    def test_dense_intervals(self):
        self.assertEqual(dense_cluster_positions([0,1,2,15], .7, 3, 3), {0,1,2})

    def make_results(self, tmp, supports):
        inv = [dict(gene_id='TEST', sha256='synthetic-sha')]
        write_tsv(tmp/'inventory.tsv', inv, ['gene_id','sha256'])
        raw = tmp/'raw'; raw.mkdir()
        events = [dict(gene='TEST', position=i*20, event_id=f'e_{i}', primary_event='yes', positional_pvalue=.001,
                       fg_support_count=f, bg_support_count=b) for i,(f,b) in enumerate(supports)]
        write_tsv(raw/'TEST.pooled.caas.events.tsv', events, ['gene','position','event_id','primary_event','positional_pvalue','fg_support_count','bg_support_count'])
        (raw/'TEST.pooled.caas.tsv').write_text('')
        receipt = dict(hypothesis_id='TEST_000', gene_id='TEST', status='success', alignment_sha256='synthetic-sha',
                       settings=read_json(ROOT/'inputs/settings.json'), strategy_id='TEST', replicate_id='000', null_family='observed',
                       settings_sha256='frozen', config_sha256='frozen', pool_sha256='frozen', tool_sha256='frozen', design_sha256='frozen',
                       coverage=dict(coverage_testable=True, coverage_testable_positions=40, present_fg=10,present_bg=10,missing_fg=[],missing_bg=[]),
                       elapsed_seconds=1,cpu_seconds=1,max_rss_bytes=1,output_bytes=1,smoke=True,
                       events_sha256=sha(raw/'TEST.pooled.caas.events.tsv'), legacy_sha256=sha(raw/'TEST.pooled.caas.tsv'))
        write_json(raw/'TEST.receipt.json',receipt)
        return raw

    def test_same_position_not_cross_position_maxima(self):
        with tempfile.TemporaryDirectory() as d:
            d=Path(d); raw=self.make_results(d,[(5,3),(3,5)])
            s=summarize('TEST_000',d/'inventory.tsv',raw,d/'summary')
            self.assertEqual(s['query_genes'],0); self.assertEqual(s['background_genes'],1)
        with tempfile.TemporaryDirectory() as d:
            d=Path(d); raw=self.make_results(d,[(4,4)])
            self.assertEqual(summarize('TEST_000',d/'inventory.tsv',raw,d/'summary')['query_genes'],1)

    def test_missing_output_blocks(self):
        with tempfile.TemporaryDirectory() as d:
            d=Path(d); raw=self.make_results(d,[])
            write_tsv(d/'inventory.tsv',[dict(gene_id='TEST',sha256='synthetic-sha'),dict(gene_id='MISSING',sha256='x')],['gene_id','sha256'])
            with self.assertRaisesRegex(ValueError,'Incomplete hypothesis'):
                summarize('TEST_000',d/'inventory.tsv',raw,d/'summary')
            self.assertEqual(read_json(d/'summary/execution_completeness.json')['missing'],['MISSING'])

    def test_zero_result_background(self):
        with tempfile.TemporaryDirectory() as d:
            d=Path(d); raw=self.make_results(d,[])
            s=summarize('TEST_000',d/'inventory.tsv',raw,d/'summary')
            self.assertEqual((s['query_genes'],s['background_genes']),(0,1))


class FunctionalTests(unittest.TestCase):
    def test_controls_size_reproducible_and_incompatible(self):
        qs={'P0_000':{'a','b'},'N2_000':{'b','c','d'},'N3_000':{'a'}}; bg=set('abcde')
        controls,status=gene_lists(qs,bg,1000,260917)
        self.assertEqual((controls,status),gene_lists(qs,bg,1000,260917))
        self.assertEqual(len(controls),2000)
        self.assertTrue(all(len(x)==2 for x in controls.values()))
        self.assertEqual(status[2]['status'],'eligible_query_too_small')
        self.assertTrue(all(x <= qs['N2_000'] for k,x in controls.items() if k.startswith('G1_')))

    def test_exact_local_enrichment(self):
        terms={'one':dict(genes=set('ab'),source='TOY',term_name='one'),'two':dict(genes=set('cd'),source='TOY',term_name='two')}
        rows=local_enrich('P0',set('ab'),set('abcd'),terms,{'threshold':.05})
        self.assertAlmostEqual(rows[0]['p_value'],1/6)
        self.assertAlmostEqual(rows[0]['p_value_adjusted'],1/3)
        self.assertEqual(rows[0]['fold_enrichment'],2)
        self.assertFalse(rows[0]['significant'])
        stat=metrics('P0',rows,set('ab'),set('abcd'),{'term_universe':['one'],'term_groups':{}})
        self.assertEqual(stat['mean_fold_enrichment'],2) # nonsignificant is not zero

    def test_conditional_tail(self):
        s=conditional_tail(5,[0]*99)
        self.assertEqual(s['conditional_tail_fraction'],.01)
        self.assertEqual(conditional_tail(5,[0]*19)['conditional_tail_fraction'],.05)
        self.assertGreater(s['monte_carlo_exceedance_ci_upper'],0)


if __name__ == '__main__': unittest.main()
