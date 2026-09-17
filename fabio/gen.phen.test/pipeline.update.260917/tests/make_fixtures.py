#!/usr/bin/env python3
"""Generate two explicitly synthetic alignments; never reuse biological outputs."""
import sys
from pathlib import Path
ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / 'scripts'))
from common import read_json, read_tsv, write_json, write_tsv, sha
from inventory_alignments import inventory


def main():
    design = ROOT / 'inputs/frozen/brain-260917-v1'
    folder = ROOT / 'tests/fixtures'
    alignments = folder / 'alignments'; alignments.mkdir(parents=True, exist_ok=True)
    species = sorted(r['species'] for r in read_tsv(design / 'sources/phenotype.tsv'))
    members = read_tsv(design / 'pool_membership.tsv')
    pools = {side: {r['species'] for r in members if r['hypothesis_id'] == 'P0_000' and r['side'] == side} for side in ('FG', 'BG')}
    for gene in ('POSITIVE', 'ZERO'):
        lines = [f'{len(species)} 20']
        for s in species:
            seq = list('A' * 20)
            if gene == 'POSITIVE':
                for pos in (0, 1, 2, 15): seq[pos] = 'A' if s in pools['FG'] else 'V' if s in pools['BG'] else 'G'
            lines.append(s + '  ' + ''.join(seq))
        (alignments / f'{gene}.phy').write_text('\n'.join(lines) + '\n')
    (folder / 'MALFORMED.phy').write_text('2 20\nUnknown_species ACDE\nOther_species VVVV\n')
    inventory(str(alignments / '*.phy'), folder / 'alignment_manifest.tsv')
    execution = read_tsv(design / 'execution_manifest.tsv')
    ids = {'P0_000','N2_000','N3_000','N4_000','N5_000','P1_000','N6_000','R0_001','R1_001'}
    write_tsv(folder / 'selection.tsv', [r for r in execution if r['hypothesis_id'] in ids], list(execution[0]))
    # Small, visibly artificial annotation universe, solely for local tests.
    ann = folder / 'toy_annotations'
    write_tsv(ann / 'term_membership.tsv', [dict(term_id='TOY:A', source='TOY', term_name='Synthetic positive', gene='POSITIVE'),
        dict(term_id='TOY:B', source='TOY', term_name='Synthetic zero', gene='ZERO')], ['term_id','source','term_name','gene'])
    write_json(ann / 'metadata.json', dict(annotation_version='toy-v1', organism='hsapiens', sources=['TOY'], toy_only=True,
                                         term_membership_sha256=sha(ann / 'term_membership.tsv')))
    cfg = read_json(ROOT / 'inputs/settings.json')
    cfg['enrichment'].update(backend='local-hypergeom-bonferroni', annotation_version='toy-v1', sources=['TOY'],
                             correction='bonferroni', implementation_approved=True)
    cfg['statistics'].update(approved=True, primary_metric='query_fraction', term_universe=['TOY:A','TOY:B'],
                             term_universe_status='synthetic_fixture_only', multiple_testing='holm', gene_list_replicates=5,
                             null_families=['R0','R1','G0'],expected_null_counts={'R0':1,'R1':1,'G0':5})
    write_json(folder / 'toy_settings.json', cfg)
    print(f'Created two synthetic {len(species)}-species alignments and a nine-hypothesis smoke selection')


if __name__ == '__main__': main()
