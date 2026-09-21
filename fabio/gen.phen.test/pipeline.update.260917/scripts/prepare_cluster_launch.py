#!/usr/bin/env python3
"""Freeze an explicit multi-strategy CAAS launch; never submit a job.

All new modes means N3/N4/N5 and P1/N6. Randomized null series are disabled.
Historical references and P2 are opt-in. Existing launch inputs are reused only
when the request and their checksums match; no silent regeneration on resume.
"""
from __future__ import annotations
import argparse
import tempfile
from collections import Counter
from pathlib import Path
from common import ROOT, checked, identifier, read_json, read_tsv, require, resolve, sha, tree_hash, verify_design, write_json, write_tsv
from inventory_alignments import inventory

MODES = {
    'deterministic': ('N3', 'N4', 'N5'),
    'paired': ('P1', 'N6'),
    'references': ('P0', 'N0', 'N1', 'N2'),
    'auxiliary-p2': ('P2',),
}
NEW_MODES = ('deterministic', 'paired')


def selected_modes(text, include_references=False, include_p2=False):
    names = list(NEW_MODES) if text == 'all' else text.split(',')
    require(names and all(n in MODES for n in names), 'Unknown --modes. Use all or: ' + ','.join(MODES))
    require(len(names) == len(set(names)), 'Repeated mode in --modes')
    if include_references and 'references' not in names: names.append('references')
    if include_p2 and 'auxiliary-p2' not in names: names.append('auxiliary-p2')
    return [name for name in MODES if name in names]


def select_rows(design, modes):
    strategies = {s for mode in modes for s in MODES[mode]}
    rows = [r for r in read_tsv(design/'execution_manifest.tsv') if r['strategy_id'] in strategies]
    require(rows and {r['strategy_id'] for r in rows} == strategies, 'Requested strategies missing from frozen design')
    require(len({r['hypothesis_id'] for r in rows}) == len(rows), 'Duplicate frozen hypothesis')
    return rows


def prepare(design, output, run_id, modes, pattern=None, inventory_path=None, resume=False):
    identifier(run_id)
    require(bool(pattern) != bool(inventory_path), 'Provide one alignment pattern OR an inventory')
    lock = verify_design(design)
    tool_hash = tree_hash(ROOT/'bin/caastools')
    require(tool_hash == lock['caastools_sha256'], 'Bundled modified CAAStools differs from the frozen version')
    rows = select_rows(design, modes)
    request = dict(run_id=run_id, modes=modes, design_sha256=sha(design/'design.lock.json'),
                   caastools_sha256=tool_hash, alignment_pattern=str(Path(pattern).resolve()) if pattern else None,
                   inventory_path=str(inventory_path.resolve()) if inventory_path else None,
                   inventory_sha256=sha(inventory_path) if inventory_path else None)
    if output.exists():
        require((output/'launch.request.json').is_file(), 'Unrecognized launch-input directory: ' + str(output))
        require(read_json(output/'launch.request.json') == request,
                'Launch modes or alignment source changed. Use a new --run-id; existing inputs are not overwritten.')
        saved = read_json(output/'launch.inputs.lock.json')
        for name, checksum in saved.items(): checked(output/name, checksum)
        return read_json(output/'launch_summary.json')
    require(not resume, 'No saved launch inputs for this run; cannot resume')
    output.parent.mkdir(parents=True, exist_ok=True)
    # Build elsewhere, then publish the complete input bundle in one rename.
    with tempfile.TemporaryDirectory(prefix='.building-', dir=output.parent) as tmp:
        tmp = Path(tmp)
        write_tsv(tmp/'selection.tsv', rows, list(rows[0]))
        if inventory_path:
            alignments = read_tsv(inventory_path)
            require(bool(alignments), 'Empty alignment inventory')
            for row in alignments:
                row['alignment_path'] = str(resolve(inventory_path.parent, row['alignment_path']))
            write_tsv(tmp/'alignment_manifest.tsv', alignments, list(alignments[0]))
        else:
            alignments = inventory(request['alignment_pattern'], tmp/'alignment_manifest.tsv')
        report = dict(run_id=run_id, modes=modes, hypotheses=len(rows), alignments=len(alignments),
                      expected_caas_tasks=len(rows)*len(alignments),
                      cycles=sum(int(r['selected_cycles']) for r in rows),
                      hypotheses_by_strategy=dict(Counter(r['strategy_id'] for r in rows)),
                      direct_discovery=True, pilot_required=False, caastools_sha256=tool_hash,
                      phylogeny_path=str((design/'sources/tree.nwk').resolve()),
                      phylogeny_sha256=sha(design/'sources/tree.nwk'),
                      phylogeny_usage='Bundled diagnostic tree; pooled CAAS does not take a gene-tree input',
                      historical_references_included='references' in modes, auxiliary_p2_included='auxiliary-p2' in modes,
                      results_imported=False, enrichment='not_requested_discovery_only')
        write_json(tmp/'launch.request.json', request)
        write_json(tmp/'launch_summary.json', report)
        write_json(tmp/'launch.inputs.lock.json', {name:sha(tmp/name) for name in
                   ('selection.tsv','alignment_manifest.tsv','launch.request.json','launch_summary.json')})
        tmp.rename(output)
    return report


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--design', type=Path, default=ROOT/'inputs/frozen/brain-260917-v1')
    p.add_argument('--run-id', required=True)
    p.add_argument('--output', type=Path)
    p.add_argument('--modes', default='all')
    p.add_argument('--include-references', action='store_true')
    p.add_argument('--include-p2', action='store_true')
    source = p.add_mutually_exclusive_group(required=True)
    source.add_argument('--alignments-pattern')
    source.add_argument('--inventory', type=Path)
    p.add_argument('--resume', action='store_true')
    a = p.parse_args()
    identifier(a.run_id)
    modes = selected_modes(a.modes, a.include_references, a.include_p2)
    output = (a.output or ROOT/'inputs/launches'/a.run_id).resolve()
    report = prepare(a.design.resolve(), output, a.run_id, modes, a.alignments_pattern,
                     a.inventory.resolve() if a.inventory else None, a.resume)
    print('Modified CAAStools verified:', report['caastools_sha256'])
    print('Selected modes:', ', '.join(modes))
    print(f"{report['hypotheses']} hypotheses x {report['alignments']} alignments = {report['expected_caas_tasks']} CAAS tasks")
    print('Cycles across hypotheses:', report['cycles'])
    print('Frozen launch inputs:', output)
    print('Bundled diagnostic phylogeny:', report.get('phylogeny_path', str(a.design.resolve()/'sources/tree.nwk')))
    print('No pilot, approval JSON, enrichment or downstream statistical prerequisites.')


if __name__ == '__main__': main()
