#!/usr/bin/env python3
"""Read-only historical inventory. Presence alone never certifies compatibility."""
import argparse
from pathlib import Path
from common import ROOT, read_tsv, sha, write_json, write_tsv
from prepare_validation import HISTORICAL


def audit(design, historical, output):
    rows=[]
    configs={r['strategy_id']:r for r in read_tsv(design/'execution_manifest.tsv') if r['historical_source']}
    for sid, old in HISTORICAL.items():
        folder=historical/old
        files=sorted(p for p in folder.rglob('*') if p.is_file()) if folder.exists() else []
        events=[p for p in files if p.name.endswith('.pooled.caas.events.tsv')]
        legacy=[p for p in files if p.name.endswith('.pooled.caas.tsv')]
        event_set, legacy_set = set(events), set(legacy)
        provenance=[p for p in files if p.suffix=='.json' or 'trace' in p.name or p.name.endswith('.log')]
        r=configs[sid]
        rows.append(dict(strategy_id=sid,historical_source=old,source_path=str(folder.resolve()),
            event_files=len(events),legacy_files=len(legacy),provenance_candidates=len(provenance),
            config_sha256=r['config_sha256'],pool_sha256=r['pool_sha256'],
            reference_status='not_imported_requires_provenance_review' if folder.exists() else 'historical_directory_absent',
            unresolved='alignment collection/checksums; complete successful-job inventory; runtime CAAStools version; full resolved discovery command',
            proposed_action='Fresh isolated reference-rerun only after approval, or provide a verified import bundle'))
        # Preserve exact candidate paths, without reading or rewriting historical results.
        write_tsv(output/f'{sid}.historical_inventory.tsv',
                  [dict(path=str(p.resolve()),bytes=p.stat().st_size,kind='events' if p in event_set else 'legacy' if p in legacy_set else 'other') for p in files],
                  ['path','bytes','kind'])
    write_tsv(output/'reference_compatibility.tsv',rows,list(rows[0]))
    write_json(output/'reference_audit.json',dict(imported_references=[],
        inventory_scope='Metadata inventory, not an unearned checksum certification of historical execution',
        original_config_checksums_preserved=True, references=rows))


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--design',type=Path,default=ROOT/'inputs/frozen/brain-260917-v1')
    p.add_argument('--historical',type=Path,default=ROOT.parent/'results/consolidated.results')
    p.add_argument('--output',type=Path,default=ROOT/'review/reference-audit')
    a=p.parse_args(); audit(a.design,a.historical,a.output)
