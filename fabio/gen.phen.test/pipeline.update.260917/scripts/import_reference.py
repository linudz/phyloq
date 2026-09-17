#!/usr/bin/env python3
"""Import a COMPLETE verified result bundle, with explicit checksum provenance.

Raw historical tables alone cannot be imported. The bundle must contain current
receipt-schema files tying each successful gene to the frozen original config,
pool, tool, settings and alignment SHA; verified reconstruction of old metadata
requires a human review, never an assumption. Original files remain read-only.
"""
import argparse
import shutil
from pathlib import Path
from common import checked, read_json, read_tsv, require, sha, write_json, write_tsv, verify_design
from summarize_hypothesis import summarize


def run(design, hypothesis, alignments, bundle, approval, output):
    verify_design(design)
    row=next(r for r in read_tsv(design/'execution_manifest.tsv') if r['hypothesis_id']==hypothesis)
    require(row['historical_source']!='','Only registered historical references can be imported')
    require(not output.exists(),'Import output exists; choose a new isolated directory')
    ap=read_json(approval)
    for key,expected in dict(bundle_manifest_sha256=sha(bundle), alignment_manifest_sha256=sha(alignments),
                             config_sha256=row['config_sha256'],pool_sha256=row['pool_sha256']).items():
        require(ap.get(key)==expected,'Reference approval mismatch: '+key)
    require(ap.get('reference_compatibility_review') is True and ap.get('approved_by'),'Reference provenance review not recorded')
    inventory=read_tsv(bundle); genes={r['gene_id'] for r in read_tsv(alignments)}
    require(len(inventory)==len(genes) and {r['gene_id'] for r in inventory}==genes,'Import bundle incomplete or duplicated')
    # Validate all evidence before copying any file.
    for r in inventory:
        receipt=read_json(checked(Path(r['receipt_path']),r['receipt_sha256']))
        for key in ('config_sha256','pool_sha256'): require(receipt[key]==row[key], 'Wrong historical '+key)
        require(receipt['hypothesis_id']==hypothesis and receipt['status']=='success','Wrong/failed receipt')
        require(receipt['gene_id']==r['gene_id'] and bool(receipt.get('command')),'Unresolved gene or command provenance')
        require(receipt['tool_sha256']==read_json(design/'design.lock.json')['caastools_sha256'],'Incompatible historical CAAStools')
        for k in ('discovery','filter'):
            require(receipt['settings'][k]==read_json(design/'settings.json')[k],'Incompatible historical '+k)
        checked(Path(r['events_path']),receipt['events_sha256']); checked(Path(r['legacy_path']),receipt['legacy_sha256'])
    for r in inventory:
        dest=output/'raw'/r['gene_id']; dest.mkdir(parents=True)
        for field,suffix in (('receipt_path','.receipt.json'),('events_path','.pooled.caas.events.tsv'),('legacy_path','.pooled.caas.tsv')):
            shutil.copyfile(r[field],dest/(r['gene_id']+suffix))
    summarize(hypothesis,alignments,output/'raw',output/hypothesis)
    write_json(output/'import.provenance.json',dict(bundle_sha256=sha(bundle),approval_sha256=sha(approval),sources=inventory,
              historical_source=row['historical_source'],originals_modified=False))


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--hypothesis',required=True)
    for n in ('design','alignments','bundle','approval','output'): p.add_argument('--'+n,type=Path,required=True)
    a=p.parse_args(); run(a.design,a.hypothesis,a.alignments,a.bundle,a.approval,a.output)
