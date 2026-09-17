#!/usr/bin/env python3
"""Create explicit proposal manifests; this script never launches anything."""
import argparse
from pathlib import Path
from common import ROOT, read_tsv, sha, verify_design, write_json, write_tsv


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--design',type=Path,default=ROOT/'inputs/frozen/brain-260917-v1')
    p.add_argument('--output',type=Path,default=ROOT/'review/launch-manifests')
    a=p.parse_args(); verify_design(a.design); allrows=read_tsv(a.design/'execution_manifest.tsv'); cols=list(allrows[0])
    selections={
        'r0-extension80':[r for r in allrows if r['strategy_id']=='R0' and 20<=int(r['replicate_id'])<=99],
        'r1-extension80':[r for r in allrows if r['strategy_id']=='R1' and 20<=int(r['replicate_id'])<=99],
        'package102':[r for r in allrows if r['strategy_id'] in ('N3','N4','N5','R0')],
        'package123':[r for r in allrows if r['strategy_id'] in ('N3','N4','N5','R0','P1','N6') or r['strategy_id']=='R1' and int(r['replicate_id'])<=19]}
    summaries=[]
    for name,rows in selections.items():
        f=a.output/(name+'.tsv'); write_tsv(f,rows,cols)
        summaries.append(dict(selection=name,hypotheses=len(rows),cycles=sum(int(r['selected_cycles']) for r in rows),sha256=sha(f)))
    write_tsv(a.output/'package_counts.tsv',summaries,list(summaries[0]))
    write_json(a.output/'launch_review.json',dict(status='not_approved_not_launched',design_sha256=sha(a.design/'design.lock.json'),
        selections=summaries, continuation='After a successful 19-replicate pilot, execute extension80 in a new run, then combine both runs through an explicit 99-hypothesis cohort. Do not rerun the first 19.'))
