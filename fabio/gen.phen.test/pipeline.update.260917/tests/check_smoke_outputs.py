#!/usr/bin/env python3
"""Assert completed synthetic end-to-end outputs, including real resume evidence."""
import argparse
from collections import Counter
from pathlib import Path
import sys
ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/'scripts'))
from common import read_json, read_tsv, require, sha, write_json


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--run-id',required=True)
    p.add_argument('--partial-resume',action='store_true')
    p.add_argument('--output',type=Path,required=True)
    a=p.parse_args(); run=ROOT/'results'/a.run_id; summaries=ROOT/'summaries'/a.run_id
    receipts=[read_json(p) for p in run.glob('*/*/raw/*/*.receipt.json')]
    require(len(receipts)==18 and len({(r['hypothesis_id'],r['gene_id']) for r in receipts})==18,'Duplicated/missing published task')
    for folder in run.glob('*/*/*_*/'):
        if not (folder/'summary.json').exists(): continue
        s=read_json(folder/'summary.json')
        require(s['completed_genes']==s['background_genes']==2 and s['status']=='complete','Incomplete/incorrect synthetic background')
        require('ZERO' in {r['gene'] for r in read_tsv(folder/'background.tsv')},'Zero-result gene lost from background')
    p0=run/'P0/000/P0_000'
    s=read_json(p0/'summary.json')
    require(s['significant_positions']==4 and s['retained_positions']==1 and s['query_genes']==1,'Unexpected positive/filter output')
    require(read_tsv(p0/'query.tsv')==[{'gene':'POSITIVE'}],'Unexpected P0 query')
    pos=read_tsv(p0/'positions.tsv')
    require({int(r['position']) for r in pos if r['cluster_filter_status']=='retained'}=={15},'Wrong retained position')
    require(all(r['fg_support_count']=='9' and r['bg_support_count']=='10' for r in pos),'Wrong pooled species support')
    status=read_json(summaries/'functional/enrichment_status.json')
    require(status['status'] in ('complete','blocked_annotation_not_frozen','not_requested_discovery_only'),'Incomplete functional stage')
    for name in ('01_query_background_fraction','02_positions_before_after_filter'):
        for ext in ('png','pdf'): require((summaries/'figures'/(name+'.'+ext)).stat().st_size>0,'Missing figure')
    traces=sorted(run.glob('nextflow.*.trace.tsv'))
    last=read_tsv(traces[-1]); caas=[r for r in last if r['name'].startswith('CAAS_POOLED')]
    counts=dict(Counter(r['status'] for r in caas))
    if a.partial_resume:
        require(len(traces)>=2 and counts.get('CACHED',0)>0 and counts.get('COMPLETED',0)>0,'No real partial resume evidence')
        require(len(caas)==18,'Unexpected resumed CAAS task count')
    write_json(a.output,dict(status='pass',run_id=a.run_id,synthetic_only=True,unique_published_tasks=18,
        hypotheses=9,alignments=2,P0_significant_positions=4,P0_retained_positions=1,P0_query_genes=1,background_genes=2,
        retained_position=15,P0_fg_support=9,P0_bg_support=10,functional_status=status,
        latest_trace_caas_statuses=counts,launch_lock_sha256=sha(run/'launch.lock.json'),
        trace_sha256={str(p.relative_to(ROOT)):sha(p) for p in traces}))
    print('PASS:',a.run_id,counts)
