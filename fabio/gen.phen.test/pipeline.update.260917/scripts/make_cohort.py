#!/usr/bin/env python3
"""Bind an explicit list of complete summaries; never glob a live result tree."""
import argparse
from pathlib import Path
from common import read_json, require, sha, write_tsv


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--summary',type=Path,action='append',required=True)
    p.add_argument('--output',type=Path,required=True)
    a=p.parse_args(); rows=[]; seen=set()
    for folder in a.summary:
        data=read_json(folder/'summary.json'); hid=data['hypothesis_id']
        require(data['status']=='complete' and hid not in seen, 'Incomplete or duplicate hypothesis: '+hid)
        require(folder.name==hid,'Summary directory basename must equal hypothesis ID')
        seen.add(hid); rows.append(dict(hypothesis_id=hid,summary_dir=str(folder.resolve()),summary_lock_sha256=sha(folder/'summary.lock.json')))
    write_tsv(a.output,sorted(rows,key=lambda r:r['hypothesis_id']),['hypothesis_id','summary_dir','summary_lock_sha256'])
