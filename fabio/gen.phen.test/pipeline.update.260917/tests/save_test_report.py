#!/usr/bin/env python3
"""Run checks and save unedited logs/exit codes outside the production outputs."""
import datetime
import os
import subprocess
import sys
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/'scripts'))
from common import tree_hash, write_json


if __name__=='__main__':
    output=ROOT/'review/tests'; output.mkdir(parents=True,exist_ok=True)
    results=[]
    for name,folder in [('validation','tests'),('bundled-caastools','bin/caastools/test')]:
        cmd=[sys.executable,'-B','-m','unittest','discover','-s',folder,'-v']
        p=subprocess.run(cmd,cwd=ROOT,stdout=subprocess.PIPE,stderr=subprocess.STDOUT,text=True,
                         env=dict(os.environ,PYTHONDONTWRITEBYTECODE='1'))
        (output/(name+'.log')).write_text(p.stdout)
        results.append(dict(suite=name,command=cmd,exit_code=p.returncode,log=name+'.log'))
        print(name,p.returncode,flush=True)
    write_json(output/'test_report.json',dict(timestamp=datetime.datetime.now(datetime.timezone.utc).isoformat(),
               scripts_sha256=tree_hash(ROOT/'scripts'),suites=results,status='pass' if all(r['exit_code']==0 for r in results) else 'failed'))
    raise SystemExit(int(any(r['exit_code'] for r in results)))
