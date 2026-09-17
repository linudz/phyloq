#!/usr/bin/env python3
"""Smoke-only fault injector, to test a real partial Nextflow resume.

Create tests/fixtures/fail_once.marker before the first run. One ZERO task
removes it and fails. All later attempts use the unchanged interpreter wrapper.
Never use this test interpreter with production data.
"""
import os
import sys
from pathlib import Path
ROOT = Path(os.environ.get('CAAS_VALIDATION_TEST_ROOT', Path(__file__).resolve().parents[1]))
args = sys.argv[1:]
marker = ROOT / 'tests/fixtures/fail_once.marker'
if args and args[0].endswith('/run_caas.py'):
    import json
    job = json.loads(Path(args[args.index('--job') + 1]).read_text())
    if not job.get('smoke'): raise SystemExit('Fault-injector forbidden for production')
    if job['gene_id'] == 'ZERO' and job['hypothesis_id'] == 'R1_001' and marker.exists():
        marker.unlink()
        raise SystemExit('Deliberate one-time synthetic failure to verify partial resume')
os.execv(sys.executable, [sys.executable, *args])
