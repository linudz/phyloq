#!/usr/bin/env python3
"""Measured pilot resources plus explicit workload counts; not a cluster price quote."""
import argparse
import statistics
from pathlib import Path
from common import read_json, read_tsv, require, write_json, write_tsv

PACKAGES={'deterministic':(3,300),'R0_pilot':(19,1900),'deterministic_R0_total99':(102,10200),
          'R1_pilot_optional':(19,1900),'paired_optional':(2,30),'combined_optional_package':(123,12130)}


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--run-root',type=Path,required=True)
    p.add_argument('--production-alignments',type=int,required=True)
    p.add_argument('--trace',type=Path,action='append',default=[])
    p.add_argument('--output',type=Path,required=True)
    a=p.parse_args(); receipts=[read_json(p) for p in a.run_root.glob('*/*/raw/*/*.receipt.json')]
    require(receipts,'No completed task receipts in this explicit run')
    require(len({(r['hypothesis_id'],r['gene_id']) for r in receipts})==len(receipts),'Duplicate task receipts')
    real=not any(r['smoke'] for r in receipts)
    rows=[dict(hypothesis_id=r['hypothesis_id'],gene=r['gene_id'],alignment_length=r['coverage']['n_positions'],
               alignment_species=r['coverage']['n_alignment_species'],elapsed_seconds=r['elapsed_seconds'],cpu_seconds=r['cpu_seconds'],
               max_rss_bytes=r['max_rss_bytes'],output_bytes=r['output_bytes']) for r in receipts]
    write_tsv(a.output/'measured_tasks.tsv',rows,list(rows[0]))
    estimates=[]
    for name,(h,c) in PACKAGES.items():
        jobs=h*a.production_alignments
        estimates.append(dict(package=name,hypotheses=h,cycles=c,alignments=a.production_alignments,discovery_tasks=jobs,
            extrapolated_cpu_hours=statistics.mean(r['cpu_seconds'] for r in rows)*jobs/3600 if real else '',
            extrapolated_raw_output_GB=statistics.mean(r['output_bytes'] for r in rows)*jobs/1e9 if real else '',
            status='rough_empirical_estimate_requires_representative_panel' if real else 'synthetic_smoke_not_extrapolated'))
    write_tsv(a.output/'workload_estimates.tsv',estimates,list(estimates[0]))
    trace=[r for f in a.trace for r in read_tsv(f)]
    failures=[r for r in trace if r.get('status') in ('FAILED','ABORTED')]
    write_json(a.output/'pilot_cost.json',dict(measured_successful_tasks=len(receipts),smoke=not real,trace_records=len(trace),failed_or_aborted_trace_records=len(failures),
        cpu_seconds=sum(r['cpu_seconds'] for r in rows),max_rss_bytes=max(r['max_rss_bytes'] for r in rows),
        raw_output_bytes=sum(r['output_bytes'] for r in rows),workloads=estimates,
        caveats=['Wall time depends on queue and concurrency; no billing tariff assumed',
                 'Include failed-task usage and driver/aggregation resources from scheduler accounting before cost approval',
                 'Inspect length-dependent runtime, scheduler delay, staging I/O and disk amplification before batching',
                 'R0-complete is 99 TOTAL, not 99 additional; reuse pilot via a checksummed cohort'],
        manual_review_required=['scheduler overhead','complete failure rate','staged work-directory disk use','actual CPU-hour tariff','representativeness of alignment panel','batching proposal']))


if __name__=='__main__': main()
