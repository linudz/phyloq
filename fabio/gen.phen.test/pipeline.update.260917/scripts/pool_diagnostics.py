#!/usr/bin/env python3
"""Document actual phylogenetic distances; genus matching is not distance matching."""
import argparse
import itertools
import statistics
from collections import defaultdict
from pathlib import Path
from Bio import Phylo
from common import ROOT, read_tsv, require, write_tsv


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--design',type=Path,default=ROOT/'inputs/frozen/brain-260917-v1')
    p.add_argument('--output',type=Path,default=ROOT/'review/pool-diagnostics')
    a=p.parse_args(); tree=Phylo.read(str(a.design/'sources/tree.nwk'),'newick')
    tips={c.name:c for c in tree.get_terminals()}; groups=defaultdict(lambda:defaultdict(list)); traits={}
    for r in read_tsv(a.design/'pool_membership.tsv'):
        groups[r['hypothesis_id']][r['side']].append(r['species']); traits[r['species']]=float(r['trait_value'])
    # Compute once across the entire union, then reuse the fixed distance matrix.
    species=sorted(traits); require(set(species)<=tips.keys(),'Selected species absent from frozen tree')
    distances={tuple(sorted((s,t))):tree.distance(tips[s],tips[t]) for s,t in itertools.combinations(species,2)}
    rows=[]; values=[]
    for hid,pools in sorted(groups.items()):
        for kind,pairs in (('within_FG',itertools.combinations(pools['FG'],2)),('within_BG',itertools.combinations(pools['BG'],2)),('FG_to_BG',itertools.product(pools['FG'],pools['BG']))):
            ds=[]
            for s,t in pairs:
                d=distances[tuple(sorted((s,t)))]; ds.append(d)
                values.append(dict(hypothesis_id=hid,comparison=kind,species1=s,species2=t,patristic_distance=d))
            rows.append(dict(hypothesis_id=hid,comparison=kind,n_pairs=len(ds),minimum=min(ds),maximum=max(ds),mean=statistics.mean(ds),median=statistics.median(ds)))
    write_tsv(a.output/'patristic_summary.tsv',rows,list(rows[0]))
    write_tsv(a.output/'patristic_distances.tsv',values,list(values[0]))
    rows=[]
    for hid,pools in sorted(groups.items()):
        for side,ss in pools.items():
            vs=[traits[s] for s in ss]; rows.append(dict(hypothesis_id=hid,side=side,n_species=len(ss),minimum=min(vs),maximum=max(vs),mean=statistics.mean(vs),median=statistics.median(vs)))
    write_tsv(a.output/'phenotype_summary.tsv',rows,list(rows[0]))


if __name__=='__main__': main()
