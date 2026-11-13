"""Validate the additive P3 bundle without rewriting the historical design."""
import math
from pathlib import Path
from common import checked, cycle_read, digest, pool_read, read_json, read_tsv, require, sha, tree_hash, ROOT

# Bound to this reviewed bundle, not an arbitrary replacement of its manifest.
P3_LOCK_SHA256 = '3f6fe38401d1d9f1b83ac1e1715d283e06652bbd09c31f3f4a705f71fb42f8d1'


def check_additional(design, additional):
    design, additional = Path(design).resolve(), Path(additional).resolve()
    checked(additional/'additional.lock.json', P3_LOCK_SHA256)
    lock = read_json(additional/'additional.lock.json')
    require(lock['base_design_sha256'] == sha(design/'design.lock.json'), 'Additional design has another base design')
    require(lock['caastools_sha256'] == tree_hash(ROOT/'bin/caastools'), 'Additional CAAStools changed')
    for rel, checksum in lock['files'].items():
        require(Path(rel).is_absolute() is False and '..' not in Path(rel).parts, 'Unsafe additional path')
        checked(additional/rel, checksum)
    rows = read_tsv(additional/'selection.tsv')
    require(len(rows) == 1 and rows[0]['hypothesis_id'] == 'P3_000' and rows[0]['strategy_id'] == 'P3', 'Expected only P3_000')
    row = rows[0]
    pools = pool_read(checked(additional/row['pool_config'], row['pool_sha256']))
    require(len(pools['FG']) == 6 and len(pools['BG']) == 13, 'Expected P3 6 FG and 13 BG')
    require(not set(pools['FG']) & set(pools['BG']), 'P3 side conflict')
    require('Homo_sapiens' not in pools['FG'] + pools['BG'], 'Human in nonhuman P3')
    traits = {r['species']:float(r['relative_brain_mass']) for r in read_tsv(design/'sources/phenotype.tsv')}
    values = read_tsv(additional/'sources/Venditti.pool.tsv')
    require(all(math.isclose(float(r['phenotype']), traits[r['species']], abs_tol=1e-10) for r in values), 'P3 phenotype differs from frozen Venditti')
    require(all(set(r['species'] for r in values if r['group'] == side) == set(pools[side]) for side in pools), 'P3 source pool differs')
    pairs = read_tsv(additional/'sources/Venditti.selected_pairs.tsv')
    require(len(pairs) == 20 and len({r['PairKey'] for r in pairs}) == 20, 'Expected 20 distinct P3 pairs')
    require(all('Homo_sapiens' not in (r['FG'],r['BG']) and traits[r['FG']] > traits[r['BG']] for r in pairs), 'P3 orientation or human exclusion changed')
    require(all({r[side] for r in pairs} == set(pools[side]) for side in pools), 'P3 endpoints differ from pool')
    cycles = cycle_read(checked(additional/row['hypotheses_config'],row['config_sha256']), pools)
    require(len(cycles) == int(row['selected_cycles']) == 100, 'Expected 100 P3 cycles')
    require(int(row['possible_cycles']) == math.comb(6,4)*math.comb(13,4) == 10725, 'Wrong P3 cycle universe')
    require(row['assignment_sha256'] == digest({s:sorted(pools[s]) for s in pools}), 'P3 assignment digest mismatch')
    return rows


def additional_identity(design, additional):
    rows = check_additional(design, additional)
    return dict(additional_id=read_json(Path(additional)/'additional.lock.json')['additional_id'],
                additional_lock_sha256=sha(Path(additional)/'additional.lock.json'),
                additional_path=str(Path(additional).resolve()),
                hypotheses=[r['hypothesis_id'] for r in rows])
