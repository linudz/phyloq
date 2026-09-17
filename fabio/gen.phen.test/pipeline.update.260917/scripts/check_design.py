#!/usr/bin/env python3
"""Independent structural checks of a frozen design and paired/random contracts."""
import argparse
import math
from collections import Counter, defaultdict
from pathlib import Path
from common import ROOT, cycle_read, digest, pool_read, read_json, read_tsv, require, verify_design, write_json


def check(design):
    lock = verify_design(design)
    expected = read_json(ROOT / "inputs/expected_pools.json")
    manifest = read_tsv(design / "execution_manifest.tsv")
    require(len({r["hypothesis_id"] for r in manifest}) == len(manifest), "Duplicate hypotheses")
    seen = defaultdict(set); pools_by_id = {}; cycles_by_id = {}
    traits = {r["species"]: float(r["relative_brain_mass"]) for r in read_tsv(design / "sources/phenotype.tsv")}
    eligible = defaultdict(set)
    for row in read_tsv(design / "eligible_species.tsv"): eligible[row["randomization"]].add(row["species"])
    require(len(eligible["R0"]) == 19 and len(eligible["R1"]) == 31, "Wrong randomization universe")
    p0 = {side: set(expected["P0"][side].split()) for side in ("FG", "BG")}
    for row in manifest:
        hid, sid = row["hypothesis_id"], row["strategy_id"]
        pools = pool_read(design / row["pool_config"])
        require(all(s in traits for items in pools.values() for s in items), f"Unknown species in {hid}")
        require(not set(pools["FG"]) & set(pools["BG"]), f"Cross-side overlap {hid}")
        pools_by_id[hid] = pools
        cc = cycle_read(design / row["hypotheses_config"], pools); cycles_by_id[hid] = cc
        require(len(cc) == int(row["selected_cycles"]), f"Wrong cycles {hid}")
        if sid in ("P0", "N3", "N4", "N5"):
            require(all(set(pools[side]) == set(expected[sid][side].split()) for side in pools), f"Unexpected {sid} pools")
        if sid in ("N5", "R0", "R1"):
            for side, index in (("FG", 1), ("BG", 2)):
                counts = Counter(s.split("_")[0] for s in pools[side])
                require(dict(counts) == {g: q[index] for g, q in expected["quotas"].items()}, f"Quota mismatch {hid} {side}")
        if sid in ("R0", "R1"):
            union = set(sum(pools.values(), []))
            require(union <= eligible[sid], f"Random species outside universe {hid}")
            if sid == "R0": require(union == eligible[sid], "R0 changed selected species")
            require(any(set(pools[side]) != p0[side] for side in pools), "P0 included in random null")
            key = digest({side: sorted(pools[side]) for side in pools})
            require(key not in seen[sid], f"Duplicate random assignment {hid}"); seen[sid].add(key)
        if sid in ("N3", "N4", "N5", "R0", "R1"):
            require(len(pools["FG"]) == 9 and len(pools["BG"]) == 10, "Wrong 9/10 sizes")
            require(int(row["possible_cycles"]) == math.comb(9,4)*math.comb(10,4) == 26460, "Wrong Cartesian cycle count")
            require(len(cc) == 100, "Not 100 unique cycles")
    n4 = pools_by_id["N4_000"]
    require(set(sum(n4.values(), [])) == eligible["R0"], "N4 union not P0")
    require(min(traits[s] for s in n4["FG"]) > max(traits[s] for s in n4["BG"]), "N4 phenotype overlap")
    pairmap = defaultdict(dict)
    for row in read_tsv(design / "pairs.tsv"):
        pairmap[row["strategy_id"]][row["linked_unit"]] = (row["foreground"], row["background"])
    subsets = {}
    for sid in ("P1", "N6"):
        pairs = pairmap[sid]
        require(len(pairs) == 6 and len(set(sum((list(x) for x in pairs.values()), []))) == 12, f"Invalid linked pairs {sid}")
        cc = cycles_by_id[sid + "_000"]
        require(len(cc) == 15, "Exactly choose(6,4)=15 cycles required")
        subsets[sid] = []
        for _, f, b in cc:
            units = {g for g, (fs, bs) in pairs.items() if fs in f and bs in b}
            require(len(units) == 4, f"Broken pair linkage {sid}")
            subsets[sid].append(sorted(units))
    require(subsets["P1"] == subsets["N6"], "P1/N6 genus subsets not corresponding")
    require(math.comb(8,3)*math.comb(5,3)*2**3 == 4480, "Wrong R0 assignment universe")
    return dict(status="pass", design_id=lock["design_id"], hypotheses=len(manifest),
                r0_assignments=4480, r0_excluded_observed=1, possible_9x10_cycles=26460,
                random_replicates=lock["random_replicates"], paired_cycles=15)


if __name__ == "__main__":
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--design", type=Path, required=True); p.add_argument("--output", type=Path, required=True)
    a = p.parse_args(); write_json(a.output, check(a.design))
