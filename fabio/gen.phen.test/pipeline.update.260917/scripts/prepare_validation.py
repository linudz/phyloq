#!/usr/bin/env python3
"""Freeze species assignments and cycles without submitting discovery jobs.

Assignments use CPython random.Random (MT19937), sequential uniform draws with
rejection of complete duplicates/P0. The prefix never depends on requested B.
Cycles follow the production SHA-256 ranking convention. Randomization uses
lexically sorted species only; phenotype is consulted solely for final display.
"""
from __future__ import annotations

import argparse
import itertools as it
import math
import platform
import random
import shutil
from collections import Counter
from pathlib import Path

from common import (ROOT, VERSION, checked, cycle_read, digest, identifier, pool_read,
                    read_json, read_tsv, require, sha, tree_hash, verify_design,
                    write_json, write_tsv)

HISTORICAL = {
    "P0": "08_pss_cercopithecidae_random_pools",
    "N0": "02_absolute_trait_tails", "N1": "01_family_extrema",
    "N2": "09_cercopithecidae_absolute_trait_tails",
    "P2": "07_pss_ranked_endpoint_disjoint_13x13",
}
SOURCES = {
    "P0": "09_cercopithecidae_pss_random_pools.tsv",
    "N0": "02_trait_distribution_tails.tsv", "N1": "01_family_trait_extrema.tsv",
    "N2": "11_cercopithecidae_absolute_trait_tails.tsv",
    "P2": "08_pss_ranked_endpoint_disjoint_13_pairs.tsv",
    "N3": "10_cercopithecidae_relative_brain_mass.tsv",
    "N4": "09_cercopithecidae_pss_random_pools.tsv",
    "N5": "10_cercopithecidae_relative_brain_mass.tsv",
    "R0": "09_cercopithecidae_pss_random_pools.tsv",
    "R1": "10_cercopithecidae_relative_brain_mass.tsv",
    "P1": "04_best_top1pct_pair_per_genus.tsv",
    "N6": "10_cercopithecidae_relative_brain_mass.tsv",
}
EXEC_COLUMNS = ["hypothesis_id", "strategy_id", "replicate_id", "null_family", "action",
                "sampling_mode", "pool_config", "hypotheses_config", "possible_cycles",
                "selected_cycles", "assignment_sha256", "pool_sha256", "config_sha256",
                "source_sha256", "source_table", "historical_source", "reference_status"]


def genus(s):
    return s.split("_", 1)[0]


def canonical(pools):
    return {side: sorted(pools[side]) for side in ("FG", "BG")}


def random_assignments(eligible, quotas, excluded, n, seed, family):
    """Every accepted complete assignment is uniform among remaining assignments.

    Each genus draw samples an ordered subset without replacement. Splitting
    it at its FG quota assigns every allowed (FG,BG) combination equal mass.
    Rejection of complete duplicate assignments gives sampling without replacement.
    """
    rng = random.Random(int(digest([seed, family]), 16))
    seen = {digest(canonical(excluded))}
    out = []
    while len(out) < n:
        pools = {"FG": [], "BG": []}
        for g, (_, nf, nb) in sorted(quotas.items()):
            candidates = sorted(s for s in eligible if genus(s) == g)
            picked = rng.sample(candidates, nf + nb)
            pools["FG"].extend(picked[:nf]); pools["BG"].extend(picked[nf:])
        key = digest(canonical(pools))
        if key not in seen:
            seen.add(key); out.append(pools)
    return out


def generate_cycles(pools, pairs, n, seed, label):
    if pairs:
        require(len(pairs) >= 4, f"Insufficient linked pairs for {label}")
        candidates = [(tuple(p[1] for p in group), tuple(p[2] for p in group),
                       tuple(p[0] for p in group)) for group in it.combinations(pairs, 4)]
    else:
        candidates = [(f, b, ()) for f in it.combinations(pools["FG"], 4)
                      for b in it.combinations(pools["BG"], 4)]
    require(len(candidates) >= n, f"Insufficient unique cycles for {label}")
    if pairs and len(pairs) == 6:
        # P1/N6 share the exact lexical genus-subset order.
        selected = candidates
    else:
        import hashlib
        def ranking(c):
            key = ",".join(sorted(c[0])) + "\t" + ",".join(sorted(c[1]))
            return hashlib.sha256(f"{seed}\t{label}\t{key}".encode()).hexdigest(), key
        selected = sorted(candidates, key=ranking)[:n]
    return len(candidates), selected


def source_inventory(pipeline, phenotype, taxonomy, tree):
    files = {"phenotype.tsv": phenotype, "taxonomy.tsv": taxonomy, "tree.nwk": tree,
             "benchmark.configs.tsv": pipeline / "inputs/benchmark.configs.tsv"}
    for name in sorted(set(SOURCES.values())):
        files["config.creation/" + name] = pipeline / "inputs/config.creation" / name
    for row in read_tsv(files["benchmark.configs.tsv"]):
        if row["approach"] in HISTORICAL.values():
            for field in ("pool_config", "hypotheses_config"):
                files[row[field]] = pipeline / row[field]
    return files


def prepare(args):
    out = args.output.resolve()
    source_lock = read_json(args.source_lock)
    if args.from_frozen:
        verify_design(args.from_frozen)
        sources = {rel:args.from_frozen/'sources'/rel for rel in source_lock['files']}
    else:
        sources = source_inventory(args.pipeline, args.phenotype, args.taxonomy, args.tree)
    for rel, path in sources.items():
        checked(path, source_lock["files"][rel])
    if out.exists():
        lock = verify_design(out)
        require(lock["seed"] == args.seed and lock["random_replicates"] == args.replicates,
                "Design already frozen with other seed/count; choose a new design ID")
        require(lock["source_lock_sha256"] == sha(args.source_lock), "Frozen source lock changed")
        print(f"Validated existing immutable design: {out}")
        return
    require(1 <= args.replicates <= 99, "Freeze 1–99 null hypotheses per family")
    identifier(args.design_id)
    traits = {r["species"]: float(r["relative_brain_mass"]) for r in read_tsv(sources['phenotype.tsv'])}
    taxonomy = {r["species"]: r["family"] for r in read_tsv(sources['taxonomy.tsv'])}
    cercos = {r["species"]: float(r["relative_brain_mass"]) for r in read_tsv(sources["config.creation/10_cercopithecidae_relative_brain_mass.tsv"])}
    require(len(cercos) == 55, "Expected 55 Cercopithecidae")
    require(all(s in traits and math.isclose(v, traits[s], abs_tol=1e-10) for s, v in cercos.items()),
            "Phenotype differs from frozen primate values")
    expected = read_json(ROOT / "inputs/expected_pools.json")
    quotas = expected["quotas"]
    ranked = lambda items: sorted(items, key=lambda s: (-traits[s], s))
    extremes = lambda species, f, b: {"FG": ranked(species)[:f], "BG": ranked(species)[-b:]}
    p0 = {side: [r["Species"] for r in read_tsv(sources["config.creation/" + SOURCES["P0"]])
                 if r["Group"] == side] for side in ("FG", "BG")}
    eligible = sorted(s for s in cercos if genus(s) in quotas)
    require(len(eligible) == 31, "Expected 31 species in five genera")
    for g, (n, nf, nb) in quotas.items():
        require(sum(genus(s) == g for s in eligible) == n, f"Eligible quota mismatch: {g}")
        require([sum(genus(s) == g for s in p0[side]) for side in ("FG", "BG")] == [nf, nb],
                f"P0 quota mismatch: {g}")
    n5 = {"FG": [], "BG": []}
    for g, (_, nf, nb) in sorted(quotas.items()):
        e = extremes([s for s in eligible if genus(s) == g], nf, nb)
        for side in n5: n5[side].extend(e[side])
    definitions = {"P0": p0, "N3": extremes(cercos, 9, 10),
                   "N4": extremes(p0["FG"] + p0["BG"], 9, 10), "N5": n5}
    for sid, pools in definitions.items():
        for side in pools:
            require(set(pools[side]) == set(expected[sid][side].split()), f"Brief/source discrepancy: {sid} {side}")
    require(set(sum(definitions["N4"].values(), [])) == set(sum(p0.values(), [])), "N4 union mismatch")
    require(min(traits[s] for s in definitions["N4"]["FG"]) > max(traits[s] for s in definitions["N4"]["BG"]), "N4 not separated")
    assignments = [(sid, "000", pools, []) for sid, pools in definitions.items() if sid != "P0"]
    for family, universe in (("R0", sum(p0.values(), [])), ("R1", eligible)):
        for index, pools in enumerate(random_assignments(universe, quotas, p0, args.replicates, args.seed, family), 1):
            assignments.append((family, f"{index:03d}", pools, []))
    paired = {"P1": [], "N6": []}
    for r in read_tsv(sources["config.creation/04_best_top1pct_pair_per_genus.tsv"]):
        if r["family"] == "Cercopithecidae":
            g = r["genus"]
            values = expected["paired"][g]
            require([r["higher_phenotype_species"], r["lower_phenotype_species"]] == values[:2], f"P1 mismatch: {g}")
            e = extremes([s for s in cercos if genus(s) == g], 1, 1)
            require([e["FG"][0], e["BG"][0]] == values[2:], f"N6 mismatch: {g}")
            paired["P1"].append((g, *values[:2])); paired["N6"].append((g, *values[2:]))
    for sid, pairs in paired.items():
        require(len(pairs) == 6, "Paired comparison requires six genera")
        pairs.sort()
        assignments.append((sid, "000", {"FG": [p[1] for p in pairs], "BG": [p[2] for p in pairs]}, pairs))
    # All source checks precede output creation. Frozen originals remain portable.
    out.mkdir(parents=True)
    for rel, path in sources.items():
        dest = out / "sources" / rel
        dest.parent.mkdir(parents=True, exist_ok=True); shutil.copyfile(path, dest)
    shutil.copyfile(args.source_lock, out / "sources.lock.json")
    shutil.copyfile(ROOT / "inputs/settings.json", out / "settings.json")
    executions, members, cycles, pair_rows = [], [], [], []
    historical_manifest = {r["approach"]: r for r in read_tsv(sources["benchmark.configs.tsv"])}
    for sid, old in HISTORICAL.items():
        row = historical_manifest[old]
        pools = pool_read(sources[row['pool_config']])
        require(all(s in traits for side in pools.values() for s in side), f"Unknown species in {sid}")
        if sid == "P0": require(canonical(pools) == canonical(p0), "Historical P0 pool mismatch")
        require(sha(sources[row['hypotheses_config']]) == row["sha256"], f"Historical config checksum mismatch: {sid}")
        assignments.append((sid, "000", pools, []))
    for sid, rid, pools, pairs in assignments:
        hid = f"{sid}_{rid}"
        require(not set(pools["FG"]) & set(pools["BG"]), f"Overlapping complete pools: {hid}")
        pools = {side: ranked(items) for side, items in pools.items()}
        folder = out / "configs" / hid; folder.mkdir(parents=True)
        pool_path, config_path = folder / "pool.cfg", folder / "cycles.cfg"
        old = HISTORICAL.get(sid, "")
        if old:
            row = historical_manifest[old]
            shutil.copyfile(sources[row['pool_config']], pool_path)
            shutil.copyfile(sources[row['hypotheses_config']], config_path)
            possible, n = int(row["possible_cycles"]), int(row["selected_cycles"])
            mode = {"N0": "within_side_unique_genus", "N1": "linked_family", "P2": "linked_pss"}.get(sid, "unrestricted")
        else:
            n = 15 if pairs else 100
            possible, selected = generate_cycles(pools, pairs, n, args.seed, hid)
            mode = "linked_genus" if pairs else "unrestricted"
            pool_path.write_text("".join(f"{s}\t{int(side == 'FG')}\n" for side in pools for s in pools[side]))
            config_path.write_text("".join(f"b_{i}\t{','.join(f)}\t{','.join(b)}\n" for i, (f, b, _) in enumerate(selected, 1)))
        observed_cycles = cycle_read(config_path, pools)
        require(len(observed_cycles) == n, f"Wrong cycle count: {hid}")
        if sid == "N0":
            require(all(len(set(map(genus, f))) == len(set(map(genus, b))) == 4 for _, f, b in observed_cycles), "N0 genus rule broken")
        if sid == "N1":
            linked = {r["maximum_species"]: r["minimum_species"] for r in read_tsv(sources["config.creation/01_family_trait_extrema.tsv"]) if r["maximum_species"] != r["minimum_species"]}
            require(all({linked[s] for s in f} == set(b) for _, f, b in observed_cycles), "N1 family linkage broken")
        for side in pools:
            for s in pools[side]:
                unit = next((p[0] for p in pairs if s in p[1:]), taxonomy[s] if sid == "N1" else "")
                members.append(dict(strategy_id=sid, replicate_id=rid, hypothesis_id=hid, side=side,
                                    species=s, genus=genus(s), family=taxonomy[s], trait_value=traits[s], linked_unit=unit))
        for cid, f, b in observed_cycles:
            units = sorted(p[0] for p in pairs if p[1] in f) if pairs else sorted({taxonomy[s] for s in f}) if sid == "N1" else []
            cycles.append(dict(hypothesis_id=hid, cycle_id=cid, foreground=",".join(f), background=",".join(b), linked_units=",".join(units)))
        for unit, f, b in pairs:
            pair_rows.append(dict(strategy_id=sid, replicate_id=rid, linked_unit=unit, foreground=f, background=b))
        executions.append(dict(hypothesis_id=hid, strategy_id=sid, replicate_id=rid,
            null_family=sid if sid in ("R0", "R1") else "observed", action="discovery",
            sampling_mode=mode, pool_config=str(pool_path.relative_to(out)), hypotheses_config=str(config_path.relative_to(out)),
            possible_cycles=possible, selected_cycles=n, assignment_sha256=digest(canonical(pools)),
            pool_sha256=sha(pool_path), config_sha256=sha(config_path), source_sha256=sha(sources["config.creation/" + SOURCES[sid]]),
            source_table="sources/config.creation/" + SOURCES[sid], historical_source=old,
            reference_status="requires_provenance_audit" if old else "new"))
    executions.sort(key=lambda r: r["hypothesis_id"])
    write_tsv(out / "execution_manifest.tsv", executions, EXEC_COLUMNS)
    replicates = [dict(r, seed=args.seed, generator=VERSION, prng="CPython-MT19937-rejection-v1", python_version=platform.python_version()) for r in executions]
    write_tsv(out / "replicate_manifest.tsv", replicates, EXEC_COLUMNS + ["seed", "generator", "prng", "python_version"])
    write_tsv(out / "pool_membership.tsv", members, ["strategy_id", "replicate_id", "hypothesis_id", "side", "species", "genus", "family", "trait_value", "linked_unit"])
    write_tsv(out / "cycles.tsv", cycles, ["hypothesis_id", "cycle_id", "foreground", "background", "linked_units"])
    write_tsv(out / "pairs.tsv", pair_rows, ["strategy_id", "replicate_id", "linked_unit", "foreground", "background"])
    eligible_rows = [dict(randomization=family, species=s, genus=genus(s), family=taxonomy[s], trait_value=traits[s], taxonomy_source="sources/taxonomy.tsv")
                     for family, universe in (("R0", sum(p0.values(), [])), ("R1", eligible)) for s in ranked(universe)]
    write_tsv(out / "eligible_species.tsv", eligible_rows, ["randomization", "species", "genus", "family", "trait_value", "taxonomy_source"])
    registry = []
    for sid in sorted(SOURCES):
        example = next(r for r in executions if r["strategy_id"] == sid)
        p = pool_read(out / example["pool_config"])
        registry.append(dict(example, fg_pool_size=len(p["FG"]), bg_pool_size=len(p["BG"]),
            hypothesis_count=sum(r["strategy_id"] == sid for r in executions)))
    write_tsv(out / "strategy_registry.tsv", registry, ["strategy_id", "historical_source", "sampling_mode", "source_table", "fg_pool_size", "bg_pool_size", "selected_cycles", "hypothesis_count", "reference_status"])
    selections = {"deterministic": {"N3", "N4", "N5"}, "paired": {"P1", "N6"},
                  "reference-rerun": {"P0", "N0", "N1", "N2"}, "auxiliary-p2": {"P2"}, "smoke": {"N3"}}
    for stage, sids in selections.items():
        write_tsv(out / "selections" / f"{stage}.tsv", [r for r in executions if r["strategy_id"] in sids], EXEC_COLUMNS)
    for family in ("R0", "R1"):
        for stage, limit in (("pilot", 19), ("complete", 99)):
            write_tsv(out / "selections" / f"{family.lower()}-{stage}.tsv",
                      [r for r in executions if r["strategy_id"] == family and int(r["replicate_id"]) <= limit], EXEC_COLUMNS)
    lock = dict(design_id=args.design_id, seed=args.seed, random_replicates=args.replicates,
                generator=VERSION, source_lock_sha256=sha(args.source_lock),
                generator_sha256=sha(Path(__file__)), common_sha256=sha(ROOT / "scripts/common.py"),
                caastools_sha256=tree_hash(ROOT / "bin/caastools"),
                files={str(p.relative_to(out)): sha(p) for p in sorted(out.rglob("*")) if p.is_file()})
    write_json(out / "design.lock.json", lock)
    verify_design(out)
    print(f"Frozen {len(executions)} hypotheses, {len(cycles)} cycles at {out}")


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--pipeline", type=Path, default=ROOT.parent / "pipeline")
    modes = p.add_mutually_exclusive_group()
    modes.add_argument('--from-frozen', type=Path, help='Portable source snapshot; default is the bundled frozen design')
    modes.add_argument('--from-project', action='store_true', help='Explicitly audit/rebuild from the original external project inputs')
    case = ROOT.parents[2] / "04.case.study.bodybrain.mammals/relative.brain.mass/trait.in.primates"
    p.add_argument("--phenotype", type=Path, default=case / "inputs/relative_brain_mass.venditti.primates.tsv")
    p.add_argument("--taxonomy", type=Path, default=case / "inputs/relative_brain_mass.primates.taxonomy.tsv")
    p.add_argument("--tree", type=Path, default=case / "results/relative_brain_mass.Kuderna_S4.analysis_tree.nwk")
    p.add_argument("--source-lock", type=Path, default=ROOT / "inputs/source.lock.json")
    p.add_argument("--lock-sources", action="store_true")
    p.add_argument("--output", type=Path)
    p.add_argument("--design-id", default="brain-260917-v1")
    p.add_argument("--seed", type=int, default=260917)
    p.add_argument("--replicates", type=int, default=99)
    args = p.parse_args()
    if not args.from_project and not args.lock_sources:
        args.from_frozen = (args.from_frozen or ROOT/'inputs/frozen/brain-260917-v1').resolve()
        args.source_lock = args.from_frozen/'sources.lock.json'
    if args.lock_sources:
        require(args.from_frozen is None, '--lock-sources requires original project inputs, not an existing frozen design')
    if args.lock_sources:
        require(not args.source_lock.exists(), "Source lock exists; review changes explicitly before making a new lock")
        sources = source_inventory(args.pipeline, args.phenotype, args.taxonomy, args.tree)
        write_json(args.source_lock, dict(files={k: sha(v) for k, v in sources.items()},
                   original_paths={k: str(v.resolve()) for k, v in sources.items()}))
    else:
        args.output = args.output or ROOT / "inputs/frozen" / args.design_id
        prepare(args)


if __name__ == "__main__":
    main()
