#!/usr/bin/env python3
"""One hypothesis × alignment; preserve raw outputs and attest successful execution."""
from __future__ import annotations
import argparse
import csv
import os
import platform
import resource
import subprocess
import sys
import time
from pathlib import Path
from common import (ROOT, checked, cycle_read, pool_read, read_json, require,
                    sha, tree_hash, write_json)


def coverage(alignment, pools, cycles, settings):
    """Pattern-independent eligibility, before allele and p-value prefilters.

    The current tool treats every character except '-' as non-gap. A gene is
    coverage-testable when at least one column passes overall gap ratio and
    at least one frozen cycle has >=3 observed non-gap species on both sides.
    Constant alignments can be testable zero-discovery genes. This is not an
    assertion that every gene has equal CAAS detection probability.
    """
    from inventory_alignments import load_alignment
    import numpy as np
    require(settings['fmt'] == 'phylip-relaxed', 'Unsupported frozen format')
    msa = load_alignment(alignment)
    species = [r.id for r in msa]
    require(len(species) == len(set(species)), "Duplicate species in alignment")
    index = {s: i for i, s in enumerate(species)}
    data = np.array([list(str(r.seq)) for r in msa]) != "-"
    valid = (1 - data.mean(axis=0)) <= float(settings["max_gaps_per_position"])
    eligible = np.zeros(msa.get_alignment_length(), dtype=bool)
    for _, f, b in cycles:
        f_idx, b_idx = [index[s] for s in f if s in index], [index[s] for s in b if s in index]
        if len(f_idx) < int(settings["min_fg_observed"]) or len(b_idx) < int(settings["min_bg_observed"]):
            continue
        eligible |= valid & (data[f_idx].sum(axis=0) >= int(settings["min_fg_observed"])) & (data[b_idx].sum(axis=0) >= int(settings["min_bg_observed"]))
    return dict(n_positions=msa.get_alignment_length(), n_alignment_species=len(species),
                coverage_testable_positions=int(eligible.sum()), coverage_testable=bool(eligible.any()),
                missing_fg=sorted(set(pools["FG"]) - set(species)), missing_bg=sorted(set(pools["BG"]) - set(species)),
                present_fg=len(set(pools["FG"]) & set(species)), present_bg=len(set(pools["BG"]) & set(species)))


def run(job_file, alignment, pool_file, config, tool_dir, output):
    job = read_json(job_file)
    checked(alignment, job["alignment_sha256"]); checked(pool_file, job["pool_sha256"])
    checked(config, job["config_sha256"])
    require(tree_hash(tool_dir) == job["tool_sha256"], "CAAStools tree checksum mismatch")
    pools = pool_read(pool_file); cycles = cycle_read(config, pools)
    require(len(cycles) == int(job["selected_cycles"]), "Runtime cycle count mismatch")
    require(Path(alignment).name.split(".")[0] == job["gene_id"], "Gene name differs from tool filename rule")
    settings = job["settings"]["discovery"]
    require(all(settings[k] == "NO" for k in ("max_fg_gaps", "max_bg_gaps", "max_gaps", "max_fg_miss", "max_bg_miss", "max_miss")),
            "Coverage contract requires frozen NO gap/missingness settings")
    cov = coverage(alignment, pools, cycles, settings)
    output.mkdir(parents=True, exist_ok=True)
    gene = job["gene_id"]
    legacy, events = output / f"{gene}.pooled.caas.tsv", output / f"{gene}.pooled.caas.events.tsv"
    cmd = [sys.executable, str(tool_dir / "ct"), "pooled-discovery", "-a", str(alignment),
           "-t", str(pool_file), "-s", str(config), "-o", str(legacy), "--event-output", str(events),
           "--hypotheses-output", "none"]
    for key, value in settings.items(): cmd += ["--" + key, str(value)]
    start = time.monotonic()
    before = resource.getrusage(resource.RUSAGE_CHILDREN)
    with (output / "caastools.log").open("w") as handle:
        process = subprocess.run(cmd, stdout=handle, stderr=subprocess.STDOUT,
                                 env=dict(os.environ, PYTHONDONTWRITEBYTECODE="1", PYTHONHASHSEED="0"))
    usage = resource.getrusage(resource.RUSAGE_CHILDREN)
    require(process.returncode == 0, f"CAAS failed with exit {process.returncode}; see caastools.log")
    require(legacy.is_file() and events.is_file(), "CAAS returned without both output files; not a zero result")
    with events.open() as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        require({"gene", "position", "primary_event", "fg_support_count", "bg_support_count", "positional_pvalue"} <= set(reader.fieldnames or []), "Malformed CAAS event header")
        count = 0
        for row in reader:
            require(row["gene"] == gene, f"Event gene does not match inventory: {row['gene']} != {gene}")
            count += 1
    receipt = dict(job, status="success", coverage=cov, raw_event_rows=count,
        command=cmd, python_version=platform.python_version(),
        elapsed_seconds=time.monotonic() - start,
        cpu_seconds=usage.ru_utime + usage.ru_stime - before.ru_utime - before.ru_stime,
        max_rss_bytes=usage.ru_maxrss if sys.platform == "darwin" else usage.ru_maxrss * 1024,
        legacy_sha256=sha(legacy), events_sha256=sha(events),
        output_bytes=legacy.stat().st_size + events.stat().st_size)
    write_json(output / f"{gene}.receipt.json", receipt)


if __name__ == "__main__":
    p = argparse.ArgumentParser(description=__doc__)
    for name in ("job", "alignment", "pool", "config", "tool-dir", "output"):
        p.add_argument("--" + name, type=Path, required=True)
    a = p.parse_args()
    run(a.job, a.alignment, a.pool, a.config, a.tool_dir, a.output)
