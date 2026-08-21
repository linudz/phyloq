#!/usr/bin/env python3

"""Materialize one deterministic CAAStools pooled-hypothesis table per run."""

import argparse
import hashlib
import sys
from pathlib import Path


def main():
    parser = argparse.ArgumentParser(
        description="Generate and document pooled CAAStools hypotheses once for a Nextflow run."
    )
    parser.add_argument("--caastools-dir", type=Path, required=True)
    parser.add_argument("--pool-file", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--metadata-output", type=Path, required=True)
    parser.add_argument("--fg-size", type=int, default=4)
    parser.add_argument("--bg-size", type=int, default=4)
    parser.add_argument("--comparisons", default="max")
    parser.add_argument("--seed", type=int, default=260811)
    arguments = parser.parse_args()

    sys.path.insert(0, str(arguments.caastools_dir.resolve()))
    from modules.init_bootstrap import pool_discovery_hypotheses

    pooled = pool_discovery_hypotheses(
        pool_file=arguments.pool_file,
        fg_size=arguments.fg_size,
        bg_size=arguments.bg_size,
        comparisons=arguments.comparisons,
        seed=arguments.seed,
    )
    pooled.print_traits(arguments.output)

    hypothesis_sha256 = hashlib.sha256(arguments.output.read_bytes()).hexdigest()
    metadata_rows = [
        ("method", "pooled-discovery"),
        ("pool_file", str(arguments.pool_file)),
        ("complete_fg_species", str(len(pooled.discovery_fg))),
        ("complete_bg_species", str(len(pooled.discovery_bg))),
        ("hypothesis_fg_size", str(pooled.pool_fg_size)),
        ("hypothesis_bg_size", str(pooled.pool_bg_size)),
        ("possible_unique_comparisons", str(pooled.pool_maximum_comparisons)),
        ("requested_comparisons", str(arguments.comparisons)),
        ("selected_comparisons", str(pooled.cycles)),
        ("seed", str(arguments.seed)),
        ("selection_method", "seeded_sha256_ranking"),
        ("hypotheses_sha256", hypothesis_sha256),
    ]
    with arguments.metadata_output.open("w") as metadata_handle:
        print("parameter\tvalue", file=metadata_handle)
        for parameter, value in metadata_rows:
            print(parameter + "\t" + value, file=metadata_handle)

    print("Pooled hypotheses:", arguments.output)
    print("Selected comparisons:", pooled.cycles)
    print("SHA-256:", hypothesis_sha256)


if __name__ == "__main__":
    main()
