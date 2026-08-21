#!/usr/bin/env python3

"""Quantify variation in fitted OU alpha among bootstrap chunks."""

import argparse
import csv
import math
import sys
from collections import defaultdict
from pathlib import Path


def parse_number(value):
    """Parse a finite numeric TSV value, returning None for missing values."""
    if value in {None, "", "NA", "NaN", "nan"}:
        return None
    number = float(value)
    return number if math.isfinite(number) else None


def parse_arguments():
    parser = argparse.ArgumentParser(
        description="Quantify fitted OU-alpha variation among bootstrap chunks."
    )
    parser.add_argument(
        "--chunks-dir",
        default="results/chunks",
        help="Directory containing *.bootstrap_fit.tsv files.",
    )
    parser.add_argument(
        "--top",
        type=int,
        default=30,
        help="Number of non-boundary groups with the largest variation to show.",
    )
    parser.add_argument(
        "--boundary",
        type=float,
        default=1e-8,
        help="Maximum alpha treated as effectively on the alpha=0 boundary.",
    )
    return parser.parse_args()


def main():
    arguments = parse_arguments()
    if arguments.top < 1:
        sys.exit("--top must be a positive integer.")
    if not math.isfinite(arguments.boundary) or arguments.boundary <= 0:
        sys.exit("--boundary must be a positive finite number.")

    files = sorted(Path(arguments.chunks_dir).glob("*.bootstrap_fit.tsv"))
    if not files:
        sys.exit(
            f"No *.bootstrap_fit.tsv files found in {arguments.chunks_dir}"
        )

    groups = defaultdict(list)
    for path in files:
        with path.open(newline="", encoding="utf-8") as handle:
            reader = csv.DictReader(handle, delimiter="\t")
            required = {"trait", "model", "alpha", "generating_model"}
            if reader.fieldnames is None or not required.issubset(reader.fieldnames):
                missing = sorted(required.difference(reader.fieldnames or []))
                sys.exit(f"{path} is missing columns: {', '.join(missing)}")

            for row in reader:
                if row["model"] != "OU":
                    continue
                alpha = parse_number(row["alpha"])
                if alpha is not None:
                    groups[(row["trait"], row["generating_model"])].append(alpha)

    rows = []
    for (trait, generating_model), values in groups.items():
        minimum = min(values)
        maximum = max(values)
        mean = sum(values) / len(values)
        absolute_range = maximum - minimum
        rows.append(
            {
                "trait": trait,
                "generating_model": generating_model,
                "chunks": len(values),
                "minimum": minimum,
                "maximum": maximum,
                "absolute_range": absolute_range,
                "relative_range": absolute_range
                / max(abs(mean), sys.float_info.epsilon),
                "near_zero": maximum < arguments.boundary,
            }
        )

    nonzero = sorted(
        (row for row in rows if not row["near_zero"]),
        key=lambda row: row["relative_range"],
        reverse=True,
    )
    near_zero = sorted(
        (row for row in rows if row["near_zero"]),
        key=lambda row: row["absolute_range"],
        reverse=True,
    )

    print(f"Fit files read: {len(files)}")
    print(f"OU trait groups: {len(rows)}")
    print(f"Boundary threshold: alpha < {arguments.boundary:g}")
    print()
    print("OU fits away from the alpha=0 boundary")
    header = (
        f"{'trait':44s} {'gen':4s} {'n':>3s} "
        f"{'alpha.min':>13s} {'alpha.max':>13s} "
        f"{'abs.range':>13s} {'rel.range':>13s}"
    )
    print(header)
    print("-" * len(header))
    for row in nonzero[: arguments.top]:
        print(
            f"{row['trait'][:44]:44s} "
            f"{row['generating_model']:4s} "
            f"{row['chunks']:3d} "
            f"{row['minimum']:13.6g} "
            f"{row['maximum']:13.6g} "
            f"{row['absolute_range']:13.6g} "
            f"{row['relative_range']:13.6g}"
        )

    print()
    print(f"OU fits effectively at alpha=0: {len(near_zero)}")
    for row in near_zero[:10]:
        print(
            f"{row['trait'][:44]:44s} "
            f"gen={row['generating_model']:2s} "
            f"n={row['chunks']:2d} "
            f"min={row['minimum']:.6g} "
            f"max={row['maximum']:.6g} "
            f"abs.range={row['absolute_range']:.6g}"
        )

    unexpected_chunks = [row for row in rows if row["chunks"] != 10]
    print()
    if unexpected_chunks:
        print("WARNING: groups not represented by exactly 10 chunks:")
        for row in sorted(unexpected_chunks, key=lambda item: item["trait"]):
            print(f"  {row['trait']}: {row['chunks']} chunks")
    else:
        print("Every OU trait group is represented by exactly 10 chunks.")


if __name__ == "__main__":
    main()
