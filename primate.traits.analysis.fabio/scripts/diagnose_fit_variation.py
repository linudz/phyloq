#!/usr/bin/env python3

"""Quantify observed-fit variation among parametric-bootstrap chunks."""

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
        description=(
            "Quantify variation in observed BM/OU fits among bootstrap chunks."
        )
    )
    parser.add_argument(
        "--chunks-dir",
        default="results/chunks",
        help="Directory containing *.bootstrap_fit.tsv files.",
    )
    parser.add_argument(
        "--top",
        type=int,
        default=20,
        help="Number of groups with the largest relative variation to show.",
    )
    return parser.parse_args()


def main():
    arguments = parse_arguments()
    if arguments.top < 1:
        sys.exit("--top must be a positive integer.")

    files = sorted(
        Path(arguments.chunks_dir).glob("*.bootstrap_fit.tsv")
    )
    if not files:
        sys.exit(
            f"No *.bootstrap_fit.tsv files found in {arguments.chunks_dir}"
        )

    groups = defaultdict(list)
    generating_models = defaultdict(set)

    for path in files:
        with path.open(newline="", encoding="utf-8") as handle:
            reader = csv.DictReader(handle, delimiter="\t")
            required = {"trait", "model", "sigma2", "generating_model"}
            if reader.fieldnames is None or not required.issubset(reader.fieldnames):
                missing = sorted(required.difference(reader.fieldnames or []))
                sys.exit(f"{path} is missing columns: {', '.join(missing)}")

            for row in reader:
                trait = row["trait"]
                model = row["model"]
                sigma2 = parse_number(row["sigma2"])
                if sigma2 is not None:
                    groups[(trait, model)].append(sigma2)
                if row["generating_model"]:
                    generating_models[trait].add(row["generating_model"])

    results = []
    for (trait, model), values in groups.items():
        minimum = min(values)
        maximum = max(values)
        mean = sum(values) / len(values)
        denominator = max(abs(mean), sys.float_info.epsilon)
        results.append(
            {
                "trait": trait,
                "model": model,
                "n_chunks": len(values),
                "sigma2_min": minimum,
                "sigma2_max": maximum,
                "relative_range": (maximum - minimum) / denominator,
            }
        )

    results.sort(key=lambda row: row["relative_range"], reverse=True)
    print(f"Fit files read: {len(files)}")
    print(f"Trait/model groups: {len(results)}")
    print()

    header = (
        f"{'trait':45s} {'model':6s} {'chunks':>7s} "
        f"{'sigma2.min':>14s} {'sigma2.max':>14s} "
        f"{'relative.range':>16s}"
    )
    print(header)
    print("-" * len(header))
    for row in results[: arguments.top]:
        print(
            f"{row['trait'][:45]:45s} "
            f"{row['model']:6s} "
            f"{row['n_chunks']:7d} "
            f"{row['sigma2_min']:14.6g} "
            f"{row['sigma2_max']:14.6g} "
            f"{row['relative_range']:16.6g}"
        )

    unstable_models = {
        trait: models
        for trait, models in generating_models.items()
        if len(models) > 1
    }
    print()
    if unstable_models:
        print("WARNING: generating model changed among chunks:")
        for trait, models in sorted(unstable_models.items()):
            print(f"  {trait}: {', '.join(sorted(models))}")
    else:
        print("Generating model was consistent among chunks for every trait.")


if __name__ == "__main__":
    main()
