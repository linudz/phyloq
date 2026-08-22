#!/usr/bin/env python3

"""Compare PSS and absolute-trait assignments in the same species set.

The input is the deduplicated Cercopithecidae pool used by benchmark 08.
The output retains exactly those species and pool sizes, but assigns the
highest relative-brain-mass values to FG and the remaining values to BG.
"""

from __future__ import annotations

import argparse
import csv
from pathlib import Path


SCRIPT_DIR = Path(__file__).resolve().parent
PROJECT_DIR = SCRIPT_DIR.parent
DEFAULT_INPUT = (
    PROJECT_DIR
    / "inputs"
    / "config.creation"
    / "09_cercopithecidae_pss_random_pools.tsv"
)
DEFAULT_OUTPUT = (
    PROJECT_DIR
    / "inputs"
    / "config.creation"
    / "10_cercopithecidae_pss_vs_absolute_trait_assignments.tsv"
)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", type=Path, default=DEFAULT_INPUT)
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    with args.input.open(newline="") as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))

    required = {"Group", "Species", "RelativeBrainMass"}
    if not rows or not required.issubset(rows[0]):
        raise ValueError(f"Input must contain: {', '.join(sorted(required))}")
    if len({row["Species"] for row in rows}) != len(rows):
        raise ValueError("Input contains duplicated species")

    foreground_count = sum(row["Group"] == "FG" for row in rows)
    background_count = sum(row["Group"] == "BG" for row in rows)
    if foreground_count + background_count != len(rows):
        raise ValueError("Every input species must be assigned to FG or BG")

    ranked = sorted(
        rows,
        key=lambda row: (-float(row["RelativeBrainMass"]), row["Species"]),
    )
    output_rows = []
    for rank, row in enumerate(ranked, start=1):
        absolute_group = "FG" if rank <= foreground_count else "BG"
        pss_group = row["Group"]
        output_rows.append(
            {
                "Rank": rank,
                "Species": row["Species"],
                "RelativeBrainMass": row["RelativeBrainMass"],
                "PSSGroup": pss_group,
                "AbsoluteTraitGroup": absolute_group,
                "AssignmentChange": (
                    "unchanged"
                    if pss_group == absolute_group
                    else f"{pss_group}_to_{absolute_group}"
                ),
            }
        )

    args.output.parent.mkdir(parents=True, exist_ok=True)
    fieldnames = [
        "Rank",
        "Species",
        "RelativeBrainMass",
        "PSSGroup",
        "AbsoluteTraitGroup",
        "AssignmentChange",
    ]
    with args.output.open("w", newline="") as handle:
        writer = csv.DictWriter(
            handle, fieldnames=fieldnames, delimiter="\t", lineterminator="\n"
        )
        writer.writeheader()
        writer.writerows(output_rows)

    print(
        f"Wrote {len(output_rows)} species to {args.output} "
        f"({foreground_count} FG, {background_count} BG by absolute trait)."
    )


if __name__ == "__main__":
    main()
