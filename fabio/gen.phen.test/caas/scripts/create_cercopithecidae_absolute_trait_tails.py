#!/usr/bin/env python3

"""Create upper and lower 10% relative-brain-mass tails in Cercopithecidae."""

from __future__ import annotations

import argparse
import csv
import math
from pathlib import Path


SCRIPT_DIR = Path(__file__).resolve().parent
PROJECT_DIR = SCRIPT_DIR.parent
DEFAULT_INPUT = (
    PROJECT_DIR
    / "inputs"
    / "config.creation"
    / "10_cercopithecidae_relative_brain_mass.tsv"
)
DEFAULT_OUTPUT = (
    PROJECT_DIR
    / "inputs"
    / "config.creation"
    / "11_cercopithecidae_absolute_trait_tails.tsv"
)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", type=Path, default=DEFAULT_INPUT)
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    parser.add_argument("--tail-fraction", type=float, default=0.10)
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    if not 0 < args.tail_fraction < 0.5:
        raise ValueError("--tail-fraction must be greater than 0 and less than 0.5")

    with args.input.open(newline="") as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    required = {"species", "relative_brain_mass"}
    if not rows or not required.issubset(rows[0]):
        raise ValueError(f"Input must contain: {', '.join(sorted(required))}")
    if len({row["species"] for row in rows}) != len(rows):
        raise ValueError("Input contains duplicated species")

    ranked = sorted(
        rows,
        key=lambda row: (-float(row["relative_brain_mass"]), row["species"]),
    )
    tail_size = math.ceil(len(ranked) * args.tail_fraction)
    foreground = ranked[:tail_size]
    background = list(reversed(ranked[-tail_size:]))
    if {row["species"] for row in foreground} & {
        row["species"] for row in background
    }:
        raise ValueError("Upper and lower tails overlap")

    output_rows = []
    for group, selected in (("FG", foreground), ("BG", background)):
        for tail_rank, row in enumerate(selected, start=1):
            species = row["species"]
            output_rows.append(
                {
                    "Group": group,
                    "Species": species,
                    "Genus": species.split("_", 1)[0],
                    "RelativeBrainMass": row["relative_brain_mass"],
                    "TailRank": tail_rank,
                    "TotalSpecies": len(ranked),
                    "TailFraction": f"{args.tail_fraction:g}",
                }
            )

    args.output.parent.mkdir(parents=True, exist_ok=True)
    fieldnames = [
        "Group",
        "Species",
        "Genus",
        "RelativeBrainMass",
        "TailRank",
        "TotalSpecies",
        "TailFraction",
    ]
    with args.output.open("w", newline="") as handle:
        writer = csv.DictWriter(
            handle, fieldnames=fieldnames, delimiter="\t", lineterminator="\n"
        )
        writer.writeheader()
        writer.writerows(output_rows)

    print(
        f"Wrote {len(foreground)} FG and {len(background)} BG species "
        f"from {len(ranked)} Cercopithecidae to {args.output}."
    )


if __name__ == "__main__":
    main()
