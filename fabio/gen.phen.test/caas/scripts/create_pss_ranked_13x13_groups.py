#!/usr/bin/env python3
"""Create 13 independent FG/BG pairs from the ranked global PSS top 1%.

Pairs are visited in decreasing PSS order. A pair is retained only when
neither endpoint has already been used by a previously retained pair. The
higher-phenotype endpoint becomes FG and the lower-phenotype endpoint becomes
BG. Selection stops after 13 pairs, yielding 13 unique FG and 13 unique BG
species while preserving the PSS link between every FG and BG member.

This is a deterministic greedy matching: changing the source ranking changes
the result, but repeated runs on the same ranked input are identical.
"""

from __future__ import annotations

import argparse
import csv
from pathlib import Path


SCRIPT_DIR = Path(__file__).resolve().parent
PIPELINE_DIR = SCRIPT_DIR.parent
DEFAULT_INPUT = (
    PIPELINE_DIR
    / "inputs"
    / "config.creation"
    / "05_all_global_top1pct_pss_pairs.tsv"
)
DEFAULT_OUTPUT_DIR = PIPELINE_DIR / "inputs" / "config.creation"


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def write_tsv(path: Path, rows: list[dict[str, object]]) -> None:
    if not rows:
        raise ValueError(f"Refusing to write an empty table: {path}")
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=list(rows[0]),
            delimiter="\t",
            lineterminator="\n",
        )
        writer.writeheader()
        writer.writerows(rows)


def epithet(taxon: str) -> str:
    return taxon.split("_", 1)[1]


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", type=Path, default=DEFAULT_INPUT)
    parser.add_argument("--output-dir", type=Path, default=DEFAULT_OUTPUT_DIR)
    parser.add_argument("--pairs", type=int, default=13)
    args = parser.parse_args()

    source = read_tsv(args.input)
    source.sort(key=lambda row: int(row["top1pct_global_rank"]))

    selected: list[dict[str, str]] = []
    used_species: set[str] = set()

    for row in source:
        fg = row["higher_phenotype_species"]
        bg = row["lower_phenotype_species"]
        # This global check prevents duplicates within a side and also prevents
        # a species from being FG in one pair and BG in another.
        if fg in used_species or bg in used_species:
            continue
        selected.append(row)
        used_species.update((fg, bg))
        if len(selected) == args.pairs:
            break

    if len(selected) != args.pairs:
        raise ValueError(
            f"Only {len(selected)} endpoint-disjoint pairs can be selected; "
            f"{args.pairs} were requested"
        )

    pair_rows: list[dict[str, object]] = []
    group_rows: list[dict[str, object]] = []
    for selected_rank, row in enumerate(selected, start=1):
        pair_rows.append({"selected_pair_rank": selected_rank, **row})

        common = {
            "selected_pair_rank": selected_rank,
            "global_top1pct_rank": row["top1pct_global_rank"],
            "pss_final_score": row["pss_final_score"],
        }
        group_rows.extend(
            [
                {
                    "group": "FG",
                    "family": row["family1"],
                    "genus": row["genus1"],
                    "species": epithet(row["higher_phenotype_species"]),
                    "taxon": row["higher_phenotype_species"],
                    "relative_brain_mass": row["higher_phenotype_value"],
                    "linked_partner": row["lower_phenotype_species"],
                    **common,
                },
                {
                    "group": "BG",
                    "family": row["family2"],
                    "genus": row["genus2"],
                    "species": epithet(row["lower_phenotype_species"]),
                    "taxon": row["lower_phenotype_species"],
                    "relative_brain_mass": row["lower_phenotype_value"],
                    "linked_partner": row["higher_phenotype_species"],
                    **common,
                },
            ]
        )

    group_rows.sort(
        key=lambda row: (
            0 if row["group"] == "FG" else 1,
            int(row["selected_pair_rank"]),
        )
    )

    pair_path = args.output_dir / "08_pss_ranked_endpoint_disjoint_13_pairs.tsv"
    group_path = args.output_dir / "08_pss_ranked_endpoint_disjoint_13x13.tsv"
    write_tsv(pair_path, pair_rows)
    write_tsv(group_path, group_rows)

    print(f"Selected pairs: {len(pair_rows)}")
    print(f"Unique FG species: {sum(r['group'] == 'FG' for r in group_rows)}")
    print(f"Unique BG species: {sum(r['group'] == 'BG' for r in group_rows)}")
    print(f"Lowest source rank reached: {selected[-1]['top1pct_global_rank']}")
    print(f"Pair table: {pair_path}")
    print(f"Group table: {group_path}")


if __name__ == "__main__":
    main()
