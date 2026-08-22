#!/usr/bin/env python3

"""Create reproducible 4-vs-4 pooled CAAStools benchmark configurations.

The strategies differ only in how biological hypotheses are assembled.
Each output row has the CAAStools resampling format:

    cycle_id<TAB>FG_species_comma_list<TAB>BG_species_comma_list

Candidate hypotheses are ranked with SHA-256 using a user-supplied seed.  This
makes the selected set reproducible across Python versions and platforms.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import itertools
from dataclasses import dataclass
from pathlib import Path
from typing import FrozenSet, Iterable, Optional, Sequence


SCRIPT_DIR = Path(__file__).resolve().parent
PROJECT_DIR = SCRIPT_DIR.parent
DEFAULT_TABLE_DIR = PROJECT_DIR / "inputs" / "config.creation"
DEFAULT_OUTPUT_DIR = PROJECT_DIR / "inputs" / "benchmark-configs"
DEFAULT_POOL_DIR = PROJECT_DIR / "inputs" / "benchmark-pools"
DEFAULT_MANIFEST = PROJECT_DIR / "inputs" / "benchmark.configs.tsv"


@dataclass(frozen=True)
class Hypothesis:
    foreground: tuple[str, ...]
    background: tuple[str, ...]

    def canonical_key(self) -> str:
        return ",".join(sorted(self.foreground)) + "\t" + ",".join(sorted(self.background))


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open(newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def seeded_rank(seed: int, label: str, hypothesis: Hypothesis) -> str:
    payload = f"{seed}\t{label}\t{hypothesis.canonical_key()}".encode()
    return hashlib.sha256(payload).hexdigest()


def select_reproducibly(
    candidates: Iterable[Hypothesis], cycles: int, seed: int, label: str
) -> list[Hypothesis]:
    unique = {candidate.canonical_key(): candidate for candidate in candidates}
    if len(unique) < cycles:
        raise ValueError(f"{label}: requested {cycles} cycles but only {len(unique)} are possible")
    ranked = sorted(unique.values(), key=lambda h: (seeded_rank(seed, label, h), h.canonical_key()))
    return ranked[:cycles]


def family_extrema(rows: Sequence[dict[str, str]], group_size: int) -> Iterable[Hypothesis]:
    """Choose families as paired units: maximum -> FG and minimum -> BG."""
    # Single-species families have the same taxon as both minimum and maximum,
    # so they do not define a directional phenotype contrast and are excluded.
    informative_rows = [
        row for row in rows if row["maximum_species"] != row["minimum_species"]
    ]
    for selected in itertools.combinations(informative_rows, group_size):
        yield Hypothesis(
            tuple(row["maximum_species"] for row in selected),
            tuple(row["minimum_species"] for row in selected),
        )


def unique_genus_combinations(
    rows: Sequence[dict[str, str]], group_size: int
) -> Iterable[tuple[str, ...]]:
    """Yield species subsets in which no two species belong to the same genus."""
    for selected in itertools.combinations(rows, group_size):
        if len({row["genus"] for row in selected}) == group_size:
            yield tuple(row["taxon"] for row in selected)


def absolute_trait_tails(rows: Sequence[dict[str, str]], group_size: int) -> Iterable[Hypothesis]:
    """Independently combine valid absolute-trait FG and BG tail subsets."""
    foreground = [row for row in rows if row["group"] == "FG"]
    background = [row for row in rows if row["group"] == "BG"]
    fg_groups = list(unique_genus_combinations(foreground, group_size))
    bg_groups = list(unique_genus_combinations(background, group_size))
    for fg, bg in itertools.product(fg_groups, bg_groups):
        yield Hypothesis(fg, bg)


def pss_pairs(rows: Sequence[dict[str, str]], group_size: int) -> Iterable[Hypothesis]:
    """Choose top-PSS genus pairs as units: higher phenotype -> FG, lower -> BG."""
    for selected in itertools.combinations(rows, group_size):
        yield Hypothesis(
            tuple(row["higher_phenotype_species"] for row in selected),
            tuple(row["lower_phenotype_species"] for row in selected),
        )


def genus(species: str) -> str:
    return species.split("_", 1)[0]


def pss_top_pair_cycles(
    rows: Sequence[dict[str, str]],
    group_size: int,
    *,
    congeneric_only: bool = False,
    allowed_genera: Optional[FrozenSet[str]] = None,
) -> Iterable[Hypothesis]:
    """Combine top-1% PSS pairs while preserving fixed pooled sides.

    Pair direction is defined only by the observed trait: the higher endpoint
    is FG and the lower endpoint is BG. Species that occur in both roles in
    the filtered pair set are removed with every pair containing them because
    species-aware pooled event reconstruction requires fixed side membership.
    Each candidate cycle must contain four distinct FG and four distinct BG
    species, with no species shared across sides.
    """
    filtered = []
    for row in rows:
        higher = row["higher_phenotype_species"]
        lower = row["lower_phenotype_species"]
        higher_genus = genus(higher)
        lower_genus = genus(lower)
        if congeneric_only and higher_genus != lower_genus:
            continue
        if allowed_genera is not None and (
            higher_genus not in allowed_genera or lower_genus not in allowed_genera
        ):
            continue
        filtered.append(row)

    higher_species = {row["higher_phenotype_species"] for row in filtered}
    lower_species = {row["lower_phenotype_species"] for row in filtered}
    ambiguous = higher_species & lower_species
    fixed_side_rows = [
        row
        for row in filtered
        if row["higher_phenotype_species"] not in ambiguous
        and row["lower_phenotype_species"] not in ambiguous
    ]

    for selected in itertools.combinations(fixed_side_rows, group_size):
        foreground = tuple(row["higher_phenotype_species"] for row in selected)
        background = tuple(row["lower_phenotype_species"] for row in selected)
        if len(set(foreground)) != group_size:
            continue
        if len(set(background)) != group_size:
            continue
        if set(foreground) & set(background):
            continue
        yield Hypothesis(foreground, background)


def all_congeneric_top_pss_pairs(
    rows: Sequence[dict[str, str]], group_size: int
) -> Iterable[Hypothesis]:
    return pss_top_pair_cycles(rows, group_size, congeneric_only=True)


def focal_three_genera_top_pss_pairs(
    rows: Sequence[dict[str, str]], group_size: int
) -> Iterable[Hypothesis]:
    return pss_top_pair_cycles(
        rows,
        group_size,
        allowed_genera=frozenset({"Macaca", "Papio", "Trachypithecus"}),
    )


def macaca_papio_top_pss_pairs(
    rows: Sequence[dict[str, str]], group_size: int
) -> Iterable[Hypothesis]:
    return pss_top_pair_cycles(
        rows,
        group_size,
        allowed_genera=frozenset({"Macaca", "Papio"}),
    )


def validate_fixed_sides(hypotheses: Sequence[Hypothesis], label: str) -> tuple[list[str], list[str]]:
    foreground = sorted({species for h in hypotheses for species in h.foreground})
    background = sorted({species for h in hypotheses for species in h.background})
    overlap = sorted(set(foreground) & set(background))
    if overlap:
        raise ValueError(f"{label}: species occur on both FG and BG sides: {', '.join(overlap)}")
    for hypothesis in hypotheses:
        if len(hypothesis.foreground) != len(set(hypothesis.foreground)):
            raise ValueError(f"{label}: duplicate species within an FG group")
        if len(hypothesis.background) != len(set(hypothesis.background)):
            raise ValueError(f"{label}: duplicate species within a BG group")
    return foreground, background


def write_hypotheses(path: Path, hypotheses: Sequence[Hypothesis]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w") as handle:
        for index, hypothesis in enumerate(hypotheses, start=1):
            handle.write(
                f"b_{index}\t{','.join(hypothesis.foreground)}\t{','.join(hypothesis.background)}\n"
            )


def write_pool(path: Path, foreground: Sequence[str], background: Sequence[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w") as handle:
        for species in foreground:
            handle.write(f"{species}\t1\n")
        for species in background:
            handle.write(f"{species}\t0\n")


def print_examples(label: str, hypotheses: Sequence[Hypothesis], number: int) -> None:
    print(f"\n{label}: first {min(number, len(hypotheses))} selected cycles")
    for index, hypothesis in enumerate(hypotheses[:number], start=1):
        print(f"  b_{index:03d} FG: {', '.join(hypothesis.foreground)}")
        print(f"        BG: {', '.join(hypothesis.background)}")


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--table-dir", type=Path, default=DEFAULT_TABLE_DIR)
    parser.add_argument("--output-dir", type=Path, default=DEFAULT_OUTPUT_DIR)
    parser.add_argument("--pool-dir", type=Path, default=DEFAULT_POOL_DIR)
    parser.add_argument("--manifest", type=Path, default=DEFAULT_MANIFEST)
    parser.add_argument("--cycles", type=int, default=100)
    parser.add_argument("--group-size", type=int, default=4)
    parser.add_argument("--seed", type=int, default=260821)
    parser.add_argument("--examples", type=int, default=3)
    args = parser.parse_args()

    definitions = [
        (
            "01_family_extrema",
            args.table_dir / "01_family_trait_extrema.tsv",
            family_extrema,
        ),
        (
            "02_absolute_trait_tails",
            args.table_dir / "02_trait_distribution_tails.tsv",
            absolute_trait_tails,
        ),
        (
            "03_pss_best_pair_per_genus",
            args.table_dir / "04_best_top1pct_pair_per_genus.tsv",
            pss_pairs,
        ),
        (
            "04_pss_all_congeneric_top1pct",
            args.table_dir / "05_all_global_top1pct_pss_pairs.tsv",
            all_congeneric_top_pss_pairs,
        ),
        (
            "05_pss_macaca_papio_trachypithecus_top1pct",
            args.table_dir / "05_all_global_top1pct_pss_pairs.tsv",
            focal_three_genera_top_pss_pairs,
        ),
        (
            "06_pss_macaca_papio_top1pct",
            args.table_dir / "05_all_global_top1pct_pss_pairs.tsv",
            macaca_papio_top_pss_pairs,
        ),
    ]

    manifest_rows = []
    for label, source, builder in definitions:
        rows = read_tsv(source)
        candidates = list(builder(rows, args.group_size))
        selected = select_reproducibly(candidates, args.cycles, args.seed, label)
        foreground, background = validate_fixed_sides(selected, label)

        config_path = args.output_dir / f"{label}.pooled.caas.cfg"
        pool_path = args.pool_dir / f"{label}.pool.caas.cfg"
        write_hypotheses(config_path, selected)
        write_pool(pool_path, foreground, background)
        digest = hashlib.sha256(config_path.read_bytes()).hexdigest()
        manifest_rows.append(
            (
                label,
                source.relative_to(PROJECT_DIR),
                config_path.relative_to(PROJECT_DIR),
                pool_path.relative_to(PROJECT_DIR),
                len(candidates),
                len(selected),
                digest,
            )
        )
        print_examples(label, selected, args.examples)
        print(
            f"        candidates={len(candidates)}; "
            f"discovery pool={len(foreground)} FG vs {len(background)} BG"
        )

    args.manifest.parent.mkdir(parents=True, exist_ok=True)
    with args.manifest.open("w", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(
            ["approach", "source_table", "hypotheses_config", "pool_config", "possible_cycles", "selected_cycles", "sha256"]
        )
        writer.writerows(manifest_rows)

    print(f"\nWrote {len(manifest_rows)} benchmark configurations to {args.output_dir}")
    print(f"Manifest: {args.manifest}")
    print(f"Seed: {args.seed}; cycles per approach: {args.cycles}; comparison: {args.group_size} vs {args.group_size}")


if __name__ == "__main__":
    main()
