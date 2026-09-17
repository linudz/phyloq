#!/usr/bin/env python3
"""Summarise CAAS discoveries by foreground and background species support.

The input is the consolidated CAAS result directory. Each retained approach is
expected to contain a ``caas-pooled-events`` directory with one event TSV per
alignment. Source files are read only; all derived tables are written to the
requested output directory.
"""

from __future__ import annotations

import argparse
import csv
import math
import os
import sys
from collections import Counter, defaultdict
from pathlib import Path
from typing import Dict, Iterable, List, Mapping, MutableMapping, Tuple


REQUIRED_COLUMNS = {
    "gene",
    "position",
    "event_id",
    "primary_event",
    "fg_support_count",
    "fg_discovery_count",
    "bg_support_count",
    "bg_discovery_count",
    "positional_pvalue",
}

APPROACH_METADATA = {
    "02_absolute_trait_tails": {
        "approach": "01_primate_wide_extremes",
        "label": "1 – Primate-wide extremes",
        "order": 1,
    },
    "01_family_extrema": {
        "approach": "02_family_level_parallel_pairs",
        "label": "2 – Family-level parallel pairs",
        "order": 2,
    },
    "08_pss_cercopithecidae_random_pools": {
        "approach": "03_cercopithecidae_pss_driven_groups",
        "label": "3 – Cercopithecidae PSS-driven groups",
        "order": 3,
    },
    "09_cercopithecidae_absolute_trait_tails": {
        "approach": "04_cercopithecidae_wide_extremes",
        "label": "4 – Cercopithecidae-wide extremes",
        "order": 4,
    },
}


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input-dir", required=True, type=Path)
    parser.add_argument("--output-dir", required=True, type=Path)
    parser.add_argument("--pvalue-threshold", type=float, default=0.05)
    parser.add_argument("--cluster-density-threshold", type=float, default=0.7)
    parser.add_argument("--cluster-min-span", type=int, default=3)
    parser.add_argument("--cluster-min-positions", type=int, default=3)
    return parser.parse_args()


def as_int(value: str, field: str, source: Path) -> int:
    try:
        return int(float(value))
    except (TypeError, ValueError) as exc:
        raise ValueError(f"invalid {field}={value!r} in {source}") from exc


def as_float(value: str, field: str, source: Path) -> float:
    try:
        result = float(value)
    except (TypeError, ValueError) as exc:
        raise ValueError(f"invalid {field}={value!r} in {source}") from exc
    if not math.isfinite(result):
        raise ValueError(f"non-finite {field}={value!r} in {source}")
    return result


def is_primary(value: str) -> bool:
    return value.strip().lower() in {"yes", "true", "1", "y"}


def write_tsv(path: Path, columns: List[str], rows: Iterable[Mapping[str, object]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=columns, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        for row in rows:
            writer.writerow({column: row.get(column, "") for column in columns})


def approach_metadata(source_approach: str) -> dict:
    if source_approach not in APPROACH_METADATA:
        raise ValueError(
            f"unrecognised approach directory {source_approach!r}; add it to APPROACH_METADATA"
        )
    return APPROACH_METADATA[source_approach]


def preferred_position_record(candidate: Mapping[str, object], current: Mapping[str, object]) -> bool:
    """Return True when candidate is the deterministic representative event."""
    candidate_rank = (
        float(candidate["positional_pvalue"]),
        -int(candidate["balanced_support_count"]),
        -int(candidate["total_support_count"]),
        str(candidate["event_id"]),
    )
    current_rank = (
        float(current["positional_pvalue"]),
        -int(current["balanced_support_count"]),
        -int(current["total_support_count"]),
        str(current["event_id"]),
    )
    return candidate_rank < current_rank


def dense_cluster_positions(
    position_values: Iterable[int],
    density_threshold: float,
    min_span: int,
    min_positions: int,
) -> set[int]:
    """Return positions belonging to at least one sufficiently dense interval.

    This is equivalent to the tested CAAS-train pruning rule, while separating
    the minimum interval span from the minimum number of hits. Every interval
    with at least ``min_positions`` distinct positions, a span of at least
    ``min_span`` residues, and ``hits / span >= density_threshold`` contributes
    all of its positions to the discarded set.
    """
    unique = sorted(set(position_values))
    n_positions = len(unique)
    if n_positions < min_positions:
        return set()

    # A difference array records the union of all qualifying index intervals
    # without repeatedly adding the same positions to a Python set.
    covered = [0] * (n_positions + 1)
    for left in range(n_positions - min_positions + 1):
        first_right = left + min_positions - 1
        for right in range(first_right, n_positions):
            span = unique[right] - unique[left] + 1
            if span < min_span:
                continue
            hit_count = right - left + 1
            if hit_count / span >= density_threshold:
                covered[left] += 1
                covered[right + 1] -= 1

    result: set[int] = set()
    coverage = 0
    for index, delta in enumerate(covered[:-1]):
        coverage += delta
        if coverage > 0:
            result.add(unique[index])
    return result


def apply_cluster_filter(
    summaries: List[dict],
    positions: List[dict],
    density_threshold: float,
    min_span: int,
    min_positions: int,
) -> None:
    """Annotate positions and approach summaries with cluster-filter status."""
    by_approach_gene: Dict[Tuple[str, str], List[dict]] = defaultdict(list)
    for row in positions:
        by_approach_gene[(str(row["approach"]), str(row["gene"]))].append(row)

    for rows in by_approach_gene.values():
        discarded = dense_cluster_positions(
            (int(row["position"]) for row in rows),
            density_threshold=density_threshold,
            min_span=min_span,
            min_positions=min_positions,
        )
        for row in rows:
            row["cluster_filter_status"] = (
                "discarded_dense_cluster" if int(row["position"]) in discarded else "retained"
            )

    positions_by_approach: Dict[str, List[dict]] = defaultdict(list)
    for row in positions:
        positions_by_approach[str(row["approach"])].append(row)

    for summary in summaries:
        rows = positions_by_approach[str(summary["approach"])]
        retained = [row for row in rows if row["cluster_filter_status"] == "retained"]
        discarded = [row for row in rows if row["cluster_filter_status"] != "retained"]
        all_genes = {str(row["gene"]) for row in rows}
        retained_genes = {str(row["gene"]) for row in retained}
        clustered_genes = {str(row["gene"]) for row in discarded}
        summary.update(
            {
                "cluster_density_threshold": density_threshold,
                "cluster_min_span": min_span,
                "cluster_min_positions": min_positions,
                "retained_significant_positions": len(retained),
                "discarded_clustered_positions": len(discarded),
                "genes_with_clustered_positions": len(clustered_genes),
                "retained_significant_genes": len(retained_genes),
                "lost_clustered_genes": len(all_genes - retained_genes),
            }
        )


def scan_approach(approach_dir: Path, pvalue_threshold: float) -> Tuple[dict, List[dict]]:
    source_approach = approach_dir.name
    metadata = approach_metadata(source_approach)
    approach = metadata["approach"]
    event_dir = approach_dir / "caas-pooled-events"
    if not event_dir.is_dir():
        raise FileNotFoundError(f"missing event directory: {event_dir}")

    metrics = Counter()
    discovery_fg_values: set[int] = set()
    discovery_bg_values: set[int] = set()
    positions: MutableMapping[Tuple[str, str], dict] = {}

    event_files = sorted(
        Path(entry.path)
        for entry in os.scandir(event_dir)
        if entry.is_file() and entry.name.endswith(".tsv")
    )
    metrics["event_files_scanned"] = len(event_files)

    for source in event_files:
        try:
            with source.open("r", encoding="utf-8", newline="") as handle:
                reader = csv.DictReader(handle, delimiter="\t")
                if reader.fieldnames is None:
                    metrics["empty_files"] += 1
                    continue
                missing = REQUIRED_COLUMNS.difference(reader.fieldnames)
                if missing:
                    raise ValueError(f"missing columns {sorted(missing)} in {source}")

                file_has_rows = False
                for row in reader:
                    file_has_rows = True
                    metrics["raw_event_rows"] += 1
                    try:
                        fg_support = as_int(row["fg_support_count"], "fg_support_count", source)
                        bg_support = as_int(row["bg_support_count"], "bg_support_count", source)
                        fg_discovery = as_int(row["fg_discovery_count"], "fg_discovery_count", source)
                        bg_discovery = as_int(row["bg_discovery_count"], "bg_discovery_count", source)
                        positional_pvalue = as_float(row["positional_pvalue"], "positional_pvalue", source)
                    except ValueError as exc:
                        metrics["malformed_rows"] += 1
                        print(f"WARNING: {exc}", file=sys.stderr)
                        continue

                    discovery_fg_values.add(fg_discovery)
                    discovery_bg_values.add(bg_discovery)

                    if not is_primary(row["primary_event"]):
                        continue
                    metrics["primary_event_rows"] += 1
                    if positional_pvalue >= pvalue_threshold:
                        continue
                    metrics["significant_primary_event_rows"] += 1

                    gene = row["gene"].strip()
                    position = row["position"].strip()
                    if not gene or not position:
                        metrics["malformed_rows"] += 1
                        continue

                    record = {
                        "approach": approach,
                        "approach_label": metadata["label"],
                        "source_approach": source_approach,
                        "gene": gene,
                        "position": position,
                        "event_id": row["event_id"].strip(),
                        "fg_support_count": fg_support,
                        "bg_support_count": bg_support,
                        "balanced_support_count": min(fg_support, bg_support),
                        "total_support_count": fg_support + bg_support,
                        "fg_discovery_species_total": fg_discovery,
                        "bg_discovery_species_total": bg_discovery,
                        "positional_pvalue": positional_pvalue,
                        "source_file": source.name,
                    }
                    key = (gene, position)
                    current = positions.get(key)
                    if current is None:
                        positions[key] = record
                    else:
                        metrics["duplicate_significant_position_rows"] += 1
                        if preferred_position_record(record, current):
                            positions[key] = record

                if file_has_rows:
                    metrics["files_with_event_rows"] += 1
        except UnicodeDecodeError as exc:
            raise ValueError(f"cannot decode {source} as UTF-8") from exc

    if not discovery_fg_values or not discovery_bg_values:
        raise ValueError(f"no event rows found for {approach}; discovery totals cannot be inferred")

    position_rows = sorted(
        positions.values(),
        key=lambda row: (str(row["gene"]), int(row["position"]), str(row["event_id"])),
    )
    summary = {
        "approach": approach,
        "approach_label": metadata["label"],
        "source_approach": source_approach,
        "event_files_scanned": metrics["event_files_scanned"],
        "files_with_event_rows": metrics["files_with_event_rows"],
        "raw_event_rows": metrics["raw_event_rows"],
        "primary_event_rows": metrics["primary_event_rows"],
        "significant_primary_event_rows": metrics["significant_primary_event_rows"],
        "unique_significant_genes": len({row["gene"] for row in position_rows}),
        "unique_significant_positions": len(position_rows),
        "fg_discovery_species_total": max(discovery_fg_values),
        "bg_discovery_species_total": max(discovery_bg_values),
        "fg_discovery_totals_observed": ",".join(map(str, sorted(discovery_fg_values))),
        "bg_discovery_totals_observed": ",".join(map(str, sorted(discovery_bg_values))),
        "max_significant_fg_support": max((int(row["fg_support_count"]) for row in position_rows), default=0),
        "max_significant_bg_support": max((int(row["bg_support_count"]) for row in position_rows), default=0),
        "duplicate_significant_position_rows": metrics["duplicate_significant_position_rows"],
        "malformed_rows": metrics["malformed_rows"],
        "pvalue_threshold": pvalue_threshold,
    }
    return summary, position_rows


def aggregate_genes(summary: Mapping[str, object], positions: List[dict]) -> List[dict]:
    by_gene: Dict[str, List[dict]] = defaultdict(list)
    for row in positions:
        by_gene[str(row["gene"])].append(row)

    output: List[dict] = []
    for gene, rows in sorted(by_gene.items()):
        retained_rows = [row for row in rows if row["cluster_filter_status"] == "retained"]
        output.append(
            {
                "approach": summary["approach"],
                "approach_label": summary["approach_label"],
                "source_approach": summary["source_approach"],
                "gene": gene,
                "significant_position_count": len(rows),
                "retained_position_count": len(retained_rows),
                "discarded_clustered_position_count": len(rows) - len(retained_rows),
                "cluster_filter_status": "retained" if retained_rows else "lost_all_positions_clustered",
                "max_fg_support_count": max(int(row["fg_support_count"]) for row in rows),
                "max_bg_support_count": max(int(row["bg_support_count"]) for row in rows),
                "max_balanced_support_same_position": max(int(row["balanced_support_count"]) for row in rows),
                "max_retained_fg_support_count": max(
                    (int(row["fg_support_count"]) for row in retained_rows), default=""
                ),
                "max_retained_bg_support_count": max(
                    (int(row["bg_support_count"]) for row in retained_rows), default=""
                ),
                "max_retained_balanced_support_same_position": max(
                    (int(row["balanced_support_count"]) for row in retained_rows), default=""
                ),
                "min_positional_pvalue": min(float(row["positional_pvalue"]) for row in rows),
                "fg_discovery_species_total": summary["fg_discovery_species_total"],
                "bg_discovery_species_total": summary["bg_discovery_species_total"],
            }
        )
    return output


def build_gene_support_distribution(summaries: List[dict], genes: List[dict]) -> List[dict]:
    genes_by_approach: Dict[str, List[dict]] = defaultdict(list)
    for row in genes:
        genes_by_approach[str(row["approach"])].append(row)

    output: List[dict] = []
    for summary in summaries:
        approach = str(summary["approach"])
        approach_genes = genes_by_approach[approach]
        for side, field, total_field in (
            ("FG", "max_fg_support_count", "fg_discovery_species_total"),
            ("BG", "max_bg_support_count", "bg_discovery_species_total"),
        ):
            counts = Counter(int(row[field]) for row in approach_genes)
            retained_counts = Counter(
                int(row[field])
                for row in approach_genes
                if row["cluster_filter_status"] == "retained"
            )
            discovery_total = int(summary[total_field])
            gene_total = len(approach_genes)
            for support in range(1, discovery_total + 1):
                count = counts[support]
                retained_count = retained_counts[support]
                output.append(
                    {
                        "approach": approach,
                        "approach_label": summary["approach_label"],
                        "source_approach": summary["source_approach"],
                        "side": side,
                        "support_species_count": support,
                        "discovery_species_total": discovery_total,
                        "gene_count": count,
                        "retained_gene_count": retained_count,
                        "lost_gene_count": count - retained_count,
                        "fraction_of_significant_genes": count / gene_total if gene_total else 0.0,
                        "retained_fraction_within_support": retained_count / count if count else 0.0,
                    }
                )
    return output


def build_joint_position_counts(summaries: List[dict], positions: List[dict]) -> List[dict]:
    positions_by_approach: Dict[str, List[dict]] = defaultdict(list)
    for row in positions:
        positions_by_approach[str(row["approach"])].append(row)

    output: List[dict] = []
    for summary in summaries:
        approach = str(summary["approach"])
        rows = positions_by_approach[approach]
        cell_positions = Counter((int(row["fg_support_count"]), int(row["bg_support_count"])) for row in rows)
        retained_cell_positions = Counter(
            (int(row["fg_support_count"]), int(row["bg_support_count"]))
            for row in rows
            if row["cluster_filter_status"] == "retained"
        )
        cell_genes: Dict[Tuple[int, int], set[str]] = defaultdict(set)
        retained_cell_genes: Dict[Tuple[int, int], set[str]] = defaultdict(set)
        for row in rows:
            key = (int(row["fg_support_count"]), int(row["bg_support_count"]))
            cell_genes[key].add(str(row["gene"]))
            if row["cluster_filter_status"] == "retained":
                retained_cell_genes[key].add(str(row["gene"]))

        for fg_support in range(1, int(summary["fg_discovery_species_total"]) + 1):
            for bg_support in range(1, int(summary["bg_discovery_species_total"]) + 1):
                key = (fg_support, bg_support)
                output.append(
                    {
                        "approach": approach,
                        "approach_label": summary["approach_label"],
                        "source_approach": summary["source_approach"],
                        "fg_support_count": fg_support,
                        "bg_support_count": bg_support,
                        "gene_count": len(cell_genes[key]),
                        "retained_gene_count": len(retained_cell_genes[key]),
                        "position_count": cell_positions[key],
                        "retained_position_count": retained_cell_positions[key],
                        "discarded_position_count": cell_positions[key] - retained_cell_positions[key],
                    }
                )
    return output


def build_balanced_threshold_counts(summaries: List[dict], positions: List[dict]) -> List[dict]:
    positions_by_approach: Dict[str, List[dict]] = defaultdict(list)
    for row in positions:
        positions_by_approach[str(row["approach"])].append(row)

    output: List[dict] = []
    for summary in summaries:
        approach = str(summary["approach"])
        rows = positions_by_approach[approach]
        largest_threshold = min(
            int(summary["fg_discovery_species_total"]),
            int(summary["bg_discovery_species_total"]),
        )
        for threshold in range(1, largest_threshold + 1):
            retained = [
                row
                for row in rows
                if int(row["fg_support_count"]) >= threshold
                and int(row["bg_support_count"]) >= threshold
            ]
            retained_after_filter = [
                row for row in retained if row["cluster_filter_status"] == "retained"
            ]
            output.append(
                {
                    "approach": approach,
                    "approach_label": summary["approach_label"],
                    "source_approach": summary["source_approach"],
                    "minimum_support_species_per_group": threshold,
                    "gene_count": len({row["gene"] for row in retained}),
                    "position_count": len(retained),
                    "retained_gene_count": len({row["gene"] for row in retained_after_filter}),
                    "retained_position_count": len(retained_after_filter),
                    "lost_gene_count": len({row["gene"] for row in retained})
                    - len({row["gene"] for row in retained_after_filter}),
                    "discarded_position_count": len(retained) - len(retained_after_filter),
                }
            )
    return output


def build_cluster_filter_summary(summaries: List[dict]) -> List[dict]:
    output: List[dict] = []
    for summary in summaries:
        position_total = int(summary["unique_significant_positions"])
        gene_total = int(summary["unique_significant_genes"])
        output.append(
            {
                "approach": summary["approach"],
                "approach_label": summary["approach_label"],
                "source_approach": summary["source_approach"],
                "significant_position_count": position_total,
                "retained_position_count": summary["retained_significant_positions"],
                "discarded_position_count": summary["discarded_clustered_positions"],
                "retained_position_fraction": int(summary["retained_significant_positions"])
                / position_total if position_total else 0.0,
                "significant_gene_count": gene_total,
                "retained_gene_count": summary["retained_significant_genes"],
                "lost_gene_count": summary["lost_clustered_genes"],
                "genes_with_clustered_positions": summary["genes_with_clustered_positions"],
                "retained_gene_fraction": int(summary["retained_significant_genes"])
                / gene_total if gene_total else 0.0,
                "cluster_density_threshold": summary["cluster_density_threshold"],
                "cluster_min_span": summary["cluster_min_span"],
                "cluster_min_positions": summary["cluster_min_positions"],
            }
        )
    return output


def main() -> None:
    args = parse_args()
    if not 0 < args.pvalue_threshold <= 1:
        raise SystemExit("--pvalue-threshold must be in (0, 1]")
    if not 0 < args.cluster_density_threshold <= 1:
        raise SystemExit("--cluster-density-threshold must be in (0, 1]")
    if args.cluster_min_span < 1:
        raise SystemExit("--cluster-min-span must be at least 1")
    if args.cluster_min_positions < 2:
        raise SystemExit("--cluster-min-positions must be at least 2")
    if not args.input_dir.is_dir():
        raise SystemExit(f"input directory does not exist: {args.input_dir}")

    approach_dirs = [
        path
        for path in args.input_dir.iterdir()
        if path.is_dir() and (path / "caas-pooled-events").is_dir()
    ]
    if not approach_dirs:
        raise SystemExit(f"no approach/caas-pooled-events directories found below {args.input_dir}")
    unknown_approaches = sorted(path.name for path in approach_dirs if path.name not in APPROACH_METADATA)
    if unknown_approaches:
        raise SystemExit(
            "unrecognised consolidated approach directories: " + ", ".join(unknown_approaches)
        )
    approach_dirs.sort(key=lambda path: APPROACH_METADATA[path.name]["order"])

    summaries: List[dict] = []
    positions: List[dict] = []
    for approach_dir in approach_dirs:
        print(f"Scanning {approach_dir.name}...", file=sys.stderr, flush=True)
        summary, approach_positions = scan_approach(approach_dir, args.pvalue_threshold)
        summaries.append(summary)
        positions.extend(approach_positions)

    print("Applying dense-position cluster filter...", file=sys.stderr, flush=True)
    apply_cluster_filter(
        summaries,
        positions,
        density_threshold=args.cluster_density_threshold,
        min_span=args.cluster_min_span,
        min_positions=args.cluster_min_positions,
    )

    genes: List[dict] = []
    for summary in summaries:
        approach_positions = [row for row in positions if row["approach"] == summary["approach"]]
        genes.extend(aggregate_genes(summary, approach_positions))

    gene_distribution = build_gene_support_distribution(summaries, genes)
    joint_position_counts = build_joint_position_counts(summaries, positions)
    balanced_counts = build_balanced_threshold_counts(summaries, positions)
    cluster_filter_summary = build_cluster_filter_summary(summaries)
    discarded_positions = [
        row for row in positions if row["cluster_filter_status"] != "retained"
    ]

    summary_columns = [
        "approach", "approach_label", "source_approach", "event_files_scanned", "files_with_event_rows",
        "raw_event_rows", "primary_event_rows", "significant_primary_event_rows",
        "unique_significant_genes", "unique_significant_positions",
        "fg_discovery_species_total", "bg_discovery_species_total",
        "fg_discovery_totals_observed", "bg_discovery_totals_observed",
        "max_significant_fg_support", "max_significant_bg_support",
        "cluster_density_threshold", "cluster_min_span", "cluster_min_positions",
        "retained_significant_positions", "discarded_clustered_positions",
        "genes_with_clustered_positions", "retained_significant_genes", "lost_clustered_genes",
        "duplicate_significant_position_rows", "malformed_rows", "pvalue_threshold",
    ]
    position_columns = [
        "approach", "approach_label", "source_approach", "gene", "position", "event_id",
        "fg_support_count", "bg_support_count", "balanced_support_count",
        "total_support_count", "fg_discovery_species_total", "bg_discovery_species_total",
        "positional_pvalue", "source_file",
        "cluster_filter_status",
    ]
    gene_columns = [
        "approach", "approach_label", "source_approach", "gene", "significant_position_count",
        "retained_position_count", "discarded_clustered_position_count", "cluster_filter_status",
        "max_fg_support_count", "max_bg_support_count",
        "max_balanced_support_same_position", "min_positional_pvalue",
        "max_retained_fg_support_count", "max_retained_bg_support_count",
        "max_retained_balanced_support_same_position",
        "fg_discovery_species_total", "bg_discovery_species_total",
    ]
    distribution_columns = [
        "approach", "approach_label", "source_approach", "side", "support_species_count",
        "discovery_species_total", "gene_count", "retained_gene_count", "lost_gene_count",
        "fraction_of_significant_genes", "retained_fraction_within_support",
    ]
    joint_columns = [
        "approach", "approach_label", "source_approach", "fg_support_count", "bg_support_count",
        "gene_count", "retained_gene_count", "position_count", "retained_position_count",
        "discarded_position_count",
    ]
    balanced_columns = [
        "approach", "approach_label", "source_approach", "minimum_support_species_per_group",
        "gene_count", "position_count", "retained_gene_count", "retained_position_count",
        "lost_gene_count", "discarded_position_count",
    ]
    cluster_summary_columns = [
        "approach", "approach_label", "source_approach",
        "significant_position_count", "retained_position_count", "discarded_position_count",
        "retained_position_fraction", "significant_gene_count", "retained_gene_count",
        "lost_gene_count", "genes_with_clustered_positions", "retained_gene_fraction",
        "cluster_density_threshold", "cluster_min_span", "cluster_min_positions",
    ]

    write_tsv(args.output_dir / "01_approach_summary.tsv", summary_columns, summaries)
    write_tsv(args.output_dir / "02_significant_positions.tsv", position_columns, positions)
    write_tsv(args.output_dir / "03_gene_maximum_support.tsv", gene_columns, genes)
    write_tsv(args.output_dir / "04_gene_counts_by_fg_bg_support.tsv", distribution_columns, gene_distribution)
    write_tsv(args.output_dir / "05_position_counts_by_joint_fg_bg_support.tsv", joint_columns, joint_position_counts)
    write_tsv(args.output_dir / "06_balanced_minimum_support_counts.tsv", balanced_columns, balanced_counts)
    write_tsv(args.output_dir / "07_cluster_filter_summary.tsv", cluster_summary_columns, cluster_filter_summary)
    write_tsv(args.output_dir / "08_cluster_discarded_positions.tsv", position_columns, discarded_positions)

    print(
        f"Wrote eight tables for {len(summaries)} approaches, "
        f"{len(genes)} approach-gene records, and {len(positions)} significant positions.",
        file=sys.stderr,
    )


if __name__ == "__main__":
    main()
