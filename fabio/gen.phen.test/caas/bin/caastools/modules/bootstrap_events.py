#                      _              _
#                     | |            |
#   ___ __ _  __ _ ___| |_ ___   ___ | |___
#  / __/ _` |/ _` / __| __/ _ \ / _ \| / __|
# | (_| (_| | (_| \__ \ || (_) | (_) | \__ \
#  \___\__,_|\__,_|___/\__\___/ \___/|_|___/

"""Species- and amino-acid-aware summaries for CAAStools bootstrap.

Bootstrap hypotheses are useful for discovering a CAAS, but overlapping
hypotheses are not independent observations.  This module collapses positive
hypotheses into amino-acid-compatible events and counts every supporting
species once per event.

The event summary is intentionally optional.  The historical bootstrap table
remains unchanged, while callers may request a companion long-format table.
"""

from collections import defaultdict
from functools import lru_cache
from math import comb
import re


EVENT_HEADER = [
    "gene",
    "position",
    "event_id",
    "primary_event",
    "fg_amino_acids",
    "fg_dominant_amino_acid",
    "bg_amino_acids",
    "bg_dominant_amino_acid",
    "hypothesis_count",
    "hypotheses",
    "event_pattern",
    "fg_support_count",
    "fg_support_species_amino_acids",
    "fg_support_amino_acid_counts",
    "fg_observed_count",
    "fg_observed_species_amino_acids",
    "fg_observed_amino_acid_counts",
    "fg_discovery_count",
    "fg_discovery_species",
    "fg_support_over_observed",
    "fg_support_over_discovery",
    "bg_support_count",
    "bg_support_species_amino_acids",
    "bg_support_amino_acid_counts",
    "bg_observed_count",
    "bg_observed_species_amino_acids",
    "bg_observed_amino_acid_counts",
    "bg_discovery_count",
    "bg_discovery_species",
    "bg_support_over_observed",
    "bg_support_over_discovery",
    "balanced_support_observed",
    "balanced_support_discovery",
    "total_support_count",
    "conflicting_species_count",
    "conflicting_species_amino_acids",
    "event_nominal_pvalue",
    "positional_pvalue",
    "trait_config",
]


def natural_key(value):
    """Return a deterministic human/numeric sorting key (b_2 before b_10)."""

    return tuple(
        int(part) if part.isdigit() else part.lower()
        for part in re.split(r"(\d+)", str(value))
    )


def infer_fixed_discovery_groups(resampled_traits):
    """Infer stable FG/BG discovery pools and validate their orientation.

    Species-aware aggregation only has a clear interpretation when a species
    never changes side across hypotheses.  Standard bootstrap output does not
    require this constraint; it is enforced only when event summaries are
    requested.
    """

    fg_species = set(getattr(resampled_traits, "discovery_fg", []))
    bg_species = set(getattr(resampled_traits, "discovery_bg", []))

    for hypothesis in set(resampled_traits.trait2fg).union(resampled_traits.trait2bg):
        fg_species.update(resampled_traits.trait2fg.get(hypothesis, []))
        bg_species.update(resampled_traits.trait2bg.get(hypothesis, []))

    if not fg_species or not bg_species:
        raise ValueError("event summaries require non-empty FG and BG discovery pools")

    overlap = fg_species.intersection(bg_species)
    if overlap:
        raise ValueError(
            "event summaries require fixed FG/BG membership across hypotheses; "
            "species found on both sides: " + ",".join(sorted(overlap, key=natural_key))
        )

    return (
        sorted(fg_species, key=natural_key),
        sorted(bg_species, key=natural_key),
    )


def amino_acids_by_species(processed_position):
    """Extract species -> amino acid from a processed alignment position."""

    return {
        species: tagged_amino_acid.split("@", 1)[0]
        for species, tagged_amino_acid in processed_position.d.items()
    }


def observed_group(species_pool, species_to_amino_acid):
    """Return non-gap discovery species observed at the position."""

    return {
        species: species_to_amino_acid[species]
        for species in species_pool
        if species in species_to_amino_acid and species_to_amino_acid[species] != "-"
    }


def signatures_compatible(first, second):
    """Whether two signatures can be represented by one disjoint CAAS event."""

    combined_fg = set(first[0]).union(second[0])
    combined_bg = set(first[1]).union(second[1])
    return combined_fg.isdisjoint(combined_bg)


def classify_event_pattern(fg_amino_acids, bg_amino_acids):
    """Classify the final merged event using the CAAStools pattern scheme."""

    fg_count = len(set(fg_amino_acids))
    bg_count = len(set(bg_amino_acids))
    if fg_count == 1 and bg_count == 1:
        return "1"
    if fg_count == 1:
        return "2"
    if bg_count == 1:
        return "3"
    return "4"


def _maximal_cliques(signatures):
    """Enumerate deterministic maximal compatible signature groups.

    Pairwise compatibility is sufficient here: if every pair is compatible,
    the union of all FG amino acids is disjoint from the union of all BG amino
    acids.  Maximal cliques therefore represent all maximal valid event
    interpretations without an order-dependent greedy merge.
    """

    nodes = list(range(len(signatures)))
    neighbours = {
        node: {
            other
            for other in nodes
            if other != node and signatures_compatible(signatures[node], signatures[other])
        }
        for node in nodes
    }
    output = []

    def bron_kerbosch(current, candidates, excluded):
        if not candidates and not excluded:
            output.append(tuple(sorted(current)))
            return

        pivot_pool = candidates.union(excluded)
        if pivot_pool:
            pivot = max(
                sorted(pivot_pool),
                key=lambda node: len(candidates.intersection(neighbours[node])),
            )
            extension = candidates.difference(neighbours[pivot])
        else:
            extension = set(candidates)

        for node in sorted(extension):
            bron_kerbosch(
                current.union({node}),
                candidates.intersection(neighbours[node]),
                excluded.intersection(neighbours[node]),
            )
            candidates.remove(node)
            excluded.add(node)

    bron_kerbosch(set(), set(nodes), set())
    return sorted(set(output))


def _format_fraction(numerator, denominator):
    if denominator == 0:
        return "NA"
    return format(numerator / float(denominator), ".12g")


def _format_amino_acids(amino_acids):
    return ",".join(sorted(amino_acids, key=natural_key))


def _format_species(species_to_amino_acid):
    return ",".join(
        species + "=" + species_to_amino_acid[species]
        for species in sorted(species_to_amino_acid, key=natural_key)
    )


def _format_counts(species_to_amino_acid):
    counts = defaultdict(int)
    for amino_acid in species_to_amino_acid.values():
        counts[amino_acid] += 1
    return ",".join(
        amino_acid + "=" + str(counts[amino_acid])
        for amino_acid in sorted(counts, key=natural_key)
    )


def format_dominant_amino_acid(species_to_amino_acid):
    """Format the dominant AA proportion among unique supporting species.

    All co-dominant amino acids are retained in deterministic order. This
    avoids silently selecting one residue when two or more have the same
    maximum support.
    """

    if not species_to_amino_acid:
        return "NA"

    counts = defaultdict(int)
    for amino_acid in species_to_amino_acid.values():
        counts[amino_acid] += 1

    maximum_count = max(counts.values())
    total = len(species_to_amino_acid)
    return ",".join(
        amino_acid + "=" + format(maximum_count / float(total), ".12g")
        for amino_acid in sorted(counts, key=natural_key)
        if counts[amino_acid] == maximum_count
    )


def exact_signature_separation_pvalue(
    fg_amino_acids,
    bg_amino_acids,
    observed_fg,
    observed_bg,
):
    """Conditional exact p-value for the event's oriented AA separation.

    The test conditions on the observed species and on the number assigned to
    FG/BG.  It enumerates every reassignment and asks how often the minimum of
    the two side-specific signature-match fractions is at least as large as
    observed.  Because the signature was selected from the same data, this is
    a *nominal* event p-value; proteome-wide multiplicity remains downstream.
    """

    n_fg = len(observed_fg)
    n_bg = len(observed_bg)
    if n_fg == 0 or n_bg == 0:
        return None

    fg_amino_acids = set(fg_amino_acids)
    bg_amino_acids = set(bg_amino_acids)

    observed_fg_matches = sum(aa in fg_amino_acids for aa in observed_fg.values())
    observed_bg_matches = sum(aa in bg_amino_acids for aa in observed_bg.values())
    observed_score = min(
        observed_fg_matches / float(n_fg),
        observed_bg_matches / float(n_bg),
    )

    all_amino_acids = list(observed_fg.values()) + list(observed_bg.values())
    fg_category_count = sum(aa in fg_amino_acids for aa in all_amino_acids)
    bg_category_count = sum(aa in bg_amino_acids for aa in all_amino_acids)
    other_category_count = len(all_amino_acids) - fg_category_count - bg_category_count

    return _exact_category_assignment_pvalue(
        fg_category_count,
        bg_category_count,
        other_category_count,
        n_fg,
        n_bg,
        observed_score,
    )


@lru_cache(maxsize=None)
def _exact_category_assignment_pvalue(
    fg_category_count,
    bg_category_count,
    other_category_count,
    n_fg,
    n_bg,
    observed_score,
):
    """Exact label enumeration collapsed to three sufficient categories."""

    favourable = 0
    total = comb(fg_category_count + bg_category_count + other_category_count, n_fg)

    # x and y are the counts of FG-signature and BG-signature species assigned
    # to the permuted FG. The remaining FG slots must come from "other".
    for x in range(fg_category_count + 1):
        for y in range(bg_category_count + 1):
            other_in_fg = n_fg - x - y
            if not 0 <= other_in_fg <= other_category_count:
                continue

            fg_matches = x
            bg_matches = bg_category_count - y
            score = min(fg_matches / float(n_fg), bg_matches / float(n_bg))
            if score + 1e-12 < observed_score:
                continue

            favourable += (
                comb(fg_category_count, x)
                * comb(bg_category_count, y)
                * comb(other_category_count, other_in_fg)
            )

    return favourable / float(total)


def summarize_position_events(
    processed_position,
    positive_hypotheses,
    discovery_fg,
    discovery_bg,
):
    """Build ranked, species-deduplicated events for one position."""

    if not positive_hypotheses:
        return []

    species_to_amino_acid = amino_acids_by_species(processed_position)
    observed_fg = observed_group(discovery_fg, species_to_amino_acid)
    observed_bg = observed_group(discovery_bg, species_to_amino_acid)

    records_by_signature = defaultdict(list)
    for record in positive_hypotheses:
        signature = (
            frozenset(record["fg_amino_acids"]),
            frozenset(record["bg_amino_acids"]),
        )
        records_by_signature[signature].append(record)

    signatures = sorted(
        records_by_signature,
        key=lambda signature: (
            tuple(sorted(signature[0], key=natural_key)),
            tuple(sorted(signature[1], key=natural_key)),
        ),
    )

    events = []
    for clique in _maximal_cliques(signatures):
        event_signatures = [signatures[index] for index in clique]
        fg_amino_acids = set().union(*(signature[0] for signature in event_signatures))
        bg_amino_acids = set().union(*(signature[1] for signature in event_signatures))
        records = [
            record
            for signature in event_signatures
            for record in records_by_signature[signature]
        ]

        hypotheses = sorted({record["hypothesis"] for record in records}, key=natural_key)
        fg_support = {}
        bg_support = {}
        for record in records:
            fg_support.update(record["fg_species_amino_acids"])
            bg_support.update(record["bg_species_amino_acids"])

        conflicting_fg = {
            species: amino_acid
            for species, amino_acid in observed_fg.items()
            if amino_acid in bg_amino_acids
        }
        conflicting_bg = {
            species: amino_acid
            for species, amino_acid in observed_bg.items()
            if amino_acid in fg_amino_acids
        }

        fg_observed_fraction = len(fg_support) / float(len(observed_fg)) if observed_fg else 0.0
        bg_observed_fraction = len(bg_support) / float(len(observed_bg)) if observed_bg else 0.0
        fg_discovery_fraction = len(fg_support) / float(len(discovery_fg))
        bg_discovery_fraction = len(bg_support) / float(len(discovery_bg))
        nominal_pvalue = exact_signature_separation_pvalue(
            fg_amino_acids,
            bg_amino_acids,
            observed_fg,
            observed_bg,
        )

        events.append({
            "fg_amino_acids": fg_amino_acids,
            "bg_amino_acids": bg_amino_acids,
            "hypotheses": hypotheses,
            "event_pattern": classify_event_pattern(fg_amino_acids, bg_amino_acids),
            "fg_support": fg_support,
            "bg_support": bg_support,
            "observed_fg": observed_fg,
            "observed_bg": observed_bg,
            "discovery_fg": list(discovery_fg),
            "discovery_bg": list(discovery_bg),
            "fg_support_over_observed": fg_observed_fraction,
            "bg_support_over_observed": bg_observed_fraction,
            "fg_support_over_discovery": fg_discovery_fraction,
            "bg_support_over_discovery": bg_discovery_fraction,
            "balanced_support_observed": min(fg_observed_fraction, bg_observed_fraction),
            "balanced_support_discovery": min(fg_discovery_fraction, bg_discovery_fraction),
            "conflicting_fg": conflicting_fg,
            "conflicting_bg": conflicting_bg,
            "event_nominal_pvalue": nominal_pvalue,
        })

    events.sort(
        key=lambda event: (
            -event["balanced_support_discovery"],
            -(len(event["fg_support"]) + len(event["bg_support"])),
            event["event_nominal_pvalue"] if event["event_nominal_pvalue"] is not None else 2.0,
            _format_amino_acids(event["fg_amino_acids"]),
            _format_amino_acids(event["bg_amino_acids"]),
        )
    )
    return events


def serialize_position_events(
    gene,
    position,
    events,
    positional_pvalue,
    trait_config,
):
    """Serialize ranked events as rows matching :data:`EVENT_HEADER`."""

    rows = []
    for index, event in enumerate(events, start=1):
        conflicting = []
        if event["conflicting_fg"]:
            conflicting.append("FG:" + _format_species(event["conflicting_fg"]))
        if event["conflicting_bg"]:
            conflicting.append("BG:" + _format_species(event["conflicting_bg"]))

        nominal_pvalue = event["event_nominal_pvalue"]
        rows.append([
            gene,
            str(position),
            "e_" + str(index),
            "yes" if index == 1 else "no",
            _format_amino_acids(event["fg_amino_acids"]),
            format_dominant_amino_acid(event["fg_support"]),
            _format_amino_acids(event["bg_amino_acids"]),
            format_dominant_amino_acid(event["bg_support"]),
            str(len(event["hypotheses"])),
            ",".join(event["hypotheses"]),
            event["event_pattern"],
            str(len(event["fg_support"])),
            _format_species(event["fg_support"]),
            _format_counts(event["fg_support"]),
            str(len(event["observed_fg"])),
            _format_species(event["observed_fg"]),
            _format_counts(event["observed_fg"]),
            str(len(event["discovery_fg"])),
            ",".join(event["discovery_fg"]),
            _format_fraction(len(event["fg_support"]), len(event["observed_fg"])),
            _format_fraction(len(event["fg_support"]), len(event["discovery_fg"])),
            str(len(event["bg_support"])),
            _format_species(event["bg_support"]),
            _format_counts(event["bg_support"]),
            str(len(event["observed_bg"])),
            _format_species(event["observed_bg"]),
            _format_counts(event["observed_bg"]),
            str(len(event["discovery_bg"])),
            ",".join(event["discovery_bg"]),
            _format_fraction(len(event["bg_support"]), len(event["observed_bg"])),
            _format_fraction(len(event["bg_support"]), len(event["discovery_bg"])),
            format(event["balanced_support_observed"], ".12g"),
            format(event["balanced_support_discovery"], ".12g"),
            str(len(event["fg_support"]) + len(event["bg_support"])),
            str(len(event["conflicting_fg"]) + len(event["conflicting_bg"])),
            ";".join(conflicting),
            "NA" if nominal_pvalue is None else format(nominal_pvalue, ".12g"),
            str(positional_pvalue),
            trait_config,
        ])
    return rows
