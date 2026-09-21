#!/usr/bin/env python3
"""Build matched enrichment inputs exclusively from declared complete summaries."""
import argparse
from pathlib import Path
from common import checked, read_json, read_tsv, require, resolve, sha, write_json, write_tsv


def build(cohort, output, summary_base=None, allow_incomplete=False):
    specs = read_tsv(cohort)
    require(specs, "Empty cohort")
    items, shared, missing_rows = {}, None, []
    for spec in specs:
        folder = summary_base / spec['hypothesis_id'] if summary_base else resolve(cohort.parent, spec["summary_dir"])
        checked(folder / "summary.lock.json", spec["summary_lock_sha256"])
        for rel, checksum in read_json(folder / "summary.lock.json").items(): checked(folder / rel, checksum)
        s = read_json(folder / "summary.json")
        require(s["status"] == "complete" or (allow_incomplete and s['status'] == 'partial'), f"Incomplete hypothesis: {folder}")
        if (folder/'non_completed_genes.tsv').exists():
            missing_rows.extend(read_tsv(folder/'non_completed_genes.tsv'))
        hid = s["hypothesis_id"]
        require(spec["hypothesis_id"] == hid and hid not in items, f"Duplicate/wrong hypothesis: {hid}")
        bg = {r["gene"] for r in read_tsv(folder / "background.tsv")}
        query = {r["gene"] for r in read_tsv(folder / "query.tsv")}
        require(query <= bg, f"Query outside testable background: {hid}")
        items[hid] = (s, bg, query)
        shared = bg if shared is None else shared & bg
    require(allow_incomplete or bool(shared), "No common coverage-testable genes")
    # All hypotheses in one matched comparison must have the same discovery
    # collection and analysis semantics. Different configs/pools are expected.
    for field in ("alignment_manifest_sha256", "discovery_filter_sha256", "downstream_sha256", "tool_sha256", "smoke"):
        require(len({item[0][field] for item in items.values()}) == 1, f"Incompatible cohort: {field}")
    output.mkdir(parents=True, exist_ok=True)
    write_tsv(output/'non_completed_genes.tsv', sorted(missing_rows, key=lambda r:(r['hypothesis_id'],r['gene'])),
              ['hypothesis_id','strategy_id','replicate_id','gene','status','reason'])
    write_tsv(output / "common_background.tsv", [dict(gene=g) for g in sorted(shared)], ["gene"])
    queries, summaries, inventory = [], [], []
    for hid, (s, bg, query) in sorted(items.items()):
        q = query & shared
        for g in sorted(q): queries.append(dict(hypothesis_id=hid, gene=g))
        for g in sorted(bg | query): inventory.append(dict(hypothesis_id=hid, gene=g, in_strategy_background=g in bg,
                                        in_common_background=g in shared, in_strategy_query=g in query, in_common_query=g in q))
        summaries.append(dict(s, common_background_genes=len(shared), common_query_genes=len(q), common_query_fraction=len(q)/len(shared) if shared else None,
                              excluded_query_genes=len(query - shared)))
    write_tsv(output / "common_queries.tsv", queries, ["hypothesis_id", "gene"])
    write_tsv(output / "background_intersections.tsv", inventory, ["hypothesis_id", "gene", "in_strategy_background", "in_common_background", "in_strategy_query", "in_common_query"])
    write_tsv(output / "strategy_comparison.tsv", summaries, list(summaries[0]))
    write_json(output / "cohort.json", dict(hypotheses=sorted(items), common_background_count=len(shared),
              execution_status='partial' if missing_rows else 'complete', non_completed_tasks=len(missing_rows),
              cohort_manifest_sha256=sha(cohort), discovery_filter_sha256=summaries[0]['discovery_filter_sha256'],
              source_settings_sha256=sorted({s['settings_sha256'] for s in summaries}),
              smoke=summaries[0]["smoke"], interpretation="Conditional on explicitly declared cohort and coverage-testable intersection"))


if __name__ == "__main__":
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--cohort", type=Path, required=True); p.add_argument("--output", type=Path, required=True)
    p.add_argument('--summary-base', type=Path)
    p.add_argument('--allow-incomplete', action='store_true')
    a = p.parse_args(); build(a.cohort, a.output, a.summary_base, a.allow_incomplete)
