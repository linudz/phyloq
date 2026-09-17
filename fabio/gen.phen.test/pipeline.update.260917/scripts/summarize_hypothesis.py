#!/usr/bin/env python3
"""Assemble every expected result before applying the frozen CAAS filters."""
from __future__ import annotations
import argparse
import math
from pathlib import Path
from common import ROOT, checked, digest, read_json, read_tsv, require, sha, write_json, write_tsv
from legacy_support import dense_cluster_positions, preferred_position_record, is_primary


def summarize(hypothesis_id, alignment_manifest, artifacts, output):
    inventory = read_tsv(alignment_manifest)
    expected = {r["gene_id"] for r in inventory}
    receipts = {}
    for path in sorted(Path(artifacts).rglob("*.receipt.json")):
        receipt = read_json(path)
        require(receipt["hypothesis_id"] == hypothesis_id, f"Wrong hypothesis receipt: {path}")
        gene = receipt["gene_id"]
        require(gene not in receipts, f"Duplicate gene output: {hypothesis_id}/{gene}")
        receipts[gene] = (path, receipt)
    output.mkdir(parents=True, exist_ok=True)
    completeness = dict(hypothesis_id=hypothesis_id, expected_count=len(expected), observed_count=len(receipts),
                        missing=sorted(expected - receipts.keys()), unexpected=sorted(receipts.keys() - expected),
                        status="complete" if expected == receipts.keys() else "incomplete")
    write_json(output / "execution_completeness.json", completeness)
    require(completeness["status"] == "complete", f"Incomplete hypothesis: {completeness}; summaries blocked")
    positions, gene_rows = [], []
    anchor = next(iter(receipts.values()))[1]
    filt = anchor["settings"]["filter"]
    for a in inventory:
        gene = a["gene_id"]
        receipt_path, receipt = receipts[gene]
        require(receipt["status"] == "success", f"Failed output for {gene}")
        require(receipt["alignment_sha256"] == a["sha256"], f"Wrong alignment version for {gene}")
        for field in ("settings_sha256", "config_sha256", "pool_sha256", "tool_sha256", "design_sha256"):
            require(receipt[field] == anchor[field], f"Mixed provenance within {hypothesis_id}: {field}")
        event_path = checked(receipt_path.parent / f"{gene}.pooled.caas.events.tsv", receipt["events_sha256"])
        checked(receipt_path.parent / f"{gene}.pooled.caas.tsv", receipt["legacy_sha256"])
        selected = {}
        for event in read_tsv(event_path):
            require(event["gene"] == gene, f"Wrong gene field in event: {event_path}")
            p = float(event["positional_pvalue"])
            require(math.isfinite(p) and 0 <= p <= 1, f"Invalid positional p value: {event_path}")
            if not is_primary(event["primary_event"]) or p >= filt["pvalue_threshold"]: continue
            pos, fg, bg = int(event["position"]), int(event["fg_support_count"]), int(event["bg_support_count"])
            require(pos >= 0 and fg >= 0 and bg >= 0, "Invalid position/support")
            row = dict(gene=gene, position=pos, event_id=event["event_id"], positional_pvalue=p,
                       fg_support_count=fg, bg_support_count=bg, balanced_support_count=min(fg, bg), total_support_count=fg + bg)
            if pos not in selected or preferred_position_record(row, selected[pos]): selected[pos] = row
        discarded = dense_cluster_positions(selected, filt["cluster_density_threshold"], filt["cluster_min_span"], filt["cluster_min_positions"])
        qualified = []
        for pos, row in sorted(selected.items()):
            row.update(hypothesis_id=hypothesis_id, strategy_id=anchor["strategy_id"], replicate_id=anchor["replicate_id"],
                       cluster_filter_status="discarded_dense_cluster" if pos in discarded else "retained")
            # Both species counts must belong to this SAME retained position.
            row["in_query"] = pos not in discarded and min(row["fg_support_count"], row["bg_support_count"]) >= filt["minimum_balanced_support"]
            if row["in_query"]: qualified.append(pos)
            positions.append(row)
        cov = receipt["coverage"]
        require(not selected or cov["coverage_testable"], f"Positive but declared untestable: {gene}")
        gene_rows.append(dict(hypothesis_id=hypothesis_id, gene=gene, completed=True,
            coverage_testable=cov["coverage_testable"], coverage_testable_positions=cov["coverage_testable_positions"],
            present_fg=cov["present_fg"], present_bg=cov["present_bg"], missing_fg=",".join(cov["missing_fg"]), missing_bg=",".join(cov["missing_bg"]),
            in_query=bool(qualified), significant_positions=len(selected), retained_positions=len(selected) - len(discarded),
            balanced_retained_positions=len(qualified), alignment_sha256=a["sha256"], events_sha256=receipt["events_sha256"],
            receipt_sha256=sha(receipt_path), elapsed_seconds=receipt["elapsed_seconds"], cpu_seconds=receipt["cpu_seconds"],
            max_rss_bytes=receipt["max_rss_bytes"], output_bytes=receipt["output_bytes"]))
    write_tsv(output / "positions.tsv", positions, ["hypothesis_id", "strategy_id", "replicate_id", "gene", "position", "event_id", "positional_pvalue", "fg_support_count", "bg_support_count", "balanced_support_count", "total_support_count", "cluster_filter_status", "in_query"])
    write_tsv(output / "gene_inventory.tsv", gene_rows, list(gene_rows[0]))
    for name, key in (("query.tsv", "in_query"), ("background.tsv", "coverage_testable"), ("completed_background.tsv", "completed")):
        write_tsv(output / name, [r for r in gene_rows if r[key]], ["gene"])
    nbg = sum(r["coverage_testable"] for r in gene_rows)
    nq = sum(r["in_query"] for r in gene_rows)
    summary = dict(hypothesis_id=hypothesis_id, strategy_id=anchor["strategy_id"], replicate_id=anchor["replicate_id"],
        null_family=anchor["null_family"], completed_genes=len(gene_rows), background_genes=nbg, query_genes=nq,
        query_fraction=nq / nbg if nbg else None, significant_positions=len(positions),
        retained_positions=sum(r["cluster_filter_status"] == "retained" for r in positions),
        elapsed_seconds=sum(r["elapsed_seconds"] for r in gene_rows), cpu_seconds=sum(r["cpu_seconds"] for r in gene_rows),
        output_bytes=sum(r["output_bytes"] for r in gene_rows), settings_sha256=anchor["settings_sha256"],
        discovery_filter_sha256=digest({k:anchor['settings'][k] for k in ('discovery','filter')}),
        downstream_sha256=digest({name:sha(ROOT/'scripts'/name) for name in ('summarize_hypothesis.py','legacy_support.py','common.py')}),
        tool_sha256=anchor["tool_sha256"], config_sha256=anchor["config_sha256"], pool_sha256=anchor["pool_sha256"],
        alignment_manifest_sha256=sha(alignment_manifest), smoke=anchor["smoke"], status="complete")
    write_json(output / "summary.json", summary)
    write_json(output / "summary.lock.json", {p.name: sha(p) for p in sorted(output.iterdir()) if p.is_file() and p.name != "summary.lock.json"})
    return summary


if __name__ == "__main__":
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--hypothesis", required=True)
    for name in ("alignments", "artifacts", "output"): p.add_argument("--" + name, type=Path, required=True)
    a = p.parse_args(); summarize(a.hypothesis, a.alignments, a.artifacts, a.output)
