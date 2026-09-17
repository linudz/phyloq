#!/usr/bin/env python3
"""Validate immutable inputs and bind a run ID to one explicit manifest/inventory."""
from __future__ import annotations
import argparse
from pathlib import Path
from common import (ROOT, VERSION, checked, cycle_read, digest, identifier, pool_read,
                    read_json, read_tsv, require, resolve, sha, tree_hash, verify_design,
                    write_json, write_tsv)


def validate(design, selection, alignments, settings, run_id, output, approvals=None, smoke=False, direct_discovery=False):
    lock = verify_design(design)
    identifier(run_id)
    settings_data = read_json(settings)
    require(settings_data["discovery"] == read_json(design / "settings.json")["discovery"],
            "Discovery settings differ from the frozen design; create a reviewed new design")
    require(settings_data["filter"] == read_json(design / "settings.json")["filter"], "Filter settings changed")
    tool_hash = tree_hash(ROOT / "bin/caastools")
    require(tool_hash == lock["caastools_sha256"], "CAAStools version changed after freezing")
    canonical = {r["hypothesis_id"]: r for r in read_tsv(design / "execution_manifest.tsv")}
    rows, inventory = read_tsv(selection), read_tsv(alignments)
    require(rows and inventory, "Empty explicit selection or alignment inventory")
    hids = [r["hypothesis_id"] for r in rows]
    require(len(hids) == len(set(hids)), "Duplicate hypothesis in selected manifest")
    genes = [r["gene_id"] for r in inventory]
    require(len(genes) == len(set(genes)), "Duplicate gene IDs in alignment manifest")
    basenames = [Path(r["alignment_path"]).name for r in inventory]
    require(len(basenames) == len(set(basenames)), "Duplicate alignment basenames")
    for a in inventory:
        identifier(a["gene_id"])
        require(a["format"] == settings_data["discovery"]["fmt"], "Alignment format mismatch")
        checked(resolve(alignments.parent, a["alignment_path"]), a["sha256"])
        require(Path(a["alignment_path"]).name.split(".")[0] == a["gene_id"], "Ambiguous gene symbol mapping")
    sources = {r["species"] for r in read_tsv(design / "sources/phenotype.tsv")}
    for r in rows:
        hid = identifier(r["hypothesis_id"])
        require(hid in canonical and r == canonical[hid], f"Selected row differs from frozen design: {hid}")
        pools = pool_read(checked(design / r["pool_config"], r["pool_sha256"]))
        require(set(sum(pools.values(), [])) <= sources, f"Unknown phenotype species in {hid}")
        require(not set(pools["FG"]) & set(pools["BG"]), f"Conflicting membership: {hid}")
        cc = cycle_read(checked(design / r["hypotheses_config"], r["config_sha256"]), pools)
        require(len(cc) == int(r["selected_cycles"]), f"Cycle count mismatch: {hid}")
    runtime_lock = dict(design_sha256=sha(design / "design.lock.json"), selection_sha256=sha(selection),
                        alignment_manifest_sha256=sha(alignments), settings_sha256=sha(settings),
                        caastools_sha256=tool_hash, implementation_sha256=tree_hash(ROOT / "scripts"),
                        run_id=run_id, hypotheses=hids, genes=genes, smoke=smoke,
                        direct_discovery=direct_discovery)
    if smoke:
        require(len(inventory) <= 3, "Smoke is limited to three synthetic alignments")
        require(all(resolve(alignments.parent, a["alignment_path"]).is_relative_to(ROOT / "tests/fixtures") for a in inventory),
                "Smoke may use only explicit local test fixtures")
    elif not direct_discovery:
        require(approvals is not None, "Production requires a reviewed approval JSON; see inputs/approval.template.json")
        approval = read_json(approvals)
        for key in ("design_sha256", "selection_sha256", "alignment_manifest_sha256", "settings_sha256"):
            require(approval.get(key) == runtime_lock[key], f"Approval missing/mismatched: {key}")
        for gate in ("pool_review", "technical_review", "production_review"):
            require(approval.get(gate) is True, f"Approval gate not recorded: {gate}")
        require(bool(approval.get("approved_by")), "Approval needs reviewer name")
        if any(r["strategy_id"] in ("R0", "R1") for r in rows):
            require(approval.get("statistical_review") is True, "Null discovery needs statistical review")
        if any(int(r["replicate_id"]) > 19 for r in rows):
            require(bool(approval.get("pilot_cost_report_sha256")), "Full null extension needs reviewed pilot cost")
    output.mkdir(parents=True, exist_ok=True)
    write_json(output / "runtime.lock.json", runtime_lock)
    # Stream the task index instead of retaining millions of job rows in RAM.
    def job_rows():
        for r in rows:
            for a in inventory:
                job_id = r["hypothesis_id"] + "__" + a["gene_id"]
                job = dict(r, gene_id=a["gene_id"], alignment_path=str(resolve(alignments.parent, a["alignment_path"])),
                           alignment_sha256=a["sha256"], alignment_format=a["format"],
                           pool_path=str(design / r["pool_config"]), config_path=str(design / r["hypotheses_config"]),
                           settings=settings_data, settings_sha256=runtime_lock['settings_sha256'], tool_sha256=tool_hash,
                           run_id=run_id, design_sha256=runtime_lock["design_sha256"],
                           selection_sha256=runtime_lock["selection_sha256"], smoke=smoke,
                           direct_discovery=direct_discovery)
                job_file = output / "jobs" / f"{job_id}.json"
                write_json(job_file, job)
                yield dict(hypothesis_id=r["hypothesis_id"], strategy_id=r["strategy_id"],
                           replicate_id=r["replicate_id"], gene_id=a["gene_id"], job_file=str(job_file.resolve()),
                           alignment_path=job["alignment_path"], pool_path=job["pool_path"], config_path=job["config_path"])
    write_tsv(output / "jobs.tsv", job_rows(), ["hypothesis_id", "strategy_id", "replicate_id", "gene_id", "job_file", "alignment_path", "pool_path", "config_path"])
    write_json(output / "input_validation.json", dict(status="pass", hypotheses=len(rows), alignments=len(inventory),
              expected_tasks=len(rows)*len(inventory), historical_results_imported=[], runtime_lock=runtime_lock))
    return runtime_lock


if __name__ == "__main__":
    p = argparse.ArgumentParser(description=__doc__)
    for name in ("design", "selection", "alignments", "settings", "output"):
        p.add_argument("--" + name, type=Path, required=True)
    p.add_argument("--run-id", required=True)
    p.add_argument("--approvals", type=Path)
    p.add_argument("--smoke", action="store_true")
    p.add_argument("--direct-discovery", action="store_true", help="Explicit CAAS-only launch without pilot/statistical approval prerequisites")
    a = p.parse_args()
    validate(a.design.resolve(), a.selection.resolve(), a.alignments.resolve(), a.settings.resolve(),
             a.run_id, a.output.resolve(), a.approvals, a.smoke, a.direct_discovery)
