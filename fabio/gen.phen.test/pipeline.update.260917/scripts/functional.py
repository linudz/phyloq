#!/usr/bin/env python3
"""Gene-list controls, frozen enrichment, and declared conditional benchmarks.

Default: manifests and an explicit blocked status (no network). Two opt-in
engines: version-checked raw g:Profiler cache, or an approved local annotation
snapshot with exact hypergeometric/BONFERRONI inference. The latter is a new
implementation and runs ALL observed and random lists anew, never mixing its
p-values with historical g:SCS. A toy-only snapshot is used in tests.
"""
from __future__ import annotations
import argparse
import math
import random
from collections import defaultdict
from pathlib import Path
from common import digest, read_json, read_tsv, require, sha, write_json, write_tsv

TERM_COLUMNS = ["hypothesis_id", "term_id", "source", "term_name", "status", "intersection_size", "term_size", "effective_query_size", "effective_domain_size", "fold_enrichment", "p_value", "p_value_adjusted", "significant", "supporting_genes"]


def gene_lists(queries, background, n, seed):
    """Exact-size independent lists; without replacement WITHIN each list."""
    require("P0_000" in queries, "P0_000 required for G0/G1 list-size controls")
    size = len(queries["P0_000"])
    require(size > 0, "P0 has an empty common query; zero-sized controls are not informative")
    controls, statuses = {}, []
    for family, universe in (("G0", background), ("G1", queries.get("N2_000")), ("G1_N3", queries.get("N3_000"))):
        if universe is None:
            statuses.append(dict(family=family, status="missing_control_query", requested_size=size)); continue
        if len(universe) < size:
            statuses.append(dict(family=family, status="eligible_query_too_small", requested_size=size)); continue
        for i in range(1, n + 1):
            rng = random.Random(int(digest([seed, family, i]), 16))
            controls[f"{family}_{i:04d}"] = set(rng.sample(sorted(universe), size))
        statuses.append(dict(family=family, status="generated", requested_size=size, replicates=n))
    return controls, statuses


def payload(query, background, cfg):
    return dict(organism=cfg["organism"], query=sorted(query), background=sorted(background), sources=cfg["sources"],
                significance_threshold_method=cfg["correction"], domain_scope=cfg["domain_scope"],
                user_threshold=cfg["threshold"], all_results=True, ordered=False, no_evidences=False)


def read_snapshot(snapshot, cfg):
    meta = read_json(snapshot / "metadata.json")
    require(meta["annotation_version"] == cfg["annotation_version"], "Annotation version mismatch")
    require(meta["organism"] == cfg["organism"], "Snapshot organism mismatch")
    require(set(meta["sources"]) == set(cfg["sources"]), "Snapshot sources mismatch")
    require(sha(snapshot / "term_membership.tsv") == meta["term_membership_sha256"], "Snapshot checksum mismatch")
    terms = {}
    for r in read_tsv(snapshot / "term_membership.tsv"):
        term = terms.setdefault(r["term_id"], dict(source=r["source"], term_name=r["term_name"], genes=set()))
        require(term["source"] == r["source"] and term["term_name"] == r["term_name"], "Conflicting term metadata")
        require(r["source"] in cfg["sources"], "Undeclared annotation source")
        term["genes"].add(r["gene"])
    require(bool(terms), "Empty snapshot")
    return terms, meta


def local_enrich(hid, query, background, terms, cfg):
    from scipy.stats import hypergeom
    annotated = set().union(*(t["genes"] for t in terms.values()))
    effective_bg = background & annotated
    effective_q = query & effective_bg
    M, N = len(effective_bg), len(effective_q)
    # Fixed family: all snapshot terms with nonzero background membership.
    active = {tid: t for tid, t in terms.items() if t["genes"] & effective_bg}
    rows = []
    for tid, term in sorted(active.items()):
        members = term["genes"] & effective_bg
        hit = members & effective_q
        p = float(hypergeom.sf(len(hit) - 1, M, len(members), N)) if N else 1.0
        adj = min(1.0, p * len(active))
        rows.append(dict(hypothesis_id=hid, term_id=tid, source=term["source"], term_name=term["term_name"],
            status="ok" if N else "empty_effective_query", intersection_size=len(hit), term_size=len(members),
            effective_query_size=N, effective_domain_size=M, fold_enrichment=(len(hit)/N)/(len(members)/M) if N else None,
            p_value=p, p_value_adjusted=adj, significant=bool(N and adj < cfg["threshold"]), supporting_genes=",".join(sorted(hit))))
    return rows


def cached_enrich(hid, request, cache, cfg):
    key = digest(dict(payload=request, annotation_version=cfg["annotation_version"]))
    path = cache / (key + ".json")
    require(path.is_file(), f"Missing frozen g:Profiler cache: {key}")
    cached = read_json(path)
    require(cached["request"] == request and cached["annotation_version"] == cfg["annotation_version"], "Cache request/version mismatch")
    raw = cached["response"]
    require(digest(raw) == cached["response_sha256"], "Cache response checksum mismatch")
    require(raw["meta"]["version"] == cfg["annotation_version"], "Returned annotation version mismatch")
    mapping = raw["meta"]["genes_metadata"]["query"]["query_1"]
    ensgs = mapping.get("ensgs", [])
    reverse = {e: s for s, ids in mapping.get("mapping", {}).items() for e in ids}
    rows = []
    for r in raw.get("result", []):
        evidence = r.get("intersections", [])
        require(len(evidence) == len(ensgs), "Cache missing auditable gene intersections")
        N, M, k, K = r["query_size"], r["effective_domain_size"], r["intersection_size"], r["term_size"]
        rows.append(dict(hypothesis_id=hid, term_id=r["native"], source=r["source"], term_name=r["name"], status="ok",
            intersection_size=k, term_size=K, effective_query_size=N, effective_domain_size=M,
            fold_enrichment=(k/N)/(K/M) if N and K and M else None,
            p_value=None, p_value_adjusted=r["p_value"], significant=bool(r["significant"]),
            supporting_genes=",".join(sorted({reverse.get(e, e) for e, hit in zip(ensgs, evidence) if hit}))))
    return rows


def metrics(hid, rows, query, background, stats):
    finite = [r for r in rows if r["status"] == "ok" and r["fold_enrichment"] is not None]
    significant = [r for r in finite if r["significant"]]
    wanted = stats["term_universe"]
    by_id = {r["term_id"]: r for r in finite}
    missing = sorted(set(wanted) - by_id.keys())
    # No substitution of zero for nonsignificance or absent service terms.
    mean_fe = sum(by_id[t]["fold_enrichment"] for t in wanted)/len(wanted) if wanted and not missing else None
    groups = {stats["term_groups"].get(r["term_id"]) for r in significant} - {None}
    return dict(hypothesis_id=hid, query_size=len(query), background_size=len(background),
        query_fraction=len(query)/len(background), significant_terms=len(significant),
        mean_fold_enrichment=mean_fe, nonredundant_significant_groups=len(groups) if stats["term_groups"] else None,
        significant_supporting_genes=len(set().union(*(set(r["supporting_genes"].split(",")) - {""} for r in significant))) if significant else 0,
        missing_declared_terms=",".join(missing), status="empty_query" if not query else "missing_declared_terms" if missing else "ok")


def conditional_tail(observed, random_values):
    from scipy.stats import beta
    import numpy as np
    B = len(random_values); require(B > 0, "No independent null hypotheses")
    k = sum(x >= observed for x in random_values)
    lo = float(beta.ppf(0.025, k, B-k+1)) if k else 0.0
    hi = float(beta.ppf(0.975, k+1, B-k)) if k < B else 1.0
    median = float(np.median(random_values))
    return dict(observed=observed, null_hypotheses=B, exceedances=k, conditional_tail_fraction=(1+k)/(B+1),
                resolution=1/(B+1), null_median=median, observed_minus_null_median=observed-median,
                null_q025=float(np.quantile(random_values,.025)), null_q975=float(np.quantile(random_values,.975)),
                monte_carlo_exceedance_ci_lower=lo, monte_carlo_exceedance_ci_upper=hi,
                mc_se_approx=math.sqrt(((k+1)/(B+2))*(1-(k+1)/(B+2))/B))


def run(matched, settings, output, snapshot=None, cache=None, discovery_only=False):
    cfg = read_json(settings); ecfg, stats = cfg["enrichment"], cfg["statistics"]
    cohort = read_json(matched / "cohort.json")
    background = {r["gene"] for r in read_tsv(matched / "common_background.tsv")}
    queries = {hid: set() for hid in cohort["hypotheses"]}
    for r in read_tsv(matched / "common_queries.tsv"): queries[r["hypothesis_id"]].add(r["gene"])
    output.mkdir(parents=True, exist_ok=True)
    if discovery_only:
        controls, control_status = {}, [dict(family="G0/G1", status="not_requested_discovery_only")]
    elif "P0_000" in queries and queries["P0_000"]:
        controls, control_status = gene_lists(queries, background, stats["gene_list_replicates"], stats["gene_list_seed"])
    else:
        controls, control_status = {}, [dict(family="G0/G1", status="P0_missing_or_empty")]
    queries.update(controls)
    write_json(output / "gene_list_control_status.json", control_status)
    list_rows = []
    for hid, query in sorted(queries.items()):
        path = output / "gene_lists" / f"{hid}.tsv"
        write_tsv(path, [dict(gene=g) for g in sorted(query)], ["gene"])
        list_rows.append(dict(hypothesis_id=hid, gene_count=len(query), genes_file=str(path.relative_to(output)), sha256=sha(path),
                             origin="gene_list_control" if hid in controls else "CAAS", seed=stats["gene_list_seed"] if hid in controls else ""))
    write_tsv(output / "gene_list_manifest.tsv", list_rows, ["hypothesis_id", "gene_count", "genes_file", "sha256", "origin", "seed"])
    write_json(output / "settings.resolved.json", cfg)
    statuses = []
    if discovery_only or ecfg["backend"] == "disabled":
        write_json(output / "enrichment_status.json", dict(
            status="not_requested_discovery_only" if discovery_only else "blocked_annotation_not_frozen", analyzed_queries=0,
            required="CAAS discovery completed separately from enrichment/statistical inference. No live requests were sent." if discovery_only else
                     "Approve frozen annotation version/engine and primary statistic. No live requests were sent."))
        for name in ("all_terms.tsv", "random_replicate_metrics.tsv", "gene_list_control_metrics.tsv", "conditional_tail_fractions.tsv"):
            write_tsv(output / name, [], ["status"])
        return
    require(bool(ecfg["annotation_version"]), "Enrichment requires an explicit annotation version")
    terms = None
    if ecfg["backend"] == "local-hypergeom-bonferroni":
        require(ecfg["correction"] == "bonferroni" and ecfg.get("implementation_approved") is True, "New local implementation needs recorded review and Bonferroni settings")
        require(snapshot is not None, "Frozen annotation snapshot required")
        terms, meta = read_snapshot(snapshot, ecfg)
        require(not meta.get("toy_only") or cohort["smoke"], "Toy annotations cannot be used for biological results")
        write_json(output / "annotation_metadata.json", dict(meta, metadata_sha256=sha(snapshot / "metadata.json"),
            consistency_check="Exact hypergeometric verified in unit tests; all observed/control lists rerun. Historical g:SCS values excluded."))
    else:
        require(ecfg["backend"] == "gprofiler-cache" and ecfg["correction"] == "g_SCS" and cache is not None, "Unknown/missing frozen enrichment backend")
    all_rows, metric_rows = [], []
    for hid, query in sorted(queries.items()):
        if not query:
            rows = []; statuses.append(dict(hypothesis_id=hid, status="empty_query"))
        elif terms is not None:
            rows = local_enrich(hid, query, background, terms, ecfg)
        else:
            request = payload(query, background, ecfg)
            key = digest(dict(payload=request, annotation_version=ecfg["annotation_version"]))
            write_json(output / "requests" / (key + ".json"), dict(hypothesis_id=hid, payload=request, annotation_version=ecfg["annotation_version"]))
            try: rows = cached_enrich(hid, request, cache, ecfg)
            except ValueError as exc:
                statuses.append(dict(hypothesis_id=hid, status="failed_or_missing_cache", reason=str(exc))); continue
        all_rows.extend(rows)
        metric_rows.append(metrics(hid, rows, query, background, stats))
    write_tsv(output / "all_terms.tsv", all_rows, TERM_COLUMNS)
    for name, subset in (("random_replicate_metrics.tsv", [r for r in metric_rows if not r["hypothesis_id"].startswith("G")]),
                         ("gene_list_control_metrics.tsv", [r for r in metric_rows if r["hypothesis_id"].startswith("G")])):
        write_tsv(output / name, subset, list(metric_rows[0]) if metric_rows else ["hypothesis_id", "status"])
    write_json(output / "query_status.json", statuses)
    failed = [s for s in statuses if s["status"] == "failed_or_missing_cache"]
    tails = []
    inferential_status = "pending_statistical_review"
    if stats["approved"]:
        require(not failed and len(metric_rows) == len(queries), "Enrichment incomplete; conditional summaries blocked")
        metric = stats["primary_metric"]
        require(metric in ("mean_fold_enrichment", "nonredundant_significant_groups", "significant_supporting_genes", "query_fraction"), "Declare a supported primary metric")
        require(stats["multiple_testing"] in ("none-single-comparison", "holm"), "Declare inferential family correction")
        by_id = {r["hypothesis_id"]: r for r in metric_rows}
        if "P0_000" in by_id:
            obs = by_id["P0_000"][metric]
            require(obs is not None, "Observed metric undefined; no threshold-based zero substitution")
            families = stats.get('null_families', [])
            require(families and len(set(families)) == len(families), 'Prespecify distinct null_families')
            for family in families:
                null = [r[metric] for hid, r in by_id.items() if hid.startswith(family + "_") and
                        (family != "G1" or not hid.startswith("G1_N3_"))]
                expected_count = stats.get('expected_null_counts', {}).get(family)
                require(expected_count is not None and len(null) == expected_count,
                        f'Incomplete declared null family {family}: {len(null)} results; expected {expected_count}')
                require(all(v is not None for v in null), f"Undefined null metric: {family}")
                tails.append(dict(null_family=family, metric=metric, **conditional_tail(obs, null)))
            if stats["multiple_testing"] == "none-single-comparison": require(len(tails) <= 1, "Multiple tests require a declared correction")
            if tails and stats["multiple_testing"] == "holm":
                ordered = sorted(range(len(tails)), key=lambda i: tails[i]["conditional_tail_fraction"])
                running = 0
                for rank, i in enumerate(ordered):
                    running = max(running, min(1, (len(tails)-rank)*tails[i]["conditional_tail_fraction"]))
                    tails[i]["holm_adjusted_tail_fraction"] = running
            inferential_status = "complete_conditional_benchmark" if tails else "no_selected_null_results"
        else: inferential_status = "P0_missing"
    write_tsv(output / "conditional_tail_fractions.tsv", tails, list(tails[0]) if tails else ["null_family", "metric", "conditional_tail_fraction"])
    write_json(output / "enrichment_status.json", dict(status="incomplete" if failed else "complete", backend=ecfg["backend"],
                annotation_version=ecfg["annotation_version"], inference=inferential_status, analyzed_queries=len(metric_rows),
                interpretation="Conditional randomization benchmark; significant terms are not independent discoveries"))
    require(not failed, "Missing/mismatched enrichment caches; see query_status.json. No null claims produced.")


if __name__ == "__main__":
    p = argparse.ArgumentParser(description=__doc__)
    for name in ("matched", "settings", "output"): p.add_argument("--" + name, type=Path, required=True)
    p.add_argument("--snapshot", type=Path); p.add_argument("--cache", type=Path)
    p.add_argument("--discovery-only", action="store_true")
    a = p.parse_args(); run(a.matched, a.settings, a.output, a.snapshot, a.cache, a.discovery_only)
