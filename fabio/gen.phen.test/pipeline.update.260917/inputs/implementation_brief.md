# CAAS strategy validation: operational Nextflow implementation brief

Date: 17 September 2026  
Scope: `score.approach/05.case.study.caas.rer.validation/CAAS`  
Companion document: `CAAS_strategy_validation_plan_20260917.odt` / `.pdf`

## Instructions to the implementation chat

Build a reproducible Nextflow DSL2 workflow to compare PSS-informed species-selection strategies with phenotype-only strategies and conditional random controls. This is an implementation brief, not a request to launch production analyses immediately.

Reuse the existing pooled CAAStools implementation and downstream filtering where possible. Keep all new inputs, outputs and manifests separate from historical production results. Deliver a working preparation stage, local smoke tests, cluster-ready launch wrappers, result aggregation and an explicit launch plan. Do not submit the full benchmark until Fabio has reviewed the generated strategy definitions, manifests and pilot cost.

The existing canonical PSS result is **P0, Cercopithecidae PSS-driven groups**, not the historical global best-pair-per-genus result. The word `random` in its historical directory name refers to sampling cycles from fixed PSS-derived pools; it does not mean that it is a species-selection null.

Do not implement RERconverge here. Do not expand this task into the deferred positional-significance audit or influential-gene/species sensitivity analyses. Preserve the current CAAS filters consistently across strategies; their use is not a new claim that the deferred audits are complete.

## 1. Inspect and reuse these project components first

Paths below are relative to the CAAS directory unless otherwise stated.

| Component | Existing location | Required use |
|---|---|---|
| Discovery workflow | `pipeline/main.nf` | Reuse the `ct pooled-discovery` invocation and event output contract. |
| Execution configuration | `pipeline/nextflow.config`, `pipeline/conf/cluster.config` | Inspect real paths, resources and environment; do not assume their defaults are portable. |
| Preparation scripts | `pipeline/scripts/create_pss_benchmark_configs.py`, `prepare_pooled_hypotheses.py`, `build_resampling.py` | Reuse deterministic cycle generation and validation where applicable. |
| Existing approach manifest | `pipeline/inputs/benchmark.configs.tsv` | Use as a schema reference, not as the default launch selection. |
| Source tables | `pipeline/inputs/config.creation/` | Authoritative phenotype values, species IDs and historical pool definitions. |
| Consolidated results | `results/consolidated.results/` | Read-only historical discovery results. |
| Downstream workflow | `results/meta-analysis/main.nf` | Reuse filtering logic, but do not run it with its default publishing destinations. |
| Downstream scripts | `results/meta-analysis/scripts/count_caas_support.py`, `run_gprofiler_enrichment.py` | Inspect for four-strategy hardcoding before adapting to replicate-aware inputs. |
| Current analysis settings | `results/meta-analysis/nextflow.config` | Freeze settings and annotation metadata in the new validation run. |

Some README command locations are historical. The preparation scripts currently live under `pipeline/scripts/`; verify paths against files, not prose alone. The current discovery workflow resolves manifest config paths against its own `projectDir`. If a new manifest lives elsewhere, explicitly support absolute paths or a declared input base directory; do not silently resolve them against the wrong directory.

Suggested new implementation directory: `CAAS/validation.workflow/`. This is a proposed layout, not an existing workflow. Keep this planning note in `CAAS/validation.planning/`.

## 2. Common experimental contract

- Use the existing adjusted/relative brain-mass values. Do not refit the phenotype separately for each strategy or random replicate.
- Use the same alignment collection, CAAStools version and discovery parameters across the benchmark.
- Each cycle has four FG and four BG species. Within an hypothesis, species have fixed FG/BG membership; the two complete pools are disjoint.
- For unrestricted pool designs, select 100 unique cycles without replacement from the Cartesian product of four-species combinations on each side.
- For linked-pair designs, select four complete pairs. Do not independently sample the two endpoints.
- Distinguish **strategy**, **replicate/hypothesis**, **cycle**, **gene/alignment** and **event** in both filenames and metadata.
- A random replicate comprises a complete hypothesis with its own pooled analysis. Its 100 overlapping cycles are not 100 independent null replicates.
- An output absent because a job failed is not a zero-result alignment. Aggregation must wait for all expected jobs or explicitly report incompleteness and block statistical summaries.
- Preserve underscore-separated species identifiers exactly as in the production inputs. Use phenotype-descending order, with species-ID lexical order as a deterministic tie-breaker.

Existing discovery settings to confirm and freeze: `phylip-relaxed`, patterns `1,2,3`, positional prefilter `0.05`, maximum gaps per position `0.5`, minimum observed non-gap species per cycle of three FG and three BG. Other gap/missingness options currently use `NO`; record the complete resolved command, not only these highlighted settings.

Retain both `*.pooled.caas.tsv` and `*.pooled.caas.events.tsv`. The event table is the primary downstream source; legacy cycle counts remain traceability information.

## 3. Strategy registry and launch stages

These IDs are validation IDs, not the numeric IDs of historical approaches.

| ID | Definition | FG / BG | Cycles per hypothesis | Default action |
|---|---|---:|---:|---|
| P0 | Current PSS-driven Cercopithecidae pools | 9 / 10 | 100 | Import existing results if compatible. |
| N0 | Primate-wide absolute extremes | 13 / 13 | 100 | Import existing reference results. |
| N1 | Family maximum/minimum linked pairs | 13 / 13 | 100 | Import existing reference results. |
| N2 | Cercopithecidae upper/lower 10% tails | 6 / 6 | 100 | Import existing reference results. |
| N3 | Cercopithecidae extremes, matching P0 pool sizes | 9 / 10 | 100 | First new deterministic control. |
| N4 | Phenotype-only reassignment of the same 19 P0 species | 9 / 10 | 100 | First new deterministic control. |
| N5 | Within-genus phenotype extremes, matching P0 genus quotas | 9 / 10 | 100 | First new deterministic control. |
| R0 | Within-genus FG/BG random reassignment of the same 19 species | 9 / 10 | 100 | Pilot: 19 hypotheses; proposed total: 99. |
| R1 | Random species and labels within the same five genera | 9 / 10 | 100 | Optional pilot of 19; optional total of 99. |
| P1 | Best top-PSS pair in each of six focal genera | 6 / 6 | 15 | Optional paired comparison with N6. |
| N6 | Maximum/minimum phenotype pair in those same six genera | 6 / 6 | 15 | Optional paired comparison with P1. |
| P2 | Global PSS-ranked endpoint-disjoint matching | 13 / 13 | 100 | Historical auxiliary design; inspect before considering rerun. |

Implement separate launch selections: `prepare`, `smoke`, `deterministic`, `r0-pilot`, `r0-complete`, `r1-pilot`, `r1-complete`, `paired`, and `summarize`. These are proposed interface names. Production selection must be explicit: never fall back to launching every historical approach.

## 4. Exact pool definitions

### P0: fixed reference

Read `pipeline/inputs/config.creation/09_cercopithecidae_pss_random_pools.tsv`. Validate the complete source pools against these IDs:

```text
FG (9): Macaca_nigra, Macaca_nemestrina, Semnopithecus_entellus, Macaca_mulatta, Papio_anubis, Colobus_polykomos, Trachypithecus_francoisi, Trachypithecus_pileatus, Trachypithecus_cristatus
BG (10): Macaca_leonina, Macaca_maura, Papio_papio, Macaca_tonkeana, Macaca_silenus, Macaca_fuscata, Semnopithecus_priam, Trachypithecus_geei, Colobus_guereza, Trachypithecus_auratus
```

For an imported reference, preserve the actual historical 100-cycle config and its checksum. Do not regenerate a different P0 config and attach historical results to it. If provenance or parameter compatibility cannot be established, list the mismatch and offer a fresh isolated reference rerun for approval.

### N3: phenotype extremes across the 55 Cercopithecidae

Sort `10_cercopithecidae_relative_brain_mass.tsv` by phenotype. Assign its nine largest values to FG and ten smallest to BG. No genus restriction on cycle sampling.

```text
FG (9): Macaca_nigra, Macaca_nemestrina, Macaca_assamensis, Cercocebus_torquatus, Theropithecus_gelada, Mandrillus_sphinx, Papio_ursinus, Semnopithecus_entellus, Macaca_mulatta
BG (10): Trachypithecus_auratus, Colobus_angolensis, Pygathrix_nigripes, Colobus_guereza, Trachypithecus_geei, Trachypithecus_phayrei, Trachypithecus_cristatus, Allenopithecus_nigroviridis, Trachypithecus_obscurus, Nasalis_larvatus
```

This matches pool size and cycle number, not taxonomic composition.

### N4: change the groups, not the species

Take the union of P0 FG and BG: exactly 19 species. Sort these species by the same phenotype. Assign the nine largest values to FG and the ten smallest to BG. No species is added or removed. No PSS value is used for this assignment.

```text
FG (9): Macaca_nigra, Macaca_nemestrina, Semnopithecus_entellus, Macaca_mulatta, Papio_anubis, Colobus_polykomos, Macaca_leonina, Macaca_maura, Papio_papio
BG (10): Trachypithecus_auratus, Colobus_guereza, Trachypithecus_geei, Trachypithecus_cristatus, Semnopithecus_priam, Macaca_fuscata, Trachypithecus_pileatus, Macaca_silenus, Macaca_tonkeana, Trachypithecus_francoisi
```

N4 tests whether the PSS-derived grouping adds information after the species set is fixed. It does not test whether PSS selected better species in the first place. Assert `union(N4) == union(P0)` and complete phenotype separation between N4 FG and BG.

The existing README mentions an aborted median-like split of these species. Do not reuse its partial outputs or obsolete approach label. N4 must have a new unambiguous validation ID and checked source/configuration hashes.

### N5: phenotype extremes with P0 genus quotas

Use all eligible species in each of P0's five focal genera from the 55-species phenotype table. For each genus, take the highest values for its FG quota and the lowest for its BG quota.

| Genus | Eligible species | FG quota | BG quota |
|---|---:|---:|---:|
| Colobus | 3 | 1 | 1 |
| Macaca | 13 | 3 | 5 |
| Papio | 5 | 1 | 1 |
| Semnopithecus | 2 | 1 | 1 |
| Trachypithecus | 8 | 3 | 2 |
| Total | 31 | 9 | 10 |

```text
FG (9): Colobus_polykomos, Macaca_nigra, Macaca_nemestrina, Macaca_assamensis, Papio_ursinus, Semnopithecus_entellus, Trachypithecus_francoisi, Trachypithecus_pileatus, Trachypithecus_germaini
BG (10): Colobus_angolensis, Macaca_fuscata, Macaca_silenus, Macaca_tonkeana, Macaca_cyclopis, Macaca_radiata, Papio_papio, Semnopithecus_priam, Trachypithecus_auratus, Trachypithecus_geei
```

Genus quotas apply to complete pools, not individual cycles. Do not add an unrequested genus restriction to 4-vs-4 sampling. Matching genera does not guarantee matched patristic distances or alignment availability; record those distributions where inputs permit and do not describe them as automatically controlled.

### R0: conditional label randomization within P0 species

Keep the same 19 species and the same per-genus FG/BG quotas. Randomly assign species to sides within each genus. Do not inspect phenotype or PSS during randomization and do not require FG to have a greater phenotype.

There are `choose(8,3) × choose(5,3) × 2 × 2 × 2 = 4,480` assignments. Exclude P0's assignment, select unique assignments without replacement, and freeze them before discovery. Each assignment receives 100 unique cycles. Use stable replicate IDs so extending from 19 to 99 preserves the first 19 byte-for-byte.

This is a conditional random-group benchmark, not automatically a phylogenetically calibrated null of no genome–phenome association. Its interpretation is restricted to alternative groupings of these selected species under the stated randomization rule.

### R1: conditional species-selection randomization

Use the 31 species in the five genera listed above. Within each genus, randomly assign distinct species to the required FG and BG quotas, leaving any remaining species unselected. Do not use phenotype or PSS. Freeze the eligible-species table and every selected assignment. Reject duplicate complete assignments across replicates, and exclude P0 itself from the random reference set.

R1 benchmarks alternatives within PSS-selected genera; it is not an independent test of why these five genera were selected from all primates. A pilot of 19 replicates is descriptive; extend to 99 only after cost review. Do not substitute an unrestricted 55-species randomization without discussing the changed biological question.

### P1 versus N6: linked pairs in the same six genera

Read `04_best_top1pct_pair_per_genus.tsv` and restrict it to Cercopithecidae for P1. Construct N6 from phenotype extrema within those exact same genera. Check against this registry:

| Genus | P1 FG | P1 BG | N6 FG | N6 BG |
|---|---|---|---|---|
| Trachypithecus | Trachypithecus_cristatus | Trachypithecus_auratus | Trachypithecus_francoisi | Trachypithecus_auratus |
| Macaca | Macaca_nigra | Macaca_tonkeana | Macaca_nigra | Macaca_fuscata |
| Semnopithecus | Semnopithecus_entellus | Semnopithecus_priam | Semnopithecus_entellus | Semnopithecus_priam |
| Colobus | Colobus_polykomos | Colobus_guereza | Colobus_polykomos | Colobus_angolensis |
| Papio | Papio_anubis | Papio_papio | Papio_ursinus | Papio_papio |
| Cercopithecus | Cercopithecus_cephus | Cercopithecus_petaurista | Cercopithecus_cephus | Cercopithecus_diana |

Each strategy has six disjoint pairs. Enumerate all `choose(6,4) = 15` cycles, preserving endpoints and using corresponding genus subsets in both strategies. Never duplicate cycles to reach 100. Compare P1 with N6 internally, not their pooled discovery volume directly with 100-cycle strategies. This comparison remains conditional on these six PSS-selected genera.

### Historical reference mapping

| Validation ID | Historical source approach | Source table |
|---|---|---|
| P0 | `08_pss_cercopithecidae_random_pools` | `09_cercopithecidae_pss_random_pools.tsv` |
| N0 | `02_absolute_trait_tails` | `02_trait_distribution_tails.tsv` |
| N1 | `01_family_extrema` | `01_family_trait_extrema.tsv` |
| N2 | `09_cercopithecidae_absolute_trait_tails` | `11_cercopithecidae_absolute_trait_tails.tsv` |
| P2 | `07_pss_ranked_endpoint_disjoint_13x13` | `08_pss_ranked_endpoint_disjoint_13x13.tsv` |

N0 has a within-side, at-most-one-species-per-genus cycle restriction. N1 samples four complete family pairs. N2 independently samples four species per side without a genus restriction. Preserve historical rules, not a universal resampling rule applied to all arms.

## 5. Input preparation and reproducibility

Implement one independently testable preparation script that outputs the registry, complete pools, pair tables, replicate assignments and pooled configs. Preparation should be possible without submitting any Nextflow discovery task.

Use a stable deterministic selection scheme, preferably the existing SHA-256 ranking convention for cycles. For random assignments, document the precise uniform sampling-without-replacement algorithm and pseudorandom generator. A seed alone is insufficient: save explicit assignment/cycle manifests, generator version and all hashes. Do not change frozen assignments when moving between machines or extending a pilot.

Required headered TSV files:

- `strategy_registry.tsv`: `strategy_id`, historical source ID if any, sampling mode, source table, pool sizes, expected cycles, reference/import status.
- `pool_membership.tsv`: `strategy_id`, `replicate_id`, `side`, `species`, `genus`, `trait_value`, pair/family ID where applicable.
- `eligible_species.tsv`: complete eligible universe for each randomization, phenotype values for documentation only, and taxonomy source.
- `replicate_manifest.tsv`: unique `hypothesis_id`, `strategy_id`, `replicate_id`, null family, assignment checksum, seed/generator metadata and pool/config paths.
- `cycles.tsv`: hypothesis ID, cycle ID, FG IDs, BG IDs, linked-unit IDs where applicable.
- `execution_manifest.tsv`: hypothesis ID, config paths, possible/selected cycle counts, source/pool/config checksums and source provenance.
- `alignment_manifest.tsv`: unique gene ID, alignment path, checksum and format; reject duplicate IDs or ambiguous gene-symbol mapping.

Existing CAAStools adapters must retain their expected formats: pooled configs are headerless `cycle_id<TAB>comma-separated-FG<TAB>comma-separated-BG`; complete pool configs are headerless species/side records using the existing numeric encoding (inspect and test the current files). Human-readable TSV manifests must not be accidentally passed to `ct` as these adapter files.

## 6. Nextflow architecture and isolation

Suggested process/module boundaries:

1. `PREPARE_VALIDATION`: generate frozen configs and manifests, or validate already frozen ones.
2. `VALIDATE_INPUTS`: check species, quotas, cycle uniqueness, reference compatibility and alignment inventory; emit a machine-readable report. Fail early on structural errors.
3. `CAAS_POOLED`: adapt the existing process; one hypothesis × one alignment task. Carry strategy and replicate IDs in channel tuples, task tags and output paths.
4. `ASSEMBLE_HYPOTHESIS`: collect all expected gene outputs by hypothesis and verify completeness.
5. `FILTER_AND_QUERY`: reuse current primary-event, positional, dense-cluster and support rules; derive exact query and background inventories.
6. `COMMON_BACKGROUND`: assemble a common documented testable gene universe and intersect queries with it.
7. `ENRICHMENT`: analyze observed and null queries with identical frozen annotation inputs.
8. `GENE_LIST_CONTROLS`: generate G0/G1 controls using the same enrichment implementation.
9. `COMPARE_STRATEGIES`: produce deterministic comparisons, null distributions and summaries with uncertainty.

Use explicit channels/declared inputs for reference files and annotation snapshots; do not discover newly published files by globbing a shared result directory while tasks are still running. Keep preparation checksums separate from runtime provenance. Changing a parameter must invalidate appropriate caches rather than reuse incompatible summaries.

Proposed isolated layout:

```text
CAAS/validation.workflow/
  main.nf
  nextflow.config
  conf/                 # local/test and cluster profiles
  modules/
  scripts/
  tests/fixtures/
  inputs/frozen/<design_id>/
  results/<run_id>/<strategy_id>/<replicate_id>/
  summaries/<run_id>/
  work/                 # or explicitly configured cluster scratch
  README.md
```

Use unique run/design IDs. Never publish into `results/consolidated.results`, historical `results/tables` or `results/figures`. Record imported references by path and checksum without editing them. Resume must require the intended run/design and explicit selected manifest. Protect against duplicate gene basenames, output collisions and automatic reruns of all historical arms.

Start with the proven gene-level discovery layout. Use the pilot to decide whether thousands of scheduler submissions or repeated alignment staging require batching. Any batching must preserve hypothesis/gene provenance, output equivalence, failure visibility and resumability. Do not silently alter statistical pooling to improve throughput.

## 7. Downstream settings and functional controls

Reuse the current definitions consistently: primary events, nominal positional `p < 0.05`, one retained record per gene/position, and current dense-cluster pruning (density ≥ 0.7, span ≥ 3, at least three positions). Query genes require a surviving position with at least four FG and four BG support at that same position. Do not combine support maxima from different positions to create this criterion.

Reconstruct backgrounds from all successfully evaluated alignment results, including genes without positive events. Verify what constitutes an effectively testable gene in the current tool; write the operational rule and coverage counts. Do not equate failed/missing output with testability. Preserve both strategy-specific backgrounds and a common background for matched enrichment comparisons.

The list-size control is important: current P0 has 1,421 query genes, whereas N2 has 11,834. For a query contained in its background, maximum fold enrichment is bounded by the inverse selected-background fraction. Therefore median significant fold enrichments alone are not a fair validation metric.

- **G0:** 1,000 random gene lists of the same size as P0 after intersection with the common background, sampled without replacement within each list from that background.
- **G1:** 1,000 same-size subsamples of the N2 query after its common-background intersection; optionally repeat for N3. If the eligible control query is too small, report the incompatibility rather than sampling with replacement.

Keep all explicit gene-list manifests. G0/G1 require no new CAAS discovery and do not replace species/group random controls: they do not reproduce gene-specific CAAS selection probabilities.

Freeze organism, annotation sources, correction and annotation version. Existing enrichment uses human annotations, g:SCS at 0.05, custom background and the sources recorded in the current meta-analysis configuration. Do not query a live service separately for thousands of replicates and assume their results form a version-matched null. Prefer a reproducible annotation snapshot and validated local enrichment, or a documented version-pinned service if available. A different enrichment implementation requires an observed/reference rerun under that implementation and a consistency check; do not mix its p-values with historical g:Profiler p-values.

If API access is needed, document the transmitted gene symbols, obtain the run decision, cache responses and metadata, and use rate limits. Live enrichment should be disabled in a default smoke test. Lack of a frozen annotation implementation must not block manifest generation or discovery pilot preparation, but it does block claims from unmatched enrichment nulls.

## 8. Prespecified comparisons and statistics

Implement metrics as explicit configuration, not values chosen after results are seen. Before full null production, Fabio must approve the primary statistic and any term universe or term-grouping algorithm.

Recommended outputs include query/background fraction, enrichment effects on a declared term universe, nonredundant functional breadth, supporting genes, and profile summaries. The current 14 P0 terms are not 14 independent discoveries. If they are frozen as targets now, label that comparison as conditional on previously discovered terms, not independent confirmation.

For a global discovery statistic, run the identical term-selection and summarization algorithm on observed and random results; do not compare each replicate's cherry-picked significant-term median. Store nonsignificant/missing-query statuses separately from failed requests. Do not assign an effect of zero solely because a term failed a significance threshold.

For a prespecified upper-tail metric, report the Monte Carlo benchmark tail fraction as `(1 + number of random values >= observed) / (B + 1)`, along with effect size, null quantiles and Monte Carlo uncertainty. With 99 replicates the minimum is 0.01; with 19 it is 0.05. A small pilot is not a broad term-by-term significance screen. If several metrics or terms are tested inferentially, declare the family and correction before production; do not silently add an exploratory correction or optimize the primary metric afterward.

These tail fractions are conditional on the randomization design. Do not present them as a universal calibrated test of genome–phenome association or proof that PSS causes better biological discovery. The fixed observed grouping and random groupings have different selection mechanisms; describe exactly what is being benchmarked.

| Comparison | Question addressed | Limit |
|---|---|---|
| P0 vs N3 | Does PSS-informed selection differ from family-wide extrema at equal pool size? | Taxonomic composition differs. |
| P0 vs N4 | Does grouping matter for the same species set? | Species selection remains PSS-derived. |
| P0 vs N5 | Does PSS differ from phenotype extrema within the same genus quotas? | Genera/quotas are selected conditionally; distances and coverage may differ. |
| P0 vs R0 | How unusual is this grouping among randomized groupings of these species? | Conditional grouping benchmark. |
| P0 vs R1 | How unusual is it among alternative pools/groupings in these genera? | Does not validate selection of genera. |
| P1 vs N6 | Does PSS pair ranking add information beyond maximum local phenotype differences? | Conditional six-genus comparison; only 15 linked cycles per arm. |
| P0 vs G0/G1 | Can query size alone explain the functional concentration? | Gene-list control, not species-selection validation. |

## 9. Pilot, budget and approval gates

Run local fixtures first. After review, start with N3/N4/N5 and a 19-hypothesis R0 pilot. A reduced alignment panel is suitable for debugging and resource estimation, not for reporting the final biological benchmark. A biological R0 pilot must use the same alignment inventory as P0.

Budget in complete hypotheses and within-hypothesis cycles:

- N3/N4/N5: three hypotheses, 300 cycles.
- R0 pilot: 19 hypotheses, 1,900 cycles.
- N3/N4/N5 plus R0 extended to 99 total: 102 hypotheses, 10,200 cycles; this includes the pilot, not 99 additional replicates.
- Optional R1 pilot: another 19 hypotheses, 1,900 cycles.
- Optional P1/N6: two hypotheses, 30 cycles.
- Combined package above: 123 hypotheses, 12,130 cycles, excluding historical reference reruns and P2.

With approximately 16,130 alignments, the 102-hypothesis package implies about 1.65 million hypothesis/alignment combinations under the current task layout. This is not a wall-time estimate. Measure CPU, memory, runtime by alignment size, failure rates, scheduler overhead, I/O and disk output in the pilot; propose resource and batching changes before the full null launch.

Approval gates:

1. Pool/strategy review: exact species, sampling rules and source checksums approved.
2. Technical review: fixtures, smoke outputs and imported-reference compatibility pass.
3. Statistical review: primary metric, annotation strategy, randomization interpretation and null size approved.
4. Production review: pilot resource estimate and explicit selected manifest approved.

If only one random family is affordable, prioritize R0 as the most tightly species-matched design, but restrict conclusions accordingly. R1 is the stronger additional benchmark for alternative species selection within these genera; do not claim that R0 alone answers that question.

## 10. Acceptance tests and required deliverables

Automated tests must cover:

- P0 exact 9/10 pools; N3/N4/N5 exact expected pools; source changes trigger a visible discrepancy, not silent regeneration.
- N4 union equality with P0 and complete phenotype separation.
- N5/R0/R1 genus quotas; no cross-side species; correct eligible universes.
- R0 assignment count of 4,480, unique selected nulls, P0 exclusion, stable pilot extension.
- P1/N6 endpoint disjointness, six pairs and exactly 15 corresponding linked cycles.
- 26,460 possible unrestricted 9/10 cycles, exactly 100 unique selected cycles per hypothesis; no repeated cycles disguised with new IDs.
- Stable assignments/configs across repeated preparation; checksums actually used in execution metadata.
- Synthetic alignments producing positive and zero-event results, plus deliberately malformed inputs that fail visibly.
- Downstream same-position support semantics, dense-cluster filtering and preservation of zero-result genes in backgrounds.
- Missing output blocks aggregation; a resumed partial run reuses compatible completed tasks without duplicate publishing.
- Duplicate gene IDs, conflicting membership, insufficient pair counts and absent species identifiers produce informative failures. Alignment-specific species absence is reported and handled by the current missingness rules, not mistaken for an unknown pool identifier.
- G0/G1 exact query sizes and reproducibility; identical enrichment settings on observed and controls.

Deliver code and tests, a README with preparation/local/cluster commands, explicit pilot/full manifests, frozen pool/cycle TSVs, a smoke-test report, and a cost-estimation template. Final validation outputs must be TSV/JSON and scientific PDF/PNG figures, never Excel/ODS spreadsheets.

Required summary files should include execution completeness, exact query/background inventories, common-background intersections, all analyzed enrichment terms and effects, random-replicate metrics, gene-list-control metrics, conditional tail fractions with uncertainty, and a strategy-level comparison table. Keep raw events and annotation metadata auditable.

End the implementation handoff by stating what was implemented, what was actually tested, what remains unrun, which references were imported and any choices requiring Fabio's review. Do not equate successful software execution with biological validation.
