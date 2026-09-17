# Software verification report

Date: 17 September 2026. Scope: isolated `pipeline.update.260917` only.

**These are synthetic software checks, not biological validation results.**
No production job, SLURM submission or live enrichment request was made.

## Preparation and regression tests

The frozen `brain-260917-v1` design contains 208 hypotheses and 20,630 cycles,
including optional strategies. R0/R1 each have 99 unique assignments; the saved
19-replicate pilot configurations are byte-identical to the first 19 of the
99-replicate design. N3/N4/N5 and P0 match the brief's exact species; genus quotas,
N4 species-set equality and complete phenotype separation, P0 exclusion, linked
P1/N6 endpoints, 15 corresponding pair cycles and 26,460 possible unrestricted
cycles are checked. Historical reference pool/cycle configs are copied exactly.

- **23 validation tests passed**, including malformed dimensions, duplicate
  genes/cycles, conflicting pools, missing phenotype identifiers, insufficient
  pair counts, source hash changes, same-position support, dense pruning,
  zero-result backgrounds, missing-output rejection, stable preparation,
  exact-sized G0/G1 lists and incompatible control-query handling.
- **18 bundled CAAStools regression tests passed**, covering pooled discovery,
  species-deduplicated event reconstruction, fixed group-size significance,
  minimum-observed coverage and invalid requests.
- Production without approval and resume without an existing run lock are
  rejected before submission. Missing frozen annotation cache is a failure,
  not an enrichment result of zero.

Unedited logs and structured status: `review/tests/validation.log`,
`bundled-caastools.log`, `test_report.json`. Regenerate with
`python3 tests/save_test_report.py`.

## Real Nextflow execution on fixtures

Two 129-species, 20-position synthetic PHYLIP alignments were generated. The
positive alignment has contrast columns at 0, 1, 2 and 15; the other is constant.
Nine hypotheses were exercised: P0, N2, N3, N4, N5, P1, N6, R0_001, R1_001.
This gives 18 real `ct pooled-discovery` tasks.

For synthetic P0, the event table contains four nominally significant positions
with 9 FG and 10 BG supporting species. Dense-cluster pruning removes 0–2 and
retains position 15. The query contains one gene; the coverage-testable background
contains both genes, including the constant zero-result alignment.

### Deliberate partial failure and resume

Run ID: `smoke-final-partial`. A smoke-only interpreter wrapper deliberately
failed R1_001/ZERO; Nextflow terminated and did not publish a completed cohort.
One additional running task was aborted. Repeating the exact run with `--resume`
reused **16 cached CAAS tasks** and completed **2 remaining tasks**. All nine
summaries, matched backgrounds, offline toy enrichment and PDF/PNG figures then
completed. The final raw inventory contains exactly 18 distinct gene/hypothesis
receipts, without duplicate publishing.

Proof: `review/tests/partial_resume.json`, named run launch logs and trace files
under `launch-plans/smoke-final-partial/` and `results/smoke-final-partial/`.

### Default annotation-disabled execution

Run ID: `smoke-default-final`. All 18 discovery tasks, filtering, common
background, gene-list-control manifests and figures completed. The enrichment
status explicitly says `blocked_annotation_not_frozen`, with zero analyzed
annotation queries and no requests sent. See `review/tests/default_smoke.json`.

### Explicit multi-run summary interface

Run ID: `summary-smoke-final`. The `summarize` stage consumed two explicitly
checksummed completed summaries from `review/smoke-cohort.tsv`, staged their
directories, and produced matched tables, functional status and figures without
any CAAS task. Different enrichment/statistical settings do not force a discovery
rerun when the discovery/filter contract is unchanged.

### Final provenance publication check

Run ID: `smoke-provenance-final`. The default frozen N3 smoke selection completed
two CAAS tasks and all downstream stages. Verified publication of
`results/smoke-provenance-final/provenance/runtime.lock.json` and
`input_validation.json`. This checks the final metadata publishing rule on a
filesystem where native symlinks are unavailable; inputs are staged by copying.

### Figures and resources

Figures were generated as PDF and 300-dpi PNG and visually inspected. Each
synthetic plot is labeled as a software test. The R script runs independently of
Nextflow. The toy hypergeometric/Bonferroni implementation was tested against an
exact small combinatorial example; it is not claimed equivalent to g:SCS.

`review/smoke-cost/` records successful task resource measurements and deliberate
failure/abort trace entries. It intentionally does **not** extrapolate synthetic
timing to a biological runtime, storage budget or cluster price. macOS Nextflow
does not provide all scheduler resource metrics; receipt-level measurements and
real SLURM accounting must be reviewed on the biological pilot.

## Historical reference compatibility

No historical reference was imported. Read-only audit found the following raw
file counts (events and legacy agree within each directory):

| Reference | Event files | Legacy files |
|---|---:|---:|
| P0 | 16,132 | 16,132 |
| N0 | 16,130 | 16,130 |
| N1 | 16,129 | 16,129 |
| N2 | 16,127 | 16,127 |

These counts are **not** a certification of unique tested genes, completeness,
or comparable backgrounds. No runtime-provenance JSON/log/trace candidates were
found within these consolidated approach directories. Original alignment
checksums, successful-job inventories, full commands and historical tool version
remain to be established. The optional P2 directory is absent from consolidated
results. Exact original config checksums and all inspected file paths are in
`review/reference-audit/`. Historical files were not modified.

## Remaining review and unrun components

1. Fabio's review of exact pools, selected manifests and interpretation.
2. Compatible historical reference evidence, or approval of isolated reruns.
3. Primary statistic, term universe/grouping, null sizes/comparison family and
   annotation snapshot/backend. No production primary metric was chosen here.
4. Actual cluster environment, approved alignment inventory and pilot resources;
   then an explicit production approval tied to those hashes.

No biological deterministic control, species-null pilot/full set, paired
comparison, production annotation analysis or real cluster costing has run.
SLURM configuration was parsed locally but not submitted. The g:Profiler-cache
backend has not been exercised against a complete real frozen service cache;
production reference-import compatibility cannot be tested without its evidence.
