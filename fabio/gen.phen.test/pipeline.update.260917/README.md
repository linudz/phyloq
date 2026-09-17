# CAAS strategy validation — isolated Nextflow DSL2 workflow

Implementation date: 17 September 2026. This directory implements the supplied
`validation.planning/CAAS_validation_nextflow_implementation_brief.md` (a copy is
retained in `inputs/implementation_brief.md`). It does **not** implement RERconverge,
the deferred significance audit, or influential-gene/species analyses.

**No production discovery or live enrichment has been launched.** Local synthetic
tests exercise the actual bundled pooled CAAStools implementation. Passing these
tests is software verification, not evidence that PSS improves biological inference.

## 0. Current cluster entry point: one launch, no mandatory pilot

Following Fabio's updated instruction, use **`launch_all_caas.sh`** to run the
CAAS validation modes together. No pilot, approval JSON, smoke test or
statistical settings decision is required to launch CAAS through this entry
point. Frozen species assignments, the modified tool version, alignment
checksums, output isolation and strict resume checks remain enforced.

From this directory on the cluster:

```bash
# Once, if not already created:
bash create_conda_environment.sh

# Reuse the previous pipeline's ../caas/inputs/alignments/*.phy:
bash launch_all_caas.sh --run-id caas-all-01
```

The default matches the old `caas/conf/cluster.config`, resolved relative to
this new sibling directory. The original research-project copy uses
`../pipeline/inputs/alignments/*.phy`; a standalone copy without either legacy
pipeline uses its own `inputs/alignments/*.phy`. Explicit alignment options
override these defaults. See `INPUT_PATHS.md` for alignment/tree provenance.

This selects **N3/N4/N5 + all 99 R0 + all 99 R1 + P1/N6**: **203 hypotheses,
20,130 cycles across hypotheses**, with one pooled CAAS task per hypothesis and
alignment. With 16,133 alignments this means 3,274,999 CAAS tasks, not 203 jobs;
the current executor limits the number of outstanding tasks to 100. All selected
modes belong to one Nextflow run, without a manual stage barrier. This does not
mean all tasks start simultaneously or that the scheduler runs them in a fixed
interleaved order. The streamed job index avoids keeping millions of job rows in
memory, but per-gene work directories and raw results still need substantial
shared storage. No alignment batching or automatic work cleanup is introduced.

Historical references **P0/N0/N1/N2** and auxiliary **P2** are excluded by default.
Add `--include-references` and/or `--include-p2` to rerun them in this new output
directory, never over the historical results. Both flags select all 208 frozen
hypotheses. `--modes deterministic,r0` selects the 102-hypothesis main package;
`--modes deterministic,r0,r1,paired` is the default `all` selection.

The launch saves its exact selection and inventory under
`inputs/launches/<run-id>/`, performs discovery, assembly, the existing position
filter and descriptive summaries/figures. Enrichment, gene-list random controls
and inferential comparisons are **not requested in this CAAS-only launch**;
they remain available separately with `summarize`. No scientific approval is
invented or set to true to make discovery run.

For an interrupted run, repeat the **same command with `--resume`**. An optional
`--plan-only` prepares/prints the plan without submitting anything. See
`RUN_ON_CLUSTER.md` for alternate input patterns, cluster resources and output
paths. The staged interface in section 6 is retained for compatibility, not as
a prerequisite for this direct entry point.

## 1. Design and scientific interpretation

P0 is the current **Cercopithecidae PSS-driven groups**, historically strategy 08;
it is not the discarded global best-pair-per-genus strategy. Its historical name
contains `random` because cycles were sampled from fixed pools, not because P0 is
a species-selection null. Phenotypes are the existing adjusted brain-mass values;
no allometric model is refitted here.

| ID | Definition | Complete FG/BG pools | Cycles | Initial action |
|---|---|---:|---:|---|
| P0 | PSS-driven Cercopithecidae groups | 9/10 | 100 | Provenance review; no automatic import |
| N0 | Primate-wide extremes | 13/13 | 100 | Preserve historical within-side genus rule |
| N1 | Family parallel pairs | 13/13 | 100 | Preserve four linked family pairs per cycle |
| N2 | Cercopithecidae upper/lower 10% tails | 6/6 | 100 | Preserve historical unrestricted cycles |
| N3 | Cercopithecidae phenotype extremes, size matched to P0 | 9/10 | 100 | First deterministic control |
| N4 | Phenotype reassignment of exactly P0's 19 species | 9/10 | 100 | First deterministic control |
| N5 | Phenotype extremes matching P0's genus quotas | 9/10 | 100 | First deterministic control |
| R0 | Random labels within genera, same 19 species | 9/10 | 100 per hypothesis | All 99 in the direct launch |
| R1 | Random species and labels within the same five genera | 9/10 | 100 per hypothesis | All 99 in the direct launch |
| P1 | Best top-1% PSS pair within six focal genera | 6/6 | **15** | Optional, compare internally with N6 |
| N6 | Phenotype extrema within those same six genera | 6/6 | **15** | Optional, compare internally with P1 |
| P2 | Historical global endpoint-disjoint PSS matching | 13/13 | 100 | Auxiliary, not selected by default |

N4 keeps species fixed but changes groups; it does not validate species selection.
R0 is conditional on PSS-selected species, R1 on PSS-selected genera/quotas. Neither
is automatically a phylogenetically calibrated null of no genome–phenome
association. Matching genus quotas does not guarantee matching patristic distances
or sequence coverage. Actual frozen-tree distances and trait summaries are in
`review/pool-diagnostics/`.

The frozen design contains **208 hypotheses / 20,630 cycles**, including all
optional arms and 99 each of R0/R1. This is an input inventory, **not a launch
selection**. The main proposed package is 102 hypotheses / 10,200 cycles; the
package adding a 19-replicate R1 pilot and P1/N6 has 123 / 12,130, excluding
reference reruns and P2.

### Reproducible sampling

Within a random genus, `random.Random` (CPython MT19937) uniformly samples an
ordered set without replacement, splits it into the required FG and BG quotas,
and leaves the remaining eligible species unused. Full assignments equal to P0
or an earlier replicate are rejected. Thus accepted assignments sample uniformly
without replacement from the admissible conditional assignment space. Inputs
are lexically sorted during randomization: phenotype/PSS never decide a draw.
The family-specific PRNG seed is the integer SHA-256 digest of `[seed, family]`.
R0 has 4,480 assignments including P0; 4,479 are eligible controls.

For each unrestricted hypothesis, rank all 26,460 possible 4-vs-4 cycles by
SHA-256 of `seed<TAB>hypothesis_id<TAB>sorted_FG<TAB>sorted_BG` and retain the first
100. The full pool is fixed and disjoint throughout the hypothesis. Genus quotas
apply to the pool, **not** to individual cycles. P1/N6 enumerate all 15 subsets of
four complete pairs, with corresponding genus-subset order in both arms. They
never receive duplicate cycles to reach 100. Display tables use descending
phenotype and lexical tie-breaking; frozen historical adapter files remain
byte-for-byte unchanged. Explicit manifests, Python version and hashes, rather
than the seed alone, define the portable analysis.

## 2. Layout and immutable inputs

```text
main.nf, nextflow.config, conf/   DSL2 and local/SLURM configuration
launch_all_caas.sh              one direct cluster launch for all new CAAS modes
run_validation.py              explicit, default-plan-only launcher
submit_validation_slurm.sh      optional SLURM driver, never auto-submitted
bin/caastools/                  copied pooled CAAStools, unmodified
scripts/                       preparation, auditing, discovery and summaries
tests/                         acceptance checks and synthetic fixtures
inputs/source.lock.json        hashes of authoritative original inputs
inputs/frozen/brain-260917-v1/  portable tables, source copies, configs and locks
review/                        proposed launch manifests, diagnostics and audit
results/<run>/<strategy>/<replicate>/raw/<gene>/
results/<run>/<strategy>/<replicate>/<hypothesis>/
summaries/<run>/{matched,functional,figures}/
work/                          Nextflow work/cache; retain for resume
launch-plans/<run>/             exact resolved commands, parameters and logs
environment.yml                independent Conda recipe, including R/Java/Nextflow
create_conda_environment.sh     create/check the dedicated runtime (no phyloq edits)
run_pipeline.sh                activate it and call the explicit launcher
```

`pool_membership.tsv`, `pairs.tsv`, `eligible_species.tsv`, `cycles.tsv`,
`strategy_registry.tsv`, `replicate_manifest.tsv` and `execution_manifest.tsv`
are headered human-readable/auditable tables. CAAStools receives only the
headerless adapters: `cycles.cfg` = cycle ID, comma-separated FG, comma-separated
BG; `pool.cfg` = species and numeric state (`1=FG`, `0=BG`).

All frozen source/config hashes are checked before execution. Preparation uses
the bundled frozen source copies by default; `--from-project` explicitly checks
the authoritative originals when working in the full research project. Changed
sources are detected, not silently regenerated. Already frozen
designs are verified rather than overwritten. A cluster only needs this complete
directory plus its explicit alignment collection: it does not need the original
project's absolute paths for discovery.

## 3. Requirements and local preparation

The bundle now includes `environment.yml` for a dedicated
**`caas-validation-260917`** environment, independent of the historical `phyloq`.
It specifies Python 3.11, Nextflow 24.04.2, Java 17, Biopython 1.84, DendroPy 4.5.2,
NumPy 1.26.4, SciPy 1.13.1, a setuptools release retaining `pkg_resources`, and
R 4.3. Figures use base R: no extra CRAN installation is needed. The initial
local testing stack recorded in `requirements-tested.txt` is a historical
record, **not** the new Conda installation recipe.

On the cluster login/setup node, from this directory:

```bash
# Load the cluster's Conda module first, if conda is not on PATH.
# Alternatively set CAAS_CONDA_SH=/path/to/etc/profile.d/conda.sh.
bash create_conda_environment.sh
bash run_pipeline.sh prepare

# OPTIONAL software smoke, not a prerequisite for launch_all_caas.sh:
# Regenerate synthetic paths on THIS machine before using this optional test:
conda run -n caas-validation-260917 python tests/make_fixtures.py
bash run_pipeline.sh smoke --run-id smoke-cluster-01 \
  --selection tests/fixtures/selection.tsv \
  --alignments tests/fixtures/alignment_manifest.tsv --execute
```

Run the local-executor smoke on a node where local computation is permitted.
Conda is discovered from PATH or `CAAS_CONDA_SH`; the known Correfoc Conda hook
is used only if it actually exists. No path into the old workflow is required.
The setup script creates the environment if absent, otherwise only checks it.
`--update` is an explicit opt-in to package changes; it never prunes or removes
environments, and refuses the base/shared `phyloq` target. Do not update an
environment while its jobs are running. `--dry-run` asks Conda to solve without
installing, and `--check` validates an existing installation.

For a custom shared-filesystem prefix:

```bash
export CAAS_CONDA_PREFIX=/ABSOLUTE/SHARED/PATH/caas-validation-260917
bash create_conda_environment.sh --prefix "$CAAS_CONDA_PREFIX"
bash run_pipeline.sh prepare
```

Keep `CAAS_CONDA_PREFIX` exported when submitting jobs. Alternatively use
`CAAS_CONDA_ENV` for a custom name. Both the environment prefix and Conda hook
must be readable on every allocated compute node. Environment creation belongs
on the setup node, **not** inside each gene task. `run_pipeline.sh` activates
the driver environment and Nextflow reactivates that same absolute prefix for
its worker processes. Hence Nextflow's *per-process environment creation*
remains disabled (`conda.enabled = false`); this does not disable our explicit
shared-environment activation.

Every executing launch checks Python imports, legacy `pkg_resources`, Java,
the pinned Nextflow version and actual headless R PDF/PNG rendering. Cluster
execution also checks SLURM commands. It saves the executable paths, package
versions and Conda package/build records in `results/<run>/environment.json`.
Changes to the recipe, installed Conda packages or runtime fingerprint block
incompatible resume. Setup reports are under `review/environments/`.

The low-level `python3 run_validation.py` interface remains available for an
already validated external environment. Prefer `bash run_pipeline.sh` on the
cluster; Docker/Singularity and automatic environment downloads per task are
not used. The setup script cannot activate an environment in its parent shell.
Before running the direct `python3` utility commands below, activate
`caas-validation-260917` in your shell, or prefix each command with
`conda run -n caas-validation-260917` (use `--prefix` for a custom prefix).

From this directory, with that environment active:

```bash
# Portable check, without discovery or Nextflow submission:
python3 run_validation.py prepare

# Verify/recreate from the portable frozen source snapshot (no old project needed):
python3 scripts/prepare_validation.py

# Optional, only in the original complete research project:
python3 scripts/prepare_validation.py --from-project

# Independent fixtures and acceptance checks:
python3 tests/make_fixtures.py
python3 -m unittest discover -s tests -v

# Proposed extension/package selections and descriptive diagnostics (portable):
python3 scripts/prepare_launch_review.py
python3 scripts/pool_diagnostics.py

# Optional read-only audit, requiring the historical result directories:
python3 scripts/audit_references.py
```

Do not run `--lock-sources` over an existing lock. A deliberately revised design
requires a separately reviewed source lock, design ID and output directory.

## 4. Explicit alignment inventory

The biological alignments remain external/ignored. On the machine that will run
the analysis, use a quoted glob covering **exactly the approved collection**:

```bash
python3 scripts/inventory_alignments.py \
  --alignments '/ABSOLUTE/CLUSTER/PATH/inputs/alignments/*.phy' \
  --output inputs/alignment_manifest.tsv
```

Adjust the extension to the real files; no extension is assumed by the workflow.
The inventory stores absolute paths, SHA-256, species and alignment dimensions.
It validates PHYLIP-relaxed contents and rejects duplicate gene IDs/basenames,
duplicate sequence IDs and header/dimension inconsistencies. The gene symbol is
the first filename token before a dot, exactly as in the bundled CAAStools:
two files resolving to the same symbol are an error. Do not treat transcript
duplicates as independent genes. Paths resolve against the inventory directory
when relative; no manifest path silently resolves against an unrelated project.

For a biological R0 pilot the inventory must be identical to P0's. A reduced
alignment panel is only for debugging/cost estimation and is not the final
biological benchmark. Absolute path inventories may be regenerated on the
cluster, with reviewed checksums; never combine unlike inventories unnoticed.

## 5. Smoke run and strict resume

```bash
bash run_pipeline.sh smoke --run-id smoke-local-01 \
  --selection tests/fixtures/selection.tsv \
  --alignments tests/fixtures/alignment_manifest.tsv --execute

# Exact same command and inputs, with --resume:
bash run_pipeline.sh smoke --run-id smoke-local-01 \
  --selection tests/fixtures/selection.tsv \
  --alignments tests/fixtures/alignment_manifest.tsv --resume --execute
```

Without `--execute`, all launch commands only save/print a plan. Smoke is limited
to at most three synthetic files inside `tests/fixtures`; it cannot bypass
production gates for real alignments. Default smoke sends no annotation requests.
For a fully synthetic local enrichment test, additionally pass
`--settings tests/fixtures/toy_settings.json --annotation tests/fixtures/toy_annotations`.
The toy backend and its statistical choices are **not production choices**.

The launcher binds the run ID to code hashes, design, settings, selection,
alignment manifest, annotation input, installed runtime, profile and output/work locations. It
rejects incompatible resume and existing directories without a matching lock.
Resume targets the named run, not Nextflow's most recent unrelated session.
Never change a selection under the same run ID. Do not delete work/cache files
if you want completed tasks reused. The pipeline terminates on task failure;
there is no `ignore` fallback converting a failed job into a zero-result gene.

`tests/fail_once_python.py` is a smoke-only fault injector used to test partial
resume. It fails one task if `tests/fixtures/fail_once.marker` exists and then
removes the marker. On filesystems that do not preserve executable bits (such as
this pCloud mount), copy the injector to a POSIX temporary directory, make that
copy executable, set `CAAS_VALIDATION_TEST_ROOT` to this pipeline directory, and
pass its absolute path with `--python`. Never use it for production.

## 6. Optional staged interface (not required for direct CAAS launch)

This section describes the original staged interface retained for compatibility.
Its approval/pilot requirements do **not** apply to `launch_all_caas.sh` / stage
`benchmark`, which implements the current request described in section 0.
Separate statistical/enrichment review does not block CAAS discovery.

**Review first:** exact pool/strategy definitions; technical tests and reference
compatibility; primary metric/term universe/annotation engine/null size; pilot
resource estimate and explicit execution manifest. `inputs/approval.template.json`
must be filled by the reviewer with matching design, selection, alignment and
settings hashes. It starts unapproved. Random-family discovery additionally
requires reviewed statistical settings. Extension beyond 19 requires an actual
pilot-cost JSON path and its matching SHA-256. Approvals do not get filled in by
the pipeline.

1. Resolve P0/N0/N1/N2 reference compatibility. Raw files alone cannot establish
   the CAAStools version, original alignment collection, commands and complete
   successful-job inventory. The audit therefore makes no automatic imports.
   Either provide a verified checksummed import bundle or approve an isolated
   `reference-rerun`. Do not attach old results to regenerated cycle configs.
2. Run `deterministic`: N3, N4, N5, 300 cycles total.
3. Run `r0-pilot`: R0_001–R0_019, 1,900 cycles total, same full alignment inventory
   as P0 for a biological comparison.
4. Review measured pilot cost and queue/I/O behavior. Only then extend R0 to
   **99 total** using the separate 80-replicate continuation manifest.
5. Optionally run R1 and P1/N6 after review. P2 is a distinct auxiliary selection.
6. Assemble an explicit cohort containing the intended references, controls and
   all intended random hypotheses; run `summarize` without new discovery.

Example plan (replace paths with real cluster locations):

```bash
bash run_pipeline.sh deterministic --run-id deterministic-01 \
  --alignments inputs/alignment_manifest.tsv --profile cluster \
  --approvals inputs/approvals/deterministic.json

bash run_pipeline.sh r0-pilot --run-id r0-pilot-01 \
  --alignments inputs/alignment_manifest.tsv --profile cluster \
  --approvals inputs/approvals/r0-pilot.json

# AFTER the first 19 have completed and the pilot cost has been reviewed:
bash run_pipeline.sh r0-complete --run-id r0-extension-01 \
  --selection review/launch-manifests/r0-extension80.tsv \
  --alignments inputs/alignment_manifest.tsv --profile cluster \
  --approvals inputs/approvals/r0-extension.json
```

Add `--execute` only after approval. A full 99-selection also exists for a fresh
run, but using it after the pilot would unnecessarily rerun the first 19.
`package102.tsv` and `package123.tsv` are review/cohort inventories, not a silent
all-strategy default. Every supported production stage has its own explicit
selection; wrong-family rows are rejected.

`conf/cluster.config` defaults to `std-cpu`, 1 CPU / 2 GB / 30 min per CAAS task,
up to 100 queued tasks and a submission rate limit. Confirm these with the real
pilot; they are inherited starting points, not measured recommendations.
Preprocessing, aggregation and plots run on the Nextflow driver. The Conda wrapper
selects Python, Rscript and Nextflow from one environment; mixed external tools
are rejected in Conda mode. Use `--cluster-config /path/to/reviewed.config` for
additional `params.task_environment`, account, memory/time and queue settings.
The driver needs enough time/memory for all submissions/aggregation; the example
SLURM wrapper is not an automatically chosen optimal resource allocation.

To submit the reviewed driver explicitly:

```bash
export CAAS_VALIDATION_ROOT='/ABSOLUTE/CLUSTER/PATH/pipeline.update.260917'
sbatch submit_validation_slurm.sh deterministic --run-id deterministic-01 \
  --alignments "$CAAS_VALIDATION_ROOT/inputs/alignment_manifest.tsv" \
  --approvals "$CAAS_VALIDATION_ROOT/inputs/approvals/deterministic.json" --execute
```

The driver activates the dedicated environment automatically; `CAAS_DRIVER_PYTHON`
is no longer required. It does not install packages or transfer alignments.
Submit from the bundle directory, or set `CAAS_VALIDATION_ROOT` as above: SLURM's
spooled script location is not assumed to be the project directory. Direct
execution on a permitted compute host is also possible. No `sbatch` was run here.

## 7. Downstream contract, backgrounds and auditability

Each task retains both legacy cycle results and pooled events, plus the full
resolved command, config/pool/alignment/tool/settings hashes, coverage, runtime,
CPU, peak RSS and raw-output bytes in its receipt. A successful exit without
both valid output files is an error. Each hypothesis checks every expected gene
before producing summaries; missing, duplicated or foreign outputs block it.

The event-table contract is unchanged: primary events; nominal positional
`p < 0.05`; deterministic selection of one record per gene/position; existing
dense-cluster pruning (density ≥0.7, span ≥3, at least three positions). Query
genes require **one surviving position with FG ≥4 AND BG ≥4**. Support maxima at
different sites must never be combined. `legacy_support.py` is a preserved copy
of the existing downstream implementation; these definitions are not a new
validation of positional significance.

Backgrounds include successfully completed genes with no discoveries. In this
workflow a **coverage-testable** gene has at least one column with whole-alignment
gap proportion ≤0.5 and at least one frozen cycle with ≥3 observed non-gap species
on each side. This eligibility check is before allele-diversity and p-value
prefilters; a constant alignment is still a valid zero-discovery background gene.
It follows the existing tool's non-gap definition (every character other than
`-`), without silently adding new missingness rules. It is a documented
pattern-independent coverage definition, not equal CAAS detection probability.

`completed_background.tsv` separately preserves **all** successful alignments;
`background.tsv` is the coverage-testable subset. Per-gene coverage counts and
absent species remain available in `gene_inventory.tsv`. Alignment-specific
absence is not an unknown phenotype identifier. The cohort's common background
is the explicit intersection of strategy backgrounds; query exclusions are
listed. Entirely incompatible inventories/settings/tool versions or an empty
common background are rejected. Compatibility compares the discovery/filter
contract and downstream implementation separately from later enrichment/statistics
settings: choosing a reviewed annotation engine does not itself require repeating
unchanged CAAS discovery. Full original settings hashes remain in provenance.
Final enrichment is conditional on that
declared cohort, not a background silently changed midway through a benchmark.

Required summaries include `execution_completeness.json`, `positions.tsv`,
`gene_inventory.tsv`, `query.tsv`, `background.tsv`, `summary.json` and hashes per
hypothesis; `matched/common_background.tsv`, `common_queries.tsv`,
`background_intersections.tsv`, `strategy_comparison.tsv`; functional terms,
gene-list manifests, null metrics and uncertainty tables. Figures are standalone
base-R, PDF and 300-dpi PNG, with their R script and Markdown description copied
to the figure directory. No Excel/ODS output is generated.

### Explicit multi-run cohort

```bash
python3 scripts/make_cohort.py \
  --summary /ABSOLUTE/RESULTS/reference-01/P0/000/P0_000 \
  --summary /ABSOLUTE/RESULTS/deterministic-01/N3/000/N3_000 \
  --output inputs/cohort.example.tsv

python3 run_validation.py summarize --run-id summary-01 \
  --cohort inputs/cohort.example.tsv --settings inputs/settings.json --execute
```

List **all** intended hypotheses explicitly, not just the illustrative two.
The cohort includes each summary lock checksum; Nextflow stages the declared
summary directories. It never glob-discovers newly published live results.
The declared inferential random-family sizes must also be complete before
conditional tail fractions can be computed.

### Verified reference imports

`scripts/import_reference.py` accepts `--design`, `--hypothesis`, `--alignments`,
`--bundle`, `--approval`, `--output`. Bundle columns are `gene_id`, `receipt_path`,
`receipt_sha256`, `events_path`, `legacy_path` (explicit absolute source paths).
Every current-schema receipt must document the compatible historical command,
tool, original config, pool, settings, alignment and both raw-output checksums.
The import approval supplies `reference_compatibility_review: true`, reviewer,
bundle/inventory/config/pool hashes. The importer checks completeness and current
compatibility, copies into a **new** isolated directory, filters, and records
source paths/checksums without modifying originals. Old tables lacking this
evidence are not sufficient; do not manufacture receipts from assumptions.

## 8. Frozen enrichment and statistical controls

Production enrichment defaults to `backend: disabled`. No live-request code is
invoked. Preparation/discovery and query inventories still work, and functional
outputs state `blocked_annotation_not_frozen` rather than inventing a result.
There are two opt-in offline paths:

- `gprofiler-cache`: exact request + annotation-version keyed raw response
  cache, validated returned version, response checksum, all-results request and
  auditable intersections. Uses the existing human/g:SCS/custom-background
  settings and source list. Incomplete caches fail visibly. The historical
  version (`e114_eg62_p19_27110d83`) is documented, **not assumed available**.
- `local-hypergeom-bonferroni`: reviewed frozen term membership, explicit organism,
  sources, version and checksum, exact upper-tail hypergeometric with Bonferroni
  over all annotated terms present in the common effective background. This is
  **not** g:SCS. It needs `implementation_approved: true`, correction `bonferroni`
  and a consistency review. It reruns **all** observed and random queries under
  one implementation; it must not mix its p-values with historical g:Profiler
  p-values. Synthetic snapshots are forbidden for biological cohorts.

Snapshot: `metadata.json` with `organism`, `sources`, `annotation_version`,
`term_membership_sha256`; `term_membership.tsv` with `term_id`, `source`,
`term_name`, `gene`. Cache files are `<sha256>.json`, SHA-256 of canonical JSON
`{payload: ..., annotation_version: ...}`; include exact `request`,
`annotation_version`, raw `response` and canonical-JSON `response_sha256`.
`functional.py` exposes the exact payload/key implementation. No network
transmission is implemented here: a future cache acquisition must document
transmitted symbols, approval, rate limits and version consistency first.

G0 generates 1,000 independent lists of P0's common-query size from the common
background, without replacement within each list. G1 uses the same size from
N2's common query; an optional N3 control is also generated. Insufficient query
size is explicitly reported, never repaired by replacement. Every gene list is
saved with seed/hash. These controls address list size, not gene-specific CAAS
selection probabilities. The bound on maximal fold enrichment depends on the
selected background fraction; significant-term medians alone are not a fair
validation metric.

Before null production, review a separate settings JSON containing:

- `statistics.approved`, a supported `primary_metric`, explicit
  `null_families` and `expected_null_counts` (e.g. R0:19 or R0:99), and
  `multiple_testing` = `none-single-comparison` or `holm` over the declared
  family of comparisons;
- a fixed `term_universe` and, if used, explicit `term_groups`; label targets
  selected from the previous 14 P0 terms as **conditional on prior discovery**;
- annotation version/backend/source/correction agreement. Nonsignificant term
  effects are retained. Missing terms and empty queries are explicit statuses,
  not effects set to zero merely for nonsignificance.

Supported metrics: query/background fraction, mean fold enrichment over the
prespecified term universe, number of significant prespecified nonredundant
groups, or distinct genes supporting significant terms. The term-group map is
not automatically optimized from results. An undefined observed/null metric or
incomplete declared null family blocks inferential summaries. No primary
production statistic has been selected by this implementation.

The upper-tail fraction is `(1 + count(null >= observed)) / (B + 1)`, with null
median/quantiles, observed-minus-median effect, approximate Monte Carlo SE and
Clopper–Pearson interval for the underlying exceedance probability. The latter
is Monte Carlo sampling uncertainty, not a biological-effect interval; for
without-replacement finite randomization it is a conservative binomial-style
description. Minimum fractions are .05 with 19 and .01 with 99 independent
hypotheses. Cycles are not null replicates; GO/pathway terms are not independent
discoveries. Optional Holm adjustment applies to the prespecified comparison
family, not an after-the-fact term screen.

## 9. Measured pilot cost and handoff

```bash
python3 scripts/cost_report.py --run-root results/r0-pilot-01 \
  --production-alignments 16130 --trace results/r0-pilot-01/nextflow.ACTUAL.trace.tsv \
  --output review/r0-pilot-cost
```

Use the actual inventory count, not the example 16,130. At 16,130 alignments,
102 hypotheses imply 1,645,260 gene/hypothesis tasks: **not** a runtime estimate.
The report distinguishes measured resources from workload counts and refuses
to extrapolate synthetic smoke timing to biological cost. Review CPU/memory by
alignment length, failures (including aborted tasks), scheduler overhead,
alignment/code staging, driver aggregation memory and work-directory disk use.
Add actual cluster billing/accounting before a monetary estimate. At millions
of tasks, scheduler and I/O overhead may justify batching, but no unvalidated
batching or altered pooling semantics has been introduced here.

See `review/COST_REVIEW_TEMPLATE.md`, `review/SMOKE_TEST_REPORT.md` and
`review/reference-audit/reference_compatibility.tsv`. Historical discovery and
consolidated result directories have not been modified. Production and online
annotation requests remain unrun pending Fabio's review.
