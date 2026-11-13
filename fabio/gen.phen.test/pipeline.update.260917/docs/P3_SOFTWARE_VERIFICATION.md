# P3 software verification on 13 November 2026

P3 adds the 6 FG and 13 BG Venditti nonhuman Upper contrast hypothesis as a
separate immutable input bundle. The original frozen design and launch inputs
are unchanged. No cluster submission or biological discovery was performed.

## Completed checks

The existing suite plus the P3 input, launcher and cache tests passed: 71 tests
with no failures, including local Nextflow integration. A subsequent recovery
suite with the added P3 receipt provenance test passed all five tests.
This covers 72 distinct test methods overall.

The full synthetic workflow integration used the actual scientific CAAS runner
on one artificial 129-species alignment. The first run executed the five old
hypotheses. A resume with the new workflow reused all five CAAS tasks and the
original validation task, ran the one new P3 task and assembled a six-strategy
cohort. A third identical resume reused all six CAAS tasks. Original raw
receipts remained byte-identical. A separate mock runner test independently
counted exactly five old calls and one additional P3 call.

Input checks verify 20 distinct oriented nonhuman source pairs, agreement with
the original frozen phenotype, 6/13 disjoint pool membership, 10725 possible
independent 4-versus-4 combinations, 100 frozen unique cycles and inclusion of
every pool species. A modified pool is rejected by its checksum.
The old PREPARE_VALIDATION, VALIDATE_INPUTS and CAAS_POOLED definitions are
compared with the pre-P3 workflow fixture from commit 3744029. The scientific
run_caas.py checksum is independently checked against its unmodified version.

Launcher tests simulate a compatible earlier run, keep its launch lock and
selection bytes intact, add P3 with a separate declaration, verify 6 analyses,
430 pooled cycles and 12 tasks for two synthetic alignments, and reject omission
of P3 after its separate run lock has been written. Existing strict tests reject
unreviewed code, settings, environment and selection changes.
Recovery tests verify the P3 receipt against its own selection digest while
retaining the historical run identity and complete/noncompleted distinction.

## Reproduce bounded tests

Generate the artificial fixtures first: their alignment inventory contains
machine-specific paths and is intentionally not committed after regeneration.
Use the installed compatible dependencies from environment.yml:

```bash
python tests/make_fixtures.py
NF_TEST_COMMAND=nextflow python -m unittest discover -s tests -v
```

Integration ran with Nextflow 24.04.2 and local execution. The temporary test
runtime used the locally available Python scientific environment and a temporary
DendroPy 4.5.2 dependency; no shared environment was modified. This is a software
routing and provenance test, not certification of the cluster environment or
a complete recreation of every pinned package build.

Keep the actual cluster environment unchanged for resume. Valid old cache also
requires the same filesystem paths, work directories and Nextflow metadata.
Missing, failed or interrupted old tasks can be executed again. The supported
code bridge does not bypass runtime or input identity checks and does not cover
the earlier per-task CAAStools staging workflow.
