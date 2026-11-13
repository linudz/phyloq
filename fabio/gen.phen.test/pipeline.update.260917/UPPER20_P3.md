# Venditti Upper contrast hypothesis P3

P3 adds the Venditti nonhuman contrast hypothesis to the five existing direct
CAAS analyses. P1 already identifies linked PSS pairs in six genera and P2 is
an archived auxiliary strategy, so the new identifier is **P3_000**. No Smaers
hypothesis is included.

## Frozen selection

1. Use the saved primate wide Venditti PSS ranking of 129 species and 8256 pairs.
2. Retain the original Upper ceil(1 percent), 83 pairs.
3. Exclude pairs containing Homo sapiens, leaving 59 pairs. Do not refit PSS.
4. Rank those pairs by absolute difference of the frozen adjusted phenotype and
   retain the first 20, with deterministic pair key tie handling.
5. Put higher endpoints in FG and lower endpoints in BG, then deduplicate.

The result has 6 FG and 13 BG, without conflicting membership. Selected
contrasts range from 0.4402262664 to 1.3638193384. This is a fixed count amplitude
filter within a PSS tail, not the historical P0 Lower tail linking rule, a
common minimum amplitude, or a measurement error calibration.

| FG | BG |
|---|---|
| Ateles paniscus | Alouatta belzebul |
| Cebus albifrons | Alouatta palliata |
| Macaca assamensis | Ateles belzebuth |
| Macaca nemestrina | Callithrix jacchus |
| Macaca nigra | Macaca fuscata |
| Semnopithecus entellus | Macaca silenus |
| | Macaca tonkeana |
| | Saguinus labiatus |
| | Saguinus midas |
| | Trachypithecus auratus |
| | Trachypithecus cristatus |
| | Trachypithecus geei |
| | Trachypithecus phayrei |

Columns list independent pools; rows do not define linked species pairs.
The original 20 species pairs define the hypothesis pool, not restrictions
on subsequent cycle combinations. Like N3 to N5, P3 samples independent 4 FG
versus 4 BG cycles from its fixed pools. The Cartesian universe has
choose(6,4) times choose(13,4) = 10725 distinct cycles. The existing SHA256
ranking convention, seed 20261113 and label P3_000 freeze 100 unique cycles.
All 19 species occur in those cycles. A cycle is internal to pooled CAAS,
not a separate scheduler task or an independent experimental replicate.

The checksummed bundle is `inputs/additional/venditti-upper20-261113-v1/`.
It contains full precision endpoint values, source pairs with original PSS
ranks, pool adapters, cycle adapters and a one-row execution selection.
The source is the research audit directory
`venditti.case/19.upper20.largest.contrasts.pools.261113/`.
Pool phenotype values match the frozen historical primate phenotype.
PSS scores remain ranking indices, not pairwise significance tests.

## Resume an existing shared tool run

First let the existing driver finish or stop that run and its own worker jobs.
Update the checkout only after its jobs have stopped. From the same directory,
with the original environment, alignment source and run ID, use:

```bash
git pull --ff-only
bash launch_all_caas.sh --run-id YOUR_EXISTING_RUN_ID --resume
```

Repeat any original alignment and resource options. `--modes all` is the
default and now appends P3. For a run originally limited to another mode,
repeat that exact mode and add `--include-upper20` explicitly. Do not switch
the frozen historical mode selection during resume.

P3 contributes **one additional task per alignment**. The default six analyses
have 430 internal cycles in total and 96798 gene/hypothesis tasks for 16133
alignments, of which 16133 are new P3 tasks. Intact completed historical CAAS
cache is reused; failed, interrupted or missing cached tasks may run again.
This is not a promise that every previous task is complete or reusable.

`--without-upper20` retains the legacy five-analysis launch before P3 has
been added to that run. After P3 is bound to the run, omitting or changing it
is rejected to prevent silently shrinking the declared cohort. P3 is automatic
only with `--modes all`; deterministic and paired launches keep their prior
meaning unless `--include-upper20` is supplied.

## Why the extension preserves the old task cache

The original frozen design, five-row selection, launch request, alignment
inventory and launch lock are not rewritten. The original PREPARE_VALIDATION,
VALIDATE_INPUTS and CAAS_POOLED process definitions and the scientific CAAS
runner are unchanged. A separate VALIDATE_ADDITIONAL_INPUTS process validates
P3 and feeds its gene jobs to the existing pooled process. Historical job
metadata and its original paths remain available through the cached original
validation process.

An exact code fingerprint bridge accepts only the published shared tool and
continuation versions while preserving the old task script fingerprint.
All other run identity and environment checks still apply. The actual current
implementation is recorded in every new launch plan. Fresh runs bind directly
to the current implementation. Older runs before shared tool staging still
require their documented new run transition and are not covered by this bridge.

The additional bundle gets its own `additional_hypotheses.lock.json` under
`results/<run-id>/`. Each new launch plan and parameter file declares P3
separately. The historical `inputs/launches/<run-id>/launch_summary.json`
continues to describe the unchanged base selection; the new launch plan is
the record of the complete six-analysis launch.

## Outputs and verification

Raw and filtered results follow the existing hierarchy under
`results/<run>/P3/000/`. The combined cohort, plots and noncompleted gene tables
include P3. Its input validation records have a `P3.` prefix in provenance.
The standalone stopped-run recovery utility also reads the additional bundle.
Failed and timed out genes remain noncompleted observations, never negative
CAAS results; tested backgrounds use completed eligible genes only.

Software tests check the exact 6/13 assignment, the 20 source pairs and
orientation, the 100 cycles and every species' inclusion, input checksum
rejection, original process definitions and cache reuse. Synthetic Nextflow
integration tests verify five old cached gene tasks plus one new P3 task,
followed by an entirely cached repeat. No biological or cluster CAAS run is
performed by adding the implementation.
