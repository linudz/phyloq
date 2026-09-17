# Pilot resource review — must be completed before full null extension

Status: **NOT APPROVED**. Reviewer/date: pending.

## Inputs and biological comparability

- Design ID and SHA-256:
- Explicit selected manifest SHA-256:
- Alignment inventory count/SHA-256; identical to P0 for biological R0 pilot:
- Settings and CAAStools hashes:
- Reduced technical panel or full biological inventory:
- Pilot output/cohort paths and completeness:

## Measured resources (not extrapolated from synthetic fixtures)

- Gene/hypothesis jobs completed, failed, aborted and retried:
- CPU and elapsed distributions by alignment length/species coverage:
- Peak resident memory and resource requests:
- Scheduler queue delay/submission overhead and concurrency:
- Alignment/code staging volume, driver memory and I/O:
- Raw outputs, total work/cache disk amplification and projected storage:
- Actual cluster charging model, if a monetary estimate is requested:
- `scripts/cost_report.py` JSON path and SHA-256:

## Proposed workload

| Package | Hypotheses | Cycles | Tasks at 16,130 alignments |
|---|---:|---:|---:|
| N3/N4/N5 | 3 | 300 | 48,390 |
| R0 pilot | 19 | 1,900 | 306,470 |
| R0 extension only, after the pilot | 80 | 8,000 | 1,290,400 |
| N3/N4/N5 + R0, 99 total | 102 | 10,200 | 1,645,260 |
| Optional R1 pilot | 19 | 1,900 | 306,470 |
| Optional P1/N6 | 2 | 30 | 32,260 |
| Combined package | 123 | 12,130 | 1,983,990 |

Replace the example alignment count with the actual approved inventory. The
combined package excludes historical reruns/P2 and does not add the R0 pilot
twice. These are task counts, not elapsed time or monetary cost.

## Decision before full launch

- Approved primary metric, term universe/grouping, annotation version/engine:
- Exact random family and total null size:
- Approved selected manifest and existing pilot cohort:
- Resource adjustments; batching proposal and equivalence tests if necessary:
- Reviewer approval of production launch:

Do not change pooling semantics for throughput. Keep failure visibility and
gene/hypothesis provenance if a future batching implementation is approved.
