# CAAS pooled validation — reduced design

The standard launch runs five analyses. R0/R1 randomized null series are disabled.

| Analysis | Selection | FG/BG pool | Pooled cycles per gene |
|---|---|---|---:|
| N3 | Cercopithecidae phenotype extremes, size matched to P0 | 9/10 | 100 |
| N4 | Same 19 PSS-selected species, reassigned by phenotype | 9/10 | 100 |
| N5 | Phenotype extremes matching P0 genus quotas | 9/10 | 100 |
| P1 | PSS-selected pairs in six genera | 6/6 | 15 |
| N6 | Phenotype-extreme pairs in those six genera | 6/6 | 15 |

One task runs one alignment and one complete pooled analysis using the bundled
modified CAAStools `pooled-discovery`. The internal cycles are NOT separate jobs.
For 16,133 alignments: **80,665 tasks**.

R0/R1 are rejected by the launcher and custom launch selections. Their frozen
files remain for provenance only. We retain controlled comparisons, not an
empirical null distribution across randomized assignments. P0 remains the
historical reference and is not rerun by default.

## Launch

**Resource fix:** the reduced run from commit 287d4ab can resume through the
reviewed resource-only compatibility bridge. CPU/memory/time flags may change
on resume; analysis inputs may not. See [RESOURCE_RESUME.md](RESOURCE_RESUME.md).

See [RUN_ON_CLUSTER.md](RUN_ON_CLUSTER.md). Nextflow processes, pooled settings
and frozen species/cycle assignments are unchanged. No pilot, enrichment or
downstream inference is required.

```bash
bash launch_all_caas.sh --run-id caas-five-01
# After an interruption, with the same inputs/code/environment:
bash launch_all_caas.sh --run-id caas-five-01 --resume
```

The old 203-analysis run cannot be resumed under this different selection/code
identity. Start a new run ID; old cache reuse is not guaranteed. Subsequent
unchanged reduced runs support normal Nextflow resume.

## Outputs

- `results/<run>/<strategy>/<replicate>/raw/<gene>/`: raw pooled outputs.
- `results/<run>/<strategy>/<replicate>/<hypothesis>/`: assembled/filtered results.
- `summaries/<run>/`: descriptive comparisons.
- `logs/<run>/`: driver/Nextflow logs; see [LOGGING.md](LOGGING.md).
- `inputs/launches/<run>/`: exact selection and alignment inventory.

Keep `work/`, `.nextflow/`, launch inputs and locks for resume.
The independent Conda specification remains in `environment.yml`.
The prior documentation in `docs/archive/` is historical, NOT current launch
instructions; its randomized series are no longer enabled.
