# Run all CAAS validation modes on the cluster

Current entry point, updated following Fabio's instruction: **one launch, no
mandatory pilot, no approval JSON, no smoke-test prerequisite**. Conda and the
automatic integrity checks are the only software setup requirements.

Copy this complete pipeline bundle to the cluster, including `bin/caastools`,
all scripts and all `inputs/frozen/` files. Keep alignments separate. Do not copy
local `work/`, `.nextflow*`, installed environments or previous development
outputs as a cluster cache. Run the commands below from `pipeline.update.260917`.

## Minimal launch

```bash
# Once; skip if the dedicated environment already exists:
bash create_conda_environment.sh

# Reuse the previous CAAS inputs; submit ONE SLURM driver.
bash launch_all_caas.sh --run-id caas-all-01
```

The script activates Conda, automatically builds the alignment inventory and
explicit selection, then launches one Nextflow workflow. The CAAS tasks activate
the same environment on worker nodes. With the repository layout, the default
is **`../caas/inputs/alignments/*.phy`**, exactly the collection configured by
the previous CAAS pipeline. Do not move or copy those alignments into the new
bundle. In the original research-project layout, the matching default is
`../pipeline/inputs/alignments/*.phy`. Only a standalone bundle without a legacy
sibling defaults to its own `inputs/alignments/*.phy`. Missing directories or
an empty default glob fail before submission, showing the expected path.

The diagnostic phylogeny is already bundled at
`inputs/frozen/brain-260917-v1/sources/tree.nwk`. CAAS pooled does not require
gene trees. The trees configured by the older bootstrap RERconverge workflow
are unrelated inputs and are not substituted here. See `INPUT_PATHS.md`.

No package installation, pilot, synthetic experiment or statistical approval is
performed inside the production launcher. Input hashes, the bundled CAAStools
version, pool membership and the frozen 4-vs-4 cycle configurations are checked
automatically. A malformed input causes a visible error.

If Conda is not on PATH, load the cluster's module or export
`CAAS_CONDA_SH=/path/to/etc/profile.d/conda.sh` before setup/launch. For a custom
installation prefix, keep `CAAS_CONDA_PREFIX` exported. Both the prefix and hook
must be accessible from all compute nodes.

## What is selected

| Mode | Strategies | Hypotheses | Cycles across hypotheses |
|---|---|---:|---:|
| deterministic | N3, N4, N5 | 3 | 300 |
| r0 | R0_001 through R0_099 | 99 | 9,900 |
| r1 | R1_001 through R1_099 | 99 | 9,900 |
| paired | P1, N6 | 2 | 30 |
| **Default `all`** | **All the above** | **203** | **20,130** |

Each hypothesis is applied to every alignment. At 16,133 alignments this is
**3,274,999 CAAS tasks**, with up to 100 outstanding tasks under the current
configuration, not 203 SLURM jobs. All modes belong to one Nextflow run; task
submission order and available resources determine actual overlap. The existing
task granularity is unchanged: budget adequate shared disk space for per-gene
work directories and raw outputs. There is no automatic cleanup.

Historical P0/N0/N1/N2 are **not silently rerun**. To include them, add
`--include-references` (207 hypotheses total). To also include P2, add
`--include-p2` (208 total). These are new results under the chosen run ID, not
imports or overwrites of historical analyses.

To select only certain modes, for example the 102-hypothesis main package:

```bash
bash launch_all_caas.sh --run-id caas-main-01 \
  --alignments-dir /ABSOLUTE/PATH/inputs/alignments \
  --modes deterministic,r0
```

## Alternate sources, resources and resume

For a directory containing non-alignment files, use a quoted pattern instead:

```bash
bash launch_all_caas.sh --run-id caas-all-01 \
  --alignments-pattern '/ABSOLUTE/PATH/inputs/alignments/*.filter2'
```

Use the real filename suffix. An existing inventory can be supplied with
`--inventory /ABSOLUTE/PATH/alignment_manifest.tsv`. Its relative alignment paths
are resolved against that inventory, not against the launch directory.

The default driver requests 1 CPU, 8 GB and 3 days. Worker defaults are `std-cpu`,
1 CPU, 2 GB, 30 minutes. These are resource limits, not runtime estimates. For
different driver resources use `sbatch` directly; it will not submit a second
driver from inside the allocation:

```bash
sbatch --time=7-00:00:00 --mem=16G launch_all_caas.sh \
  --run-id caas-all-01 --alignments-dir /ABSOLUTE/PATH/inputs/alignments
```

Use a time/memory request permitted by your cluster. Worker resources and queues
can be overridden using `--cluster-config /ABSOLUTE/PATH/custom.config` with
Nextflow parameters such as `params.partition`, `params.account`,
`params.task_cpus`, `params.task_memory`, `params.task_time` and executor settings.
When submitting from outside this directory, export `CAAS_VALIDATION_ROOT` to
the bundle's absolute path.

For an interruption, repeat the exact original command with **`--resume`**:

```bash
bash launch_all_caas.sh --run-id caas-all-01 \
  --alignments-dir /ABSOLUTE/PATH/inputs/alignments --resume
```

Keep the same modes, input source, code, environment and cluster configuration.
Keep `work/` and `.nextflow/`; do not submit a second driver for a run that is
still active. Changed inputs/settings require a new run ID. `--plan-only` is an
optional no-submission preview, not a compulsory step.

## Outputs and analysis scope

- `inputs/launches/<run>/`: exact selection, alignment inventory, checksums and
  a readable launch summary.
- `results/<run>/<strategy>/<replicate>/raw/<gene>/`: pooled CAAS tables,
  species-level events, logs and execution receipts.
- `results/<run>/<strategy>/<replicate>/<hypothesis>/`: assembled, filtered
  positions and per-hypothesis gene/background tables.
- `summaries/<run>/`: descriptive comparisons, gene lists and figures.
- `slurm-<job-id>.out`: driver output; Nextflow traces/logs remain in their
  run-specific directories.

The direct launcher runs **CAAS discovery, filtering and descriptive summaries**.
It does not run enrichment, gene-list resampling controls or inferential null
comparisons. Those can use the finished results later through `summarize`;
their undecided settings do not block discovery. The functional status is
`not_requested_discovery_only`, not an enrichment failure.

The old staged commands remain available but are not required for this launch.
