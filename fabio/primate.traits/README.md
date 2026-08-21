# Primate PSS bootstrap — Correfoc export

This is a self-contained export of **only the parametric-bootstrap analysis**.
It does not depend on the original project path and does not rerun trait
filtering, empirical PSS estimation, or taxonomic classification.

The frozen input contains 153 retained traits (redundancy below 10% and more
than 20 observed species), the corresponding classified empirical PSS tables,
the S4 tree, and the taxonomy used by the analysis.

## What it runs

For every trait and simulation, the workflow performs both versions:

1. **Conditional bootstrap:** simulate under the observed selected model and
   fitted parameters, then calculate PSS with those parameters held fixed.
2. **Full bootstrap:** simulate under the observed selected model, refit BM and
   OU, select the model using the current AIC rule, and recalculate PSS.

For both top and bottom 1% tails, it records family/superfamily/order excess,
median depth percentile, and the distributions of S, Delta, T, and PSS.

## Configure

Edit `params.yaml`. The main controls are at the top of that file:

- `bootstrap_replicates`: total replications per trait (default 1,000).
- `replicates_per_task`: replications in each Slurm job (default 100).
- `max_parallel_jobs`: throttle for simultaneous Slurm jobs.
- `slurm_partition`, `slurm_account`, and `slurm_qos`: fill only if required.
- per-job memory and time limits.

At the defaults, the workflow creates 10 chunks per trait: **1,530 independent
bootstrap jobs**, followed by one validation and aggregation job. Simulation
IDs and random seeds are global, so changing chunk size does not change the
simulated datasets.

## Install once on Correfoc

From the export directory:

```bash
conda env create -f environment.yml
conda activate primate-pss-bootstrap
```

Nextflow also uses the same YAML for every Slurm task. Its cached task
environment is stored in `.nextflow-conda`, which must be visible to compute
nodes. Keep the export directory on a shared filesystem.

## Run on Slurm

```bash
nextflow run main.nf -profile slurm -params-file params.yaml -resume
```

The same command can be safely issued again after an interruption. Nextflow
will reuse completed chunks and submit only missing or failed work.

## Small local validation

Before submitting the complete run, this command tests two simulations of one
trait using the already active R environment:

```bash
nextflow run main.nf -profile local -params-file params.yaml \
  --trait_regex '^Body_mass$' \
  --expected_traits 1 \
  --bootstrap_replicates 2 \
  --replicates_per_task 1 \
  --output_dir results.smoke
```

## Outputs

Partial chunk files are written to `results/chunks/`. Final validated tables
are written directly under `results/`:

- `primate.traits.bootstrap.simulations.tsv`
- `primate.traits.bootstrap.observed.tsv`
- `primate.traits.bootstrap.fits.tsv`
- `primate.traits.bootstrap.summary.tsv`
- `primate.traits.bootstrap.model.stability.tsv`

The aggregator stops with an error if a trait, simulation number, bootstrap
type, or tail is missing or duplicated. Execution report, timeline, and trace
files are also saved under `results/`.

## Frozen-input provenance

- `inputs/traits/`: filtered trait tables used in the current score approach.
- `inputs/classified/`: observed PSS results with taxonomic levels.
- `inputs/tree/`: S4 primate phylogeny.
- `inputs/taxonomy/`: species-to-family and primate-group mapping.
- `inputs/metadata/`: domain mapping and filtering audit, included for context.

No input from the old `statistical.test.approach` is used.

## Continuous phylogenetic-depth profiles

The independent workflow under `depth_profiles/` recreates the same
deterministic bootstrap while retaining the complete depth distribution of
top and bottom PSS tails. It produces cumulative curves, pointwise and
simultaneous null envelopes, calibrated global curve statistics, excursion
intervals and BH-adjusted trait-level probabilities.

Its outputs are written to the separate, Git-ignored
`depth_profile_results/` directory. See `depth_profiles/README.md` for the
Correfoc command, parameters, outputs and validation rules.
