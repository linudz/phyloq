# Continuous phylogenetic-depth profiles

This independent Nextflow workflow recreates the validated 1,000-replicate
parametric bootstrap for all 153 retained primate traits and retains the
information needed to analyse the complete phylogenetic-depth distribution of
the top and bottom PSS tails.

The simulation seed formula, BM/OU fitting, AIC model selection, PSS
calculation, tail rule, tree and frozen trait inputs are the same as in the
completed bootstrap. Empirical BM/OU parameters are frozen from that completed
run in `inputs/observed_model_fits.tsv`, eliminating the tiny optimizer-level
variation that previously occurred when each chunk refitted the observed
trait. The run is deterministic from this point forward. Because the earlier
run did not serialize full-precision in-memory covariance objects, the new
replicates are calibrated redraws from the same fitted null rather than a
promise of bit-for-bit identity with its discarded pairwise tables.

## What is retained

For every conditional and full-bootstrap replicate, and for both top and
bottom 1% tails, each chunk RDS stores:

- the legacy taxonomic and median PSS-component statistics;
- model selection and AIC values;
- the opportunity-relative depth percentile of every selected pair.

The depth percentile is the average rank of a pair's patristic distance among
all pairs available for that trait, minus 0.5 and divided by the number of
pairs. Selected-depth vectors are more compact than complete pairwise PSS
tables and allow the cumulative curves to be regenerated on any grid without
rerunning the bootstrap.

## Configuration

Edit `params.yaml`. Defaults are the intended full analysis:

- 153 traits;
- 1,000 simulations per trait;
- 100 simulations per Slurm task;
- top and bottom 1% tails;
- conditional and full bootstraps;
- 101 depth-grid points;
- at most 200 simultaneous jobs on `std-cpu`.

The output directory is `depth_profiles/results/`. It contains a local
`.gitignore`, so generated results can change on Correfoc without affecting
pulls. The workflow's `work/`, Nextflow cache and Conda cache are also ignored.

## Run on Correfoc

Create the existing repository environment once if needed:

```bash
conda env create -f environment.yml
```

Submit the driver **from repository root**, where the `logs/` directory is
already present:

```bash
sbatch depth_profiles/submit_depth_profiles.sh
```

The wrapper activates `primate-pss-bootstrap` and runs:

```bash
cd depth_profiles
nextflow run main.nf -c nextflow.config \
  -profile slurm -params-file params.yaml -resume
```

`-resume` reuses completed tasks from this depth-profile workflow. It cannot
reuse the previous aggregate-only bootstrap tasks, because those tasks did not
retain selected depth distributions.

## Outputs

The ignored `depth_profiles/results/` directory contains:

- `chunks/*.depth_profiles.rds`: compact selected-depth vectors and replicate
  statistics for each trait and chunk;
- `primate.traits.depth_profiles.replicates.tsv.gz`: per-replicate legacy
  statistics retained for quality control;
- `primate.traits.depth_profiles.observed.tsv`: observed tail summaries;
- `primate.traits.depth_profiles.fits.tsv`: empirical BM/OU fits;
- `primate.traits.depth_profiles.curves.tsv`: observed curve, null median,
  pointwise interval, simultaneous band and deviation at every grid point;
- `primate.traits.depth_profiles.summary.tsv`: signed/absolute areas, maximum
  deviations, empirical probabilities, BH-adjusted probabilities, model
  stability, resolution flags and excursion intervals;
- Nextflow report, timeline and trace.

The full bootstrap is the primary inferential result. The conditional
bootstrap remains a diagnostic.

## Aggregation and inference

For each trait, bootstrap type and tail, the aggregator calculates:

- the median null cumulative curve;
- pointwise 2.5% and 97.5% limits;
- a 95% simultaneous unstudentized maximum-deviation envelope;
- signed and absolute area between observed and null curves;
- maximum absolute, deep and shallow deviations;
- plus-one empirical probabilities;
- BH correction across 153 traits within bootstrap type and tail;
- the number, position, length and direction of simultaneous-band excursions.

Resolution is flagged as `very_low` for fewer than five selected pairs,
`cautious` for five to nine, and `shape_suitable` for at least ten.

## Small local validation

From `depth_profiles/`, with R, ape and geiger available:

```bash
nextflow run main.nf -c nextflow.config \
  -profile local -params-file params.yaml \
  --trait_regex '^Body_mass$' \
  --expected_traits 1 \
  --bootstrap_replicates 2 \
  --replicates_per_task 1 \
  --output_dir results.smoke
```

The aggregator validates complete simulation IDs, unique keys, selected-depth
lengths, medians, nondecreasing curves, endpoints and consistency of repeated
observed fits before writing final tables.
