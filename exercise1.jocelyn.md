# Exercise 1: BodyMass_kg PhyloQ Bootstrap Check

## Before Starting: Update The Repository

Before running the exercise, make sure that your local copy of the repository is
up to date.

Open a terminal and move into the `phyloq` repository:

```bash
cd /path/to/phyloq
```

Then update the repository from GitHub:

```bash
git pull
```

This step is important because the exercise depends on recent changes to the
Nextflow pipeline, the alignment files, and the configuration files.

After `git pull`, move into the bootstrap directory:

```bash
cd bootstrap
```

All commands in the rest of this exercise should be run from inside the
`bootstrap` directory.

## Background

We tested the R implementation of PhyloQ on the primate phenotype `BodyMass_kg`.
The analysis was run as a standalone R workflow and used the available trait table
and phylogenetic tree to estimate pairwise PhyloQ statistics among species.

The workflow performed the following steps:

1. Read the primate trait table and phylogenetic tree.
2. Selected the phenotype `BodyMass_kg`.
3. Pruned the tree to species with available `BodyMass_kg` values.
4. Fitted a Brownian Motion model with `geiger::fitContinuous()`.
5. Simulated null trait datasets under the fitted model.
6. Standardized simulated traits by z-score.
7. Computed pairwise absolute phenotypic distances.
8. Estimated `Nu` values for all species pairs.

The main output table contains all pairwise comparisons and includes observed
phenotypic differences, phylogenetic distances, empirical statistics, and the
smooth PhyloQ statistic `Nu`.

## BodyMass_kg Results

The `BodyMass_kg` run produced pairwise PhyloQ results for primate species.
We then annotated the pairwise results with family-level information and isolated
within-family contrasts.

For this exercise, the most important downstream result is the set of
Cercopithecidae configuration files:

- `bootstrap/cfg.files/cfg/div/`
- `bootstrap/cfg.files/cfg/ctrl/`

The `div` configuration files represent divergent within-family contrasts.
The `ctrl` configuration files represent control within-family pairs.
Each configuration file contains eight Cercopithecidae species:

- four species assigned to group `1`
- four species assigned to group `0`

These files are used as input for `caastools`.

## Nextflow Pipeline

The Nextflow pipeline is defined in:

```text
bootstrap/launch_pipeline.nf
```

The pipeline runs `caastools discovery` by combining:

- the alignments in `bootstrap/alignments/`
- the configuration files in `bootstrap/cfg.files/`

The pipeline was updated so that configuration files can be selected with an
explicit glob pattern. This makes it possible to run the divergent and control
groups separately.

## Exercise Goal

The goal of this exercise is to check whether the Nextflow pipeline runs
correctly on the available alignments using the `BodyMass_kg` configuration
files.

The exercise is not to interpret biological results yet. The goal is only to
verify that:

1. Nextflow starts correctly.
2. The available alignments are detected.
3. The selected configuration files are detected.
4. `caastools discovery` runs for every alignment/configuration combination.
5. Output `.caas` files are written without errors.

## Commands To Run

From the `bootstrap` directory, run the divergent configuration files:

```bash
nextflow run launch_pipeline.nf \
  --configs 'cfg.files/cfg/div/*.cfg' \
  --outdir 'pipeline.results/bodymasskg_div'
```

Then run the control configuration files:

```bash
nextflow run launch_pipeline.nf \
  --configs 'cfg.files/cfg/ctrl/*.cfg' \
  --outdir 'pipeline.results/bodymasskg_ctrl'
```

## Expected Output

If there are four alignments and 100 configuration files in each group, each run
should produce 400 tasks:

```text
4 alignments x 100 cfg files = 400 caastools runs
```

The expected output folders are:

```text
bootstrap/pipeline.results/bodymasskg_div/
bootstrap/pipeline.results/bodymasskg_ctrl/
```

Each folder should contain `.caas` output files grouped by alignment.

## Basic Verification

After the runs finish, check the number of `.caas` files:

```bash
find pipeline.results/bodymasskg_div -type f -name '*.caas' | wc -l
find pipeline.results/bodymasskg_ctrl -type f -name '*.caas' | wc -l
```

With four alignments and 100 configuration files, both commands should return:

```text
400
```

## Notes

The current exercise is a technical validation step. If the pipeline completes
successfully, the next step is to inspect the `.caas` outputs and decide how to
summarize the results across divergent and control groups.
