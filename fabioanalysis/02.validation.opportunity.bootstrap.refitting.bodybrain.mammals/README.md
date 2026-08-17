# Mammalian brain/body PSS validation

SLURM/Nextflow workflow for opportunity-depth validation and parametric
bootstrap/refitting of three mammalian traits:

- `log_body_mass`
- `log_brain_mass`
- `relative_brain_mass`

The directory is self-contained for cluster execution. Frozen observed-state
files are under `inputs/states/`; the exact PSS implementation used to create
them is under `inputs/core/`. The workflow uses the repository-level Conda
environment at `../../environment.yml` (`phyloq`).

## Before the first run on the cluster

From the root of the cloned `phyloq` repository:

```bash
conda env update -n phyloq -f environment.yml --prune
conda activate phyloq
```

## Pilot

The pilot runs 3 conditional simulations and 1 full-refit simulation per
trait:

```bash
cd fabioanalysis/02.validation.opportunity.bootstrap.refitting.bodybrain.mammals
bash submit_pilot.sh
```

Inspect `results/pilot/execution_trace.tsv` and the SLURM logs before changing
the production resource requests.

## Production

The current production plan is 1,000 conditional and 100 full-refit
simulations per trait:

```bash
bash submit_production.sh
```

Nextflow writes `work/` beside the workflow so that the driver and all SLURM
worker nodes see the same task directories on shared SamanthaFS. The Conda
cache uses `$SCRATCH` when that variable is defined and otherwise stays beside
the workflow. Published outputs remain under `results/pilot/` or `results/production/`.
Both results and logs are ignored by Git, except for their `.gitkeep` files.

## Returning results

Copy the completed `results/` directory to the local project folder:

`score.approach/04.case.study.bodybrain.mammals/02.validation.opportunity.bootstrap.refitting.bodybrain.mammals.fromcluster/`

Input checksums are recorded in `MANIFEST.sha256`.
