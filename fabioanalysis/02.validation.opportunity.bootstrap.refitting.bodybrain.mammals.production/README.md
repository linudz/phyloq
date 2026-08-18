# Mammalian brain/body PSS production validation

Production-only SLURM/Nextflow workflow for opportunity-depth validation and
parametric bootstrap/refitting of three mammalian traits:

- `log_body_mass`
- `log_brain_mass`
- `relative_brain_mass`

This directory is the production sibling of the completed pilot workflow in
`../02.validation.opportunity.bootstrap.refitting.bodybrain.mammals/`.
It contains frozen observed states and the exact PSS implementation required
by the analysis.

## Production design

For each trait the workflow runs:

- 1,000 conditional bootstrap simulations, grouped 25 per SLURM job;
- 100 full bootstrap/refitting simulations, one per SLURM job.

Across the three traits this produces 3,300 simulations in 420 worker jobs.
At most 100 worker jobs run concurrently. Aggregation begins only after every
worker result has completed successfully.

Resource requests were calibrated from the completed pilot: conditional jobs
request 2 GB RAM and full-refit jobs request 4 GB. Correfoc `std-cpu` billing
therefore remains at one computing unit per running job. The pilot-based
expected internal UPF compute charge is approximately EUR 3, subject to actual
runtime, retries and the applicable project rate.

## Cluster setup

From the root of the cloned `phyloq` repository:

```bash
git pull
conda env update -n phyloq -f environment.yml --prune
conda activate phyloq
cd fabioanalysis/02.validation.opportunity.bootstrap.refitting.bodybrain.mammals.production
sha256sum -c MANIFEST.sha256
```

## Launch

```bash
bash submit.sh
```

The command prints the SLURM identifier of the Nextflow driver. Monitor the
driver and worker jobs with the usual SLURM queue commands. Re-running
`bash submit.sh` resumes completed Nextflow tasks because the driver invokes
Nextflow with `-resume` and `work/` is on shared SamanthaFS.

## Outputs

Final outputs and execution reports are written under `results/`. Individual
bootstrap chunks are in `results/chunks/`; the three aggregate tables are
written directly under `results/`. SLURM driver logs are written under `logs/`.

Copy the completed `results/` directory back to the local case-study results
folder for interpretation. Results, logs, work directories and Conda caches
are intentionally ignored by Git.
