# Cluster quickstart — current direct launch

From `pipeline.update.260917` on the cluster:

```bash
# Once, if the dedicated environment does not yet exist:
bash create_conda_environment.sh

# Replace the alignment path; this submits one SLURM driver.
bash launch_all_caas.sh --run-id caas-all-01 \
  --alignments-dir /ABSOLUTE/PATH/inputs/alignments
```

The script inventories the alignments and runs N3/N4/N5, all 99 R0 hypotheses,
all 99 R1 hypotheses and P1/N6 in one Nextflow run. **No pilot, approval JSON or
smoke test is required.** Integrity checks run automatically.

Default: 203 hypotheses / 20,130 cycles across hypotheses. Historical
P0/N0/N1/N2 and P2 are excluded unless explicitly requested using
`--include-references` and `--include-p2`.

For an interrupted run, repeat the same command with `--resume`.
Do not delete `work/` or `.nextflow/`; do not resume an active driver.

See [RUN_ON_CLUSTER.md](RUN_ON_CLUSTER.md) for the exact workload, all options,
alternative input patterns, resource settings, outputs and Conda setup.
This replaces the earlier quickstart that required staged approvals/pilots.
