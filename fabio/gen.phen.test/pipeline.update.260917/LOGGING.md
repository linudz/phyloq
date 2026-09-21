# Execution logs

Recommended entry point (from the bundle):

```bash
bash launch_all_caas.sh --run-id caas-all-01
```

Logs are ignored by Git and grouped in `logs/<run-id>/`:

- `submission-<UTC timestamp>-<PID>.log`: submission or plan-only output.
- `slurm-<job-id>.out` and `.err`: scheduler stdout/stderr for automatic submission.
- `driver-<job-id>.log`: combined driver stdout/stderr, including Conda initialization errors before Nextflow starts.
- `<launch-id>.nextflow.log`: Nextflow engine log.
- `nextflow.<launch-id>.trace.tsv`: task execution trace.

Existing logs are not moved or deleted. New submissions use unique IDs; a requeued driver appends to its existing log.

For direct `sbatch`, submit from the bundle directory. The tracked `logs/.gitkeep` ensures that `logs/` exists before SLURM opens its output files. Direct submissions place the scheduler's bootstrap output in `logs/slurm-<job-id>.out/.err`; the driver still writes to the run-specific subdirectory. If submitting from elsewhere, supply absolute `sbatch --output` and `--error` paths to an existing directory and export `CAAS_VALIDATION_ROOT`.

Inspect a run:

```bash
ls -lt logs/caas-all-01/
tail -n 100 logs/caas-all-01/driver-<job-id>.log
```

Nextflow's native `.command.*` task logs remain in `work/`, and published per-gene `caastools.log` files remain alongside their results. These are not relocated or duplicated: keeping task directories intact preserves Nextflow's normal diagnostics and resume behavior. Launch plans and parameter JSON files remain in `launch-plans/`.
