# Recover summaries without rerunning CAAS

Stop the original run's driver and ensure no tasks from that run are still
active. Do not remove work/, .nextflow/, published results or launch inputs.
Run from this bundle, with the original environment available:

```bash
bash recover_results.sh \
  --run-id caas-five-shared-01 \
  --recovery-id recovery-01 \
  --log logs/caas-five-shared-01/20260921T141233519309Z.nextflow.log \
  --confirm-stopped
```

This submits ONE `CAAS_RECOVER` driver on std-cpu (1 CPU, 8 GB, 1 day).
It never starts Nextflow or CAAStools, nor submits per-gene jobs.
The original results and metadata are read-only. Existing recovery IDs are
refused; use a new recovery ID after a failed recovery.

Outputs: `summaries/caas-five-shared-01/recovery/recovery-01/`:

- `hypotheses/`: assembled and filtered summaries per hypothesis.
- `matched/strategy_comparison.tsv`: comparison across the five analyses.
- `matched/non_completed_genes.tsv`: missing gene/hypothesis results, with
  `failed` or `submission_failed` from the source log (not inferred timeouts).
- `figures/`: PDF/PNG summaries and their R script.
- `source_task_audit.tsv`: original task states and source receipt checksums.
- `recovery.complete.json`: written only after successful tables and plots;
  separates recovery completion from partial biological completeness.

Driver logs: `logs/<run-id>/recovery-<recovery-id>/driver-<job-id>.log`.

The recovery binds the source log to its explicit launch params, verifies the
selection/inventory/settings/design against the original launch lock, and
requires terminal log records for every expected task. It uses only successful
or cached task outputs. Published receipts are preferred; if unavailable, the
recorded work directory/cache hash is resolved unambiguously. Receipt run IDs,
species/cycle config hashes, alignment hashes, settings, tool version and output
checksums must match. Missing/corrupt outputs for a logged success cause a
visible failure, not silent exclusion. A changed source log also fails recovery.

No global glob of other runs, directory copying, or CAAS rerun occurs. Existing
filter and summary functions are reused. No enrichment is performed.

Advanced local invocation (already inside the correct environment):

```bash
python recovery/recover_completed.py --run-id RUN --recovery-id RECOVERY \
  --log SOURCE_LOG --tables-only
```

`--tables-only` intentionally omits figures and records that in the completion
receipt; the cluster wrapper always generates figures.
