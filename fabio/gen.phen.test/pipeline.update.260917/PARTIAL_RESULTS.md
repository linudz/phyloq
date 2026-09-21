# Continue pooled CAAS after task failures

In the direct benchmark launched by `launch_all_caas.sh`, CAAS task failures
(including SLURM timeouts) use `errorStrategy 'ignore'`, without automatic retries.
Other tasks continue. The 30-minute limit, 2 GB memory and shared CAAStools path
remain unchanged. Input validation, checksum/provenance violations during
assembly and downstream script failures still terminate the workflow.

## Reports

The consolidated list is:

`summaries/<run-id>/matched/non_completed_genes.tsv`

Columns: `hypothesis_id`, `strategy_id`, `replicate_id`, `gene`, `status`, `reason`.
Each row represents one expected gene/hypothesis with no successful declared
task output. Status is `not_completed`, NOT a negative CAAS result. This report
does not infer a specific cause: consult the Nextflow trace/log for timeout,
exit status, scheduler job ID and task work directory. A gene that failed in
two strategies appears twice. The table exists with a header even if empty.

Each hypothesis also has `non_completed_genes.tsv` and
`execution_completeness.json` under:

`results/<run-id>/<strategy>/<replicate>/<hypothesis>/`

The comparison table includes expected/completed/non-completed counts and
`status=partial` when applicable. `matched/cohort.json` records the overall
execution status and number of missing gene/hypothesis results. Even a
hypothesis with no successful tasks remains in the reports.

Gene inventories, queries and backgrounds include successful results only.
Failed genes are excluded, not counted as zero discoveries or untestable genes.
The common background is the intersection of successfully completed,
coverage-testable genes. If it is empty, query fractions are undefined and the
figure explicitly says they cannot be estimated. No null inference/enrichment
is introduced by this change.

## Resume the current shared-tool run

After confirming the previous driver and its tasks have stopped, update and run:

```bash
bash launch_all_caas.sh --run-id caas-five-shared-01 --resume
```

The reviewed compatibility bridge accepts the exact shared-tool version
0c55eff (implementation fingerprint 0557b68e...). It retains valid completed
CAAS cache: inputs, run_caas.py, the CAAS command and tool are unchanged.
Assembly and comparison tasks are recomputed for the new partial-result policy.
Other code/input/runtime changes are rejected. Keep work/, .nextflow/, inputs
and results. This is not a bridge from the earlier per-task-tool-copy workflow.

An ignored failed task may be attempted again on a future explicit resume;
there is no automatic retry loop within a run. Updated downstream summaries
replace earlier summaries, avoiding stale missing-gene lists after recovery.
Do not run concurrent drivers against the same run ID.

A successful driver now means the workflow finished, not that every gene
succeeded. Always inspect the missing-gene table and completeness status.
