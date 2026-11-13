# Six analysis pooled CAAS launch

P3 adds the Venditti nonhuman Upper contrast pool with 6 FG and 13 BG.
Its 100 pooled cycles are defined in [UPPER20_P3.md](UPPER20_P3.md).
P1 retains its original linked genus comparison; it is not renamed.

## Resume the existing five analysis shared tool run

Wait for the driver to finish or stop that run and its own workers before
updating. Do not edit files under active jobs or start concurrent drivers.
From the same checkout and working directory, keep the original run ID:

```bash
git pull --ff-only
bash launch_all_caas.sh --run-id YOUR_EXISTING_RUN_ID --resume
```

Repeat any original alignment source, mode and resource options. Default
`--modes all` keeps the frozen N3/N4/N5/P1/N6 selection and appends P3 through
a separate validated input channel. Completed intact CAAS cache can be reused;
P3 contributes one new task per alignment. Old failed or missing tasks may
be retried. Downstream combined summaries are refreshed for the expanded cohort.

Keep `work/`, `.nextflow/`, `inputs/launches/`, results and all launch locks.
The original launch inputs and lock are not rewritten. A separate additional
hypothesis lock binds P3 to the run. Code and runtime must match the supported
compatibility contract; no general relaxation of input checks is introduced.

For a previous `--modes deterministic` or `--modes paired` run, repeat that
exact mode and add `--include-upper20`. Once P3 is bound, it cannot be omitted
on subsequent resumes. `--without-upper20` is available for unchanged legacy
runs to which P3 has not been added.

## Fresh launch

Run from `phyloq/fabio/gen.phen.test/pipeline.update.260917`:

```bash
bash launch_all_caas.sh --run-id caas-six-01
```

One driver is submitted to `std-cpu`; no nohup is needed. The existing shared
Conda environment and its overrides are unchanged. For fresh environment
installation see `environment.yml` and `conf/conda_helpers.sh`.
The default alignment source remains `../caas/inputs/alignments/*.phy` in
phyloq or `../pipeline/inputs/alignments/*.phy` in the research copy. Explicit
`--alignments-dir` and `--inventory` remain supported.

The default six analyses have 430 internal cycles and 96798 tasks for 16133
genes. Internal cycles are not individual scheduler tasks. R0/R1 remain
disabled. Historical P0/N0/N1/N2 and auxiliary P2 remain opt in.
The task resource defaults remain 30 minutes, 2 GB and 1 CPU; repeat nondefault
resource flags on resume. The launcher calls Nextflow named `-resume`.

The old 203-analysis launch and the workflow before shared-tool staging are
not covered by this P3 extension bridge. Their transition still requires a
new run ID; do not force their selection or locks into a shared-tool run.

## Monitor and interpret completion

```bash
squeue -u "$USER"
ls -lt logs/YOUR_EXISTING_RUN_ID/
tail -f logs/YOUR_EXISTING_RUN_ID/driver-<job-id>.log
```

Results are under `results/<run-id>/`, comparisons under `summaries/<run-id>/`.
Inspect `matched/non_completed_genes.tsv` and each hypothesis execution
completeness report. Successful driver completion does not mean every gene
succeeded. See [PARTIAL_RESULTS.md](PARTIAL_RESULTS.md), [LOGGING.md](LOGGING.md)
and [RECOVER_RESULTS.md](RECOVER_RESULTS.md).
