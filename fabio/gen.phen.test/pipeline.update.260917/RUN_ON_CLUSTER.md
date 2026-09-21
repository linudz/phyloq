# Five-analysis pooled CAAS launch

Run from `phyloq/fabio/gen.phen.test/pipeline.update.260917`.
The existing shared Conda environment and its hook/prefix overrides are unchanged.
For fresh installation, see `environment.yml` and `conf/conda_helpers.sh`.

## Transition from the large run

First stop the old driver and any remaining CAAS jobs belonging to that run.
Inspect `squeue -u "$USER"`; do not cancel unrelated jobs. Preserve existing
results and Nextflow metadata. Do not edit a live run or run two drivers against
the same output directory.

After updating the checkout, use a NEW run ID:

```bash
bash launch_all_caas.sh --run-id caas-five-01
```

One driver is submitted to `std-cpu`; no nohup is needed. The default source
remains `../caas/inputs/alignments/*.phy` in phyloq, or
`../pipeline/inputs/alignments/*.phy` in the research copy.
Explicit `--alignments-dir` and `--inventory` remain supported.

Default `--modes all` selects N3/N4/N5/P1/N6: five analyses, 330 internal pooled
cycles across analyses, 80,665 tasks for 16,133 genes. R0/R1 are rejected.
You can select `--modes deterministic` or `--modes paired` separately.
Historical references/P2 are opt-in only.

## Resume

```bash
bash launch_all_caas.sh --run-id caas-five-01 --resume
```

Repeat the original options and keep code, inputs, settings and environment
unchanged. Keep `work/`, `.nextflow/`, `inputs/launches/` and run locks.
The launcher invokes Nextflow's named resume.

Do NOT reuse the old 203-analysis run ID: strict locks reject the changed
selection/code. Cache reuse from that older run is not guaranteed; the first
reduced launch is a new run. Its subsequent resumes are supported.

## Monitor

```bash
squeue -u "$USER"
ls -lt logs/caas-five-01/
tail -f logs/caas-five-01/driver-<job-id>.log
```

Results: `results/caas-five-01/`; descriptive summaries:
`summaries/caas-five-01/`. See [README.md](README.md) and [LOGGING.md](LOGGING.md).
