# Independent environment verification

Date: 17 September 2026.

## Scope

The bundle now carries its own `environment.yml`, environment creation/check
script, activation wrapper and cluster quickstart. It no longer needs the
historical `phyloq` environment or original project paths for preparation and
discovery. Biological alignments remain external and must be inventoried on the
execution machine. The optional historical-result audit still requires its
explicit historical inputs.

Environment activation applies to the Nextflow driver and its worker jobs.
Installation is performed once on a setup node, not inside individual tasks.
Existing environments are not updated without `--update`, and the setup refuses
base/shared `phyloq` targets. Runtime checks cover Python dependencies, Java,
Nextflow and actual headless R PDF/PNG output; cluster execution additionally
requires the SLURM commands. Each executing run records its runtime fingerprint
and Conda package records and checks these on resume.

## Verification completed

- **51 automated tests passed:** 33 workflow tests, including nine environment
  wrapper tests and a standalone frozen-source preparation test, plus 18 tests
  for the bundled, unmodified CAAStools. Logs and machine-readable status are in
  `review/tests/test_report.json` and the associated log files.
- **Shell syntax checks passed** for the setup, activation, SLURM submission and
  shared Conda helper scripts.
- **Linux dependency resolution passed**, using Conda 22.11.1 with strict channel
  priority, `CONDA_SUBDIR=linux-64` and `CONDA_OVERRIDE_GLIBC=2.17`. A `conda env
  create --dry-run --json` invocation against this `environment.yml` exited with
  status 0. Its package cache and proposed prefix were isolated in a temporary
  directory; no environment was installed. Among the resolved packages were
  Python 3.11.9, Nextflow 24.04.2, OpenJDK 17.0.11, Biopython 1.84, DendroPy 4.5.2,
  NumPy 1.26.4, SciPy 1.13.1, setuptools 80.10.2 and R 4.3.3. This is a dependency
  compatibility check, not a Linux execution test or a complete package lock.
- **Actual local synthetic workflow passed**, run ID
  `standalone-environment-check`: nine hypotheses times two alignments, giving
  18 unique completed CAAS tasks, assembly/filtering, summaries and PNG/PDF
  figures. The positive P0 fixture recovered four significant positions, one
  retained position and one query gene; the background correctly retained both
  the positive and zero-result genes. See
  `review/tests/standalone_environment_smoke.json`.
- **Exact-runtime resume passed**: all 18 CAAS tasks and all other workflow
  processes were reused from cache; no discovery task was repeated. The final
  smoke report records both the original and resumed trace hashes.
- **Portable preparation passed** in a temporary isolated copy with nonexistent
  original-project paths. Regenerated 19-replicate pilot pool/cycle files agreed
  byte-for-byte with the frozen snapshot.

The execution test used the already available local macOS runtime, whose exact
versions are recorded in `results/standalone-environment-check/environment.json`;
it did not install or execute the proposed Linux Conda environment. A changed
Python executable was correctly rejected by the resume identity guard.

## Remaining execution-machine checks

Create the dedicated environment on the cluster using
`bash create_conda_environment.sh`. This also validates the actual installed
runtime. Regenerate fixture paths there and execute the documented synthetic
smoke on a permitted node before production. Review storage visibility, SLURM
partition/account settings and the real alignment inventory.

No cluster job, biological discovery or live enrichment was launched for this
environment update. The default synthetic enrichment stage deliberately reports
`blocked_annotation_not_frozen`; production annotation/statistical approvals
remain necessary. Scientific assignments, frozen input tables, sampling rules
and historical results were not changed.
