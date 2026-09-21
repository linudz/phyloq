# Resource limits and safe resume

## Shared-tool staging update

The current version uses the shared absolute CAAStools directory as a value
input, never as a staged directory. Workers verify its frozen checksum and
disable Python bytecode writes. Alignment, job and cycle files still use copy
staging. No old work/results directories are removed automatically.

This changes main.nf and therefore deliberately does NOT qualify for the
resource-only bridge described below. Stop an old run before updating the
checkout and start a NEW ID, e.g. `caas-five-shared-01`. Subsequent resumes of
this new run use the same ID and can adjust resource flags. Keep the shared
tool path mounted and unchanged for the entire run.

## Historical resource-only transition (before shared-tool update)

The cluster selector explicitly sets CPU, memory and time. Nextflow 24.04.2
replaced the base selector when loading the cluster profile, previously losing
these directives. Defaults are now 1 CPU, 2 GB and 30 minutes per CAAS task.
Driver walltime remains 3 days. A task timeout still fails the workflow; it is
not silently treated as a negative result or skipped.

Stop the existing driver and check for remaining jobs belonging to that run
before updating. Do not remove `work/`, `.nextflow/`, launch inputs or results.
After updating, resume the SAME five-analysis run ID with the original options:

```bash
bash launch_all_caas.sh --run-id caas-five-01 --resume
```

Replace the example ID with the actual five-analysis run ID. This is NOT a
migration from the old 203-analysis selection. No new run ID is needed for this
resource fix when resuming the reduced version published as commit 287d4ab.

Resources can subsequently be changed without editing tracked configuration:

```bash
bash launch_all_caas.sh --run-id caas-five-01 --resume --task-time 30m --task-memory '2 GB' --task-cpus 1
```

Resource defaults are applied each launch: repeat any nondefault resource flags
when resuming. Every launch plan and params JSON records the requested resources.
Custom configuration files remain strictly checksummed; editing them is NOT
automatically accepted as a resource-only change.

The compatibility bridge accepts only the exact reviewed old/new implementation
fingerprints while requiring all other identity fields and the installed runtime
to match. It preserves the old task-script fingerprint for Nextflow cache reuse,
and records the actual new implementation separately in every launch plan.
It never overwrites the original launch lock. Analysis scripts, main.nf, pooled
settings and CAAStools are unchanged. Unknown code/input/environment changes fail.

Completed tasks with intact valid cache can be reused. Interrupted/failed tasks
are rerun. Already submitted SLURM jobs are not modified retroactively.
Verify new jobs with `squeue`/`scontrol`: expected time limit `00:30:00`, memory
2 GB and 1 CPU. The definitive submission directives are in `.command.run`.
