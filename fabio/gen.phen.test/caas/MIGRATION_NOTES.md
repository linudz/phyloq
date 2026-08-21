# Migration notes

## Source

The workflow was copied on 2026-08-20 from:

`/Users/fabio/active.research/aging.update/agingprimates/05.gen.phen/lq.table2.bootstrap.nextflow`

The bundled `bin/caastools` is the modified pooled-discovery implementation,
including its tests and event-pooling logic.

## Copied

- Nextflow workflow and configuration;
- local pooled CAAStools implementation and tests;
- pooled-hypothesis preparation scripts;
- Conda environment definition and creation helper;
- local and SLURM launch scripts;
- cluster configuration and original workflow documentation.

## Intentionally not copied

- input configurations and resampling tables, because they will be generated
  by step 4 of the PSS case study;
- previous `.nextflow` state and cache;
- Nextflow logs;
- previous cluster results and archives;
- Python bytecode caches.

## Adaptation still required

The copied files still contain longevity-specific defaults, filenames and
paths. Before execution, adapt the workflow to the PSS-informed pools and the
alignment location that will be used in `phyloq`, then validate locally before
transfer to the cluster.
