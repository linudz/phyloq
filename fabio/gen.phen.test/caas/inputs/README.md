# Benchmark inputs

- `config.creation/` contains the source selections used to define the nine
  benchmark strategies.
- `benchmark-configs/` contains the nine generated 100-cycle, 4-vs-4 pooled
  CAAStools configurations.
- `benchmark-pools/` contains the complete fixed FG/BG pools required for
  pooled event reconstruction.
- `benchmark.configs.tsv` is the Nextflow manifest connecting each approach to
  its source, pooled configuration and complete pool.
- `benchmark.configs.pss-ranked-13x13.tsv` contains only the hierarchical
  13-vs-13 PSS strategy for a focused Nextflow launch.
- `benchmark.configs.pss-cercopithecidae-random-pools.tsv` contains only the
  cercopithecid PSS-defined random-pool strategy.
- `benchmark.configs.cercopithecidae-absolute-trait-tails.tsv` contains only
  the cercopithecid upper/lower 10% trait-tail control and should be used with
  `-resume` to avoid selecting the complete benchmark manifest.
- `alignments/` is populated directly on Correfoc and is excluded from Git.

Regenerate all derived inputs with:

```bash
python3 scripts/create_pss_benchmark_configs.py
```

See the project-level `README.md` for the biological rationale, selection
rules, exclusions, examples and execution procedure.
