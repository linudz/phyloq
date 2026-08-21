# Benchmark inputs

- `config.creation/` contains the three source selections used to define the
  benchmark strategies.
- `benchmark-configs/` contains the three generated 100-cycle, 4-vs-4 pooled
  CAAStools configurations.
- `benchmark-pools/` contains the complete fixed FG/BG pools required for
  pooled event reconstruction.
- `benchmark.configs.tsv` is the Nextflow manifest connecting each approach to
  its source, pooled configuration and complete pool.

Regenerate all derived inputs with:

```bash
python3 scripts/create_pss_benchmark_configs.py
```

See the project-level `README.md` for the biological rationale, selection
rules, exclusions, examples and execution procedure.
