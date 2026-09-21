# Cluster quickstart

Follow [RUN_ON_CLUSTER.md](RUN_ON_CLUSTER.md).

```bash
bash launch_all_caas.sh --run-id caas-five-01
# Later, after interruption of this reduced run:
bash launch_all_caas.sh --run-id caas-five-01 --resume
```

Runs N3/N4/N5/P1/N6 only. Stop the old 203-analysis run first; do not reuse its
run ID. Conda configuration and alignment locations are unchanged.
