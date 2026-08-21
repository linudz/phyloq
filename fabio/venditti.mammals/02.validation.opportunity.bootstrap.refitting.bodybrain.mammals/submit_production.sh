#!/usr/bin/env bash
set -euo pipefail
cd "$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
mkdir -p logs
exec sbatch slurm_driver.sh params.production.yaml "$@"

