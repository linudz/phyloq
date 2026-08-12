#!/usr/bin/env bash

#SBATCH --job-name=phyloq-driver
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --time=7-00:00:00
#SBATCH --output=slurm.logs/%x-%j.out
#SBATCH --error=slurm.logs/%x-%j.err

set -euo pipefail

script_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
cd "${script_dir}"

echo "SLURM driver job: ${SLURM_JOB_ID:-unknown}"
echo "Host: $(hostname)"
echo "Started: $(date --iso-8601=seconds)"

export PHYLOQ_NEXTFLOW_CONFIG="${script_dir}/slurm.config"
exec bash run_pipeline.sh "$@"
