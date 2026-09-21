#!/usr/bin/env bash
# Submit THIS driver explicitly with sbatch after reviewing a dry-run launch plan.
#SBATCH --job-name=CAAS_VALIDATION_DRIVER
#SBATCH --cpus-per-task=1
#SBATCH --partition=std-cpu
#SBATCH --mem=8G
#SBATCH --time=3-00:00:00
set -euo pipefail
# SLURM copies this script to its spool, so BASH_SOURCE alone is not the project
# location. Run sbatch from the bundle, or export CAAS_VALIDATION_ROOT explicitly.
caas_root="${CAAS_VALIDATION_ROOT:-${SLURM_SUBMIT_DIR:-$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)}}"
[[ -f "$caas_root/run_pipeline.sh" && -f "$caas_root/main.nf" ]] || {
    echo 'Run sbatch from the pipeline directory or set CAAS_VALIDATION_ROOT to its absolute path.' >&2
    exit 1
}
cd "$caas_root"
# Activates the independent environment on the driver; Nextflow also activates
# the same absolute prefix in worker jobs. Does not install on compute nodes.
exec bash "$caas_root/run_pipeline.sh" "$@" --profile cluster
