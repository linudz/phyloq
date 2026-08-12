#!/usr/bin/env bash

set -euo pipefail

script_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
log_dir="${script_dir}/slurm.logs"
driver="${script_dir}/slurm_driver.sh"

if ! command -v sbatch >/dev/null 2>&1; then
    echo "Cannot find sbatch. Load SLURM or run this command on a SLURM login node." >&2
    exit 1
fi

mkdir -p "${log_dir}"
cd "${script_dir}"

sbatch_options=(--parsable)

if [[ -n "${PHYLOQ_SLURM_ACCOUNT:-}" ]]; then
    sbatch_options+=("--account=${PHYLOQ_SLURM_ACCOUNT}")
fi

if [[ -n "${PHYLOQ_SLURM_PARTITION:-}" ]]; then
    sbatch_options+=("--partition=${PHYLOQ_SLURM_PARTITION}")
fi

if [[ -n "${PHYLOQ_SLURM_QOS:-}" ]]; then
    sbatch_options+=("--qos=${PHYLOQ_SLURM_QOS}")
fi

job_id="$(sbatch "${sbatch_options[@]}" "${driver}" "$@")"

echo "Submitted PhyloQ SLURM driver job: ${job_id}"
echo "Monitor it with: squeue -j ${job_id}"
echo "Driver logs: ${log_dir}/phyloq-driver-${job_id}.out"
echo "Errors: ${log_dir}/phyloq-driver-${job_id}.err"
