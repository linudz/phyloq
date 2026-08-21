#!/usr/bin/env bash

#SBATCH --job-name=genphen-pooled-driver
#SBATCH --cpus-per-task=1
#SBATCH --mem=16G
#SBATCH --time=7-00:00:00
#SBATCH --partition=std-cpu
#SBATCH --output=logs/genphen-pooled-driver-%j.out
#SBATCH --error=logs/genphen-pooled-driver-%j.err

# Uncomment and edit the directives required by the cluster.
##SBATCH --account=YOUR_ACCOUNT
##SBATCH --qos=YOUR_QOS
##SBATCH --constraint=YOUR_CONSTRAINT

set -euo pipefail

if [[ -z "${SLURM_JOB_ID:-}" ]]; then
    echo "This is a SLURM batch script. Submit it with:" >&2
    echo "  sbatch submit_pipeline_slurm.sh" >&2
    exit 1
fi

# Slurm runs a spooled copy of this file. Return to the directory from which
# sbatch was invoked, where the complete pipeline checkout must be located.
pipeline_dir="${SLURM_SUBMIT_DIR:?SLURM_SUBMIT_DIR is not set}"
if [[ ! -f "${pipeline_dir}/run_pipeline.sh" || ! -f "${pipeline_dir}/main.nf" ]]; then
    echo "Submit this script from the lq.table2.bootstrap.nextflow directory." >&2
    echo "SLURM submission directory: ${pipeline_dir}" >&2
    exit 1
fi
cd "${pipeline_dir}"

export GENPHEN_NEXTFLOW_CONFIG="${GENPHEN_NEXTFLOW_CONFIG:-${pipeline_dir}/conf/cluster.config}"

# Initialize the Correfoc Miniconda installation and activate the environment
# shared by the Nextflow driver and all submitted pooled-discovery tasks.
conda_init_script="${GENPHEN_CONDA_SH:-/homes/aplic/noarch/software/Miniconda3/23.9.0-0/etc/profile.d/conda.sh}"
conda_environment="${GENPHEN_CONDA_ENV:-phyloq}"
if [[ ! -r "${conda_init_script}" ]]; then
    echo "Cannot read Conda initialization script: ${conda_init_script}" >&2
    exit 1
fi
source "${conda_init_script}"
conda activate "${conda_environment}"

# Keep the Nextflow JVM comfortably below the driver's SLURM memory limit.
export NXF_OPTS="${NXF_OPTS:--Xms512m -Xmx8g}"

echo "GenPhen CAAStools pooled-discovery driver"
echo "SLURM job: ${SLURM_JOB_ID}"
echo "Host: $(hostname)"
echo "Started: $(date --iso-8601=seconds)"
echo "Nextflow config: ${GENPHEN_NEXTFLOW_CONFIG}"
echo "Conda environment: ${CONDA_DEFAULT_ENV:-${conda_environment}}"
echo "Nextflow JVM options: ${NXF_OPTS}"

exec bash run_pipeline.sh "$@"
