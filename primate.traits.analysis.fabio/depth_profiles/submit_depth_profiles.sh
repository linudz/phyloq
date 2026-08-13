#!/bin/bash
#SBATCH --job-name="primate_depth_profiles"
#SBATCH --output=logs/depth_profiles_driver_%j.out
#SBATCH --error=logs/depth_profiles_driver_%j.err
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --time=7-00:00:00
#SBATCH --partition=std-cpu

set -euo pipefail

# Initialize Conda directly: Correfoc's Lmod hierarchy is not stable.
CORREFOC_CONDA_SH=/homes/aplic/noarch/software/Miniconda3/23.9.0-0/etc/profile.d/conda.sh
if [[ ! -r "${CORREFOC_CONDA_SH}" ]]; then
    echo "Cannot read Conda initialization script: ${CORREFOC_CONDA_SH}" >&2
    exit 1
fi
source "${CORREFOC_CONDA_SH}"
conda activate primate-pss-bootstrap

# Submit from repository root so Slurm can open logs/ before the script starts.
if [[ ! -f "${SLURM_SUBMIT_DIR}/depth_profiles/main.nf" ]]; then
  echo "Submit this wrapper from primate.traits.analysis.fabio/." >&2
  exit 1
fi
pipeline_dir="${SLURM_SUBMIT_DIR}/depth_profiles"

repository_dir="$(dirname "${pipeline_dir}")"
mkdir -p \
  "${repository_dir}/logs" \
  "${pipeline_dir}/results" \
  "${pipeline_dir}/work"

cd "${pipeline_dir}"

nextflow run main.nf \
  -c nextflow.config \
  -profile slurm \
  -params-file params.yaml \
  -resume

exit 0
