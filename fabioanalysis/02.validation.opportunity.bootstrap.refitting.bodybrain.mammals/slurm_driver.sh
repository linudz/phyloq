#!/usr/bin/env bash
#SBATCH --job-name=mammal_pss_validation
#SBATCH --output=logs/mammal_pss_validation_driver_%j.out
#SBATCH --error=logs/mammal_pss_validation_driver_%j.err
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --time=7-00:00:00
#SBATCH --partition=std-cpu
set -euo pipefail
PARAMS_FILE="${1:?Pass pilot or production params file}"
shift
CASE_DIR="${SLURM_SUBMIT_DIR}"
CONDA_SH=/homes/aplic/noarch/software/Miniconda3/23.9.0-0/etc/profile.d/conda.sh
[[ -r "$CONDA_SH" ]] && source "$CONDA_SH" || source "$(conda info --base)/etc/profile.d/conda.sh"
conda activate phyloq
cd "$CASE_DIR"
nextflow run main.nf -c nextflow.config -profile slurm \
  -params-file "$PARAMS_FILE" -resume "$@"
