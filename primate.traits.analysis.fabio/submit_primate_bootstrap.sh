#!/bin/bash
#SBATCH --job-name="primate_traits_bootstrap"
#SBATCH --output=logs/bootstrap_driver_%j.out
#SBATCH --error=logs/bootstrap_driver_%j.err
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --time=7-00:00:00
#SBATCH --partition=std-cpu

# Uncomment and complete these directives if required by Correfoc.
##SBATCH --account=ACCOUNT_NAME
##SBATCH --qos=QOS_NAME

set -euo pipefail

module load Miniconda3
source /homes/aplic/noarch/software/Miniconda3/23.9.0-0/etc/profile.d/conda.sh
conda activate primate-pss-bootstrap

# Slurm executes a spooled copy of this script. Return to the directory from
# which sbatch was invoked, where main.nf and params.yaml are located.
cd "${SLURM_SUBMIT_DIR:?SLURM_SUBMIT_DIR is not set}"
mkdir -p .nextflow logs results

nextflow run main.nf \
  -profile slurm \
  -params-file params.yaml \
  -resume

exit 0
