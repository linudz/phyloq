#!/usr/bin/env bash

set -euo pipefail

script_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
conda_init_script="${GENPHEN_CONDA_SH:-/homes/aplic/noarch/software/Miniconda3/23.9.0-0/etc/profile.d/conda.sh}"
conda_environment="${GENPHEN_CONDA_ENV:-phyloq}"

if [[ ! -r "${conda_init_script}" ]]; then
    echo "Cannot read Conda initialization script: ${conda_init_script}" >&2
    exit 1
fi

source "${conda_init_script}"

if conda run -n "${conda_environment}" true >/dev/null 2>&1; then
    echo "Updating Conda environment: ${conda_environment}"
    # Do not prune: the shared phyloq environment may also contain RERconverge.
    conda env update \
        --name "${conda_environment}" \
        --file "${script_dir}/environment.yml"
else
    echo "Creating Conda environment: ${conda_environment}"
    conda env create \
        --name "${conda_environment}" \
        --file "${script_dir}/environment.yml"
fi

conda activate "${conda_environment}"

python3 -c 'import Bio, dendropy, numpy, scipy; print("CAAStools Python dependencies: OK")'
nextflow -version

echo "Conda environment ${conda_environment} is ready."
