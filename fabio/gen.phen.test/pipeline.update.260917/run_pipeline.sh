#!/usr/bin/env bash
# Portable entry point: activate the dedicated environment for the whole driver.
set -euo pipefail
caas_root="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "$caas_root/conf/conda_helpers.sh"
caas_init_conda
caas_activate_conda
exec "$CONDA_PREFIX/bin/python" "$caas_root/run_validation.py" "$@" \
    --conda-prefix "$CONDA_PREFIX" --conda-init "$CAAS_CONDA_SH"
