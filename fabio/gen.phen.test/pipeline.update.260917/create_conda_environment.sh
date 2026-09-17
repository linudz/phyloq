#!/usr/bin/env bash
# Create the dedicated environment once on a login/setup node, never per CAAS job.
set -euo pipefail
caas_root="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "$caas_root/conf/conda_helpers.sh"
caas_mode=create
while [[ $# -gt 0 ]]; do
    case "$1" in
        --name) export CAAS_CONDA_ENV="${2:?Missing environment name}"; unset CAAS_CONDA_PREFIX; shift 2 ;;
        --prefix) export CAAS_CONDA_PREFIX="${2:?Missing absolute prefix}"; shift 2 ;;
        --update|--dry-run|--check) caas_mode="${1#--}"; shift ;;
        -h|--help)
            echo 'Usage: bash create_conda_environment.sh [--name NAME | --prefix ABSOLUTE_PATH] [--dry-run | --check | --update]'
            echo 'Default: caas-validation-260917. Existing environments are checked, not silently updated.'
            exit 0 ;;
        *) echo "Unknown argument: $1" >&2; exit 2 ;;
    esac
done
caas_init_conda
if [[ -n "${CAAS_CONDA_PREFIX:-}" ]]; then
    [[ "$CAAS_CONDA_PREFIX" == /* && "$CAAS_CONDA_PREFIX" != / ]] || { echo 'Use an absolute, dedicated environment prefix.' >&2; exit 2; }
    caas_selector=(--prefix "$CAAS_CONDA_PREFIX")
else
    caas_name="${CAAS_CONDA_ENV:-caas-validation-260917}"
    [[ "$caas_name" != base && "$caas_name" != root && "$caas_name" != phyloq && "$caas_name" != */* && -n "$caas_name" ]] || { echo 'Use a dedicated environment name, not base/phyloq.' >&2; exit 2; }
    caas_selector=(--name "$caas_name")
fi
export CONDA_CHANNEL_PRIORITY=strict
caas_validate_conda_target
if [[ "$caas_mode" == dry-run ]]; then
    conda env create "${caas_selector[@]}" --file "$caas_root/environment.yml" --dry-run
    exit 0
fi
if [[ "$caas_mode" == check ]]; then
    caas_activate_conda
elif [[ "$caas_mode" == update ]]; then
    # Explicit opt-in only. Never use --force, remove an environment, or prune.
    conda env update "${caas_selector[@]}" --file "$caas_root/environment.yml"
    caas_activate_conda
elif caas_activate_conda >/dev/null 2>&1; then
    echo "Checking existing environment: $CONDA_PREFIX (no packages changed)."
else
    conda env create "${caas_selector[@]}" --file "$caas_root/environment.yml"
    caas_activate_conda
fi
caas_stamp="$(date -u +%Y%m%dT%H%M%SZ)"
"$CONDA_PREFIX/bin/python" "$caas_root/scripts/check_environment.py" --conda-prefix "$CONDA_PREFIX" \
    --output "$caas_root/review/environments/environment.$caas_stamp.json"
echo 'Environment ready. Test the bundle with: bash run_pipeline.sh prepare'
