#!/usr/bin/env bash
# Source from this bundle; no reference to the old workflow or phyloq required.
caas_init_conda() {
    if [[ -n "${CAAS_CONDA_SH:-}" ]]; then
        [[ -r "$CAAS_CONDA_SH" ]] || { echo "Cannot read CAAS_CONDA_SH: $CAAS_CONDA_SH" >&2; return 1; }
    elif command -v conda >/dev/null 2>&1; then
        local caas_conda_base
        caas_conda_base="$(conda info --base)" || return 1
        CAAS_CONDA_SH="$caas_conda_base/etc/profile.d/conda.sh"
    elif [[ -r /homes/aplic/noarch/software/Miniconda3/23.9.0-0/etc/profile.d/conda.sh ]]; then
        # Known Correfoc installation, used only if it actually exists.
        CAAS_CONDA_SH=/homes/aplic/noarch/software/Miniconda3/23.9.0-0/etc/profile.d/conda.sh
    else
        echo 'Conda not found. Load your cluster Conda module or set CAAS_CONDA_SH to its conda.sh.' >&2
        return 1
    fi
    [[ -r "$CAAS_CONDA_SH" ]] || { echo "Cannot read Conda initialization: $CAAS_CONDA_SH" >&2; return 1; }
    export CAAS_CONDA_SH
    # Some older Conda hooks reference unset variables; restore caller options.
    local caas_had_nounset=0
    [[ $- == *u* ]] && caas_had_nounset=1
    set +u
    source "$CAAS_CONDA_SH"
    local caas_status=$?
    [[ $caas_had_nounset == 1 ]] && set -u
    return "$caas_status"
}

caas_validate_conda_target() {
    local caas_target="${CAAS_CONDA_PREFIX:-${CAAS_CONDA_ENV:-caas-validation-260917}}"
    if [[ "$caas_target" == base || "$caas_target" == root || "${caas_target##*/}" == phyloq ]]; then
        echo 'Use a dedicated CAAS validation environment, not base or the shared phyloq environment.' >&2
        return 1
    fi
    if [[ -n "${CAAS_CONDA_PREFIX:-}" && -d "$CAAS_CONDA_PREFIX" ]]; then
        local caas_base_path caas_target_real caas_base_real
        caas_base_path="$(conda info --base)" || return 1
        caas_base_real="$(cd "$caas_base_path" && pwd -P)" || return 1
        caas_target_real="$(cd "$CAAS_CONDA_PREFIX" && pwd -P)" || return 1
        [[ "$caas_base_real" != "$caas_target_real" ]] || { echo 'Refusing to use/modify the base Conda prefix.' >&2; return 1; }
    fi
}

caas_activate_conda() {
    caas_validate_conda_target || return 1
    local caas_target="${CAAS_CONDA_PREFIX:-${CAAS_CONDA_ENV:-caas-validation-260917}}"
    local caas_had_nounset=0
    [[ $- == *u* ]] && caas_had_nounset=1
    set +u
    if conda activate "$caas_target"; then
        [[ $caas_had_nounset == 1 ]] && set -u
        return 0
    else
        [[ $caas_had_nounset == 1 ]] && set -u
        echo "Environment unavailable: $caas_target. First run: bash create_conda_environment.sh" >&2
        return 1
    fi
}
