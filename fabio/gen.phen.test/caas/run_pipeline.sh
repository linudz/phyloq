#!/usr/bin/env bash

set -euo pipefail

script_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
state_file="${script_dir}/.last_run_id"
config_file="${GENPHEN_NEXTFLOW_CONFIG:-${script_dir}/conf/cluster.config}"
resume=false

for argument in "$@"; do
    if [[ "${argument}" == "-resume" ]]; then
        resume=true
        break
    fi
done

if [[ ! -f "${config_file}" ]]; then
    echo "Nextflow configuration not found: ${config_file}" >&2
    exit 1
fi

if [[ "${resume}" == true ]]; then
    if [[ ! -s "${state_file}" ]]; then
        echo "Cannot resume: ${state_file} does not contain a previous run ID." >&2
        exit 1
    fi
    run_id="$(<"${state_file}")"
else
    run_id="run-$(date +%y%m%d%H%M%S)"
    printf '%s\n' "${run_id}" > "${state_file}"
fi

mkdir -p "${script_dir}/logs"

echo "GenPhen CAAStools pooled-discovery run: ${run_id}"
echo "Nextflow config: ${config_file}"

cd "${script_dir}"
exec nextflow -c "${config_file}" run main.nf \
    --run_id "${run_id}" \
    "$@"
