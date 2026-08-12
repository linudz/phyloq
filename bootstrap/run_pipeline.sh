#!/usr/bin/env bash

set -euo pipefail

script_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
state_file="${script_dir}/.last_run_id"
resume=false

for argument in "$@"; do
    if [[ "${argument}" == "-resume" ]]; then
        resume=true
        break
    fi
done

if [[ "${resume}" == true ]]; then
    if [[ ! -s "${state_file}" ]]; then
        echo "Cannot resume: ${state_file} does not contain a previous run ID." >&2
        exit 1
    fi
    run_id="$(<"${state_file}")"
else
    run_id="run-$(date +%y%m%d%H%M)"
    printf '%s\n' "${run_id}" > "${state_file}"
fi

results_dir="${script_dir}/../results/${run_id}"
if [[ "${resume}" == false && -e "${results_dir}" ]]; then
    echo "Cannot create a new run: ${results_dir} already exists." >&2
    echo "Wait one minute or use -resume for the previous run." >&2
    exit 1
fi
mkdir -p "${results_dir}"

echo "PhyloQ run: ${run_id}"
echo "Results: ${results_dir}"

cd "${script_dir}"
exec nextflow run launch_pipeline.nf --run_id "${run_id}" "$@"
