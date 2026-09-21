#!/usr/bin/env bash
# One driver, one Nextflow run, all explicitly selected CAAS modes.
#SBATCH --job-name=CAAS_VALIDATION_ALL
#SBATCH --cpus-per-task=1
#SBATCH --partition=std-cpu
#SBATCH --mem=8G
#SBATCH --time=3-00:00:00
#SBATCH --output=logs/slurm-%j.out
#SBATCH --error=logs/slurm-%j.err
#SBATCH --open-mode=append
set -euo pipefail
caas_script_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
caas_root="${CAAS_VALIDATION_ROOT:-$caas_script_dir}"
if [[ ! -f "$caas_root/main.nf" && -n "${SLURM_SUBMIT_DIR:-}" ]]; then
    caas_root="$SLURM_SUBMIT_DIR"
fi
caas_original_args=("$@")
caas_run_id=''
caas_modes=all
caas_source=()
caas_extras=()
caas_launch_extras=()
caas_resume=()
caas_plan_only=0
while [[ $# -gt 0 ]]; do
    case "$1" in
        --run-id) caas_run_id="${2:?Missing run ID}"; shift 2 ;;
        --modes) caas_modes="${2:?Missing modes}"; shift 2 ;;
        --alignments-dir)
            [[ ${#caas_source[@]} == 0 ]] || { echo 'Use only one alignment source.' >&2; exit 2; }
            [[ -d "${2:?Missing alignment directory}" ]] || { echo "Alignment directory not found: $2" >&2; exit 2; }
            caas_dir="$(cd "$2" && pwd)"
            caas_source=(--alignments-pattern "$caas_dir/*"); shift 2 ;;
        --alignments-pattern|--inventory)
            [[ ${#caas_source[@]} == 0 ]] || { echo 'Use only one alignment source.' >&2; exit 2; }
            caas_source=("$1" "${2:?Missing alignment source}"); shift 2 ;;
        --include-references|--include-p2) caas_extras+=("$1"); shift ;;
        --resume) caas_resume=(--resume); shift ;;
        --plan-only) caas_plan_only=1; shift ;;
        --cluster-config) caas_launch_extras+=(--cluster-config "${2:?Missing config path}"); shift 2 ;;
        -h|--help)
            echo 'Usage: bash launch_all_caas.sh --run-id NAME [--alignments-dir DIR] [options]'
            echo 'One SLURM driver runs N3/N4/N5 + all 99 R0 + all 99 R1 + P1/N6 (203 hypotheses).'
            echo 'No pilot or approval JSON. Uses the bundled modified pooled CAAStools.'
            echo 'Sources: --alignments-dir DIR, --alignments-pattern QUOTED_GLOB, or --inventory TSV'
            echo 'Default source in phyloq: ../caas/inputs/alignments/*.phy (the existing pipeline inputs).'
            echo 'In the research project: ../pipeline/inputs/alignments/*.phy; standalone: inputs/alignments/*.phy.'
            echo 'Options: --modes all|deterministic,r0,r1,paired --include-references --include-p2'
            echo '         --resume --plan-only --cluster-config FILE'
            echo 'Historical references/P2 are NOT rerun unless explicitly included.'
            echo 'Outside SLURM, bash submits with sbatch. Direct sbatch use is also supported.'
            exit 0 ;;
        *) echo "Unknown argument: $1" >&2; exit 2 ;;
    esac
done
[[ -f "$caas_root/main.nf" ]] || { echo 'Set CAAS_VALIDATION_ROOT to the pipeline directory.' >&2; exit 2; }
[[ "$caas_run_id" =~ ^[A-Za-z0-9][A-Za-z0-9_-]*$ ]] || { echo 'Supply --run-id using letters, digits, underscores or hyphens.' >&2; exit 2; }
source "$caas_root/conf/logging_helpers.sh"
caas_start_logging
if [[ ${#caas_source[@]} == 0 ]]; then
    # Preserve the location and *.phy glob from the previous CAAS cluster.config.
    # The new bundle is a sibling, not the owner of another alignment collection.
    if [[ -f "$caas_root/../caas/conf/cluster.config" ]]; then
        caas_alignment_dir="$caas_root/../caas/inputs/alignments"
    elif [[ -f "$caas_root/../pipeline/conf/cluster.config" ]]; then
        caas_alignment_dir="$caas_root/../pipeline/inputs/alignments"
    else
        caas_alignment_dir="$caas_root/inputs/alignments"
    fi
    [[ -d "$caas_alignment_dir" ]] || {
        echo "Expected alignment directory: $caas_alignment_dir" >&2
        echo 'Keep the old CAAS inputs in place, or supply an explicit --alignments-dir / --alignments-pattern / --inventory.' >&2
        exit 2
    }
    caas_alignment_dir="$(cd "$caas_alignment_dir" && pwd)"
    caas_pattern="$caas_alignment_dir/*.phy"
    compgen -G "$caas_pattern" >/dev/null || { echo "No alignments match the old pipeline pattern: $caas_pattern" >&2; exit 2; }
    caas_source=(--alignments-pattern "$caas_pattern")
    echo "Using existing CAAS alignments: $caas_pattern"
fi
export CAAS_VALIDATION_ROOT="$caas_root"
# Resolve any relative source/config paths against the caller's working directory.
# Preserve that location when SLURM runs its spooled script.
export CAAS_LAUNCH_CALLER_DIR="${CAAS_LAUNCH_CALLER_DIR:-$PWD}"
if [[ -z "${SLURM_JOB_ID:-}" && "$caas_plan_only" == 0 ]]; then
    command -v sbatch >/dev/null || { echo 'sbatch not found: run this command on the cluster.' >&2; exit 1; }
    echo "Submitting one CAAS driver; selected modes: $caas_modes. No pilot."
    exec sbatch --export=ALL --output="$caas_log_dir/slurm-%j.out" --error="$caas_log_dir/slurm-%j.err" --open-mode=append "$caas_root/launch_all_caas.sh" "${caas_original_args[@]}"
fi
cd "$CAAS_LAUNCH_CALLER_DIR"
source "$caas_root/conf/conda_helpers.sh"
caas_init_conda
caas_activate_conda
export PYTHONDONTWRITEBYTECODE=1 PYTHONNOUSERSITE=1
caas_input_dir="$caas_root/inputs/launches/$caas_run_id"
caas_prepare_command=("$CONDA_PREFIX/bin/python" "$caas_root/scripts/prepare_cluster_launch.py" \
    --run-id "$caas_run_id" --modes "$caas_modes" --output "$caas_input_dir" "${caas_source[@]}")
# Conditional appends also work with Bash 3/4.3 and nounset + empty arrays.
if [[ ${#caas_extras[@]} -gt 0 ]]; then caas_prepare_command+=("${caas_extras[@]}"); fi
if [[ ${#caas_resume[@]} -gt 0 ]]; then caas_prepare_command+=("${caas_resume[@]}"); fi
"${caas_prepare_command[@]}"
caas_run_command=("$CONDA_PREFIX/bin/python" "$caas_root/run_validation.py" benchmark \
    --run-id "$caas_run_id" --selection "$caas_input_dir/selection.tsv" \
    --alignments "$caas_input_dir/alignment_manifest.tsv" --profile cluster \
    --conda-prefix "$CONDA_PREFIX" --conda-init "$CAAS_CONDA_SH")
if [[ ${#caas_launch_extras[@]} -gt 0 ]]; then caas_run_command+=("${caas_launch_extras[@]}"); fi
if [[ ${#caas_resume[@]} -gt 0 ]]; then caas_run_command+=("${caas_resume[@]}"); fi
if [[ "$caas_plan_only" == 0 ]]; then caas_run_command+=(--execute); fi
exec "${caas_run_command[@]}"
