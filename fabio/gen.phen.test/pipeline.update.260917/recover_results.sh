#!/usr/bin/env bash
# One recovery driver, zero CAAS jobs. Run from the bundle directory.
#SBATCH --job-name=CAAS_RECOVER
#SBATCH --partition=std-cpu
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --time=1-00:00:00
#SBATCH --output=logs/recovery-%j.out
#SBATCH --error=logs/recovery-%j.err
set -euo pipefail
recovery_root="${CAAS_RECOVERY_ROOT:-${SLURM_SUBMIT_DIR:-$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)}}"
recovery_original=("$@")
recovery_run=''; recovery_id=''; recovery_log=''; recovery_stopped=0
while [[ $# -gt 0 ]]; do
    case "$1" in
        --run-id) recovery_run="${2:?Missing run ID}"; shift 2 ;;
        --recovery-id) recovery_id="${2:?Missing recovery ID}"; shift 2 ;;
        --log) recovery_log="${2:?Missing source log}"; shift 2 ;;
        --confirm-stopped) recovery_stopped=1; shift ;;
        *) echo 'Usage: bash recover_results.sh --run-id NAME --recovery-id NAME --log FILE --confirm-stopped' >&2; exit 2 ;;
    esac
done
[[ "$recovery_run" =~ ^[A-Za-z0-9][A-Za-z0-9_-]*$ && "$recovery_id" =~ ^[A-Za-z0-9][A-Za-z0-9_-]*$ ]] || { echo 'Valid run and recovery IDs required.' >&2; exit 2; }
[[ "$recovery_stopped" == 1 ]] || { echo 'Stop the original driver and its tasks, then pass --confirm-stopped.' >&2; exit 2; }
cd "$recovery_root"
[[ -f "$recovery_log" ]] || { echo "Source log missing: $recovery_log" >&2; exit 2; }
export CAAS_RECOVERY_ROOT="$PWD"
recovery_logs="$PWD/logs/$recovery_run/recovery-$recovery_id"
mkdir -p "$recovery_logs"
if [[ -z "${SLURM_JOB_ID:-}" ]]; then
    exec sbatch --export=ALL --output="$recovery_logs/driver-%j.log" --error="$recovery_logs/driver-%j.log" \
        "$PWD/recover_results.sh" "${recovery_original[@]}"
fi
echo "Recovering $recovery_run into recovery/$recovery_id; no CAAS or Nextflow execution."
source "$PWD/conf/conda_helpers.sh"
caas_init_conda
caas_activate_conda
export PYTHONDONTWRITEBYTECODE=1 PYTHONNOUSERSITE=1
exec "$CONDA_PREFIX/bin/python" -u "$PWD/recovery/recover_completed.py" \
    --run-id "$recovery_run" --recovery-id "$recovery_id" --log "$recovery_log" \
    --rscript "$CONDA_PREFIX/bin/Rscript"
