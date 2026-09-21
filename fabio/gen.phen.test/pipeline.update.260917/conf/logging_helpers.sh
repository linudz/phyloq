# Source after resolving caas_root and validating caas_run_id, before Conda.
caas_start_logging() {
    caas_root="$(cd "$caas_root" && pwd)"
    caas_log_dir="$caas_root/logs/$caas_run_id"
    mkdir -p "$caas_log_dir"
    if [[ -n "${SLURM_JOB_ID:-}" ]]; then
        caas_log_file="$caas_log_dir/driver-${SLURM_JOB_ID}.log"
        echo "Driver log: $caas_log_file"
        exec >>"$caas_log_file" 2>&1
    else
        caas_log_file="$caas_log_dir/submission-$(date -u +%Y%m%dT%H%M%SZ)-$$.log"
        echo "Launch log: $caas_log_file"
        exec > >(tee -a "$caas_log_file") 2>&1
    fi
    echo "[$(date -u +%Y-%m-%dT%H:%M:%SZ)] Run: $caas_run_id; job: ${SLURM_JOB_ID:-not-submitted}"
}
