nextflow.enable.dsl = 2

process CAASTOOLS_POOLED_BENCHMARK {
    tag "${approach}:${gene_id}"

    stageInMode 'copy'
    publishDir path: { "${params.results_root}/${params.run_id}/${approach}/caas-pooled" },
               pattern: "*.pooled.caas.tsv", mode: 'copy', overwrite: true
    publishDir path: { "${params.results_root}/${params.run_id}/${approach}/caas-pooled-events" },
               pattern: "*.pooled.caas.events.tsv", mode: 'copy', overwrite: true

    input:
    tuple val(approach), val(gene_id), path(alignment), path(pool_config), path(pooled_hypotheses)
    path caastools_dir

    output:
    tuple val(approach), val(gene_id), path("${gene_id}.pooled.caas.tsv"), emit: legacy_results
    tuple val(approach), val(gene_id), path("${gene_id}.pooled.caas.events.tsv"), emit: event_results

    script:
    """
    "${params.python_command}" "${caastools_dir}/ct" pooled-discovery \
        -a "${alignment}" \
        -t "${pool_config}" \
        -s "${pooled_hypotheses}" \
        -o "${gene_id}.pooled.caas.tsv" \
        --event-output "${gene_id}.pooled.caas.events.tsv" \
        --hypotheses-output none \
        --fmt "${params.caas_alignment_format}" \
        --filter_significant "${params.caas_filter_significant}" \
        --patterns "${params.caas_patterns}" \
        --max_bg_gaps "${params.caas_max_bg_gaps}" \
        --max_fg_gaps "${params.caas_max_fg_gaps}" \
        --max_gaps "${params.caas_max_gaps}" \
        --max_gaps_per_position "${params.caas_max_gaps_per_position}" \
        --max_bg_miss "${params.caas_max_bg_miss}" \
        --max_fg_miss "${params.caas_max_fg_miss}" \
        --max_miss "${params.caas_max_miss}" \
        --min_fg_observed "${params.caas_min_fg_observed}" \
        --min_bg_observed "${params.caas_min_bg_observed}"
    """
}

workflow {
    if (!params.run_id) {
        error "Missing --run_id. Use run_pipeline.sh or submit_pipeline_slurm.sh."
    }
    if (!params.alignments) {
        error "params.alignments is empty. Edit conf/cluster.config."
    }
    if (!params.benchmark_manifest) {
        error "params.benchmark_manifest is empty. Edit conf/cluster.config."
    }

    caastools_dir = file(params.caastools_dir, checkIfExists: true)

    benchmark_ch = Channel
        .fromPath(params.benchmark_manifest, checkIfExists: true)
        .splitCsv(header: true, sep: '\t')
        .map { row ->
            tuple(
                row.approach,
                file("${projectDir}/${row.pool_config}", checkIfExists: true),
                file("${projectDir}/${row.hypotheses_config}", checkIfExists: true)
            )
        }

    alignments_ch = Channel
        .fromPath(params.alignments, checkIfExists: true)
        .filter { alignment -> !alignment.isDirectory() }
        .map { alignment ->
            tuple(alignment.name.replaceFirst(/\.[^.]+$/, ''), alignment)
        }

    benchmark_jobs_ch = benchmark_ch
        .combine(alignments_ch)
        .map { approach, pool_config, pooled_hypotheses, gene_id, alignment ->
            tuple(approach, gene_id, alignment, pool_config, pooled_hypotheses)
        }

    CAASTOOLS_POOLED_BENCHMARK(benchmark_jobs_ch, caastools_dir)
}
