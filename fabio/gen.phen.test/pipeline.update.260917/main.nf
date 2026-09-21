nextflow.enable.dsl = 2

process PREPARE_VALIDATION {
    tag "validate frozen species assignments"
    executor 'local'
    stageInMode 'copy'
    input:
    path design
    output:
    path 'design_check.json'
    script:
    """
    # implementation ${params.code_fingerprint}
    "${params.python_command}" "${projectDir}/scripts/check_design.py" --design "${design}" --output design_check.json
    """
}

process VALIDATE_INPUTS {
    tag "bind selected hypotheses to alignment inventory"
    executor 'local'
    stageInMode 'copy'
    publishDir "${params.results_root}/${params.run_id}/provenance", mode: 'copy', pattern: 'validation/*.json',
               saveAs: { filename -> filename.toString().tokenize('/').last() }
    input:
    path design
    path selection
    path alignments
    path settings
    path ready
    output:
    path 'validation/jobs.tsv', emit: jobs
    path 'validation/runtime.lock.json', emit: lock
    path 'validation/input_validation.json', emit: report
    script:
    def approvalArg = params.approvals ? "--approvals '${params.approvals}'" : ''
    def smokeArg = params.stage == 'smoke' || params.synthetic_inputs ? '--smoke' : ''
    def directArg = params.direct_discovery ? '--direct-discovery' : ''
    """
    # implementation ${params.code_fingerprint}
    "${params.python_command}" "${projectDir}/scripts/validate_inputs.py" \
      --design "${design}" --selection "${selection}" --alignments "${alignments}" \
      --settings "${settings}" --run-id "${params.run_id}" --output validation ${approvalArg} ${smokeArg} ${directArg}
    """
}

process CAAS_POOLED {
    // Failed discovery tasks are missing observations, never negative results.
    errorStrategy { params.direct_discovery ? 'ignore' : 'terminate' }
    tag "${strategy}/${replicate}:${gene}"
    stageInMode 'copy'
    publishDir path: { "${params.results_root}/${params.run_id}/${strategy}/${replicate}/raw" }, mode: 'copy', overwrite: false
    input:
    tuple val(hypothesis), val(strategy), val(replicate), val(gene), path(job), path(alignment), path(pool), path(cycles)
    // Shared read-only installation: a value, not a staged directory input.
    // run_caas.py verifies its frozen checksum before every invocation.
    val caastools
    output:
    tuple val(hypothesis), val(strategy), val(replicate), path("${gene}")
    script:
    """
    # implementation ${params.code_fingerprint}
    "${params.python_command}" "${projectDir}/scripts/run_caas.py" \
      --job "${job}" --alignment "${alignment}" --pool "${pool}" --config "${cycles}" \
      --tool-dir "${caastools}" --output "${gene}"
    """
}

process ASSEMBLE_FILTER_AND_QUERY {
    tag "${hypothesis}: verify completeness and filter"
    executor 'local'
    stageInMode 'copy'
    publishDir path: { "${params.results_root}/${params.run_id}/${strategy}/${replicate}" }, mode: 'copy', overwrite: true
    input:
    tuple val(hypothesis), val(strategy), val(replicate), path(artifacts, stageAs: 'raw/*')
    path alignments
    output:
    tuple val(hypothesis), path("${hypothesis}")
    script:
    def partialArg = params.direct_discovery ? '--allow-incomplete' : ''
    """
    # implementation ${params.code_fingerprint}
    "${params.python_command}" "${projectDir}/scripts/summarize_hypothesis.py" \
      --hypothesis "${hypothesis}" --alignments "${alignments}" --artifacts raw --output "${hypothesis}" ${partialArg}
    """
}

process COMMON_BACKGROUND {
    tag "common background for the declared cohort"
    executor 'local'
    stageInMode 'copy'
    publishDir "${params.summaries_root}/${params.run_id}", mode: 'copy', overwrite: true
    input:
    path cohort
    path summaries, stageAs: 'cohort_summaries/*'
    output:
    path 'matched'
    script:
    """
    # implementation ${params.code_fingerprint}
    "${params.python_command}" "${projectDir}/scripts/common_background.py" --cohort "${cohort}" --summary-base cohort_summaries --output matched ${params.direct_discovery ? '--allow-incomplete' : ''}
    """
}

process FUNCTIONAL_CONTROLS_AND_COMPARE {
    tag "frozen enrichment and conditional comparisons"
    executor 'local'
    stageInMode 'copy'
    publishDir "${params.summaries_root}/${params.run_id}", mode: 'copy', overwrite: true
    input:
    path matched
    path settings
    path annotation
    output:
    path 'functional'
    script:
    def cacheArg = params.annotation_backend == 'gprofiler-cache' ? "--cache '${annotation}'" : ''
    def annotationArg = params.annotation_backend == 'local-hypergeom-bonferroni' ? "--snapshot '${annotation}'" : ''
    def directArg = params.direct_discovery ? '--discovery-only' : ''
    """
    # implementation ${params.code_fingerprint}
    "${params.python_command}" "${projectDir}/scripts/functional.py" \
      --matched "${matched}" --settings "${settings}" --output functional ${cacheArg} ${annotationArg} ${directArg}
    """
}

process PLOT_SUMMARIES {
    executor 'local'
    stageInMode 'copy'
    publishDir "${params.summaries_root}/${params.run_id}", mode: 'copy', overwrite: true
    input:
    path matched
    path functional
    output:
    path 'figures'
    script:
    """
    # implementation ${params.code_fingerprint}
    "${params.r_command}" "${projectDir}/scripts/plot_validation.R" matched/strategy_comparison.tsv figures functional
    cp "${projectDir}/scripts/plot_validation.R" figures/
    cp "${projectDir}/scripts/figures.md" figures/README.md
    """
}

workflow {
    def stages = ['prepare', 'smoke', 'benchmark', 'deterministic', 'r0-pilot', 'r0-complete', 'r1-pilot', 'r1-complete', 'paired', 'summarize', 'reference-rerun', 'auxiliary-p2']
    if (!(params.stage in stages)) error "Choose explicit --stage: ${stages.join(', ')}"
    if (!params.code_fingerprint) error 'Use run_validation.py, which binds implementation hashes and run identity.'
    if (params.stage == 'prepare') {
        PREPARE_VALIDATION(file(params.design, checkIfExists: true))
    } else if (params.stage == 'summarize') {
        if (!params.cohort_manifest || !params.run_id) error 'summarize requires a cohort manifest and run ID'
        cohortFile = file(params.cohort_manifest, checkIfExists: true)
        summaries = Channel.fromPath(params.cohort_manifest, checkIfExists: true).splitCsv(header: true, sep: '\t').map { row ->
            def p = java.nio.file.Paths.get(row.summary_dir)
            file(p.isAbsolute() ? p : cohortFile.parent.resolve(p), checkIfExists: true)
        }.collect()
        COMMON_BACKGROUND(cohortFile, summaries)
        FUNCTIONAL_CONTROLS_AND_COMPARE(COMMON_BACKGROUND.out, file(params.settings, checkIfExists: true), file(params.annotation_input, checkIfExists: true))
        PLOT_SUMMARIES(COMMON_BACKGROUND.out, FUNCTIONAL_CONTROLS_AND_COMPARE.out)
    } else {
        if (!params.selection_manifest || !params.alignment_manifest || !params.run_id) error 'Explicit run, selection and alignment manifests are required'
        design = file(params.design, checkIfExists: true)
        alignments = file(params.alignment_manifest, checkIfExists: true)
        PREPARE_VALIDATION(design)
        VALIDATE_INPUTS(design, file(params.selection_manifest, checkIfExists: true), alignments,
                        file(params.settings, checkIfExists: true), PREPARE_VALIDATION.out)
        jobs = VALIDATE_INPUTS.out.jobs.splitCsv(header: true, sep: '\t').map { row ->
            tuple(row.hypothesis_id, row.strategy_id, row.replicate_id, row.gene_id,
                  file(row.job_file, checkIfExists: true), file(row.alignment_path, checkIfExists: true),
                  file(row.pool_path, checkIfExists: true), file(row.config_path, checkIfExists: true))
        }
        CAAS_POOLED(jobs, file(params.caastools_dir, checkIfExists: true).toAbsolutePath().toString())
        // Seed every hypothesis with one validated job metadata file. Even a
        // hypothesis with zero successful tasks must be reported, not disappear.
        seeds = jobs.unique { row -> row[0] }.map { row -> tuple(row[0], row[1], row[2], row[4]) }
        grouped = CAAS_POOLED.out.mix(seeds).groupTuple(by: [0,1,2])
        ASSEMBLE_FILTER_AND_QUERY(grouped, alignments)
        // The cohort comes from declared outputs; no glob of a publishing directory.
        cohort = ASSEMBLE_FILTER_AND_QUERY.out.collect(flat: false).map { items ->
            def rows = items.sort { a,b -> a[0] <=> b[0] }.collect { item ->
                def lockFile = item[1].resolve('summary.lock.json')
                def checksum = java.security.MessageDigest.getInstance('SHA-256').digest(lockFile.bytes).encodeHex().toString()
                "${item[0]}\t${item[1].toAbsolutePath()}\t${checksum}"
            }
            "hypothesis_id\tsummary_dir\tsummary_lock_sha256\n" + rows.join('\n') + '\n'
        }.collectFile(name: 'cohort.tsv', newLine: false)
        summaries = ASSEMBLE_FILTER_AND_QUERY.out.map { hid, folder -> folder }.collect()
        COMMON_BACKGROUND(cohort, summaries)
        FUNCTIONAL_CONTROLS_AND_COMPARE(COMMON_BACKGROUND.out, file(params.settings, checkIfExists: true), file(params.annotation_input, checkIfExists: true))
        PLOT_SUMMARIES(COMMON_BACKGROUND.out, FUNCTIONAL_CONTROLS_AND_COMPARE.out)
    }
}
