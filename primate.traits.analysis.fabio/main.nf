/*
 * Independent parametric-bootstrap workflow for the primate PSS analysis.
 * The retained empirical traits and classified PSS tables are frozen inputs;
 * this workflow runs only conditional and full parametric bootstraps.
 */

params {
    traits_dir: String = 'inputs/traits'
    classified_dir: String = 'inputs/classified'
    tree: String = 'inputs/tree/science.abn7829_data_s4.nex.tree'
    taxonomy: String = 'inputs/taxonomy/species_family_primate_group.tsv'
    pss_core: String = 'scripts/pss.core.R'
    bootstrap_script: String = 'scripts/parametric_bootstrap_chunk.R'
    summary_script: String = 'scripts/summarize_parametric_bootstrap.R'
    output_dir: String = 'results'
    trait_regex: String = '.*'
    bootstrap_replicates: Integer = 1000
    replicates_per_task: Integer = 100
    bootstrap_seed: Integer = 20260812
    tail_proportion: Float = 0.01
    bootstrap_max_attempts: Integer = 10
    expected_traits: Integer = 153
}

def project_path(value) {
    def text = value.toString()
    java.nio.file.Paths.get(text).isAbsolute() ? file(text) : file("${projectDir}/${text}")
}

process BOOTSTRAP_CHUNK {
    tag "${trait_id}: simulations ${simulation_start}-${simulation_start + chunk_replicates - 1}"
    label 'bootstrap'

    publishDir "${params.output_dir}/chunks", mode: 'copy', overwrite: true

    input:
    tuple val(trait_id), path(trait_file), path(classified_file),
          val(chunk_id), val(simulation_start), val(chunk_replicates)
    path tree_file
    path taxonomy_file
    path pss_core
    path bootstrap_script

    output:
    tuple val(trait_id), val(chunk_id),
          path("${trait_id}.chunk_${chunk_id}.bootstrap_simulations.tsv"),
          emit: simulations
    tuple val(trait_id), val(chunk_id),
          path("${trait_id}.chunk_${chunk_id}.bootstrap_observed.tsv"),
          emit: observed
    tuple val(trait_id), val(chunk_id),
          path("${trait_id}.chunk_${chunk_id}.bootstrap_fit.tsv"),
          emit: fits

    script:
    """
    Rscript '${bootstrap_script}' \
      --trait '${trait_id}' \
      --trait-file '${trait_file}' \
      --tree '${tree_file}' \
      --taxonomy '${taxonomy_file}' \
      --pss-core '${pss_core}' \
      --observed-classified '${classified_file}' \
      --simulation-start '${simulation_start}' \
      --replicates '${chunk_replicates}' \
      --tail-proportion '${params.tail_proportion}' \
      --seed '${params.bootstrap_seed}' \
      --max-attempts '${params.bootstrap_max_attempts}' \
      --simulations-output '${trait_id}.chunk_${chunk_id}.bootstrap_simulations.tsv' \
      --observed-output '${trait_id}.chunk_${chunk_id}.bootstrap_observed.tsv' \
      --fit-output '${trait_id}.chunk_${chunk_id}.bootstrap_fit.tsv'
    """
}

process AGGREGATE_BOOTSTRAP {
    tag 'validate and aggregate bootstrap results'
    label 'aggregate'

    publishDir params.output_dir, mode: 'copy', overwrite: true

    input:
    path simulation_files
    path observed_files
    path fit_files
    path summary_script

    output:
    path 'primate.traits.bootstrap.simulations.tsv', emit: simulations
    path 'primate.traits.bootstrap.observed.tsv', emit: observed
    path 'primate.traits.bootstrap.fits.tsv', emit: fits
    path 'primate.traits.bootstrap.summary.tsv', emit: summary
    path 'primate.traits.bootstrap.model.stability.tsv', emit: stability

    script:
    """
    Rscript '${summary_script}' \
      --search-dir . \
      --expected-replicates '${params.bootstrap_replicates}' \
      --expected-traits '${params.expected_traits}' \
      --simulations-output primate.traits.bootstrap.simulations.tsv \
      --observed-output primate.traits.bootstrap.observed.tsv \
      --fits-output primate.traits.bootstrap.fits.tsv \
      --summary-output primate.traits.bootstrap.summary.tsv \
      --stability-output primate.traits.bootstrap.model.stability.tsv
    """
}

workflow {
    def traits_dir = project_path(params.traits_dir)
    def classified_dir = project_path(params.classified_dir)
    def trait_pattern = ~params.trait_regex.toString()

    if (!traits_dir.isDirectory()) {
        error "Trait directory not found: ${traits_dir}"
    }
    if (!classified_dir.isDirectory()) {
        error "Classified-results directory not found: ${classified_dir}"
    }

    def trait_files = traits_dir.toFile().listFiles()
        .findAll { entry -> entry.isFile() && entry.name.endsWith('.tsv') }
        .collect { entry -> file(entry) }
        .findAll { entry -> entry.baseName ==~ trait_pattern }
    def classified_files = classified_dir.toFile().listFiles()
        .findAll { entry -> entry.isFile() && entry.name.endsWith('.score_results.classified.tsv') }
        .collect { entry -> file(entry) }
        .findAll { entry ->
            entry.name.replaceFirst(/\.score_results\.classified\.tsv$/, '') ==~ trait_pattern
        }

    def trait_map = trait_files.collectEntries { entry -> [(entry.baseName): entry] }
    def classified_map = classified_files.collectEntries { entry ->
        [(entry.name.replaceFirst(/\.score_results\.classified\.tsv$/, '')): entry]
    }
    def trait_ids = trait_map.keySet().sort()
    def classified_ids = classified_map.keySet().sort()

    if (trait_ids != classified_ids) {
        error "Trait/classified inputs do not match. Trait-only: ${trait_ids - classified_ids}; classified-only: ${classified_ids - trait_ids}"
    }
    if (trait_ids.size() != params.expected_traits as int) {
        error "Expected ${params.expected_traits} matched traits but found ${trait_ids.size()}."
    }

    def total_replicates = params.bootstrap_replicates as int
    def per_task = params.replicates_per_task as int
    if (total_replicates < 1 || per_task < 1) {
        error 'bootstrap_replicates and replicates_per_task must be positive integers.'
    }
    def number_of_chunks = Math.ceil(total_replicates / per_task.toDouble()) as int
    def chunks = (1..number_of_chunks).collect { chunk_id ->
        def start = (chunk_id - 1) * per_task + 1
        def count = Math.min(per_task, total_replicates - start + 1)
        tuple(chunk_id, start, count)
    }
    def tasks = trait_ids.collectMany { trait_id ->
        chunks.collect { chunk ->
            tuple(
                trait_id, trait_map[trait_id], classified_map[trait_id],
                chunk[0], chunk[1], chunk[2]
            )
        }
    }

    log.info "Bootstrap export: ${trait_ids.size()} traits, ${total_replicates} replicates, ${number_of_chunks} chunks per trait (${tasks.size()} tasks)."

    def task_channel = channel.fromList(tasks)
    BOOTSTRAP_CHUNK(
        task_channel,
        project_path(params.tree),
        project_path(params.taxonomy),
        project_path(params.pss_core),
        project_path(params.bootstrap_script)
    )

    def simulation_files = BOOTSTRAP_CHUNK.out.simulations.map { _trait, _chunk, path -> path }.collect()
    def observed_files = BOOTSTRAP_CHUNK.out.observed.map { _trait, _chunk, path -> path }.collect()
    def fit_files = BOOTSTRAP_CHUNK.out.fits.map { _trait, _chunk, path -> path }.collect()

    AGGREGATE_BOOTSTRAP(
        simulation_files,
        observed_files,
        fit_files,
        project_path(params.summary_script)
    )
}
