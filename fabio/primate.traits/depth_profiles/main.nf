/*
 * Independent continuous phylogenetic-depth profile workflow.
 *
 * It recreates the validated conditional and full parametric bootstraps with
 * the same deterministic seeds, retaining the top/bottom-tail depth
 * percentiles needed for cumulative curves and global curve inference.
 */

params {
    traits_dir: String = '../inputs/traits'
    classified_dir: String = '../inputs/classified'
    tree: String = '../inputs/tree/science.abn7829_data_s4.nex.tree'
    taxonomy: String = '../inputs/taxonomy/species_family_primate_group.tsv'
    pss_core: String = '../scripts/pss.core.R'
    observed_fits: String = 'inputs/observed_model_fits.tsv'
    bootstrap_script: String = 'scripts/parametric_bootstrap_depth_chunk.R'
    summary_script: String = 'scripts/summarize_depth_profiles.R'
    output_dir: String = 'results'
    trait_regex: String = '.*'
    bootstrap_replicates: Integer = 1000
    replicates_per_task: Integer = 100
    bootstrap_seed: Integer = 20260812
    tail_proportion: Float = 0.01
    bootstrap_max_attempts: Integer = 10
    depth_grid_points: Integer = 101
    expected_traits: Integer = 153
}

def project_path(value) {
    def text = value.toString()
    java.nio.file.Paths.get(text).isAbsolute() ? file(text) : file("${projectDir}/${text}")
}

process DEPTH_PROFILE_CHUNK {
    tag "${trait_id}: depth profiles ${simulation_start}-${simulation_start + chunk_replicates - 1}"
    label 'depth_bootstrap'

    publishDir "${projectDir}/${params.output_dir}/chunks", mode: 'copy', overwrite: true

    input:
    tuple val(trait_id), path(trait_file), path(classified_file),
          val(chunk_id), val(simulation_start), val(chunk_replicates)
    path tree_file
    path taxonomy_file
    path pss_core
    path observed_fits
    path bootstrap_script

    output:
    tuple val(trait_id), val(chunk_id),
          path("${trait_id}.chunk_${chunk_id}.depth_profiles.rds"),
          emit: profiles

    script:
    """
    Rscript '${bootstrap_script}' \
      --trait '${trait_id}' \
      --trait-file '${trait_file}' \
      --tree '${tree_file}' \
      --taxonomy '${taxonomy_file}' \
      --pss-core '${pss_core}' \
      --observed-fits '${observed_fits}' \
      --observed-classified '${classified_file}' \
      --simulation-start '${simulation_start}' \
      --replicates '${chunk_replicates}' \
      --tail-proportion '${params.tail_proportion}' \
      --seed '${params.bootstrap_seed}' \
      --max-attempts '${params.bootstrap_max_attempts}' \
      --profile-output '${trait_id}.chunk_${chunk_id}.depth_profiles.rds'
    """
}

process AGGREGATE_DEPTH_PROFILES {
    tag 'validate and aggregate continuous depth profiles'
    label 'depth_aggregate'

    publishDir "${projectDir}/${params.output_dir}", mode: 'copy', overwrite: true

    input:
    path profile_files
    path summary_script

    output:
    path 'primate.traits.depth_profiles.replicates.tsv.gz', emit: replicates
    path 'primate.traits.depth_profiles.observed.tsv', emit: observed
    path 'primate.traits.depth_profiles.fits.tsv', emit: fits
    path 'primate.traits.depth_profiles.curves.tsv', emit: curves
    path 'primate.traits.depth_profiles.summary.tsv', emit: summary

    script:
    """
    Rscript '${summary_script}' \
      --search-dir . \
      --expected-replicates '${params.bootstrap_replicates}' \
      --expected-traits '${params.expected_traits}' \
      --tail-proportion '${params.tail_proportion}' \
      --grid-points '${params.depth_grid_points}' \
      --replicates-output primate.traits.depth_profiles.replicates.tsv.gz \
      --observed-output primate.traits.depth_profiles.observed.tsv \
      --fits-output primate.traits.depth_profiles.fits.tsv \
      --curves-output primate.traits.depth_profiles.curves.tsv \
      --summary-output primate.traits.depth_profiles.summary.tsv
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

    log.info "Depth-profile bootstrap: ${trait_ids.size()} traits, ${total_replicates} replicates, ${number_of_chunks} chunks per trait (${tasks.size()} tasks)."

    DEPTH_PROFILE_CHUNK(
        channel.fromList(tasks),
        project_path(params.tree),
        project_path(params.taxonomy),
        project_path(params.pss_core),
        project_path(params.observed_fits),
        project_path(params.bootstrap_script)
    )

    def profile_files = DEPTH_PROFILE_CHUNK.out.profiles
        .map { _trait, _chunk, path -> path }
        .collect()

    AGGREGATE_DEPTH_PROFILES(
        profile_files,
        project_path(params.summary_script)
    )
}
