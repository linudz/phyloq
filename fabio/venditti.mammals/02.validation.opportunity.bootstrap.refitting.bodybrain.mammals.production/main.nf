nextflow.enable.dsl = 2

params.states_dir = 'inputs/states'
params.pss_core = 'inputs/core/pss.core.R'
params.output_dir = 'results'
params.conditional_replicates = 1000
params.conditional_replicates_per_task = 25
params.full_replicates = 100
params.full_replicates_per_task = 1
params.bootstrap_seed = 20260817
params.bootstrap_max_attempts = 10
params.depth_grid_points = 101
params.expected_traits = 3

def projectPath(value) {
    def text = value.toString()
    java.nio.file.Paths.get(text).isAbsolute() ? file(text) : file("${projectDir}/${text}")
}

process BOOTSTRAP_CHUNK {
    tag "${trait_id} ${bootstrap_type}: simulations ${simulation_start}-${simulation_start + chunk_replicates - 1}"
    label 'bootstrap'
    memory { bootstrap_type == 'full' ? params.full_memory : params.conditional_memory }
    time { bootstrap_type == 'full' ? params.full_time : params.conditional_time }
    publishDir "${projectDir}/${params.output_dir}/chunks", mode: 'copy', overwrite: true
    input:
    tuple val(trait_id), path(state_file), val(bootstrap_type), val(chunk_id),
          val(simulation_start), val(chunk_replicates)
    path pss_core
    path script
    output:
    tuple val(trait_id), val(bootstrap_type), val(chunk_id),
          path("${trait_id}.${bootstrap_type}.chunk_${chunk_id}.bootstrap.rds"), emit: chunks
    script:
    """
    \$CONDA_PREFIX/bin/Rscript '${script}' --state '${state_file}' --pss-core '${pss_core}' \
      --simulation-start '${simulation_start}' --replicates '${chunk_replicates}' \
      --bootstrap-types '${bootstrap_type}' --seed '${params.bootstrap_seed}' \
      --max-attempts '${params.bootstrap_max_attempts}' \
      --output '${trait_id}.${bootstrap_type}.chunk_${chunk_id}.bootstrap.rds'
    """
}

process AGGREGATE_BOOTSTRAP {
    tag 'aggregate opportunity-depth validation'
    label 'aggregate'
    publishDir "${projectDir}/${params.output_dir}", mode: 'copy', overwrite: true
    input:
    path chunks
    path script
    output:
    path 'mammal.brain_body.bootstrap.replicates.tsv.gz'
    path 'mammal.brain_body.depth_curves.tsv'
    path 'mammal.brain_body.depth_summary.tsv'
    script:
    """
    \$CONDA_PREFIX/bin/Rscript '${script}' --search-dir . \
      --conditional-replicates '${params.conditional_replicates}' \
      --full-replicates '${params.full_replicates}' \
      --grid-points '${params.depth_grid_points}' \
      --replicates-output mammal.brain_body.bootstrap.replicates.tsv.gz \
      --curves-output mammal.brain_body.depth_curves.tsv \
      --summary-output mammal.brain_body.depth_summary.tsv
    """
}

workflow {
    def statesDir = projectPath(params.states_dir)
    if (!statesDir.isDirectory()) error "Observed state directory not found: ${statesDir}."
    def stateFiles = statesDir.toFile().listFiles()
        .findAll { it.isFile() && it.name.endsWith('.observed.state.rds') }
        .collect { file(it) }.sort { it.name }
    if (stateFiles.size() != params.expected_traits as int) {
        error "Expected ${params.expected_traits} observed states but found ${stateFiles.size()}."
    }
    def stateMap = stateFiles.collectEntries { f ->
        [(f.name.replaceFirst(/\.observed\.state\.rds$/, '')): f]
    }
    def specifications = [
        [type: 'conditional', total: params.conditional_replicates as int,
         perTask: params.conditional_replicates_per_task as int],
        [type: 'full', total: params.full_replicates as int,
         perTask: params.full_replicates_per_task as int]
    ]
    specifications.each { spec ->
        if (spec.total < 1 || spec.perTask < 1) error "Invalid replicate settings for ${spec.type}."
    }
    def tasks = stateMap.keySet().sort().collectMany { trait ->
        specifications.collectMany { spec ->
            def nChunks = Math.ceil(spec.total / spec.perTask.toDouble()) as int
            (1..nChunks).collect { id ->
                def start = (id - 1) * spec.perTask + 1
                tuple(trait, stateMap[trait], spec.type, id, start,
                      Math.min(spec.perTask, spec.total - start + 1))
            }
        }
    }
    log.info "Validation tasks: ${tasks.size()} (${params.conditional_replicates} conditional and ${params.full_replicates} full replicates per trait)."
    BOOTSTRAP_CHUNK(channel.fromList(tasks), projectPath(params.pss_core),
                    projectPath('scripts/bootstrap_mammal_pss_chunk.R'))
    AGGREGATE_BOOTSTRAP(
        BOOTSTRAP_CHUNK.out.chunks.map { trait, type, id, f -> f }.collect(),
        projectPath('scripts/aggregate_mammal_bootstrap.R')
    )
}
