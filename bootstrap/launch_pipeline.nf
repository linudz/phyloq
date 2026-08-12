nextflow.enable.dsl=2

process CAASTOOLS_DISCOVERY {

    tag "${gene_id} | ${cfg_id}"

    stageInMode 'copy'

    publishDir path: { "${params.outdir}/${gene_id}_results" }, mode: 'copy'

    input:
    tuple val(gene_id), path(alignment), val(cfg_id), path(config)

    output:
    path "${gene_id}.${cfg_id}.caas"

    script:
    """
    python "${params.caastools}" discovery \\
        -a "${alignment}" \\
        -t "${config}" \\
        -o "${gene_id}.${cfg_id}.caas" \\
        --fmt ${params.fmt}

    if [ ! -f "${gene_id}.${cfg_id}.caas" ]; then
        touch "${gene_id}.${cfg_id}.caas"
    fi
    """
}

process INSTALL_RERCONVERGE {

    tag "check/install RERconverge"
    cache false

    output:
    path "rerconverge.ready", emit: ready

    script:
    """
    Rscript -e '
      if (!requireNamespace("remotes", quietly = TRUE)) {
        stop("The R package remotes is required to install RERconverge.")
      }

      if (!requireNamespace("RERconverge", quietly = TRUE)) {
        message("RERconverge is not installed; installing it from GitHub...")
        remotes::install_github(
          "nclark-lab/RERconverge",
          ref = "master",
          dependencies = TRUE,
          upgrade = "never"
        )
      }

      suppressPackageStartupMessages(library(RERconverge))
      cat("RERconverge ", as.character(packageVersion("RERconverge")), " is ready.\\n", sep = "")
    '

    touch rerconverge.ready
    """
}

process RERCONVERGE {

    tag "${cfg_id} | ${params.rerconverge_max_trees} genes"
    stageInMode 'copy'
    publishDir path: params.rerconverge_outdir, mode: 'copy', overwrite: true

    input:
    tuple val(cfg_id), path(config), path(tree_manifest), path(master_tree), path(rer_script), path(rerconverge_ready)

    output:
    path "rerconverge.${cfg_id}"

    script:
    """
    Rscript "${rer_script}" \\
        --trees "${tree_manifest}" \\
        --master-tree "${master_tree}" \\
        --phenotype "${config}" \\
        --outdir "rerconverge.${cfg_id}" \\
        --max-trees ${params.rerconverge_max_trees} \\
        --min-trees ${params.rerconverge_min_trees} \\
        --min-valid ${params.rerconverge_min_valid} \\
        --min-species ${params.rerconverge_min_species} \\
        --min-foreground ${params.rerconverge_min_foreground}
    """
}

workflow RERCONVERGE_FAST {
    config_glob = "${params.config_dir}/${params.config_pattern}"

    rerconverge_configs_ch = Channel
        .fromPath(config_glob, checkIfExists: true)
        .filter { file -> !file.isDirectory() }
        .map { config ->
            tuple(
                config.baseName,
                config,
                file(params.rerconverge_trees, checkIfExists: true),
                file(params.rerconverge_master, checkIfExists: true),
                file(params.rerconverge_script, checkIfExists: true)
            )
        }

    INSTALL_RERCONVERGE()

    rerconverge_ready_ch = rerconverge_configs_ch
        .combine(INSTALL_RERCONVERGE.out.ready)

    RERCONVERGE(rerconverge_ready_ch)
}

workflow {
    config_glob = "${params.config_dir}/${params.config_pattern}"

    alignments_ch = Channel
        .fromPath(params.alignments, checkIfExists: true)
        .map { file ->
            tuple(file.baseName, file)
        }

    configs_ch = Channel
        .fromPath(config_glob, checkIfExists: true)
        .filter { file -> !file.isDirectory() }
        .map { file ->
            tuple(file.baseName, file)
        }

    rerconverge_configs_ch = Channel
        .fromPath(config_glob, checkIfExists: true)
        .filter { file -> !file.isDirectory() }
        .map { config ->
            tuple(
                config.baseName,
                config,
                file(params.rerconverge_trees, checkIfExists: true),
                file(params.rerconverge_master, checkIfExists: true),
                file(params.rerconverge_script, checkIfExists: true)
            )
        }

    alignments_ch
        .combine(configs_ch)
        .map { gene_id, alignment, cfg_id, config ->
            tuple(gene_id, alignment, cfg_id, config)
        }
        | CAASTOOLS_DISCOVERY

    INSTALL_RERCONVERGE()

    rerconverge_ready_ch = rerconverge_configs_ch
        .combine(INSTALL_RERCONVERGE.out.ready)

    RERCONVERGE(rerconverge_ready_ch)
}
