nextflow.enable.dsl=2

params.caastools  = "${baseDir}/caastools/ct"
params.alignments = "${baseDir}/alignments/*"
params.configs    = "${baseDir}/cfg.files/*"
params.outdir     = "${baseDir}/pipeline.results"
params.fmt        = "fasta"

process CAASTOOLS_DISCOVERY {

    tag "${gene_id} | ${cfg_id}"

    stageInMode 'copy'

    publishDir(
        "${params.outdir}/${gene_id}_results",
        mode: 'copy'
    )

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

workflow {

    alignments_ch = Channel
        .fromPath(params.alignments, checkIfExists: true)
        .map { file ->
            tuple(file.baseName, file)
        }

    configs_ch = Channel
        .fromPath(params.configs, checkIfExists: true)
        .map { file ->
            tuple(file.baseName, file)
        }

    alignments_ch
        .combine(configs_ch)
        .map { gene_id, alignment, cfg_id, config ->
            tuple(gene_id, alignment, cfg_id, config)
        }
        | CAASTOOLS_DISCOVERY
}