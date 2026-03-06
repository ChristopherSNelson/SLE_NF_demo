process FASTQC {
    label 'process_low'
    conda "${projectDir}/envs/fastqc.yml"
    publishDir "${params.outdir}/fastqc", mode: 'copy'

    input:
    tuple val(sample_id), path(reads)

    output:
    tuple val(sample_id), path("*.html"), emit: html
    tuple val(sample_id), path("*.zip"),  emit: zip

    script:
    """
    fastqc --threads ${task.cpus} ${reads}
    """
}
