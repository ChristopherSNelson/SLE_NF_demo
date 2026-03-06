process TRIM_GALORE {
    label 'process_low'
    conda "${projectDir}/envs/trim_galore.yml"
    publishDir "${params.outdir}/trim_galore", mode: 'copy'

    input:
    tuple val(sample_id), path(reads)

    output:
    tuple val(sample_id), path("*_val_{1,2}.fq.gz"), emit: trimmed_reads
    tuple val(sample_id), path("*_trimming_report.txt"), emit: reports
    tuple val(sample_id), path("*_val_{1,2}_fastqc.{html,zip}"), emit: fastqc

    script:
    """
    trim_galore \\
        --paired \\
        --fastqc \\
        --cores ${task.cpus} \\
        ${reads[0]} ${reads[1]}
    """
}
