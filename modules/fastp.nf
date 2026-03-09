process FASTP {
    label 'process_medium'
    conda "${projectDir}/envs/fastp.yml"
    publishDir "${params.outdir}/fastp", mode: 'copy', pattern: '*.{html,json}'

    input:
    tuple val(sample_id), path(reads)

    output:
    tuple val(sample_id), path("${sample_id}_trimmed_{1,2}.fastq.gz"), emit: trimmed_reads
    tuple val(sample_id), path("${sample_id}_fastp.html"),             emit: html
    tuple val(sample_id), path("${sample_id}_fastp.json"),             emit: json

    script:
    """
    fastp \\
        --in1 ${reads[0]} \\
        --in2 ${reads[1]} \\
        --out1 ${sample_id}_trimmed_1.fastq.gz \\
        --out2 ${sample_id}_trimmed_2.fastq.gz \\
        --html ${sample_id}_fastp.html \\
        --json ${sample_id}_fastp.json \\
        --thread ${task.cpus} \\
        --detect_adapter_for_pe
    """
}
