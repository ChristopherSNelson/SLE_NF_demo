process BWAMETH_ALIGN {
    label 'process_high'
    conda "${projectDir}/envs/bwameth.yml"
    publishDir "${params.outdir}/bwameth_align", mode: 'copy'

    input:
    tuple val(sample_id), path(reads)
    path genome
    path index

    output:
    tuple val(sample_id), path("${sample_id}.sorted.bam"),     emit: bam
    tuple val(sample_id), path("${sample_id}.sorted.bam.bai"), emit: bai
    tuple val(sample_id), path("${sample_id}.flagstat.txt"),    emit: flagstat

    script:
    """
    bwameth.py \\
        --threads ${task.cpus} \\
        --reference ${genome} \\
        ${reads[0]} ${reads[1]} \\
    | samtools sort -@ ${task.cpus} -o ${sample_id}.sorted.bam

    samtools index ${sample_id}.sorted.bam
    samtools flagstat ${sample_id}.sorted.bam > ${sample_id}.flagstat.txt
    """
}
