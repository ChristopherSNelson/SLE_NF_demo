process MARK_DUPLICATES {
    label 'process_medium'
    conda "${projectDir}/envs/picard.yml"
    publishDir "${params.outdir}/mark_duplicates", mode: 'copy', pattern: '*.dup_metrics.txt'

    input:
    tuple val(sample_id), path(bam)

    output:
    tuple val(sample_id), path("${sample_id}.dedup.bam"),              emit: bam
    tuple val(sample_id), path("${sample_id}.dedup.bam.bai"),          emit: bai
    tuple val(sample_id), path("${sample_id}.dup_metrics.txt"),        emit: metrics

    script:
    """
    picard MarkDuplicates \\
        INPUT=${bam} \\
        OUTPUT=${sample_id}.dedup.bam \\
        METRICS_FILE=${sample_id}.dup_metrics.txt \\
        REMOVE_DUPLICATES=false \\
        CREATE_INDEX=false \\
        VALIDATION_STRINGENCY=LENIENT

    samtools index ${sample_id}.dedup.bam
    """
}
