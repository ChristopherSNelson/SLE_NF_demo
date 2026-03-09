process BAM_TO_CRAM {
    label 'process_medium'
    conda "${projectDir}/envs/picard.yml"
    publishDir "${params.outdir}/cram", mode: 'copy'

    input:
    tuple val(sample_id), path(bam), path(bai)
    path genome
    path fasta_fai

    output:
    tuple val(sample_id), path("${sample_id}.cram"),     emit: cram
    tuple val(sample_id), path("${sample_id}.cram.crai"), emit: crai

    script:
    """
    samtools view -C -T ${genome} -o ${sample_id}.cram ${bam}
    samtools index ${sample_id}.cram
    """
}
