process FETCH_SRA {
    label 'process_medium'
    conda "${projectDir}/envs/sra_tools.yml"
    publishDir "${params.outdir}/fastqs", mode: 'copy'

    input:
    tuple val(sample_id), val(srr_accession)

    output:
    tuple val(sample_id), path("${srr_accession}_{1,2}.fastq.gz"), emit: reads

    script:
    """
    fasterq-dump \\
        --split-files \\
        --threads ${task.cpus} \\
        --temp . \\
        --disk-limit-tmp 0 \\
        ${srr_accession}

    gzip ${srr_accession}_1.fastq ${srr_accession}_2.fastq
    """
}
