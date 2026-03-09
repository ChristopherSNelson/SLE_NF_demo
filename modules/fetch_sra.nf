process FETCH_SRA {
    label 'process_medium'
    conda "${projectDir}/envs/sra_tools.yml"
    storeDir "${projectDir}/fastq_cache"

    input:
    tuple val(sample_id), val(srr_accession)

    output:
    tuple val(sample_id), path("${srr_accession}_{1,2}.fastq.gz"), emit: reads

    script:
    """
    # prefetch first to avoid fasterq-dump APFS disk-limit bug (fs_type=unexpected)
    prefetch ${srr_accession} --max-size 100G -o ${srr_accession}.sra

    fasterq-dump \\
        --split-files \\
        --threads ${task.cpus} \\
        --temp . \\
        ${srr_accession}.sra

    rm -f ${srr_accession}.sra
    gzip ${srr_accession}_1.fastq ${srr_accession}_2.fastq
    """
}
