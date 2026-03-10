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
    # prefetch first to download .sra locally
    prefetch ${srr_accession} --max-size 100G -o ${srr_accession}.sra

    # fastq-dump streams to gzip with no temp space (fasterq-dump needs ~3x SRA size)
    fastq-dump --split-files --gzip ${srr_accession}.sra
    rm -f ${srr_accession}.sra
    """
}
