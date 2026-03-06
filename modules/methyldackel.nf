process METHYLDACKEL {
    label 'process_medium'
    conda "${projectDir}/envs/methyldackel.yml"
    publishDir "${params.outdir}/methyldackel", mode: 'copy'

    input:
    tuple val(sample_id), path(bam), path(bai)
    path genome
    path fasta_fai

    output:
    tuple val(sample_id), path("${sample_id}_CpG.bedGraph"), emit: bedgraph
    tuple val(sample_id), path("${sample_id}_mbias.txt"),    emit: mbias

    script:
    """
    MethylDackel mbias ${genome} ${bam} ${sample_id}_mbias

    MethylDackel extract \\
        --minDepth ${params.min_depth} \\
        -o ${sample_id} \\
        ${genome} \\
        ${bam}

    # Rename default output to expected name
    if [ -f "${sample_id}_CpG.bedGraph" ]; then
        echo "bedGraph found"
    elif [ -f "${sample_id}.bedGraph" ]; then
        mv ${sample_id}.bedGraph ${sample_id}_CpG.bedGraph
    fi

    # Collect mbias outputs into single file
    if [ ! -f "${sample_id}_mbias.txt" ]; then
        cat ${sample_id}_mbias*.txt > ${sample_id}_mbias.txt 2>/dev/null || touch ${sample_id}_mbias.txt
    fi
    """
}
