process METHYLDACKEL {
    label 'process_medium'
    conda "${projectDir}/envs/methyldackel.yml"
    publishDir "${params.outdir}/methyldackel", mode: 'copy', pattern: "*.bedGraph"
    publishDir "${params.outdir}/methyldackel/mbias", mode: 'copy', pattern: "*_mbias_*.png"

    input:
    tuple val(sample_id), path(aln), path(aln_idx)
    path genome
    path fasta_fai

    output:
    tuple val(sample_id), path("${sample_id}_CpG.bedGraph"),  emit: bedgraph
    tuple val(sample_id), path("${sample_id}_mbias_OT.png"),  emit: mbias_ot
    tuple val(sample_id), path("${sample_id}_mbias_OB.png"),  emit: mbias_ob

    script:
    """
    MethylDackel mbias ${genome} ${aln} ${sample_id}_mbias

    MethylDackel extract \\
        --minDepth ${params.min_depth} \\
        -o ${sample_id} \\
        ${genome} \\
        ${aln}

    # Rename default output to expected name
    if [ -f "${sample_id}_CpG.bedGraph" ]; then
        echo "bedGraph found"
    elif [ -f "${sample_id}.bedGraph" ]; then
        mv ${sample_id}.bedGraph ${sample_id}_CpG.bedGraph
    fi

    # Convert mbias SVGs to PNG
    rsvg-convert -f png -o ${sample_id}_mbias_OT.png ${sample_id}_mbias_OT.svg
    rsvg-convert -f png -o ${sample_id}_mbias_OB.png ${sample_id}_mbias_OB.svg
    """
}
