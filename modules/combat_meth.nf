process COMBAT_METH {
    label 'process_medium'
    conda "${projectDir}/envs/r_methylation.yml"
    publishDir "${params.outdir}/combat_meth", mode: 'copy'

    input:
    path bedgraphs
    path sample_sheet

    output:
    path "mvalues_raw.tsv",        emit: mvalues_raw
    path "mvalues_corrected.tsv",  emit: mvalues_corrected
    path "beta_matrix.rds",        emit: beta_matrix
    path "cpg_manifest.tsv",       emit: cpg_manifest

    script:
    """
    combat_meth.R \\
        --bedgraph_dir . \\
        --sample_sheet ${sample_sheet} \\
        --outdir .
    """
}
