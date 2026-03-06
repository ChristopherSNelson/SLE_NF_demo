process HOUSEMAN_DECONV {
    label 'process_low'
    conda "${projectDir}/envs/r_methylation.yml"
    publishDir "${params.outdir}/houseman", mode: 'copy'

    input:
    path beta_matrix
    path sample_sheet

    output:
    path "cell_fractions.tsv",  emit: fractions
    path "cell_fractions.png",  emit: plot

    script:
    """
    houseman_deconv.R \\
        --beta_matrix ${beta_matrix} \\
        --sample_sheet ${sample_sheet} \\
        --outdir .
    """
}
