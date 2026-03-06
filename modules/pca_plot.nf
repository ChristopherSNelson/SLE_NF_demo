process PCA_PLOT {
    label 'process_low'
    conda "${projectDir}/envs/r_methylation.yml"
    publishDir "${params.outdir}/pca", mode: 'copy'

    input:
    path mvalues_raw
    path mvalues_corrected
    path sample_sheet

    output:
    path "*.png",              emit: plots
    path "pca_variance.tsv",   emit: variance

    script:
    """
    pca_plot.R \\
        --mvalues_raw ${mvalues_raw} \\
        --mvalues_corrected ${mvalues_corrected} \\
        --sample_sheet ${sample_sheet} \\
        --outdir .
    """
}
