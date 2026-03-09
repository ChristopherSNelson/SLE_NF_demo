process NMF_STRATIFY {
    label 'process_medium'
    conda "${projectDir}/envs/python_ml.yml"
    publishDir "${params.outdir}/nmf", mode: 'copy'

    input:
    path mvalues_corrected
    path sample_sheet
    path cell_fractions
    path dmr_bed

    output:
    path "nmf_clusters.tsv",      emit: clusters
    path "W_matrix.tsv",          emit: w_matrix
    path "H_matrix.tsv",          emit: h_matrix
    path "rank_selection.png",    emit: rank_plot
    path "nmf_umap.png",         emit: umap_plot
    path "stability_loo.tsv",    emit: loo_stability
    path "clinical_correlations.tsv", emit: clinical_corr, optional: true

    script:
    def cell_arg = cell_fractions.name != 'NO_CELL_FRACTIONS' ? "--cell_fractions ${cell_fractions}" : ''
    def dmr_arg = dmr_bed.name != 'NO_DMR_BED' ? "--dmr_bed ${dmr_bed}" : ''
    def clinical_arg = params.clinical_metadata ? "--clinical_metadata ${params.clinical_metadata}" : ''
    """
    nmf_stratify.py \\
        --mvalues ${mvalues_corrected} \\
        --sample_sheet ${sample_sheet} \\
        --k_min ${params.nmf_kmin} \\
        --k_max ${params.nmf_kmax} \\
        --n_runs ${params.nmf_nruns} \\
        ${cell_arg} \\
        ${dmr_arg} \\
        ${clinical_arg} \\
        --outdir .
    """
}
