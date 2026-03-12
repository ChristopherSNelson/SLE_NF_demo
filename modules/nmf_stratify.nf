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
    path "nmf_clusters.tsv",           optional: true, emit: clusters
    path "W_matrix.tsv",               optional: true, emit: w_matrix
    path "H_matrix.tsv",               optional: true, emit: h_matrix
    path "rank_selection.png",         optional: true, emit: rank_plot
    path "H_matrix_nmf_ordered.png",   optional: true, emit: heatmap_ordered
    path "H_matrix_clustered.png",     optional: true, emit: heatmap_clustered
    path "H_matrix_unclustered.png",   optional: true, emit: heatmap_unclustered
    path "stability_loo.tsv",          optional: true, emit: loo_stability
    path "clinical_correlations.tsv",  optional: true, emit: clinical_corr

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
    if [ -f nmf_clusters.tsv ] && [ -f H_matrix.tsv ]; then
        generate_nmf_heatmaps.py --h_matrix H_matrix.tsv --clusters nmf_clusters.tsv --outdir .
    fi
    """
}
