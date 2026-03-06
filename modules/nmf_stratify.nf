process NMF_STRATIFY {
    label 'process_medium'
    conda "${projectDir}/envs/python_ml.yml"
    publishDir "${params.outdir}/nmf", mode: 'copy'

    input:
    path mvalues_corrected
    path sample_sheet

    output:
    path "nmf_clusters.tsv",      emit: clusters
    path "W_matrix.tsv",          emit: w_matrix
    path "H_matrix.tsv",          emit: h_matrix
    path "rank_selection.png",    emit: rank_plot
    path "nmf_umap.png",         emit: umap_plot

    script:
    """
    nmf_stratify.py \\
        --mvalues ${mvalues_corrected} \\
        --sample_sheet ${sample_sheet} \\
        --k_min ${params.nmf_kmin} \\
        --k_max ${params.nmf_kmax} \\
        --n_runs ${params.nmf_nruns} \\
        --outdir .
    """
}
