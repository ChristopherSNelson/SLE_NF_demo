process REGION_DETECT {
    label 'process_medium'
    conda "${projectDir}/envs/python_ml.yml"
    publishDir "${params.outdir}/region_detect", mode: 'copy'

    input:
    path mvalues_corrected
    path sample_sheet

    output:
    path "candidate_dmrs.bed",    emit: dmrs
    path "window_results.tsv",    emit: windows
    path "dmr_manhattan.png",     emit: plot

    script:
    """
    region_detect.py \\
        --mvalues ${mvalues_corrected} \\
        --sample_sheet ${sample_sheet} \\
        --window_size ${params.window_size} \\
        --step_size ${params.step_size} \\
        --min_cpgs ${params.min_cpgs} \\
        --outdir .
    """
}
