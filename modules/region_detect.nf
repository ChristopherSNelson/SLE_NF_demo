process REGION_DETECT {
    label 'process_medium'
    conda "${projectDir}/envs/r_methylation.yml"
    publishDir "${params.outdir}/region_detect", mode: 'copy'

    input:
    path bedgraphs
    path sample_sheet

    output:
    path "candidate_dmrs.bed",    emit: dmrs
    path "window_results.tsv",    emit: windows
    path "dmr_manhattan.png",     emit: plot

    script:
    def bg_list = bedgraphs.collect { it.name }.join(',')
    """
    region_detect.R \\
        --bedgraphs ${bg_list} \\
        --sample_sheet ${sample_sheet} \\
        --outdir .
    """
}
