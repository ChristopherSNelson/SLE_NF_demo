process BWAMETH_INDEX {
    label 'process_high'
    conda "${projectDir}/envs/bwameth.yml"
    publishDir "${params.outdir}/bwameth_index", mode: 'copy'

    input:
    path genome

    output:
    path "${genome}*", emit: index

    script:
    """
    bwameth.py index-mem2 ${genome} || bwameth.py index ${genome}
    """
}
