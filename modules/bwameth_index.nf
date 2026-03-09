process BWAMETH_INDEX {
    label 'process_high'
    conda "${projectDir}/envs/bwameth.yml"
    storeDir "${projectDir}/genome_index"

    input:
    path genome

    output:
    path "${genome}*", emit: index

    script:
    """
    bwameth.py index-mem2 ${genome} || bwameth.py index ${genome}
    """
}
