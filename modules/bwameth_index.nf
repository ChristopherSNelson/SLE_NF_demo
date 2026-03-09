process BWAMETH_INDEX {
    label 'process_high'
    conda "${projectDir}/envs/bwameth.yml"
    // Index is intermediate — no publishDir (stays in work/ for resume)

    input:
    path genome

    output:
    path "${genome}*", emit: index

    script:
    """
    bwameth.py index-mem2 ${genome} || bwameth.py index ${genome}
    """
}
