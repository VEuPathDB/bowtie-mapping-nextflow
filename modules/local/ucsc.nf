#!/usr/bin/env nextflow
nextflow.enable.dsl=2


process bedgraph2bigwig {
    container 'quay.io/biocontainers/ucsc-bedgraphtobigwig:469--h9b8f530_0'

    publishDir params.outputDir, mode: 'copy'

    input:
    tuple val(meta), path(bed)
    path sizes

    output:
    path "${meta.id}.bw"

    script:
    """
    bedGraphToBigWig $bed $sizes ${meta.id}.bw
    """

}

process wigToBigWig {
    container = "quay.io/biocontainers/ucsc-wigtobigwig:472--h664eb37_2"

    publishDir params.outputDir, mode: 'copy'

    input:
    tuple val(meta), path(wig)
    path(chromSizes)

    output:
    path "${meta.id}.bw"

    script:
    """
    wigToBigWig -clip $wig $chromSizes ${meta.id}.bw
    """

}
