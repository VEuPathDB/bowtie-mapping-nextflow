#!/usr/bin/env nextflow
nextflow.enable.dsl=2



process chromSizes {
    container = 'veupathdb/bowtiemapping:1.0.0'

    input:
    path(genome)

    output:
    path("chrom.sizes")

    script:
    """
    samtools faidx $genome
    cut -f 1,2 ${genome}.fai >chrom.sizes
    """
}


