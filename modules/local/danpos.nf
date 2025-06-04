#!/usr/bin/env nextflow
nextflow.enable.dsl=2



process danpos {
    container = 'biocontainers/danpos:v2.2.2_cv3'

    input:
    tuple val(meta), path(bam), path(bai)

    output:
    tuple val(meta), path("out/pooled/output.smooth.wig")


    script:
    """
    mkdir out
    danpos.py dpos $bam -o out
    """
}
