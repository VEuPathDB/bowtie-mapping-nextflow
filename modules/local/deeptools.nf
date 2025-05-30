#!/usr/bin/env nextflow
nextflow.enable.dsl=2


process bamCompare {
    container "quay.io/biocontainers/deeptools:3.5.6--pyhdfd78af_0"
    publishDir "$params.outputDir", pattern: "*.bw", mode: "copy", saveAs: { filename -> "${meta.id}_vs_${refMeta.id}.bw" }

    input:
    tuple val(meta), path(bam), path(bai), val(refMeta), path(refBam, name: 'reference.bam'), path(refBai, name: 'reference.bam.bai')

    output:
    tuple val(meta), path('ratio.bw')


    script:
    """
    bamCompare -b1 $bam -b2 reference.bam --operation ratio -o ratio.bw
    """
    

}
