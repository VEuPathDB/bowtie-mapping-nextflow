#!/usr/bin/env nextflow
nextflow.enable.dsl=2


process homer2Bed {
    //container "quay.io/biocontainers/homer:5.1--pl5262h9948957_0"
    container "jbrestel/bowtie"

    input:
    tuple val(meta), path(bam), path(bai)

    output:
    tuple val(meta), path("output.bedgraph")


    script:
    def fragLengthArg = ""
    if(task.ext.fragLength) {
        fragLengthArg = "fragLength {task.ext.fragLength}";
    }
    """
    makeTagDirectory $fragLengthArg homerTagDir $bam
    makeUCSCfile $fragLengthArg homerTagDir/ -o output.bedgraph
    gunzip output.bedgraph.gz
    """
}
