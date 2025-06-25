#!/usr/bin/env nextflow
nextflow.enable.dsl=2


process homerMakeTagDir {
    container = 'veupathdb/bowtiemapping:1.0.0'

    input:
    tuple val(meta), path(bam), path(bai)

    output:
    tuple val(meta), path("homerTagDir")


    script:
    def fragLengthArg = ""
    if(task.ext.fragLength) {
        fragLengthArg = "-fragLength {task.ext.fragLength}";
    }
    """
    makeTagDirectory $fragLengthArg homerTagDir $bam
    """
}

process homer2Bed {
    container = 'veupathdb/bowtiemapping:1.0.0'

    input:
    tuple val(meta), path(homerTagDir)

    output:
    tuple val(meta), path("output.bedgraph")


    script:
    def fragLengthArg = ""
    if(task.ext.fragLength) {
        fragLengthArg = "-fragLength {task.ext.fragLength}";
    }
    """
    makeUCSCfile $fragLengthArg $homerTagDir -o output.bedgraph
    gunzip output.bedgraph.gz
    """
}


process homerFindPeaks_histonemod {
    container = 'veupathdb/bowtiemapping:1.0.0'

    publishDir params.outputDir, mode: 'copy', pattern: "*_peaks.txt"
    
    input:
    tuple val(meta), path(homerTagDir), val(inputMeta), path(inputHomerTagDir, name: "inputTagDir")

    output:
    tuple val(meta), path("${meta.id}_peaks.gff"), path("${meta.id}_peaks.txt"), path("${meta.id}.config")

    script:
    def fragLengthArg = ""
    if(task.ext.fragLength) {
        fragLengthArg = "-fragLength {task.ext.fragLength}";
    }
    """
    findPeaks $homerTagDir -style histone -o auto -i $inputHomerTagDir $fragLengthArg
    parseHomerRegions.pl --sample ${meta.id} --inputFile ${homerTagDir}/regions.txt --outputGFF ${meta.id}_peaks.gff --outputTab ${meta.id}_peaks.txt --outputConfig ${meta.id}.config --profileSetName "${params.profileSetName}"
    """
}
