#!/usr/bin/env nextflow
nextflow.enable.dsl=2


process createIndex {
    container = 'veupathdb/bowtiemapping:1.0.0'

    input:
    path fasta

    output:
    path 'index.*'

    script:
    """
    mkdir bowtie2
    bowtie2-build $fasta index
    """
}



process bowtie2 {
    container = 'veupathdb/bowtiemapping:1.0.0'

    input:
    path index
    tuple val(meta), path(reads)

    output:
    tuple val(meta), path('tmpOut.bam')


    script:
    def fileTypeArg = params.inputFileType == "fastq" ? "-q" : "-f"
    def readsArg = ""
    if (params.hasPairedReads) {
        readsArg = "-1 ${reads[0]} -2 ${reads[1]}"

    } else {
        readsArg = "-U ${reads}"

    }
    
    """
    bowtie2 ${fileTypeArg} -x index ${readsArg} -S tmpOut.sam ${task.ext.args}
    samtools view -buS tmpOut.sam | samtools sort -o tmpOut.bam
    """
}


process removePCRDuplicates {
    container = 'veupathdb/bowtiemapping:1.0.0'

//    publishDir "$params.outputDir", pattern: "*.bam", mode: "copy", saveAs: { filename -> "${sampleName}.bam" }
//    publishDir "$params.outputDir", pattern: "*.bam.bai", mode: "copy", saveAs: { filename -> "${sampleName}.bam.bai" }
//    publishDir "$params.outputDir", pattern: "*.bed", mode: "copy", saveAs: { filename -> "${sampleName}.bed" }

    input:
    tuple val(meta), path(bam)

    output:
    tuple val(meta), path('rmdup.bam')

    script:
    """
    samtools rmdup -S $bam rmdup.bam
    """
}


process indexBam {
    container = 'veupathdb/bowtiemapping:1.0.0'

    publishDir "$params.outputDir", pattern: "*.bam", mode: 'copy',saveAs: { filename -> params.saveAlignments ? "${meta.id}.bam" : null }
    publishDir "$params.outputDir", pattern: "*.bai", mode: 'copy',saveAs: { filename -> params.saveAlignments ? "${meta.id}.bam.bai" : null }

    input:
    tuple val(meta), path(bam, name: "output.bam")

    output:
    tuple val(meta), path("output.bam"), path("output.bam.bai")

    script:
    """
    samtools index output.bam
    """
    
    
}



workflow bowtieMapping {

    main:
    samples = Channel.fromPath(params.input + "/" + params.samplesheetFileName )
        .splitCsv( skip:1)

    sample_ch = samples.map { row ->

        def fasta1;
        if (row[1].startsWith("/")) {
            fasta1 = file(row[1]);
        }
        else {
            fasta1 = file(params.input + "/" + row[1])
        }

        if(params.hasPairedReads && row[2]) {
            def fasta2;
            if (row[2].startsWith("/")) {
                fasta2 = file(row[2]);
            }
            else {
                fasta2 = file(params.input + "/" + row[2])
            }

            return [ [id: row[0], ref: row[3] ], [fasta1, fasta2] ]
        }
        return [ [id: row[0], ref: row[3] ], [fasta1] ]
    }

    index = createIndex(params.genome)
    alignments = bowtie2(index, sample_ch)


    if(params.removePCRDuplicates) {
        output = removePCRDuplicates(alignments)
    }
    else {
        output = alignments
    }

    indexBam(output)

    emit:
    indexBam.out
    
    

}
