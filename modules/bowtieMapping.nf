#!/usr/bin/env nextflow
nextflow.enable.dsl=2


process createIndex {
  container = 'veupathdb/bowtiemapping:1.0.0'
  
  input:
    path databaseFasta
    
  output:
    path 'index.*'

  script:
    template 'createIndex.bash'
}


process downloadFiles {
  container = 'veupathdb/bowtiemapping:1.0.0'

  input:
    val id
    
  output:
    tuple val(id), path("${id}**.fastq")

  script:
    template 'downloadFiles.bash'
}


process bowtieFromLocalSingle {
  container = 'veupathdb/bowtiemapping:1.0.0'

  input:
    path indexfiles
    path mateA
    val index
    val hasPairedReads
    val extraBowtieParams

  output:
    path '*.bam'

  script:
    template 'bowtieFromLocalSingle.bash'
}

process bowtieFromLocalPaired {
  container = 'veupathdb/bowtiemapping:1.0.0'

  input:
    path indexfiles
    path mateA
    path mateB
    val index
    val hasPairedReads
    val extraBowtieParams

  output:
    path '*.bam'

  script:
    template 'bowtieFromLocalPaired.bash'
}


process bowtieFromSra {
  container = 'veupathdb/bowtiemapping:1.0.0'

  input:
    path indexfiles
    tuple val(sample), path(readsFastq)
    val index
    val hasPairedReads
    val extraBowtieParams

  output:
    path '*.bam'

  script:
    template 'bowtieFromSra.bash'
}


process PCRDuplicates {
  container = 'veupathdb/bowtiemapping:1.0.0'

  publishDir "$params.outputDir", pattern: "*.bam", mode: "copy", saveAs: { filename -> "${sampleName}.bam" }
  publishDir "$params.outputDir", pattern: "*.bam.bai", mode: "copy", saveAs: { filename -> "${sampleName}.bam.bai" }
  publishDir "$params.outputDir", pattern: "*.bed", mode: "copy", saveAs: { filename -> "${sampleName}.bed" }

  input:
    path bamfile
    val removePCRDuplicates
    val writeBedFile
    val sampleName

  output:
    path 'out.*'

  script:
    template 'pcrDuplicates.bash'
}


workflow makeIndex {
  main:
    if(params.preconfiguredDatabase) {
      indexfiles = file( params.databaseFileDir + "/*.bt*" )
      indexFileBasename = params.indexFileBasename
    }
    else {
      indexfiles = createIndex(params.databaseFasta)
      indexFileBasename = "index"
    }

  emit:
    indexfiles
    indexFileBasename
}

workflow sra {
  take:
    indexfiles
    indexFileBasename
    accessions

  main:
    ids = Channel.fromList( accessions )
    files = downloadFiles( ids )
    bowtieResults = bowtieFromSra( indexfiles, files, indexFileBasename, params.hasPairedReads, params.extraBowtieParams ) 
    PCRDuplicates(bowtieResults, params.removePCRDuplicates, params.writeBedFile, params.sampleName)
}

workflow local {
  take:
    indexfiles
    indexFileBasename
    file

  main:
    if (params.hasPairedReads) {
      bowtieResults = bowtieFromLocalPaired( indexfiles, file, params.mateB, indexFileBasename, params.hasPairedReads, params.extraBowtieParams )
    }
    else {
      bowtieResults = bowtieFromLocalSingle( indexfiles, file, indexFileBasename, params.hasPairedReads, params.extraBowtieParams )
    }
    PCRDuplicates(bowtieResults, params.removePCRDuplicates, params.writeBedFile, params.sampleName)
}

workflow bowtieMapping {
  take:
    accessions

  main:
    makeIndex()
    if(params.downloadMethod.toLowerCase() == 'sra') {
      sra(makeIndex.out, accessions)
    }
    else if(params.downloadMethod.toLowerCase() == 'local') {
      local(makeIndex.out, accessions)
    }
}