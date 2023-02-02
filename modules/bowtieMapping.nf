#!/usr/bin/env nextflow
nextflow.enable.dsl=2


process createIndex {
  input:
    path databaseFasta
    val isColorSpace
    
  output:
    path 'index.*'

  script:
    template 'createIndex.bash'
}


process downloadFiles {
  input:
    val id
    
  output:
    tuple val(id), path("${id}**.fastq")

  script:
    template 'downloadFiles.bash'
}


process bowtie {
  input:
    path indexfiles
    tuple val(sample), path(readsFastq)
    val index
    val isColorSpace
    val isSingleEnd

  output:
    path '*.bam'

  script:
    template 'bowtie.bash'
}


process PCRDuplicates {
  publishDir params.outputDir, mode: "copy"

  input:
    path bamfile
    val removePCRDuplicates
    val writeBedFile

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
      indexfiles = createIndex(params.databaseFasta, params.isColorSpace)
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
    bowtieResults = bowtie( indexfiles, files, indexFileBasename, params.isColorSpace, params.isSingleEnd ) 
    PCRDuplicates(bowtieResults, params.removePCRDuplicates, params.writeBedFile)
}

workflow local {
  take:
    indexfiles
    indexFileBasename
    files

  main:
    bowtieResults = bowtie( indexfiles, files, indexFileBasename, params.isColorSpace, params.isSingleEnd )
    PCRDuplicates(bowtieResults, params.removePCRDuplicates, params.writeBedFile)
}

workflow bowtieMapping {
  take:
    accessions

  main:
    makeIndex()
    if(params.downloadMethod.toLowerCase() == 'sra') {
      sra(makeIndex.out, accessions)
    }
    else if(params.downLoadMethod.toLowerCase() == 'local') {
      local(makeIndex.out, accessions)
    }
}