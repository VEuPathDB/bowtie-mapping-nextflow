#!/usr/bin/env nextflow
nextflow.enable.dsl=2


process createIndex {
  input:
    path databaseFasta

  output:
    path 'index.*'

  script:
    if(params.isColorspace)
      template 'createIndexColor.bash'
    else
      template 'createIndexNoColor.bash'
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

  output:
    path '*.bam'

  script:
    if(params.isColorspace && params.isSingleEnd)
      template 'bowColorSingle.bash'
    else if(!params.isColorspace && params.isSingleEnd)
      template 'bowNoColorSingle.bash'
    else if(params.isColorspace && !params.isSingleEnd)
      template 'bowColorPaired.bash'
    else if(!params.isColorspace && !params.isSingleEnd)
      template 'bowNoColorPaired.bash'
}


process PCRDuplicates {
  publishDir params.outputDir, mode: "copy"

  input:
    path bamfile

  output:
    path 'out.*'

  script:
    if(params.removePCRDuplicates  && params.writeBedFile)
      template 'pcrRemBed.bash'
    else if(params.removePCRDuplicates && !params.writeBedFile)
      template 'pcrRemNoBed.bash'
    else if(!params.removePCRDuplicates && params.writeBedFile)
      template 'pcrNoRemBed.bash'
    else
      template 'pcrNoRemNoBed.bash'
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
    bowtie( indexfiles, files, indexFileBasename ) | PCRDuplicates
}

workflow local {
  take:
    indexfiles
    indexFileBasename
    files

  main:
    bowtie( indexfiles, files, indexFileBasename ) | PCRDuplicates
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