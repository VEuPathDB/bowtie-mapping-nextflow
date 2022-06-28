nextflow.enable.dsl=2

process createIndex {

  output:
  path 'index.*'

  """
  createIndex.pl --isColorspace $params.isColorspace --bowtieIndex $params.bowtieIndex --sampleName index
  """
}

process bowtieMappingPaired {
  input:
  path 'mateA.fq' 
  path 'mateB.fq'
  // Retaining names from the input into db_vch
  path indexfiles
  
  output:
  path 'runBowtieMapping.log' 
  path '*.bam'
  path '*.sam' optional true

  """
  runBowtieMapping.pl --ma mateA.fq --mb mateB.fq --bowtie2 /usr/bin/bowtie2 --bowtieIndex index --sampleName "sample" --isColorspace $params.isColorspace --removePCRDuplicates $params.removePCRDuplicates --extraBowtieParams $params.extraBowtieParams --deleteIntermediateFiles $params.deleteIntermediateFiles
  """
}

process bowtieMappingSingle {
  input:
  path 'mateA.fq'
  path 'mateB.fq'
  path 'mateA.fq.qual'
  // Retaining names from the input into db_vch
  path indexfiles
  
  output:
  path 'runBowtieMapping.log' 
  path '*.bam' 
  path '*.sam' 

  """
  runBowtieMapping.pl --ma mateA.fq --mb mateB.fq --bowtie2 /usr/bin/bowtie2 --bowtieIndex index --sampleName $params.sampleName --isColorspace $params.isColorspace --removePCRDuplicates $params.removePCRDuplicates --extraBowtieParams $params.extraBowtieParams --deleteIntermediateFiles $params.deleteIntermediateFiles
  """
}

workflow {
  if (params.preformattedDatabase == "false") {
  index_files = createIndex()
    if(params.singleEnd == "false" ) {
      mateA = channel.fromPath(params.mateA).splitFasta( by:1, file:true  )
      mateB = Channel.fromPath(params.mateB).splitFasta( by:1, file:true  )
      results = bowtieMappingPaired(mateA, mateB, index_files)
    }
    if(params.singleEnd == "true" ) {
      mateA = channel.fromPath(params.mateA).splitFasta( by:1, file:true  )
      mateB = channel.fromPath(params.mateA).splitFasta( by:1, file:true  )
      results = bowtieMappingSingle(mateA, mateB, params.mateAQual, index_files)
    }
  } else if (params.preformattedDatabase == "true") {
    index_files = file(params.databaseFileDir + "/*.bt2")
    if(params.singleEnd == "false" ) {
      mateA = channel.fromPath(params.mateA).splitFasta( by:1, file:true  )
      mateB = Channel.fromPath(params.mateB).splitFasta( by:1, file:true  )
      results = bowtieMappingPaired(mateA, mateB, index_files)
    }
    if(params.singleEnd == "true" ) {
      mateA = channel.fromPath(params.mateA).splitFasta( by:1, file:true  )
      mateB = channel.fromPath(params.mateA).splitFasta( by:1, file:true  )
      results = bowtieMappingSingle(mateA, mateB, params.mateAQual, index_files)
    }
  }
  results[0] | collectFile(storeDir: params.outputDir, name:params.logFile)
  results[1] | collectFile(storeDir: params.outputDir, name:params.bamFile)
  results[2] | collectFile(storeDir: params.outputDir, name:params.samFile)  
}
