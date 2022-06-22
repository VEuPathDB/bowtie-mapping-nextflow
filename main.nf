nextflow.enable.dsl=2

process createIndex {
    output:
    path 'index.*'
    """
    createIndex.pl --isColorspace $params.isColorspace --bowtieIndex $params.bowtieIndex --sampleName index
    """
}

process bowtieMappingPaired {
   publishDir params.outputDir, mode: 'copy'
   input:
   path 'mateA.fq' 
   path 'mateB.fq'
   path indexfiles
   // Retaining names from the input into db_vch
   output:
   path 'runBowtieMapping.log' 
   path '*.bam'
   path '*.sam' optional true
   """
   runBowtieMapping.pl --ma mateA.fq --mb mateB.fq --bowtie2 /usr/bin/bowtie2 --bowtieIndex index --sampleName $params.sampleName --isColorspace $params.isColorspace --removePCRDuplicates $params.removePCRDuplicates --extraBowtieParams $params.extraBowtieParams --deleteIntermediateFiles $params.deleteIntermediateFiles
   """
}

process bowtieMappingSingle {
   publishDir params.outputDir, mode: 'copy'
   input:
   path 'mateA.fq'
   path 'mateB.fq'
   path 'mateA.fq.qual'
   path indexfiles
   // Retaining names from the input into db_vch
   output:
   path 'runBowtieMapping.log' 
   path '*.bam' 
   path '*.sam' 
   """
   runBowtieMapping.pl --ma mateA.fq --mb mateB.fq --bowtie2 /usr/bin/bowtie2 --bowtieIndex index --sampleName $params.sampleName --isColorspace $params.isColorspace --removePCRDuplicates $params.removePCRDuplicates --extraBowtieParams $params.extraBowtieParams --deleteIntermediateFiles $params.deleteIntermediateFiles
   """
}

workflow {
  db_vch = Channel.value()
  index_files = createIndex()
  if(params.singleEnd == "false" ) {
    mateA = channel.fromPath(params.mateA).splitFasta( by:1, file:true  )
    mateB = Channel.fromPath(params.mateB).splitFasta( by:1, file:true  )
    bowtieMappingPaired(mateA, mateB, index_files)
  }
  if(params.singleEnd == "true" ) {
    mateA = channel.fromPath(params.mateA).splitFasta( by:1, file:true  )
    mateB = channel.fromPath(params.mateA).splitFasta( by:1, file:true  )
    bowtieMappingSingle(mateA, mateB, params.mateAQual, index_files)
}
}
