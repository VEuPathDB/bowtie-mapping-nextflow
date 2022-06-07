nextflow.enable.dsl=1
mateA_qch = Channel.fromPath(params.mateA).splitFasta( by:1, file:true  )
mateB_qch = Channel.fromPath(params.mateB).splitFasta( by:1, file:true  )
db_vch = Channel.value()
process createDatabase {
    output:
    path 'index.*' into db_vch
    """
    createIndex.pl --isColorspace $params.isColorspace --bowtieIndex $params.bowtieIndex --sampleName index
    """
}

process blastSimilarity {
   input:
   path 'mateA.fq' from mateA_qch
   path 'mateB.fq' from mateB_qch
   // Retaining names from the input into db_vch
   path indexFiles from db_vch
   output:
   path 'runBowtieMapping.log' into log_qch
   path '*.bam' into bam_qch
   path '*.sam' optional true into sam_qch
   """
   runBowtieMapping.pl --ma mateA.fq --mb mateB.fq --bowtie2 /usr/bin/bowtie2 --bowtieIndex index --sampleName $params.sampleName --isColorspace $params.isColorspace --removePCRDuplicates $params.removePCRDuplicates --extraBowtieParams $params.extraBowtieParams --deleteIntermediateFiles $params.deleteIntermediateFiles
   """
}

bamResults = bam_qch.collectFile(storeDir: params.outputDir, name: params.bamFile)
logResults = log_qch.collectFile(storeDir: params.outputDir, name: params.logFile)
samResults = sam_qch.collectFile(storeDir: params.outputDir, name: params.samFile)

