#!/usr/bin/env nextflow
nextflow.enable.dsl=2

include { bowtieMapping } from './workflows/bowtieMapping.nf'

include { bamCompare } from './modules/local/deeptools.nf'


workflow {

  bam = bowtieMapping()

    
  if(params.experimentType == "originsOfReplication") {
    reference = bam.filter { v -> v[0].id == params.referenceSample }
    nonReference = bam.filter { v -> v[0].id != params.referenceSample }

    bamCompare(nonReference, reference.first())
  }

  if(params.experimentType == "chipSeq") {
  }

  if(params.experimentType == "splicedLeader") {
  }

  if(params.experimentType == "polyA") {
  }



  
  
}
