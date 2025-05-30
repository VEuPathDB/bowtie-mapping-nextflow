#!/usr/bin/env nextflow
nextflow.enable.dsl=2

include { bowtieMapping } from './workflows/bowtieMapping.nf'

include { coverage } from './workflows/coverage.nf'

workflow {

  bam = bowtieMapping()

  if(params.saveCoverage) {
    coverage(bam)
  }
  
  
}
