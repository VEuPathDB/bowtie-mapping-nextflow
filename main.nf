#!/usr/bin/env nextflow
nextflow.enable.dsl=2

include { bowtieMapping } from './workflows/bowtieMapping.nf'

include { coverageAndPeaks } from './workflows/coverageAndPeaks.nf'

workflow {

  bam = bowtieMapping()

  coverageAndPeaks(bam)
  
}
