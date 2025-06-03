#!/usr/bin/env nextflow
nextflow.enable.dsl=2



process chromSizes {
    container = 'veupathdb/bowtiemapping:1.0.0'

    input:
    path(genome)

    output:
    path("chrom.sizes")

    script:
    """
    samtools faidx $genome
    cut -f 1,2 ${genome}.fai >chrom.sizes
    """
}


process indexGff {
  container = 'biocontainers/tabix:v1.9-11-deb_cv1'

  publishDir params.outputDir, mode: 'copy'

  input:
    path gff

  output:
    path '*.gz'
    path '*.gz.tbi'

  script:
  """
  sort -k1,1 -k4,4n $gff > ${params.gffFileName}
  bgzip ${params.gffFileName}
  tabix -p gff ${params.gffFileName}.gz
  """
}
