#!/usr/bin/env nextflow
nextflow.enable.dsl=2

// pipeline input parameters
params.genomeFasta = "${projectDir}/data/genome.fasta"
params.csv = "${projectDir}/data/samplesheet.csv"
params.paired = false
params.index = "index"
// params.fq1 = "${projectDir}/data/*_1.fastq.gz"
// params.fq2 = "${projectDir}/data/*_2.fastq.gz"
params.reads = '${projectDir}/data/*_{1,2}.fastq.gz'

process bowtie2Index {
    container "quay.io/biocontainers/bowtie2:2.5.4--he96a11b_5"
    publishDir "$projectDir/data", mode: 'copy'

    input:
    path genomeFasta
    val  index

    output:
    path "index.*"

    script:
    """
    bowtie2-build $genomeFasta $index
    """
}


process alignBowtie2PE {
  container "quay.io/biocontainers/bowtie2:2.5.4--he96a11b_5"
  publishDir "$projectDir/results", mode: 'copy'

  input:
  tuple val(sample_id), path(reads)
  path bowtieDB
  val dbName

  output:
  path "${sample_id}.sam"
 

  script:
  """
  bowtie2 -x "${bowtieDB}/${dbName}" -p 8  -1 ${reads[0]} -2 ${reads[1]} -S "${sample_id}.sam}"
  """
}


workflow {

    // Run bowtie2-build on the genomeFasta file
    bowtieDB = bowtie2Index(params.genomeFasta, params.index)

//  Read samplesheet rows into samplename,fq_1,fq_2  

    reads = Channel.fromFilePairs(params.reads, flat: true)

//    reads = Channel.fromPath(params.csv)
//        | splitCsv(header: ['group','fastq_1,'fastq_2'])
//        | map { row-> tuple(row.sampleId, file(row.fastq_1), file(row.fastq_2)) }


    // Run bowtie2 on the read pairs
    samFilePE = alignBowtie2PE (reads, "$projectDir/data", params.index)

}


