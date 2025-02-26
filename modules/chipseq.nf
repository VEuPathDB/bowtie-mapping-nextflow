#!/usr/bin/env nextflow
nextflow.enable.dsl=2

// pipeline input parameters
params.genomeFasta = "${projectDir}/data/genome.fasta"
params.csv = "${projectDir}/data/samplesheet.csv"
params.paired = false
params.index = "index"

process samtools {
    publishDir "${projectDir}/results", mode: 'copy'
    container "quay.io/biocontainers/samtools:1.20--h50ea8bc_0"
 
    input:
    path sampleSam
    tuple val(sample_id)

    output:
    path "${sample_id}.bam"
    path "${sample_id}.bam.bai"

    script:
    """
    samtools view -bS $sampleSam | samtools sort - -o "${sample_id}.bam"
    samtools index "${sample_id}.bam"
    """
}

  
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
  tuple val(sample_id), path(reads_1), path(reads_2)
  path bowtieDB
  val dbName

  output:
  val $sample_id
  path "${sample_id}.sam"
 

  script:
  """
  bowtie2 -x "${bowtieDB}/${dbName}"  -1 $reads_1 -2 $reads_2 -S "${sample_id}.sam}"
  """
}


process alignBowtie2SE {
  container "quay.io/biocontainers/bowtie2:2.5.4--he96a11b_5"  
  publishDir "$projectDir/results", mode: 'copy'
    
  input:  
  tuple val(sample_id), path(reads_1)
  path bowtieDB
  val dbName
    
  output:  
  val $sample_id
  path "${sample_id}.sam"  
    
  script:
  """
  bowtie2 -x "${bowtieDB}/${dbName}"  -U $reads_1 -S "${sample_id}.sam}"
  """
}


//process bamCompare {
//    container "quay.io/biocontainers/deeptools"
//    publishDir "$projectDir/results", mode: 'copy'
//    input: 
//	  path bamFile
//	  path controlBamFile
//	  val sample_id
//    output: 
//	  "${sample_id}_compare.bw"
//    	script:
//	"""
//	  deeptools bamCompare -b1 $bamFile -b2 $controlBamFile -o "${sample_id}_compare.bw"
//	"""

//process bamCoverage {   
//    container "quay.io/biocontainers/deeptools"
//    publishDir "$projectDir/results", mode: 'copy'
//    input:
//      path bamFile
//      val sample_id
//    output:
//      "${sample_id}.bw"
//    script:
//	"""
//      deeptools bamCoverage -b $bamFile  -o "${sample_id}.bw"
//	"""

//process callPeaks {
//    container "quay.io/biocontainers/macs3"
//    publishDir "$projectDir/results", mode: 'copy'
//    input:
//      path bamFile
//      path controlBamFile
//      val sample_id
//    output:
//      "${sample_id}_*.*
//    script:
//	"""
//      macs3 callpeak -t $bamFile -c $controlBamFile -n $sample_id
//	"""

workflow {

    // Run bowtie2-build on the genomeFasta file
    bowtieDB = bowtie2Index(params.genomeFasta, params.index)

    if (params.paired) {

    // Define the input channel from samplesheet.csv file, sample_id in reads list
      reads = Channel.fromPath(params.csv)
          .splitCsv(header: true)
          .map { row -> [row.sample, file(row.fastq_1), file(row.fastq_2),row.replicate,row.antibody,file(row.control),row.control_replicate] }

    samFilePE = alignBowtie2PE (reads, "$projectDir/data", params.index)
    bamFile = samtools(samFilePE, reads)
    }

    else {
// single end data, input channel with only one fastq and run alignBowtie2SE
      reads = Channel.fromPath(params.csv)
	.splitCsv(header: true)
        .map { row -> [row.sample, file(row.fastq_1)] }
    
    samFileSE = alignBowtie2SE (reads, "$projectDir/data", params.index)
    bamFile = samtools(samFileSE, reads)
    }

//    Peaks = callPeaks(bamFile, ControlBamFile, sample_id)

//    CompBigWig = bamCompare(bamFile, ControlBamFile, sample_id)

//    Coverage = bamCoverage(bamFile, sample_id)

}
