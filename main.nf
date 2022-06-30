nextflow.enable.dsl=2

process createIndex {

  input:
  path 'databaseFasta'

  output:
  path 'index.*'

  script:
  if(params.isColorspace)
    """
    bowtie-build databaseFasta index
    """
  else
    """
    bowtie2-build databaseFasta index
    """
}

process prepSRASingle {

  input:
  tuple val(genomeName), path(genomeReads) 

  output:
  path '*1.fastq'

  """
  gzip -d --force $genomeReads 
  """ 
}

process prepSRAPaired {

  input:
  tuple val(genomeName), path(genomeReads) 

  output:
  path '*1.fastq'
  path '*2.fastq'

  """
  gzip -d --force ${genomeReads[0]} 
  gzip -d --force ${genomeReads[1]} 
  """
}

process bowtieSingle {

  input:
  path indexfiles
  path '1.fastq'
  val index
  output:
  path '*.bam'

  script:
  if(params.isColorspace)
      """
      bowtie -f -C -a -S -n 3 --best --strata --sam-RG 'ID:EuP' --sam-RG 'SM:TU114' --sam-RG 'PL:Illumina' -x $index -1 1.fastq -Q $params.mateAQual > tmpOut.sam
      samtools view -buS tmpOut.sam | samtools sort -o tmpOut.bam
      """
  else
      """
      bowtie2 --rg-id EuP --rg 'SM:TU114' --rg 'PL:Illumina' -x $index -U $params.mateA -S tmpOut.sam
      samtools view -buS tmpOut.sam | samtools sort -o tmpOut.bam
      """
}

process bowtiePaired {

  input:
  path indexfiles
  path 'mateA'
  path 'mateB'
  val index
  output:
  path '*.bam'

  script:
  if(params.isColorspace)
      """
      bowtie -f -C -a -S -n 3 --best --strata --sam-RG 'ID:EuP' --sam-RG 'SM:TU114' --sam-RG 'PL:Illumina' -x $index -1 mateA --Q1 $params.mateAQual -2 mateB --Q2 $params.mateBQual > tmpOut.sam
      samtools view -buS tmpOut.sam | samtools sort -o tmpOut.bam
      """
  else
      """
      bowtie2 --rg-id EuP --rg 'SM:TU114' --rg 'PL:Illumina' -x $index -1 mateA -2 mateB -S tmpOut.sam
      samtools view -buS tmpOut.sam | samtools sort -o tmpOut.bam
      """
}

process PCRDuplicates {
  publishDir params.outputDir, mode: "copy"

  input:
  path 'bamfile'

  output:
  path 'out.*'

  script:
  if(params.removePCRDuplicates  && params.writeBedFile)
      """
      samtools rmdup -S bamfile out.bam
      samtools index out.bam
      bedtools bamutobed -i out.bam > output.bed
      """
  else if(params.removePCRDuplicates && !params.writeBedFile)
      """
      samtools rmdup -S bamfile out.bam
      """
  else if(!params.removePCRDuplicates && params.writeBedFile)
      """
      mv bamfile out.bam
      samtools index out.bam
      bedtools bamtobed -i out.bam > out.bed
      """
  else
      """
      mv bamfile out.bam
      """
}

workflow makeIndex {
  main:
    if(params.preconfiguredDatabase) {
      indexfiles = file(params.databaseFileDir + "/*.bt*")
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

workflow processing {
  take:
    indexfiles
    indexFileBasename
  main:
    if(params.isSingleEnd && params.fromSRA) {
      files = channel.fromSRA( params.sraID, apiKey: params.apiKey, protocol: "http" )
      seqs = prepSRASingle(files)
      bowtieSingle(indexfiles,seqs, indexFileBasename) | PCRDuplicates
    }
    else if(!params.isSingleEnd && params.fromSRA) {
      files = channel.fromSRA( params.sraID, apiKey: params.apiKey, protocol: "http" )
      seqs = prepSRAPaired(files)
      bowtiePaired(indexfiles, seqs[0], seqs[1], indexFileBasename) | PCRDuplicates
    }
    else if(params.isSingleEnd && !params.fromSRA) {
      bowtieSingle(indexfiles, params.mateA, indexFileBasename) | PCRDuplicates
    }
    else {
      bowtiePaired(indexfiles, params.mateA, params.mateB, indexFileBasename) | PCRDuplicates
    }
}

workflow { 
  main:
    makeIndex()
    processing(makeIndex.out)
}