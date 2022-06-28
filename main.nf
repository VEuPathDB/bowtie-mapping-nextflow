nextflow.enable.dsl=2

process createIndex {
  output:
  path 'index.*'
  """
  createIndex.pl --isColorspace $params.isColorspace --bowtieIndex $params.databaseFasta --sampleName index
  """
}

process bowtie {
  input:
  path indexfiles
  output:
  path '*.bam'
  script:
  if(params.isColorspace == "true" && params.isSingleEnd == "true" )
      """
      bowtie -f -C -a -S -n 3 --best --strata --sam-RG 'ID:EuP' --sam-RG 'SM:TU114' --sam-RG 'PL:Illumina' -x index -1 $params.mateA -Q $params.mateAQual > tmpOut.sam
      samtools view -buS tmpOut.sam | samtools sort -o tmpOut.bam
      """
  else if(params.isColorspace == "true" && params.isSingleEnd == "false")
      """
      bowtie -f -C -a -S -n 3 --best --strata --sam-RG 'ID:EuP' --sam-RG 'SM:TU114' --sam-RG 'PL:Illumina' -x index -1 $params.mateA --Q1 $params.mateAQual -2 $params.mateB --Q2 $params.mateBQual > tmpOut.sam
      samtools view -buS tmpOut.sam | samtools sort -o tmpOut.bam
      """
  else if(params.isColorspace == "false" && params.isSingleEnd == "true")
      """
      bowtie2 --rg-id EuP --rg 'SM:TU114' --rg 'PL:Illumina' -x index -U $params.mateA -S tmpOut.sam
      samtools view -buS tmpOut.sam | samtools sort -o tmpOut.bam
      """
  else if(params.isColorspace == "false" && params.isSingleEnd == "false")
      """
      bowtie2 --rg-id EuP --rg 'SM:TU114' --rg 'PL:Illumina' -x index -1 $params.mateA -2 $params.mateB -S tmpOut.sam
      samtools view -buS tmpOut.sam | samtools sort -o tmpOut.bam
      """
}

process PCRDuplicates {
  publishDir params.outputDir, mode: "copy", saveAs: { filename -> params.bamFile }
  input:
  path 'bamfile'
  output:
  path 'out.bam'
  script:
  if(params.removePCRDuplicates == "true")
      """
      samtools rmdup -S bamfile out.bam
      """
  else if(params.removePCRDuplicates == "false")
      """
      mv bamfile out.bam
      """
}

workflow {
  if(params.preconfiguredDatabase == "true") {
    indexfiles = file(params.databaseFileDir + "/*.bt*")
    bowtie(indexfiles) | PCRDuplicates
  }
  else if(params.preconfiguredDatabase == "false") {
    indexfiles = createIndex()
    bowtie(indexfiles) | PCRDuplicates
  }
}