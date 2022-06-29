nextflow.enable.dsl=2

process createIndex {

  output:
  path 'index.*'

  """
  createIndex.pl --isColorspace $params.isColorspace --bowtieIndex $params.databaseFasta --sampleName index
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

  output:
  path '*.bam'

  script:
  if(params.isColorspace == "true" && params.isSingleEnd == "true" )
      """
      bowtie -f -C -a -S -n 3 --best --strata --sam-RG 'ID:EuP' --sam-RG 'SM:TU114' --sam-RG 'PL:Illumina' -x index -1 1.fastq -Q $params.mateAQual > tmpOut.sam
      samtools view -buS tmpOut.sam | samtools sort -o tmpOut.bam
      """
  else if(params.isColorspace == "false" && params.isSingleEnd == "true")
      """
      bowtie2 --rg-id EuP --rg 'SM:TU114' --rg 'PL:Illumina' -x index -U $params.mateA -S tmpOut.sam
      samtools view -buS tmpOut.sam | samtools sort -o tmpOut.bam
      """
}

process bowtiePaired {

  input:
  path indexfiles
  path 'mateA'
  path 'mateB'

  output:
  path '*.bam'

  script:
  if(params.isColorspace == "true" && params.isSingleEnd == "false")
      """
      bowtie -f -C -a -S -n 3 --best --strata --sam-RG 'ID:EuP' --sam-RG 'SM:TU114' --sam-RG 'PL:Illumina' -x index -1 mateA --Q1 $params.mateAQual -2 mateB --Q2 $params.mateBQual > tmpOut.sam
      samtools view -buS tmpOut.sam | samtools sort -o tmpOut.bam
      """
  else if(params.isColorspace == "false" && params.isSingleEnd == "false")
      """
      bowtie2 --rg-id EuP --rg 'SM:TU114' --rg 'PL:Illumina' -x index -1 mateA -2 mateB -S tmpOut.sam
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

workflow makeIndex {
  main:
    if(params.preconfiguredDatabase == "true") {
      indexfiles = file(params.databaseFileDir + "/*.bt*")
    }
    else if(params.preconfiguredDatabase == "false") {
      indexfiles = createIndex()
    }
  emit:
      indexfiles
}

workflow processing {
  take: indexfiles
  main:
    if(params.isSingleEnd == "true" && params.fromSRA == "true") {
      files = channel.fromSRA( params.sraID, apiKey: params.apiKey, protocol: "http" )
      seqs = prepSRASingle(files)
      bowtieSingle(indexfiles,seqs) | PCRDuplicates
    }
    else if(params.isSingleEnd == "false" && params.fromSRA == "true") {
      files = channel.fromSRA( params.sraID, apiKey: params.apiKey, protocol: "http" )
      seqs = prepSRAPaired(files)
      bowtiePaired(indexfiles, seqs[0], seqs[1]) | PCRDuplicates
    }
    else if(params.isSingleEnd == "true" && params.fromSRA == "false") {
      bowtieSingle(indexfiles, params.mateA, params.mateB) | PCRDuplicates
    }
    else if(params.isSingleEnd == "false" && params.fromSRA == "false") {
      bowtiePaired(indexfiles, params.mateA, params.mateB) | PCRDuplicates
    }
}

workflow { 
  main:
    makeIndex()
    processing(makeIndex.out)
}