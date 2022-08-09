nextflow.enable.dsl=2

import nextflow.splitter.CsvSplitter


def fetchRunAccessions( tsv ) {
    def splitter = new CsvSplitter().options( header:true, sep:'\t' )
    def reader = new BufferedReader( new FileReader( tsv ) )
    splitter.parseHeader( reader )
    List<String> run_accessions = []
    Map<String,String> row
    while( row = splitter.fetchRecord( reader ) ) {
       run_accessions.add( row['run_accession'] )
    }
    return run_accessions
}


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


process downloadFiles {
  input:
    val id
    
  output:
    tuple val(id), path("${id}**.fastq")

  """
  fasterq-dump --split-3 ${id}
  """ 
}


process bowtie {
  input:
    path indexfiles
    tuple val(sample), path(readsFastq)
    val index

  output:
    path '*.bam'

  script:
    if(params.isColorspace && params.isSingleEnd)
      """
      bowtie \
        -f -C -a -S -n 3 \
        --best \
        --strata \
        --sam-RG 'ID:EuP' \
        --sam-RG 'SM:TU114' \
        --sam-RG 'PL:Illumina' \
        -x $index \
        -1 ${readsFastQ} \
        -Q $params.mateAQual > tmpOut.sam
      samtools view -buS tmpOut.sam | samtools sort -o tmpOut.bam
      """
    else if(!params.isColorspace && params.isSingleEnd)
      """
      bowtie2 \
        --rg-id EuP \
        --rg 'SM:TU114' \
        --rg 'PL:Illumina' \
        -x $index \
        -U ${readsFastq} \
        -S tmpOut.sam
      samtools view -buS tmpOut.sam | samtools sort -o tmpOut.bam
      """
    else if(params.isColorspace && !params.isSingleEnd)
      """
      bowtie \
        -f -C -a -S -n 3 \
        --best \
        --strata \
        --sam-RG 'ID:EuP' \
        --sam-RG 'SM:TU114' \
        --sam-RG 'PL:Illumina' \
        -x $index \
        -1 ${sample}_1.fastq \
        --Q1 $params.mateAQual \
        -2 ${sample}_2.fastq \
        --Q2 $params.mateBQual > tmpOut.sam
      samtools view -buS tmpOut.sam | samtools sort -o tmpOut.bam
      """
    else if(!params.isColorspace && !params.isSingleEnd)
      """
      bowtie2 \
      --rg-id EuP \
      --rg 'SM:TU114' \
      --rg 'PL:Illumina' \
      -x $index \
      -1 ${sample}_1.fastq \
      -2 ${sample}_2.fastq \
      -S tmpOut.sam
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
    if(params.downloadMethod == 'sra') {
      accessions = fetchRunAccessions( params.studyIdFile )
      ids = Channel.fromList(accessions)
      files = downloadFiles(ids)
    }
    else if(params.downloadMethod == 'local') {
      files = Channel.fromFilePairs(params.localFileLocation + "/*_{1,2}.fa*")
    }
    bowtie(indexfiles, files, indexFileBasename) | PCRDuplicates
}

workflow { 
  main:
    makeIndex()
    processing(makeIndex.out)
}