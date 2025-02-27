#!/usr/bin/env nextflow
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

//--------------------------------------------------------------------------
// Param Checking
//--------------------------------------------------------------------------

if(params.preconfiguredDatabase) {

  if (!params.databaseFileDir) {
    throw new Exception("Missing params.databaseFileDir")
  }
  else if (!params.indexFileBasename) {
    throw new Exception("Missing params.indexFileBasename")
  }
}

if(!params.preconfiguredDatabase && !params.databaseFasta) {
  throw new Exception("Missing params.databaseFasta")
}

if(params.downloadMethod.toLowerCase() == 'sra') {
  accessions = fetchRunAccessions( params.input )
}
else if(params.downloadMethod.toLowerCase() == 'local' && !params.hasPairedReads) {
  file = channel.fromPath([params.inputDir + '*_1.fastq'])
}
else if(params.downloadMethod.toLowerCase() == 'local' && params.hasPairedReads) {
  file = Channel.fromFilePairs([params.inputDir + '*_{1,2}.fastq'])
}
else {
  throw new Exception("Invalid value for params.downloadMethod")
}

//--------------------------------------------------------------------------
// Includes
//--------------------------------------------------------------------------

include { bowtieMapping } from './modules/bowtieMapping.nf'

//--------------------------------------------------------------------------
// Main Workflow
//--------------------------------------------------------------------------

workflow {
  if(params.downloadMethod.toLowerCase() == 'sra')
    bowtieMapping(accessions)
  else {
    bowtieMapping(file)
  }
}