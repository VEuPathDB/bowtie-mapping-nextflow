#!/usr/bin/env nextflow
nextflow.enable.dsl=2

include { bamCompare } from '../modules/local/deeptools.nf'

include { homer2Bed; homerMakeTagDir; homerFindPeaks_histonemod } from '../modules/local/chipSeq.nf'

include { bedgraph2bigwig; wigToBigWig } from '../modules/local/ucsc.nf'

include { chromSizes; indexGff } from '../modules/local/util.nf'

include { danpos } from '../modules/local/danpos.nf'

workflow coverageAndPeaks {
    take:
    bam

    main:

    //  make the cross product of pairs and then filter by the ref field
    pairs = bam.combine(bam).filter { v -> v[0].ref && v[0].ref == v[3].id }

    bamCompare(pairs)

    chromSizesRes = chromSizes(params.genome)

    if(params.datasetType == "chipSeq") {

        if(params.experimentType ==  "mnase") {
            danpos(bam)

            if(params.saveCoverage) {
                wigToBigWig(danpos.out, chromSizesRes)
            }
        }
        else {

            homerTagDir = homerMakeTagDir(bam)
            if(params.saveCoverage) {
                bed = homer2Bed(homerTagDir)
                bedgraph2bigwig(bed, chromSizesRes)
            }

            // only call peaks for histonemod experiments and samples with input
            if(params.experimentType == "histonemod") {
                //homerTagDir.combine(homerTagDir).view()
                homerTagDirPairs = homerTagDir.combine(homerTagDir).filter { v -> v[0].ref && v[0].ref == v[2].id }
                homerPeaks = homerFindPeaks_histonemod(homerTagDirPairs)

                // collect gff files and index
                indexGff(homerPeaks.map { v -> v[1] }.collectFile())

                // collect the config file for the loader
                homerPeaks.map { v -> v[3] }.collectFile(name: "insert_study_results", storeDir: params.outputDir, keepHeader: true, skip: 1)
            }
        }
    }
}
