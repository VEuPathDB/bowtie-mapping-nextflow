#!/usr/bin/env nextflow
nextflow.enable.dsl=2

include { bamCompare } from '../modules/local/deeptools.nf'

include { homer2Bed } from '../modules/local/chipSeq.nf'

include { bedgraph2bigwig; wigToBigWig } from '../modules/local/ucsc.nf'

include { chromSizes } from '../modules/local/util.nf'

include { danpos } from '../modules/local/danpos.nf'

workflow coverage {
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
            wigToBigWig(danpos.out, chromSizesRes)
        }
        else {
            bed = homer2Bed(bam)
            bedgraph2bigwig(bed, chromSizesRes)
        }
    }
}
