#!/usr/bin/env bash

set -euo pipefail

if [ "$removePCRDuplicates" = true ] && [ "$writeBedFile" = true ]; then

    samtools rmdup -S $bamfile out.bam
    samtools index out.bam
    bedtools bamutobed -i out.bam > output.bed    

elif [ "$removePCRDuplicates" = true ] && [ "$writeBedFile" != true ]; then

    samtools rmdup -S $bamfile out.bam    

elif [ "$removePCRDuplicates" != true ] && [ "$writeBedFile" = true ]; then

    mv $bamfile out.bam
    samtools index out.bam
    bedtools bamtobed -i out.bam > out.bed
    
else
    
    mv $bamfile out.bam
    
fi
