#!/usr/bin/env bash

set -euo pipefail

if [ "$hasPairedReads" != true ]; then

    bowtie2 --rg-id EuP --rg 'SM:TU114' --rg 'PL:Illumina' -x $index -U ${readsFastq} -S tmpOut.sam $extraBowtieParams
    samtools view -buS tmpOut.sam | samtools sort -o tmpOut.bam

elif [ "$hasPairedReads" = true ]; then

    bowtie2 --rg-id EuP --rg 'SM:TU114' --rg 'PL:Illumina' -x $index -1 ${sample}_1.fastq -2 ${sample}_2.fastq -S tmpOut.sam $extraBowtieParams
    samtools view -buS tmpOut.sam | samtools sort -o tmpOut.bam
    
fi
