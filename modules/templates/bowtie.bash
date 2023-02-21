#!/usr/bin/env bash

set -euo pipefail

if [ "$isColorSpace" = true ] && [ "$hasPairedReads" != true ]; then

    bowtie -f -C -a -S -n 3 --best --strata --sam-RG 'ID:EuP' --sam-RG 'SM:TU114' --sam-RG 'PL:Illumina' -x $index -1 ${readsFastq} -Q $params.mateAQual > tmpOut.sam
    samtools view -buS tmpOut.sam | samtools sort -o tmpOut.bam

elif [ "$isColorSpace" != true ] && [ "$hasPairedReads" != true ]; then

    bowtie2 --rg-id EuP --rg 'SM:TU114' --rg 'PL:Illumina' -x $index -U ${readsFastq} -S tmpOut.sam
    samtools view -buS tmpOut.sam | samtools sort -o tmpOut.bam

elif [ "$isColorSpace" = true ] && [ "$hasPairedReads" = true ]; then

    bowtie -f -C -a -S -n 3 --best --strata --sam-RG 'ID:EuP' --sam-RG 'SM:TU114' --sam-RG 'PL:Illumina' \
        -x $index -1 ${sample}_1.fastq --Q1 $params.mateAQual -2 ${sample}_2.fastq --Q2 $params.mateBQual > tmpOut.sam
    samtools view -buS tmpOut.sam | samtools sort -o tmpOut.bam
    
elif [ "$isColorSpace" != true ] && [ "$hasPairedReads" = true ]; then

    bowtie2 --rg-id EuP --rg 'SM:TU114' --rg 'PL:Illumina' -x $index -1 ${sample}_1.fastq -2 ${sample}_2.fastq -S tmpOut.sam
    samtools view -buS tmpOut.sam | samtools sort -o tmpOut.bam
    
fi
