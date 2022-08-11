#!/usr/bin/env bash

bowtie2 \
  --rg-id EuP \
  --rg 'SM:TU114' \
  --rg 'PL:Illumina' \
  -x $index \
  -1 ${sample}_1.fastq \
  -2 ${sample}_2.fastq \
  -S tmpOut.sam
samtools view -buS tmpOut.sam | samtools sort -o tmpOut.bam
