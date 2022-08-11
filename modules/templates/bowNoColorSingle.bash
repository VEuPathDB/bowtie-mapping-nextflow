#!/usr/bin/env bash

bowtie2 \
  --rg-id EuP \
  --rg 'SM:TU114' \
  --rg 'PL:Illumina' \
  -x $index \
  -U ${readsFastq} \
  -S tmpOut.sam
samtools view -buS tmpOut.sam | samtools sort -o tmpOut.bam
