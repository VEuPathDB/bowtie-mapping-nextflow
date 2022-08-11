#!/usr/bin/env bash

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
