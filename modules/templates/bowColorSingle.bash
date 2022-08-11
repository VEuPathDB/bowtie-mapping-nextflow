#!/usr/bin/env bash

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
