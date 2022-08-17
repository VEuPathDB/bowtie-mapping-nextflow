#!/usr/bin/env bash

set -euo pipefail
bowtie2 \
  --rg-id EuP \
  --rg 'SM:TU114' \
  --rg 'PL:Illumina' \
  -x $index \
  -U ${readsFastq} \
  -S tmpOut.sam
samtools view -buS tmpOut.sam | samtools sort -o tmpOut.bam
