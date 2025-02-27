#!/usr/bin/env bash

set -euo pipefail

bowtie2 --rg-id EuP --rg 'SM:TU114' --rg 'PL:Illumina' -x $index -1 ${sampleName}_1.fastq -2 ${sampleName}_2.fastq -S tmpOut.sam $extraBowtieParams
samtools view -buS tmpOut.sam | samtools sort -o tmpOut.bam

