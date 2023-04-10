#!/usr/bin/env bash

set -euo pipefail

bowtie2 --rg-id EuP --rg 'SM:TU114' --rg 'PL:Illumina' -x $index -1 $mateA -2 $mateB -S tmpOut.sam $extraBowtieParams
samtools view -buS tmpOut.sam | samtools sort -o tmpOut.bam

