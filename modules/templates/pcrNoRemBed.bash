#!/usr/bin/env bash

set -euo pipefail
mv bamfile out.bam
samtools index out.bam
bedtools bamtobed -i out.bam > out.bed
