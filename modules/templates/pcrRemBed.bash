#!/usr/bin/env bash

set -euo pipefail
samtools rmdup -S $bamfile out.bam
samtools index out.bam
bedtools bamutobed -i out.bam > output.bed
