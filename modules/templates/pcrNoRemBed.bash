#!/usr/bin/env bash

mv bamfile out.bam
samtools index out.bam
bedtools bamtobed -i out.bam > out.bed
