#!/usr/bin/env bash

samtools rmdup -S bamfile out.bam
samtools index out.bam
bedtools bamutobed -i out.bam > output.bed
