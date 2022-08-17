#!/usr/bin/env bash

set -euo pipefail
samtools rmdup -S bamfile out.bam
