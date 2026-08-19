# bowtie-mapping-nextflow

A Nextflow pipeline that aligns short reads to a reference genome with Bowtie2 and derives coverage tracks and, for ChIP-seq histone modification data, called peaks.

## Overview

This pipeline covers the common short-read mapping and signal-generation steps used across several VEuPathDB genomic data types: origins-of-replication profiling, splice-site read mapping, and ChIP-seq (histone modification, transcription factor binding, MNase, FAIRE, and DNase experiments). Samples are aligned to a reference genome with Bowtie2, optionally deduplicated, and then compared against a paired reference/input sample with `deeptools bamCompare` to produce ratio coverage tracks. For ChIP-seq datasets, the pipeline additionally builds HOMER tag directories (or, for MNase data, runs DANPOS) to generate coverage bigWigs, and — for histone modification experiments with a matched input sample — calls peaks with HOMER and indexes the results as a tabix-indexed GFF for loading into VEuPathDB's genome browser and study infrastructure.

## Requirements

- [Nextflow](https://www.nextflow.io/) (DSL2)
- Docker or Singularity/Apptainer — processes run in the `veupathdb/bowtiemapping:1.0.0` image plus `quay.io/biocontainers/deeptools`, `quay.io/biocontainers/ucsc-bedgraphtobigwig`, `quay.io/biocontainers/ucsc-wigtobigwig`, `biocontainers/danpos`, and `biocontainers/tabix` (select the engine via the `docker` or `singularity` profile/config in `conf/`)

## Usage

```
nextflow run VEuPathDB/bowtie-mapping-nextflow -r main \
  --input /path/to/input \
  --genome /path/to/genome.fasta \
  --datasetType chipSeq \
  --experimentType histonemod \
  --inputFileType fastq \
  --hasPairedReads true \
  --removePCRDuplicates true \
  --outputDir /path/to/output \
  -C conf/docker.config \
  -resume
```

The pipeline has a single (default) entry point, made up of two chained sub-workflows:

1. **`bowtieMapping`** — reads `params.input/<samplesheetFileName>`, builds a Bowtie2 index from `params.genome`, aligns each sample's reads with `bowtie2`, optionally removes PCR duplicates (`samtools rmdup`), and indexes the resulting BAM.
2. **`coverageAndPeaks`** — pairs each sample with the reference/input sample named in its `ref` field, runs `deeptools bamCompare` to produce a ratio bigWig per pair, and then — when `datasetType` is `chipSeq` — either runs DANPOS (for `experimentType = mnase`) or builds HOMER tag directories and bedGraph/bigWig coverage tracks (for `histonemod`, `tfbinding`, `faire`, `dnase`). For `experimentType = histonemod`, peaks are additionally called with HOMER `findPeaks` against the paired input sample and indexed as a tabix-compressed GFF.

### Samplesheet

`params.input/<samplesheetFileName>` is a CSV (header row skipped) with columns: sample ID, read file 1, read file 2 (if paired), and reference/input sample ID (the `ref` column used to pair a sample with its control for `bamCompare` and peak calling). Read file paths are resolved relative to `params.input` unless given as absolute paths.

## Key Parameters

| Parameter | Default | Description |
|---|---|---|
| `input` | `$launchDir/input` | Directory containing the samplesheet and (relative) read files |
| `samplesheetFileName` | `samplesheet.csv` | Samplesheet filename, resolved under `input` |
| `genome` | `$launchDir/input/genome.fasta` | Reference genome FASTA used for the Bowtie2 index and chromosome sizes |
| `datasetType` | — | Data type being processed: `originsOfReplication`, `spliceSites`, or `chipSeq` |
| `experimentType` | — | For `chipSeq` datasets: `histonemod`, `tfbinding`, `mnase`, `faire`, or `dnase` |
| `inputFileType` | — | `fastq` or `fasta`, passed to `bowtie2` as `-q`/`-f` |
| `hasPairedReads` | `true` | Whether samples have paired-end (`-1`/`-2`) or single-end (`-U`) reads |
| `removePCRDuplicates` | `false` | Whether to run `samtools rmdup` on alignments before downstream analysis |
| `saveAlignments` | `false` | Publish indexed BAM/BAI files to `outputDir` |
| `saveCoverage` | `false` | Publish bigWig coverage tracks to `outputDir` |
| `profileSetName` | — | Profile set name written into browser track/peak config output for the study loader |
| `gffFileName` | `output.gff` | Filename used for the sorted, tabix-indexed peaks GFF |
| `outputDir` | `$launchDir/output` | Directory published output files are written to |

## Output

Depending on `datasetType`/`experimentType` and the `save*` flags, published to `outputDir`:

- Indexed BAM alignments (`<sample>.bam`, `<sample>.bam.bai`) when `saveAlignments` is `true`
- Ratio coverage bigWigs from `bamCompare` (`<sample>_vs_<ref>.bw`) and their browser track config, collected into `metadata_ratios`
- For ChIP-seq data with `saveCoverage` enabled: unlogged coverage bigWigs (HOMER-derived bedGraph→bigWig for most experiment types, or DANPOS smoothed signal→bigWig for MNase), with browser track configs collected into `metadata_unlogged`
- For `histonemod` experiments: called peaks per sample (`<sample>_peaks.txt`), a merged tabix-indexed peaks GFF (`<gffFileName>.gz` / `.gz.tbi`), and a collected `insert_study_results` config file for loading peak results
