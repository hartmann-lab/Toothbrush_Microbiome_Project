#!/bin/bash
#<add job submission details>

###
### Downloading the HMP-oral metagenome subset sequencing files
###

cd <path to Toothbrush_Microbiome_Project>
mkdir HMP_oral_subset
cd HMP_oral_subset

## load SRA toolkit
module load sratoolkit/2.8.1

## download raw read files for "centroid nearest neighbors" from the NCBI SRA (accession numbers match Table S3)
fastq-dump SRR1804378 --split-3 --gzip --skip-technical
fastq-dump SRR514849 --split-3 --gzip --skip-technical
fastq-dump SRR062389 --split-3 --gzip --skip-technical 
fastq-dump SRR1804274 --split-3 --gzip --skip-technical
fastq-dump SRR059875 --split-3 --gzip --skip-technical
fastq-dump SRR513813 --split-3 --gzip --skip-technical
fastq-dump SRR1804260 --split-3 --gzip --skip-technical
fastq-dump SRR1952532 --split-3 --gzip --skip-technical
fastq-dump SRR1031677 --split-3 --gzip --skip-technical
fastq-dump SRR1804718 --split-3 --gzip --skip-technical
fastq-dump SRR628262 --split-3 --gzip --skip-technical
fastq-dump SRR513802 --split-3 --gzip --skip-technical
fastq-dump SRR1565744 --split-3 --gzip --skip-technical
fastq-dump SRR514244 --split-3 --gzip --skip-technical
fastq-dump SRR2175813 --split-3 --gzip --skip-technical
fastq-dump SRR062389 --split-3 --gzip --skip-technical
fastq-dump SRR1564324 --split-3 --gzip --skip-technical
fastq-dump SRR1804818 --split-3 --gzip --skip-technical
fastq-dump SRR513376 --split-3 --gzip --skip-technical
fastq-dump SRR059905 --split-3 --gzip --skip-technical
fastq-dump SRR061249 --split-3 --gzip --skip-technical
fastq-dump SRR1952578 --split-3 --gzip --skip-technical
fastq-dump SRR060389 --split-3 --gzip --skip-technical
fastq-dump SRR1952524 --split-3 --gzip --skip-technical
fastq-dump SRR061257 --split-3 --gzip --skip-technical
fastq-dump SRR1803859 --split-3 --gzip --skip-technical
fastq-dump SRR1803853 --split-3 --gzip --skip-technical
fastq-dump SRR1952557 --split-3 --gzip --skip-technical
fastq-dump SRR533811 --split-3 --gzip --skip-technical
fastq-dump SRR1804292 --split-3 --gzip --skip-technical
fastq-dump SRR061251 --split-3 --gzip --skip-technical
fastq-dump SRR346700 --split-3 --gzip --skip-technical
fastq-dump SRR060041 --split-3 --gzip --skip-technical
fastq-dump SRR514204 --split-3 --gzip --skip-technical

## NOTE: repeat the kneaddata, metaphlan, and shortbred analysis described for the toothbrush samples on this dataset to generate outputs for the comparison analysis
