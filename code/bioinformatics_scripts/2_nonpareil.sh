#!/bin/bash
#<add job submission details>

###
### NONPAREIL for metagenome sequence coverage
###

## navigate to cleaned reads directory; unzip files
cd <path to Toothbrush_Microbiome_Project/tb_clean_reads>
gunzip *fastq.gz
mkdir NP_coverage

## load work environment
module load anaconda3
source activate nonpareil

## use loop to run the redundancy analysis
files_knead=`ls | grep "_1.fastq"`
for reads in $files_knead 
	do out=`echo $reads | cut -d "_" -f1`_np
	nonpareil -s $reads -T kmer -f fastq -b NP_coverage/$out
done

## zip files
gzip *fastq

## NOTE: the output 'npo' files are used for downstream analysis in R
