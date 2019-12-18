#!/bin/bash
#<add job submission details>

###
### KNEADDATA for raw read quality control
###

### Begin by downloading the toothbrush metagenomes raw seq files from NCBI SRA
### Set these files set to Toothbrush_Microbiome_Project/tb_raw_reads

## unzip files
cd <path to Toothbrush_Microbiome_Project>
gunzip tb_raw_reads/*fastq.gz
mkdir kneaddata

## load modules for accessing biobakery tools (at NU, this is run on singularity container)
module load singularity

## use loop to run kneaddata on multiple samples
# first, you must create reference databases (https://bitbucket.org/biobakery/kneaddata/wiki/Home#markdown-header-how-to-run)
# the databases used here were generated from 
# (1) human genome seq data downloaded from the above link and
# (2) from the merged negative control sequence files
# move generated 'human_reference_db' and 'NC_reference_db' to Toothbrush_Microbiome_Project/kneaddata

files_knead=`ls tb_raw_reads | grep "_R1_001.fastq"`
for read1 in $files_knead 
	do read2=`echo $read1 | cut -d R -f1`R2_001.fastq
	out=`echo $read1 | cut -d R -f1`out
	singularity exec <path to biobakery-diamondv0822.simg> kneaddata -t 24 --input tb_raw_reads/$read1 --input tb_raw_reads/$read2 -db kneaddata/human_reference_db -db kneaddata/NC_reference_db --output kneaddata/knead_out/$out
done

## zip input and output files
gzip tb_raw_reads/*fastq
gzip kneaddata/knead_out/*/*fastq

## create clean reads directory; copy in kneaddata output files (use output 1) 
mkdir tb_clean_reads
cd tb_clean_reads
cp <path to tb_raw_reads/kneaddata/knead_out>/*/*kneaddata_paired_1.fastq.gz .
