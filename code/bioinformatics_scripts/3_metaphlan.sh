#!/bin/bash
#<add job submission details>

###
### METAPHLAN for taxonomic analysis
###

cd <path to Toothbrush_Microbiome_Project>
mkdir metaphlan

## load modules for accessing biobakery tools (at NU, this is run on singularity container)
module load singularity

## use loop to run metaphlan on multiple samples
knead_seq=`ls tb_clean_reads | grep fastq.gz""`
for input in $knead_seq
	do output=`echo $input | cut -d "_" -f1`.txt
	singularity exec <path to biobakery-diamondv0822.simg> metaphlan2.py tb_clean_reads/$input --input_type fastq --nproc 4 > metaphlan/$output
done

## clean by removing generated bowtie files
rm tb_clean_reads/*bowtie2out.txt

## merge the taxon abundance table
singularity exec <path to biobakery-diamondv0822.simg> merge_metaphlan_tables.py metaphlan/*.txt > metaphlan/merged_abundance_table.txt
