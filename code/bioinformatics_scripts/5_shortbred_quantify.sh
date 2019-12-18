#!/bin/bash
#<add job submission details>

###
### SHORTBRED_QUANTIFY to determine ARG protein families metagenome samples
###

cd <path to Toothbrush_Microbiome_Project>

## load modules for accessing biobakery tools (at NU, this is run on singularity container)
module load singularity

## run loop for shortbred quantify on all samples
knead_seq=`ls tb_clean_reads | grep fastq.gz""`
for input in $knead_seq
	do output=`echo $input | cut -d "_" -f1`_sb.txt
	tmpdir=`echo $input | cut -d "_" -f1`_tmp_u90
	singularity exec <path to biobakery-diamondv0822.simg> shortbred_quantify.py --threads 8 --markers shortbred/card_markers_u90 --wgs tb_clean_reads/$input --results shortbred/output/$output --tmp shortbred/tmpdir/$tmpdir
done


###
### merging output files
###

## navigate to output directory
cd <path to Toothbrush_Microbiome_Project/shortbred/output>

## rename columns in output files to have sample names
# loop for all samples
name=`ls | grep txt""`
for type in $name
	do fam=`echo $type | cut -d "_" -f1`_Family
	count=`echo $type | cut -d "_" -f1`_count
	hits=`echo $type | cut -d "_" -f1`_Hits
	tml=`echo $type | cut -d "_" -f1`_TotMarkerLength
	sed -i "1s/.*/$fam $count $hits $tml/" $type
done

## pull 'counts' columns
mkdir counts
# use loop for all samples
name=`ls | grep txt""`
for input in $name
	do output=`echo $input | cut -d "_" -f1`counts.txt
	awk '{print $2}' $input > counts/$output
done

### combine counts from all samples into single file
cd counts
paste *txt > sb_all_CARD.txt
