#!/bin/bash
#<add job submission details>

###
### SHORTBRED_IDENTIFY to build marker database from CARD and uniprot as reference
###

cd <path to Toothbrush_Microbiome_Project>
mkdir shortbred

## move to shortbred directory; add directory for the marker gene files
cd shortbred
mkdir sb_identify

## download and input card-data files (https://card.mcmaster.ca/download) and uniref90 (https://www.uniprot.org/downloads)
mkdir sb_identify/card-data # add card fasta files
mkdir sb_identify/uniref_db # add uniref.fasta.gz

## merge the card protein files (4 total)
cat sb_identify/card-data/protein_fasta_protein_homolog_model.fasta sb_identify/card-data/protein_fasta_protein_knockout_model.fasta sb_identify/card-data/protein_fasta_protein_overexpression_model.fasta sb_identify/card-data/protein_fasta_protein_variant_model.fasta > sb_identify/card-data/CARD_proteins.fasta

## load modules for accessing biobakery tools (at NU, this is run on singularity container)
module load singularity

## run shortbred identify to build CARD marker set with reference to uniref90
gunzip sb_identify/uniref_db/uniref90.fasta.gz
singularity exec <path to biobakery-diamondv0822.simg> shortbred_identify.py --threads 24 --goi sb_identify/card-data/CARD_proteins_adj.fasta --ref sb_identify/uniref_db/uniref90.fasta --markers CARD_markers_u90.faa --tmp tmpdir
gzip sb_identify/uniref_db/uniref90.fasta
