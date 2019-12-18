#!/bin/bash
#<add job submission details>

###
### HALLA for correlations between metadata features and species/ARGS in toothbrush samples
###

cd <path to Toothbrush_Microbiome_Project>
mkdir halla
cd halla

## load modules for accessing biobakery tools (at NU, this is run on singularity container)
module load singularity

## halla with species table
singularity exec <path to biobakery-diamondv0822.simg> halla -X species_table.txt -Y metadata_halla.txt --output output_taxa --diagnostics-plot

## halla with ARGs table
singularity exec <path to biobakery-diamondv0822.simg> halla -X ARGs_table.txt -Y metadata_halla.txt --output output_ARG --diagnostics-plot

## NOTE: similarity tables and associations tables from the output are used for the creation of the 'hallagrams' (Fig. 6)
