###
### SOURCETRACKER (Knights et al. 2011)
###

## activate qiime (Caporaso et al. 2010; sourcetracker uses qiime)
source activate qiime1

## navigate to access data
cd <path to Toothbrush_Microbiome_Project>

## define env sourcetracker path, and run the program (default parameters)
env SOURCETRACKER_PATH='<path to Toothbrush_Microbiome_Project>/sourcetracker-1.0.1/' R --slave --vanilla --args -i data/taxa_tables/genus/genus_merge_tb_hmp_oreg_gerl_1000x_edit.txt -m data/map_tables/mapping_tb_oreg_gerl_hmp_clean_edit.txt -o sourcetracker_out < $SOURCETRACKER_PATH/sourcetracker_for_qiime.r


## NOTE: genus and map table were set to only include TMP samples as sink data and HMP samples as source data (i.e., subset from the respective data/taxa_tables and data/map_tables files)

## deactivate qiime
source deactivate
