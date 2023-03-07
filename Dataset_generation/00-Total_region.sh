#!/bin/bash

wave_1="2000-01-01 2020-07-20"
wave_2="2020-07-20 2020-12-31"
wave_both="2000-01-01 2500-01-01"
# unwanted_1="./data/wave_1_unwanted.tsv"
# unwanted_2="./data/wave_2_unwanted.tsv"
# wanted_regions_1="./data/wanted_regions_1.tsv" # both wanted_regions are the same, ARA PAQ IDF
# wanted_regions_2="./data/wanted_regions_2.tsv"
# output="output-region-wave_1-100_trees"
output="20230208-region-wave_2-100"

selection_1="./data/abn_select-region-wave_1-500.csv"
selection_2="./data/abn_select-region-wave_2-500.csv"


# For proof of french abysmal sequencing
# selection_2="./data/abn_select-French_LowCov-regions2022.csv"
# output="output-region-LowCov"

n_samples=100

# Rscript --vanilla './script/01-panel_selection_region.R' ${output} ${wave_both} ${selection_1} ${n_samples}
Rscript --vanilla './script/01-panel_selection_region.R' ${output} ${wave_both} ${selection_2} ${n_samples}


Rscript --vanilla './script/02-genome_extraction.R' ${output}

Rscript --vanilla './script/03-msa.R' ${output}
Rscript --vanilla './script/033-msa_masking.R' ${output}
rm -r './${output}/02-genomes/*'

Rscript --vanilla './script/04-tree.R' ${output}
rm './${output}/04-tree/*uniqueseq.phy'
rm './${output}/04-tree/*mldist'

Rscript --vanilla './script/05-dating.R' ${output}