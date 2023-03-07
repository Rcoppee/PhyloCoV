#!/bin/bash
# output="output-world-wave_1-1000"
# output="output-world-wave_2-spike_1-1500"
# output="output-euro-wave_1-1000"
output="20230113-world-wave_both_2000-25"
wave_1="2000-01-01 2020-07-20"
wave_2="2020-07-21 2020-12-31"
wave_both="2019-01-01 2020-12-31"

selec_1="./data/abn_select-continent-wave_1-1000.csv"
selec_2="./data/abn_select-continent-wave_2-1000.csv"
selec_both="./data/abn_select-continent-wave_1_1000-wave_2_1000.csv"
selec_eur_1="./data/abn_select-euro-wave_1-1000.csv"
selec_eur_2="./data/abn_select-euro-wave_2-1000.csv"

n_sample=100

# Rscript --vanilla './script/01-panel_selection_world.R' ${output} ${wave_both} ${n_sample} ${selec_1}
# Rscript --vanilla './script/01-panel_selection_world.R' ${output} ${wave_both} ${n_sample} ${selec_2}
Rscript --vanilla './script/01-panel_selection_world.R' ${output} ${wave_both} ${n_sample} ${selec_both}
# Rscript --vanilla './script/01-panel_selection_world.R' ${output} ${wave_both} ${n_sample} ${selec_eur_1}
# Rscript --vanilla './script/01-panel_selection_world.R' ${output} ${wave_both} ${n_sample} ${selec_eur_2}

Rscript --vanilla './script/02-genome_extraction.R' ${output}

Rscript --vanilla './script/03-msa.R' ${output}
Rscript --vanilla './script/033-msa_masking.R' ${output}

rm ./${output}/02-genomes/*

Rscript --vanilla './script/04-tree.R' ${output}

rm './output/03-msa/*'
rm './output/04-tree/*uniqueseq.phy'
rm './output/04-tree/*mldist'

Rscript --vanilla './script/05-dating.R' ${output} "FALSE"
# 
# Rscript --vanilla './script/06-finalize_trees.R' ${output}

