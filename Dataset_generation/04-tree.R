library(glue)
frise <- paste0(rep('*', 80), collapse = '')
message(paste0(frise, '\n***** 04 Tree\n', frise), collapse = '')

out_path <- commandArgs(trailingOnly = TRUE)[1]
selec_dir <- glue::glue('./{out_path}/01-selections/')
align_dir <- glue::glue('./{out_path}/03-msa/')
tree_dir <- glue::glue('./{out_path}/04-tree/')


iqtree <- '/2020-SARS_CoV_2-seq/00_docker/tools/iqtree-2.2.0-Linux/bin/iqtree2'
fasttree <- '/2020-SARS_CoV_2-seq/00_docker/tools/fasttree/FastTree'

outgroup <- 'EPI_ISL_402125'
outgroup_DDP <- lubridate::ymd('2019-12-26')


if (! dir.exists(tree_dir)) dir.create(tree_dir)

all_msa <- list.files(path = align_dir, pattern = '*.msa', full.names = T) # *.msa_masked


for (i in 1:length(all_msa)) {
  message(paste0(frise, '\n Tree number ', i), collapse = '')
  curr_msa <- all_msa[i]
  curr_name <- stringr::str_replace(string = curr_msa, pattern = '.*/(.*)\\.msa', replacement = '\\1')
  curr_dir <- glue('{tree_dir}{curr_name}')
  if (! dir.exists(curr_dir)) dir.create(curr_dir)
  curr_prefix <- glue('{tree_dir}{curr_name}/{curr_name}')
  curr_selection <- sub(x = curr_name, pattern = '^(.*)_masked$', replacement = '\\1.tsv')
  curr_dates <- 
    readr::read_tsv(file = glue('{selec_dir}{curr_selection}'), 
                    col_names = TRUE)[, c('Accession ID', 'Collection date')]
  curr_dates_path <- glue('{curr_dir}/{curr_name}.dates')
  write.table(curr_dates, file = curr_dates_path, 
              row.names = FALSE, col.names = FALSE,
              sep = '\t')  
  
  
  iqtree_command <- glue('{iqtree} ')
  iqtree_command <- glue('{iqtree_command} -s {curr_msa}')
  iqtree_command <- glue('{iqtree_command} --prefix {curr_prefix}')
  iqtree_command <- glue('{iqtree_command} -T AUTO')
  iqtree_command <- glue('{iqtree_command} --threads-max 32')
  iqtree_command <- glue('{iqtree_command} -m "GTR+I+G" ')
  iqtree_command <- glue('{iqtree_command} -o {outgroup} ')
  iqtree_command <- glue('{iqtree_command} --ufboot 1000 ')
  # iqtree_command <- glue('{iqtree_command} --date {curr_dir}/{curr_name}.dates ')
 
  
  # GTR + G, best compromise https://github.com/roblanf/sarscov2phylo/blob/master/tree_estimation.md
  fasttree_command <- glue('{fasttree} -nt ')
  fasttree_command <- glue('{fasttree_command} -gtr ')
  fasttree_command <- glue('{fasttree_command} -gamma ')
  fasttree_command <- glue('{fasttree_command} -nosupport ')
  # fasttree_command <- glue('{fasttree_command} -boot 100 ')
  fasttree_command <- glue('{fasttree_command} {curr_msa} ')
  fasttree_command <- glue('{fasttree_command} > {curr_dir}.nwk')

print(fasttree_command)
system(fasttree_command)
  # print(iqtree_command)
  # system(iqtree_command)
}

