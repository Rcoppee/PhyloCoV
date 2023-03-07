#
# msa
#

frise <- paste0(rep('*', 80), collapse = '')
message(paste0(frise, '\n***** 03 MSA\n', frise), collapse = '')


# Packages Loading ----


# Path ----
out_path <- commandArgs(trailingOnly = TRUE)[1]
fasta_dir <- glue::glue('./{out_path}/02-genomes')
ref <- './data/'

# Data Loading ----
all_sel <- 
  list.files(path = fasta_dir, pattern = '*\\.fasta', full.names = T)

# msa ----
for (i in 1:length(all_sel)) {
  curr_sel <- all_sel[i]
  curr_name <- 
    stringr::str_replace(string = curr_sel, 
                                    pattern = '.*/(.*).fasta', 
                                    replacement = '\\1')
  
  msa_command <- 
    glue::glue('bash ./script/msa_progressive.sh {curr_sel} {curr_name} {out_path}')
  print(glue::glue('##### Launching msa for selection: \n##### {curr_sel} in {out_path}'))
  system(msa_command)
}
