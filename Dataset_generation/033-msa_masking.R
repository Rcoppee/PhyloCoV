##
## Masking the problematic sites in SARS-CoV2 genomes 
## 

frise <- paste0(rep('*', 80), collapse = '')
message(paste0(frise, '\n***** 033 MSA masking\n', frise), collapse = '')

library(tidyverse)

out_path <- 
  commandArgs(trailingOnly = TRUE)[1]

in_align_path <- 
  glue::glue('{out_path}/03-msa/')

message(in_align_path)

prob_sites_path <- 
  './data/problematic_sites_sarsCov2.vcf'
  
prob_sites <- 
  read_tsv(prob_sites_path, comment = '#', 
           col_names = c('CHROM', 'POS', 'ID', 'REF', 
                         'ALT', 'QUAL', 'FILTER', 'INFO'))
pos_to_mask <- 
  prob_sites$POS

all_msa <- 
  list.files(path = in_align_path, pattern = '*.msa$', full.names = TRUE)

mask_line <- function(align, index, site_list) {
  line <- align[index]
  if (index %% 2 > 0) { 
    return(line) 
  }
  line <- str_split(line, pattern = '')[[1]]
  line[site_list] <- 'n'
  line <- paste0(line, collapse = '')
  return(line)
}



for (i in 1:length(all_msa)) {
  
  # curr_out <- as.character(glue::glue(all_msa[i], '_masked'))
  curr_out <- sub(x = all_msa[i], pattern = '(.*)\\.msa', replacement = '\\1_masked.msa')
  curr_align <- read_lines(all_msa[i])
  l_align <- length(curr_align)
  curr_masked <- character(length = l_align)
  cat('Starting ', all_msa[i], ' ',i,'/', length(all_msa),'\n')
  for (j in 1:l_align) {
    curr_masked[j] <-
      mask_line(align = curr_align, index = j, site_list = pos_to_mask)
  }
  message(curr_out)
  write_lines(x = curr_masked, path = curr_out, sep = '\n', append = FALSE)
  file.remove(all_msa[i])
}



