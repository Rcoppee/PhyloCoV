

#
# Genomes extraction
#

frise <- paste0(rep('*', 80), collapse = '')
message(paste0(frise, '\n***** 02 Genome Extraction\n', frise), collapse = '')

# Packages Loading ----
to_load <- c('tidyverse', 'lubridate', 'glue')
for (packus in to_load) {
  print(packus)
  suppressPackageStartupMessages(library(packus,
                                         character.only = T))
}
rm(list = c('packus', 'to_load'))

# tools & Path ----
out_path <- commandArgs(trailingOnly = TRUE)[1]
seqtk <- '/tools/seqtk/seqtk'

sel_dir <- glue('./{out_path}/01-selections')
out_dir <- glue('./{out_path}/02-genomes/')

if (!dir.exists(out_dir))
  dir.create(out_dir, recursive = TRUE)

# Data Loading ----
msa_epi_path <-
  './data/msa_0527-2019_2020.fasta'

wanted_name <- '2021_10_28-sRu5FXAtOb'

all_sel <-
  list.files(
    path = sel_dir,
    #pattern = glue('{wanted_name}.*\\.tsv'),
    pattern = glue('.*\\.tsv'),
    full.names = T
  )

glue('Launching {length(all_sel)} genome subsets ')

if (length(all_sel) < 1) {
  message("Empty input !")
  quit(save = 'no', status = 1)
}
# Extraction ----
system.time({
  for (i in 1:length(all_sel)) {
    system.time({
      curr_sel <- all_sel[i]
      if(is.na(curr_sel)) next()
      curr_out_name <- 
        str_replace(string = curr_sel,
                                   pattern = '.*/(.*)\\.tsv',
                                   replacement = '\\1.fasta')
      sub_command <-
        glue(
          'bash -c "{seqtk} subseq {msa_epi_path} <(cat {curr_sel} | cut -f 3) > {out_dir}{curr_out_name}"')
      print(glue('Launching iteration number {i}:\n {sub_command}'))
      system(sub_command)
    })
  }# internal system.time
}) # whole loop system.time