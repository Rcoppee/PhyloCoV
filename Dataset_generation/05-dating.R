#
# Tree Dating
#

frise <- paste0(rep('*', 80), collapse = '')
message(paste0(frise, '\n***** 05 TREE DATING\n', frise), collapse = '')

library(tidyverse)
# Tools & Paths ----
lsd2_path <- '../00_docker/tools/lsd2_unix'
out_path <- commandArgs(trailingOnly = TRUE)[1]
oldest_as_outgroup <- as.logical(commandArgs(trailingOnly = TRUE)[2])
if (is.na(oldest_as_outgroup)) {
  print("No second arg")
  oldest_as_outgroup <- FALSE
}
# metadata <- './data/2022_06_09-metadata_filt_final.tsv'
# gis_dates <- './data/gis_dates_sorted.tsv'
# date_command <- 
# glue::glue('cat {metadata} | awk \'{{print $3"\\t"$4}}\' | sort > {gis_dates}')
# print(date_command)
# system(date_command)
selec_dir <- glue::glue('./{out_path}/01-selections/')
msa_dir <- glue::glue('./{out_path}/03-msa/')
align_size <- 29703
rate <- 5e-4
input_dir <- glue::glue('./{out_path}/04-tree/')
output_dir <- glue::glue('./{out_path}/05-dating/')

if (! dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)

print(commandArgs(trailingOnly = TRUE))
print(oldest_as_outgroup)
# Loading ----
all_trees <-
  list.files(
    path = input_dir,
    pattern = '.*nwk$', # '.*treefile' for ML tree '.*contree' for consensus and bootstrap '.*nwk' for fasttree
    full.names = TRUE,
    recursive = TRUE
  )

print(glue::glue('Found {length(all_trees)} trees ...'))

# dates_df <- 
  # read_tsv(gis_dates, col_names = F) # BADDDDD get dates from 01_panel_selection

# Dating ----

## create the outgroup file ----
system('bash -c "echo -e \'1\nEPI_ISL_402125\' > temp.outgroup"')

for (i in 1:length(all_trees)) {
# for (i in 1) {
  curr_tree <- all_trees[i]
  curr_name <- str_remove(string = basename(curr_tree), pattern = '\\..*')
  curr_out_dir <- glue::glue('{output_dir}{curr_name}')
  curr_selection <- str_replace(string = curr_name, 
                                pattern = '^(.*)_masked$', 
                                replacement = '\\1')
  curr_selection <- paste0(selec_dir, curr_selection,'.tsv', collapse = '')
  
  curr_dates <- 
    read_tsv(file = curr_selection, col_names = TRUE)
  
  curr_dates %>%
    select(`Accession ID`, `Collection date`) %>%
    rbind(c('EPI_ISL_402125', '2019-12-26')) %>%
    mutate(`Collection date` = as.character(`Collection date`)) %>%
    identity() -> curr_dates
  names(curr_dates) <- c('X1', 'X2')
  
  message(paste0('\n\nStarting iteration ', i, ' named ', curr_name), collapse = '')
  #if (! dir.exists(curr_out_dir))  dir.create(curr_out_dir)
  curr_msa <-
    list.files(path = msa_dir,
               pattern = curr_name,
               full.names = TRUE)
  print(curr_msa)
  if (length(curr_msa) == 0){
    message('MSA empty !')
    next()
  } 
  
  ## Create the current ID file ----
  id_command <- 
    glue::glue('cat {curr_msa} | gawk "/^>EPI/ {{print}}" | tr -d ">" | sort > temp_{i}.ID')
  print(id_command)
  system(id_command)
  id_df <- 
    read_tsv(glue::glue('temp_{i}.ID'), col_names = F)

  id_df %>% 
    left_join(curr_dates) %>% 
    # mutate(X2 = lubridate::decimal_date(lubridate::ymd(X2))) %>% # test with decimal dates
    identity() -> id_df
  
  ## Use internal outgroup ----
  if (oldest_as_outgroup) {
    print('REROOTING ON OLDEST TIP')
    ### remove Hu-1/EPI_ISL_402125 from ID file ----
    id_df %>%
      filter(X1 != 'EPI_ISL_402125') %>% 
      identity() -> id_df
    write_tsv(x = id_df, glue::glue('temp_{i}.ID'),col_names = FALSE)
    
    ### prune it from the tree ----
    full_tree <- 
      ape::read.tree(file = curr_tree)
    pruned_tree <- 
      ape::drop.tip(phy = full_tree,
                    tip = 'EPI_ISL_402125',
                    trim.internal = TRUE)
    ### reroot the tree ----
    oldest <- 
       sort(lubridate::ymd(id_df$X2))[1]
    
    oldest_tip <-
      id_df$X1[id_df$X2 == oldest][1]
    
    pruned_tree <-
      ape::root.phylo(phy = pruned_tree, 
                    outgroup = oldest_tip, 
                    resolve.root = TRUE)
    new_tree_path <-
      glue::glue('{curr_tree}_reroot')
    ape::write.tree(phy = pruned_tree, file = new_tree_path)
    curr_tree <- new_tree_path
    ### new outgroup file ----
    curr_outgroup <- glue::glue('1\n{oldest_tip}')
    write(curr_outgroup, 'temp.outgroup')
  }
  
  
  ## Create the current ID \t dates file ----
    date_command_header <- 
      glue::glue('wc -l temp_{i}.ID | cut -d \' \' -f 1 > temp_{i}.date')
    system(date_command_header)
    # write_tsv(x = id_df, glue::glue('temp_{i}.date'), append = TRUE)
    write.table(x = id_df, file = glue::glue('temp_{i}.date'), quote = FALSE, 
                append = TRUE, col.names = FALSE, row.names = FALSE)
    
  
  ## create the rates file ----
    write(x = rate, glue::glue('temp_{i}.rate'))
    
  # lsd2 ----
  lsd_command <-
  glue::glue('{lsd2_path} ')
  lsd_command <-
    glue::glue(lsd_command, '-i {curr_tree} ')
  lsd_command <-
    glue::glue(lsd_command, '-d temp_{i}.date ')
  lsd_command <-
    glue::glue(lsd_command, '-s {align_size} ')
  # lsd_command <-
  # glue::glue(lsd_command, '-b {1e-6} ') # variance param
  # lsd_command <-
  #   glue::glue(lsd_command, '-v 1 ') # variance
  lsd_command <-
    glue::glue(lsd_command, '-r k ') # l to estimate around, a for everywhere, k on outgroup branch
  lsd_command <-
    glue::glue(lsd_command, '-e 2 ') # Zscore outlier
  lsd_command <-
    glue::glue(lsd_command, '-l {0.0/align_size} ')   # nullBlen default 0.5/align_size   
  lsd_command <-
    glue::glue(lsd_command, '-R 365 ') # round to day
  # lsd_command <-
    # glue::glue(lsd_command, '-F ')  # removes the constrains constraints that the date of every node is equal or smaller then the
                                    # dates of its descendants
  lsd_command <-
    glue::glue(lsd_command, '-u 0 ') # minBlen forces every branch of the time scaled tree to be >= minBlen (option e) or >= 0 (option 0)
  # lsd_command <-
  #   glue::glue(lsd_command, '-E ') # estimate internal and external branches separately (takes all eternity to run)
  # lsd_command <-
  #   glue::glue(lsd_command, '-S 0.7 ') # min support
  lsd_command <-
    glue::glue(lsd_command, '-g temp.outgroup ')
  # lsd_command <-
  #   glue::glue(lsd_command, '-w temp_{i}.rate ')
  lsd_command <-
    glue::glue(lsd_command, '-o {curr_out_dir}')

  print(lsd_command)        
  system(lsd_command)
}

# cleanup ----
system('rm temp*')
