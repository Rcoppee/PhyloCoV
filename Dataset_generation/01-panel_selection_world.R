#
# Panel Selection
#

frise <- paste0(rep('*', 80), collapse = '')
message(paste0(frise, '\n***** 01 Panel Selection\n', frise), collapse = '')

# Packages Loading ----
to_load <- c('tidyverse', 'lubridate', 'glue')
for (packus in to_load) {
  print(packus)
  suppressPackageStartupMessages(
    library(packus, 
            character.only = T)
  )
}

rm(list = c('packus', 'to_load'))

# Function ----
pick_from_float <-
  function(vec) {
    l_vec <- length(vec)
    # out_vec <- rep(0, l_vec)
    r_vec <- runif(n = l_vec, min = 0, 1)
    out_vec <- r_vec <= vec 
    return(out_vec)
}

# Args parsing ----
out_path <- commandArgs(trailingOnly = TRUE)[1]
date_begin <- (commandArgs(trailingOnly = TRUE)[2])
date_end <- (commandArgs(trailingOnly = TRUE)[3])
# date_begin <- '2020-07-21'
# date_end <- '2020-12-31'
n_selec <- commandArgs(trailingOnly = TRUE)[4]
abn_select_path <- commandArgs(trailingOnly = TRUE)[5]

# Data Loading ----
gis_meta_path <- './data/2022_06_09-metadata_filt_final.tsv'
system.time({
gis_meta <- 
  read_tsv(gis_meta_path, col_names = T, col_types = cols())
})

abn_select <- 
  read_csv(abn_select_path)

rm(list = c('gis_meta_path', 'abn_select_path'))

# Data Massaging & Filtering ----

min_length <- 29000 # minimum genome length
max_N <- 0.05 # max percentage of unresolved positions

gis_meta %>%
  mutate(date = ymd(`Collection date`)) %>%
  filter(date >= date_begin & date <= date_end ) %>%
  mutate(epi_week = epiweek(date)) %>% 
  mutate(iso_year = isoyear(date)) %>%
  filter(iso_year == 2020) %>%
  filter(Host == 'Human') %>%
  filter(`Sequence length` >= min_length) %>%
  mutate(`N-Content` = if_else(condition = is.na(`N-Content`), 
                               true = 0, 
                               false = `N-Content`, 
                               missing = 0)) %>% # GISAID is strangely missing lots of data about N-content
  filter(`N-Content` <= max_N) %>%
  filter(`Is complete?` == TRUE) %>%
  filter(`Is high coverage?` == TRUE) %>%
  identity() -> gis_meta_filt

# second round of filtering based on the number of NS mutations
max_percentile <- 0.95
NS_df <- 
  tibble(ID = gis_meta_filt$`Accession ID`, 
         count_NS = str_count(gis_meta_filt$`AA Substitutions`, ','))

my_cdf <- ecdf(NS_df$count_NS)

NS_df$quantile <- my_cdf(NS_df$count_NS)
NS_df$outlier <- NS_df$quantile > max_percentile

gis_meta_filt <- gis_meta_filt[!NS_df$outlier, ]

print(glue('Elimininated {sum(NS_df$outlier)} outliers with too much mutations'))
#


gis_meta_filt %>% 
  mutate(Location = str_replace_all(string = Location, pattern = ' ', replacement = '-')) %>%
  mutate(Location = str_replace_all(string = Location, pattern = '-/-', replacement = '/')) %>%
  separate(col = Location, into = c('Continent', 'Country', 'Region', 'City', 'Hood'), sep = '/') %>%
  identity() -> gis_meta_filt

# Stupid country names corrections
abn_select %>% 
  mutate(Country = str_replace_all(string = Country, pattern = ' ', replacement = '-')) %>%
  mutate(Country  = if_else(Country == 'Congo', 'Republic-of-the-Congo', Country)) %>%
  mutate(Country  = if_else(Country == 'Democratic-Republic-of-Congo', 'Democratic-Republic-of-the-Congo', Country)) %>%
  mutate(Country  = if_else(Country == 'United-States', 'USA', Country)) %>%
  identity() -> abn_select

rm(gis_meta)

# Selections ----

max_week <- max(gis_meta_filt$epi_week)
min_week <- min(gis_meta_filt$epi_week)
dayday <- format(Sys.time(), "%Y_%m_%d")
rand_name <- 
  paste0(sample(x = c(letters, LETTERS, 1:9), 
        size = 10, replace = T), 
        collapse = '')
out_dir <- glue('./{out_path}/01-selections/')
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

n_zeroes <- nchar(n_selec)
system.time({
for (i in 1:n_selec) {
  message(paste0('*** selection ', i, collapse = ''))
  zeroes_to_add <-
    n_zeroes - nchar(i)
  if (zeroes_to_add < 0) zeroes_to_add <- 0
i_ <- paste0(paste0(rep(0, each = zeroes_to_add), 
                      collapse = ''),i , collapse = '')
    
  curr_name <- glue('{dayday}-{rand_name}-{i_}-selection')
  curr_out <- glue('{out_dir}/{curr_name}.tsv')
  curr_log <- glue('{out_dir}/{curr_name}.log')
  
  epi_list <- vector('list', length = max_week * length(unique(abn_select$Country)))
  epi_n <- 1
  
  log_df <- tibble(
                  country = NA_character_, 
                  week = NA_integer_, 
                  wanted = NA_integer_, 
                  available = NA_integer_)
  abn_select %>% 
    mutate(infection_proxy_normalised = n_integ + pick_from_float(vec = n_float)) %>% 
    identity() -> abn_select
  
    for (j in min_week:max_week) {
      curr_week_abn <- 
        abn_select %>% 
        filter(week == j) %>% 
        filter(infection_proxy_normalised > 0) %>% 
        identity()
      if (nrow(curr_week_abn) == 0) next()
      
      curr_week_gis <-
        gis_meta_filt %>% 
        filter(epi_week == j) %>% 
        identity()
      
      for (k in 1:length(unique(curr_week_abn$Country))) {
        curr_country <- unique(curr_week_abn$Country)[k]
        n_to_sample <- 
          curr_week_abn %>% 
          filter(Country == curr_country) %>% 
          pull(infection_proxy_normalised)
        
        curr_country_gis <-
          curr_week_gis %>% 
          filter(Country == curr_country) %>% 
          identity()
        
        if (length(n_to_sample) > 1) {
          n_to_sample <- n_to_sample[2]
          message("WARNING, multiple lines")
        }
        
        if (n_to_sample > nrow(curr_week_gis)) {
          warning_message <-
            glue('!!!\nWARNING, not enough genomes !\n
                 Iteration {i}\t Week {j}\t Country {k}')
          message(warning_message)
        }
        
        log_df[nrow(log_df) + 1, ] <-
          list(country = curr_country, 
               week = j, 
               wanted = n_to_sample, 
               available = nrow(curr_country_gis))
        
        if (n_to_sample == 0) next()
        if (nrow(curr_country_gis) == 0) next()
        
        epi_list[[epi_n]] <-
          curr_country_gis %>%
          slice_sample(n = n_to_sample, replace = F) %>%
          identity()
        epi_n <- epi_n + 1
      }
    }
  epi_list <- epi_list[1:epi_n - 1]
  gis_meta_sub <- do.call(rbind, epi_list)
  rm(epi_list)
  write_tsv(gis_meta_sub, curr_out)
  write_tsv(log_df, curr_log)
  rm(list = ls()[startsWith(x = ls(), prefix = 'curr')])
}

}) # system.time



