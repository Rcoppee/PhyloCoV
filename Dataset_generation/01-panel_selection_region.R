#
# Panel Selection
#
Sys.setlocale(locale ="fr_FR.UTF-8")
# Packages Loading ----
to_load <- c('tidyverse', 'lubridate', 'glue')
for (packus in to_load) {
  print(packus)
  suppressPackageStartupMessages(
    library(packus,
            character.only = T)
  )
}
rm(to_load)

# sourcing and defining functions ---- 
source('./script/region_translate.R')

pick_from_float <-
  function(vec) {
    l_vec <- length(vec)
    # out_vec <- rep(0, l_vec)
    r_vec <- runif(n = l_vec, min = 0, 1)
    out_vec <- r_vec <= vec 
    return(out_vec)
}

# Data Loading ----

out_path <- commandArgs(trailingOnly = TRUE)[1]
date_begin <- (commandArgs(trailingOnly = TRUE)[2])

date_end <- (commandArgs(trailingOnly = TRUE)[3])

fb_select_path <- commandArgs(trailingOnly = TRUE)[4]

n_selec <- commandArgs(trailingOnly = TRUE)[5]

if (is.na(n_selec)) {
  stop("No selection asked for...")
  
}
# wanted_regions_path <- commandArgs(trailingOnly = TRUE)[5]
# wanted_regions <- read_tsv(wanted_regions_path)

gis_meta_path <- 
  './data/2022_06_09-metadata_filt_final.tsv'

gis_meta <- 
  read_tsv(gis_meta_path, col_names = T, col_types = cols(), )

# gis_problem_path <- 
#   './data/2019_2020-metadata_2021-11-08_17-08_NotInMSA.tsv'
# gis_problem <- read_tsv(gis_problem_path)


#unwanted_region <- read_tsv(file = unwanted_path)
# fb_select_path <- 
#   './data/ABN-infection_per_region_week-FB.csv'

fb_select <- 
  read_csv(fb_select_path)

rm(list = c('gis_meta_path', 'fb_select_path'))

# Data Massaging & Filtering ----
fb_select %>%
  rename(week = week_shifted) %>%
  mutate(infection_proxy_normalised = round(infection_proxy_normalised * 1)) %>%
  mutate(region = lib_reg) %>% # the translation to ISO code is already done in the subset step as well as the filtering
  #filter( region %in% wanted_regions$region) %>% 
  mutate(Country = region) %>% # the rest of the code expect a variable named Country ...
  identity() -> fb_select

min_length <- 29000 # minimum genome length
max_N <- 0.05 # max percentage of unresolved positions

gis_meta %>%
  mutate(date = ymd(`Collection date`)) %>%
  # filter(date >= date_begin & date <= date_end ) %>%
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

dim(gis_meta_filt)

# second round of filtering based on the number of NS mutations
max_percentile <- 0.95
NS_df <- 
  tibble(ID = gis_meta_filt$`Accession ID`, 
         count_NS = str_count(gis_meta_filt$`AA Substitutions`, ','))

my_cdf <- ecdf(NS_df$count_NS)

NS_df$quantile <- my_cdf(NS_df$count_NS)
NS_df$outlier <- NS_df$quantile > max_percentile

gis_meta_filt <- gis_meta_filt[!NS_df$outlier, ]
dim(gis_meta_filt)
print(glue('Elimininated {sum(NS_df$outlier)} outliers with too much mutations'))


# Region Translation
gis_meta_filt %>%
  mutate(Location = str_replace_all(string = Location, pattern = ' ', replacement = '-')) %>%
  mutate(Location = str_replace_all(string = Location, pattern = '-/-', replacement = '/')) %>%
  separate(col = Location, into = c('Continent', 'Country', 'Region', 'City', 'Hood'), sep = '/') %>%
  filter(Country == 'France') %>%
  filter(! is.na(Region)) %>% 
  mutate(Region = region_translate_short(Region)) %>%
  mutate(Country = Region) %>% # the rest of the code expect a variable named Country ...
  identity() -> gis_meta_filt

dim(gis_meta_filt)


cat('Metadata filtered !\n')

# Check regions ----
bad_regions <- 
sum(! unique(fb_select$region) %in% unique(gis_meta_filt$Region)) # 1 could be region "Unknown" in gis_meta_filt and is ok
glue::glue('Bad regions =  {bad_regions}')
# rm(gis_problem)
rm(gis_meta)

# Selections ----
w_in_y <- max(gis_meta_filt$epi_week)
dayday <- format(Sys.time(), "%Y_%m_%d")
rand_name <- 
  paste0(sample(x = c(letters, LETTERS, 1:9), 
                size = 10, replace = T), 
         collapse = '')

out_dir <- glue('{out_path}/01-selections/')
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

n_zeroes <- nchar(n_selec)
cat('Starting selections !\n')

system.time({
  for (i in 1:n_selec) {
    message(glue('Selection {i}\n'))
    zeroes_to_add <-
      n_zeroes - nchar(i)
    if (zeroes_to_add < 0) zeroes_to_add <- 0
    i_ <- paste0(paste0(rep(0, each = zeroes_to_add), 
                        collapse = ''),i , collapse = '')
    
    curr_name <- glue('{dayday}-{rand_name}-{i_}-selection')
    curr_out <- glue('{out_dir}/{curr_name}.tsv')
    curr_log <- glue('{out_dir}/{curr_name}.log')
    
    epi_list <- vector('list', nrow(gis_meta_filt))
    epi_n <- 1
    
    log_df <- tibble(
      country = NA_character_, 
      week = NA_integer_, 
      wanted = NA_integer_, 
      available = NA_integer_)
    
    fb_select %>% 
      mutate(infection_proxy_normalised = n_integ + pick_from_float(vec = n_float)) %>% 
      identity() -> fb_select
    
    for (j in 1:w_in_y) {
      curr_week_fb <- 
        fb_select %>% 
        filter(week == j) %>% 
        filter(infection_proxy_normalised > 0) %>% 
        identity()
      if (nrow(curr_week_fb) == 0) next()
      
      curr_week_gis <-
        gis_meta_filt %>% 
        filter(epi_week == j) %>% 
        identity()
      for (k in 1:length(unique(curr_week_fb$Country))) {
        
        curr_country <- unique(curr_week_fb$Country)[k]
        
        n_to_sample <- 
          curr_week_fb %>% 
          filter(Country == curr_country) %>% 
          pull(infection_proxy_normalised)
        curr_cntr_gis <-
          curr_week_gis %>% 
          filter(Country == curr_country) %>% 
          identity()
        
        if (length(n_to_sample) > 1) {
          n_to_sample <- n_to_sample[2]
          print("WARNING")
        }
        
        log_df[nrow(log_df) + 1, ] <-
          list(country = curr_country, 
               week = j, 
               wanted = n_to_sample, 
               available = nrow(curr_cntr_gis))
        
        if (n_to_sample == 0) next()
        if (nrow(curr_cntr_gis) == 0) next()
        
        epi_list[[epi_n]] <-
          curr_week_gis %>%
          filter(Country == curr_country) %>%
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
  
})