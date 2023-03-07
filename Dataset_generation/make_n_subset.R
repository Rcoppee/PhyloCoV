# Get n representative SARS-CoV-2 genomes in 2020
# antoine.bridier-nahmias@inserm.fr
# 2022-01-26

# Packages Loading ----
to_load <- 
  c('tidyverse', 'lubridate', 'glue')
for (packus in to_load) {
  print(packus)
  suppressPackageStartupMessages(
    library(packus, 
            character.only = T)
  )
}

rm(list = c('packus', 'to_load'))

# Data Loading ----
death_ori <- 
  read_csv('data/20220121-covid_death_verif.csv')

# Variables ----
n_week <- 53
death_lag <- 2
wave_limit_week <- 30

n_tot <- 1e3
wave_wanted <- 2

# Massaging ---
death_ori %>% 
  select(c(1:4,9)) %>%
  filter(date < '2021-02-01') %>% 
  filter(! startsWith(iso_code, 'OWID')) %>% 
  mutate(death_noshift = if_else(is.na(new_deaths), 0, new_deaths)) %>%
  mutate(death_noshift = if_else(death_noshift < 0, 0, death_noshift)) %>% 
  mutate(year = year(date)) %>% 
  mutate(week = isoweek(date)) %>% 
  mutate(week = if_else(year == 2021, week + n_week, week)) %>%
  filter(week <= 55) %>% 
  group_by(iso_code, week) %>% 
  mutate(death_noshift = sum(death_noshift)) %>% 
  slice(1) %>%
  select(c(-new_deaths)) %>% 
  ungroup() %>% 
  identity() -> death_df

my_shifter <- function(vec, by) {
  if (by == 0) return(vec)
  c(tail(vec, -by), head(vec, by))
}


death_df %>% 
  group_by(iso_code) %>%
  arrange(week) %>% 
  mutate(death = my_shifter(vec = death_noshift, by = death_lag)) %>%
  filter(year == 2020) %>% 
  ungroup() %>% 
  identity() -> death_df
  
  
  

# Selection ----

# Continents ----
death_df %>%
  mutate(wave = if_else(week <= wave_limit_week, 1, 2)) %>%
  filter(wave == wave_wanted) %>% 
  # mutate(continent = if_else(iso_code == 'FRA', "France", continent)) %>% 
  group_by(location, week) %>% 
  mutate(death = sum(death)) %>%
  slice(1) %>% 
  ungroup() %>% 
  mutate(prop = death / sum(death)) %>%
  mutate(n_to_sample = prop * n_tot) %>% 
  mutate(n_integ = floor(n_to_sample)) %>% 
  mutate(n_float = n_to_sample - n_integ) %>%
  mutate(continent = str_replace_all(string = continent, 
                                     pattern = ' ', 
                                     replacement = '-')) %>%
  rename('Country' = 'location') %>% 
  # select(-iso_code, -location, -date) %>% 
  identity() -> death_sample_cont

death_sample_cont %>% 
  group_by(continent) %>% 
  summarize(n = sum(n_to_sample))

death_sample_cont %>% 
  group_by(Country) %>% 
  summarize(n = sum(n_to_sample)) %>% arrange(n) |> print(n = Inf)


death_sample_cont %>% 
  ggplot() +
    # geom_smooth(aes(x = week, y = n_to_sample, colour = continent), n = 100, method = 'gam', se = FALSE) +
    geom_line(aes(x = week, y = n_to_sample, colour = continent)) +
    geom_point(aes(x = week, y = n_to_sample, colour = continent)) +
    geom_vline(xintercept = 30) +
    # scale_colour_viridis_d() +
    theme_bw()

death_sample_cont %>%
  group_by(continent) %>% 
  arrange(week) %>% 
  mutate(death_cum = cumsum(death)) %>% 
  ggplot() +
  # geom_smooth(aes(x = week, y = n_to_sample, colour = continent), n = 100, method = 'gam', se = FALSE) +
  geom_line(aes(x = week, y = death_cum, colour = continent)) +
  geom_point(aes(x = week, y = death_cum, colour = continent)) +
  geom_vline(xintercept = 32) +
  # xlim(c(0,35)) +
  # ylim(c(0,3e5)) +
  # scale_colour_viridis_d() +
  theme_bw()


write_csv(death_sample_cont,
          file = glue('./data/abn_select-continent-wave_{wave_wanted}-{n_tot}.csv'),
          col_names = TRUE)
# write_csv(death_sample_cont,
#           file = glue('./data/abn_select-continent-wave_{wave_wanted}-{n_tot}.csv'),
#           col_names = c(TRUE, FALSE)[wave_wanted])

# Europe ----

# Quick death study
death_df |> 
  filter(continent == 'Europe') |> 
  group_by(location) |> 
  mutate(n_death = sum(death)) |> 
  slice(1) |> 
  select(iso_code, location, n_death) |> 
  identity() -> death_euro_sum

euro_death_cutoff <- 1e4

ggplot(death_euro_sum) +
  # geom_histogram(aes(x = n_death), binwidth = 1e3) +
  geom_col(aes(y = reorder(iso_code, n_death), x = n_death)) +
  geom_vline(xintercept = euro_death_cutoff) +
  
  theme_bw()



unwanted_euro <- 
  c('UKR', 'CZE', 'HUN')

death_euro_sum |> 
  filter(n_death > euro_death_cutoff) |>
  filter(! iso_code %in% unwanted_euro) |> 
  pull(iso_code) -> wanted_euro_list

print(wanted_euro_list)

n_wanted_europe <- c(0, 1e3)

death_df |> 
  filter(continent == 'Europe') |> 
  filter(iso_code %in% wanted_euro_list) |> 
  mutate(prop = death / sum(death)) |>
  mutate(wave = if_else(condition = week <=30, 
                        true = 1, 
                        false = 2)) |>
  group_by(wave) %>% 
  mutate(prop = death / sum(death)) %>%
  ungroup() |> 
  mutate(n_sample = prop * n_wanted_europe[wave]) %>%
  mutate(n_integ = floor(n_sample)) %>%
  mutate(n_float = n_sample - n_integ) %>% 
  mutate(infection_proxy_normalised = 0) %>%
  rename('Country' = 'location') |> 
  identity() -> death_df_euro

sum(death_df_euro$n_sample)

# write_csv(death_df_euro, file = './data/abn_select-euro-wave_1-1000.csv')
write_csv(death_df_euro,
          file = glue('./data/abn_select-euro-wave_{wave_wanted}-{sum(n_wanted_europe)}.csv'))
# OLD
# death_sample_cont %>% 
#   filter(continent == 'Europe') %>% 
#   group_by(Country) %>%
#   summarize(n = sum(n_to_sample)) %>%
#   ungroup() %>%
#   arrange(desc(n)) %>% 
#   mutate(cum = cumsum(n)) %>% 
#   identity() %>% 
#   print(n = Inf)
# 
# death_sample_cont %>% 
#   filter(continent == 'Europe') %>% 
#   identity() %>% 
#   write_csv(file = './data/abn_select_europe-wave_1-500.csv')
# 
# read_csv(file = './data/abn_select_europe-wave_1-500.csv') %>% 
#   group_by(Country) %>% 
#   arrange(week) %>% 
#   mutate(samp_cum = cumsum(n_to_sample)) %>% 
#   ggplot() +
#   geom_line(aes(x = week, y = samp_cum, colour = Country)) +
#   theme_bw() +
#   theme(legend.position = "none") -> p
# 
# plotly::ggplotly(p)


# Regions ----
## Modify FB_selection for regions

source('./script/region_translate.R')

fb_total <- 
  read_csv('./data/infection_per_region_week-FB.csv')

wanted_regions_path <- "./data/wanted_regions.tsv"
# wanted_regions_2="./data/wanted_regions_2.tsv"
wanted_regions <- read_tsv(wanted_regions_path, comment = '#')

fb_total %>% 
  mutate(wave = if_else(condition = week_shifted <=30, 
                 true = 1, 
                 false = 2)) %>%
  mutate(lib_reg = region_translate_short(lib_reg)) %>% 
  filter(lib_reg %in% wanted_regions$region) %>% 
  group_by(wave) %>% 
  mutate(prop = infection_proxy / sum(infection_proxy)) %>%
  identity() -> fb_total_


n_wanted <- 
  c(0.5e3, 0.5e3)

fb_total_ %>% 
  mutate(n_sample = prop * n_wanted[wave]) %>%
  mutate(n_integ = floor(n_sample)) %>%
  mutate(n_float = n_sample - n_integ) %>% 
  mutate(infection_proxy_normalised = 0) %>% 
  identity() -> fb_total_out

sum(fb_total_out$n_sample)
sum(fb_total_out$n_sample[fb_total_out$wave == 1])
sum(fb_total_out$n_sample[fb_total_out$wave == 2])


output_wave_1 <- glue('./data/abn_select-region-wave_1-{n_wanted[1]}.csv')
output_wave_2 <- glue('./data/abn_select-region-wave_2-{n_wanted[2]}.csv')

write_csv(x = fb_total_out[fb_total_out$wave == 1, ], 
          output_wave_1, col_names = TRUE) # col_names TRUE for wave_1, FALSE for wave_2
write_csv(x = fb_total_out[fb_total_out$wave == 2, ], output_wave_2, 
          col_names = TRUE) # col_names TRUE for wave_1, FALSE for wave_2

# for wave 2 spike wave 1
output_wave_2_spike_1 <- 
  './data/ABN-infection_per_region_week-wave_2_spike_1-FB.csv'

concat_command <-
  glue('cat {output_wave_1} {output_wave_2} > {output_wave_2_spike_1}')

system(concat_command)

