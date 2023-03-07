library(ape)
library(tidyverse)
library(cowplot)

# Simul truth ----
n_samples = 1e2
n_tip <- 1e6
n_tip_sub <- n_tip / 1e3
n_states <- 5
set.seed(1234)
tree <- 
  rtree(n = n_tip, rooted = TRUE, )

median_length <- median(diag(vcv.phylo(tree))) # root to tip distances

my_rates <- sample(1:4, replace = TRUE, 
                   size = n_states*(n_states-1)) * 1/median_length

trait <- 
  rTraitDisc(phy = tree, model = 'ARD', k = n_states, 
             rate = my_rates, 
             ancestor = TRUE)
trait_colours <- 
  c('tomato3', 'steelblue', 'grey50', 'sienna2', 'yellow4')[as.numeric(trait)]

num_tip <- 1:n_tip
num_nodes <- (num_tip[-n_tip]) + (tree$Nnode + 1)

# plot(tree, node.color = trait_colours, show.tip.label = FALSE)
# tiplabels(text = tree$tip.label, frame = 'none', adj = -4)
# tiplabels(adj = -2, col = trait_colours[num_tip], frame = 'none')
# tiplabels(text = trait[num_tip], adj = 0, col = trait_colours[num_tip], frame = 'none')
# nodelabels(col = trait_colours[num_nodes], frame = 'none', adj = -1)
# nodelabels(text = trait[num_nodes], col = trait_colours[num_nodes], frame = 'none', adj = 0)

full_ancestry <-
  data.frame(node = tree$edge[,2], 
             node_trait = trait[tree$edge[,2]], 
             parent = tree$edge[,1], 
             parent_trait = trait[tree$edge[,1]], 
             row.names = NULL)
full_trait_change <- full_ancestry$parent_trait != full_ancestry$node_trait
full_ancestry <- full_ancestry[full_trait_change, ]
true_trait_change <- as.data.frame(table(paste(full_ancestry$parent_trait, full_ancestry$node_trait)))
names(true_trait_change)[2] <- 'True_Freq'

true_trait_change$True_prop <- 
  true_trait_change$True_Freq / sum(true_trait_change$True_Freq)
#-----Sampling---------------------------------

all_sample_ancestry <- list()
system.time({
for (i in 1:n_samples) {
  title <- c("Sample ", i)
  cat(title, date(), '\n')
  curr_tips <- 
    sample(x = num_tip, size = n_tip_sub, replace = FALSE)
  
  curr_tree <- keep.tip(phy = tree, tip = curr_tips)
  
  curr_tip_trait <- 
    trait[curr_tips]
  curr_tip_trait <- 
    curr_tip_trait[curr_tree$tip.label]
  
  curr_ace <- 
    ace(x = curr_tip_trait, phy = curr_tree, type = 'discrete', model = 'ARD')
  
  curr_lik.anc <- Re(curr_ace$lik.anc)
  
  curr_node_trait <- 
    apply(X = curr_lik.anc, MARGIN = 1, 
          function(x) factor(levels(trait))[which(x == max(x))])
  curr_trait <- c(curr_tip_trait, curr_node_trait)
  
  curr_num_tip <- 1:Ntip(curr_tree)
  curr_num_nodes <- (curr_num_tip[-Ntip(curr_tree)]) + (curr_tree$Nnode + 1)
  
  curr_trait_colours <-
    c('tomato3', 'steelblue',
      'grey50', 'sienna2', 'yellow4')[as.numeric(curr_trait)]

  # plot(curr_tree, show.tip.label = FALSE)
  # title(title)
  # tiplabels(text = curr_trait[curr_num_tip], adj = 0, col = curr_trait_colours[curr_num_tip], frame = 'none')
  # nodelabels(col = curr_trait_colours[curr_num_nodes], frame = 'none', adj = -1)
  # nodelabels(text = curr_trait[curr_num_nodes], col = curr_trait_colours[curr_num_nodes], frame = 'none', adj = 0)

  curr_ancestry <-
    data.frame(node = curr_tree$edge[,2], 
               node_trait = curr_trait[curr_tree$edge[,2]], 
               parent = curr_tree$edge[,1], 
               parent_trait = curr_trait[curr_tree$edge[,1]], 
               iter = i,
               row.names = NULL)
  curr_trait_change <- curr_ancestry$parent_trait != curr_ancestry$node_trait
  curr_ancestry <- curr_ancestry[curr_trait_change, ]
  curr_ancestry <- table(paste(curr_ancestry$parent_trait, curr_ancestry$node_trait))
  curr_ancestry <- 
    data.frame(curr_ancestry, iter = i)
  all_sample_ancestry[[i]] <- curr_ancestry
}
})
all_trait_change <- do.call(rbind, all_sample_ancestry)

# delta calc ----
out_delta <-
  list()
for (i in 1:n_samples) {
  all_trait_change |> 
    mutate(N_sub = i) |> 
    filter(iter <= i) |>
    group_by(iter) |> # calculate proportion of Freq by iter and Var1
    mutate(Prop = Freq/sum(Freq)) |>
    ungroup() |>
    left_join(true_trait_change, by = c('Var1' = 'Var1')) |>
    mutate(delta = abs(Prop - True_prop)) |> 
    select(N_sub, delta) |> 
    identity() -> out_delta[[i]]
}

out_delta <- 
  as_tibble(do.call(what = rbind, args = out_delta))

file_out <- glue::glue('./{format(Sys.Date(), format = "%y%M%d")}-out_delta.tsv')
write_tsv(x = out_delta, 
          file = file_out)
# Plotting ----

out_delta_ <- read_tsv(file = '230027-out_delta.tsv')

ggplot(out_delta_) +
  stat_summary(fun.data = 'mean_cl_boot', aes(x = N_sub, y = delta)) +
  scale_x_continuous(name = 'Total subsamples') +
  scale_y_continuous(name = 'Mean error proportion ', 
                     labels = scales::percent_format(scale = 100)) +
  theme_cowplot() +
  geom_blank() -> plot
plot

plot_out <- glue::glue('./{format(Sys.Date(), format = "%y%M%d")}-simul_plot.pdf')
save_plot(plot, filename = plot_out, base_width = 8, base_height = 5)

summary(out_delta_$delta[out_delta_$N_sub == 100])
boxplot(out_delta_$delta[out_delta_$N_sub == 100])
