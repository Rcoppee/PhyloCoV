# clear workspace
rm(list = ls())

# set working directory
#setwd("C:/Users/romain.coppee/Documents/PhyloCoV/World/")

# import all packages
list.of.packages <- c("adephylo","tidyverse", "lubridate", "glue",
                      "ape", "phytools", "tidytree", "ggplot2",
                      "treeio", "data.table", "dplyr",
                      "parallel", "foreach", "doParallel", "sqldf",
                      "ggtree", "zoo")

new.packages <- list.of.packages[!(list.of.packages %in%
                                     installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages)

lapply(list.of.packages, require, character.only = TRUE)

n_replicates <- 100

#initialize a vector of samples
list_sample = c()

#import metadata
#here is for continent, change the variable "continent_wave" depending on the geographic scale investigated
metadata <- readxl::read_xlsx("metadata_file4.xlsx", sheet = "Feuil1")

size_set=c()

for(replicat in 1:n_replicates){
	print(replicat)
	#import the tree
  my_tree = read.nexus(file = gsub(" ", "", paste(replicat, ".nexus")))
  #convert to binary tree
	my_tree = multi2di(my_tree)
	#get tips
	my_tips = my_tree$tip.label
    nb_sample = 0
	#we count the number of samples
	nb_sample = length(my_tips)
	#concatenation of the number of samples in the replicate and the previous replicates
    size_set = c(size_set, nb_sample)
	#add the identifiers
	list_sample = c(list_sample, my_tips)
}

#get the mean number of sequences sampled in a subsample
mean(size_set)

#get the total number of different sequences
nrow(metadata[metadata$gisaid_id %in% list_sample,])

