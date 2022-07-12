# clear workspace
rm(list = ls())

# set working directory
setwd("C:/Users/romain.coppee/Documents/PhyloCoV/World/")

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
metadata <- readxl::read_xlsx("metadata_file4.xlsx", sheet = "Feuil1") %>%
	mutate(geo_loc=continent_wave)

#if wave 2, add:
#metadata = metadata[metadata$`Collection date`>="2020-07-20",]

size_set=c()

for(replicat in 1:n_replicates){
	print(replicat)
	#import the tree
	my_tree = read.nexus(file = gsub(" ", "", paste(replicat,"_tree_dated.nexus")))
	#convert to binary tree
	my_tree = multi2di(my_tree)
	#get tips
	my_tips = my_tree$tip.label
    nb_sample = 0
	#we count the number of samples
	for(i in my_tips){
		if(i %in% metadata$gisaid_id){
			nb_sample = nb_sample + 1
		}
    }
	#concatenation of the number of samples in the replicate and the previous replicates
    size_set = c(size_set, nb_sample)
	#add the identifiers
	list_sample = c(list_sample, my_tips)
}

#get the metadata of the different samples across the replicates
metadata2 = metadata[metadata$gisaid_id %in% list_sample,]
#write table if needed
#write.table(metadata2, "world_w1.txt", sep="\t", col.names = T, row.names = F)