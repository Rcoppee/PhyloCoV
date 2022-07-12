#here is the script for getting dates from world datasets
#the phylogenis have a specific name, and a number corresponding to the replicate

# clear workspace
rm(list = ls())

# set working directory
setwd("C:/Users/romain.coppee/Documents/PhyloCoV/World/")

# import all packages
list.of.packages <- c("adephylo","tidyverse", "lubridate", "glue", "cowplot",
                      "ape", "phytools", "tidytree", "ggplot2", "ggtree",
                      "treeio", "data.table", "xlsx", "RColorBrewer", "dplyr",
                      "ggridges","ggmuller","circlize","compare",
                      "parallel", "foreach", "doParallel", "sqldf", "ggthemes",
                      "viridis","hrbrthemes")
options(knitr.table.format = "html")

new.packages <- list.of.packages[!(list.of.packages %in%
                                     installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages)

lapply(list.of.packages, require, character.only = TRUE)

#import metadata. Replace continent_wave by the variable, depending on the geographic scale investigated
metadata <- readxl::read_xlsx("metadata_file4.xlsx", sheet = "Feuil1") %>%
	mutate(geo_loc = continent_wave)

#number of replicates
n_replicat = 100

for(replicat in 1:n_replicat){
	print(replicat)
	#import the tree
	my_tree = read.tree(file = gsub(" ", "", paste("2022_06_24-ZueYqbus3U-",replicat,"-selection.nwk")))
	tree_df = as_tibble(my_tree) 
	
	#get tips
	my_tips = my_tree$tip.label
	metadata_sub = metadata[metadata$gisaid_id %in% my_tips,]
	
	#get only the gisaid id and the date, corresponding to the input for treetime
	metadata_sub = metadata_sub[,c(1,2)]

	#rename columns for treetime
	names(metadata_sub)[1] <- "name"
	names(metadata_sub)[2] <- "date"
	
	#write the table
	write.table(metadata_sub, paste(replicat, "_date.txt", sep=""), sep=",", col.names = T, row.names = F, quote = F)
}