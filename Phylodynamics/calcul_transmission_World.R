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

#import metadata
metadata <- readxl::read_xlsx("metadata_file4.xlsx", sheet = "Feuil1") %>%
	mutate(geo_loc = continent)

#fix number of replicates
n_replicates <- 100		

for(replicat in 1:n_replicates){
	print(replicat)
	#import tree
	my_tree = read.nexus(file = gsub(" ", "", paste(replicat,"_tree_dated.nexus")))
	my_tree = multi2di(my_tree)
	
	my_tips = my_tree$tip.label
	#keep only data from the tree
	metadata_sub = metadata[metadata$gisaid_id %in% my_tips,]
	
	#add France in continent category
	for(i in 1:nrow(metadata_sub)){
		if(metadata_sub$country[i]=="France"){
			metadata_sub$geo_loc[i]="France"
		}
	}
	
	#generate a tibble
	tree_table <- as_tibble(my_tree) %>%
		left_join(select(metadata_sub, c('gisaid_id', 'geo_loc', 'Clade')), by = c('label' = 'gisaid_id'))

	as.treedata(tree_table) -> tree_data
  
	#calculate distances
	root_node <- max(tree_table$node[tree_table$node <= Ntip(tree_data)])+1
	distances = dist.nodes(as.phylo(tree_table))
	distance_root = distances[,root_node]
	distance_root = as.data.frame(distance_root)
  
	#define the date origin of the tree (manual)
	tmrca = Date2decimal(as.Date("2020-01-07"))
	#if second wave, replace the last line to fix the origin by
	#tmrca = Date2decimal(as.Date("2020-02-01"))
  
	for (i in 1:nrow(distance_root)) {
		distance_root[i,2] = (decimal2Date(tmrca+distance_root[i,1]))
	}
	distance_root$gisaid_id = tree_table$label
	colnames(distance_root)[2] <- "old_dates"
  
	get_branch_length = my_tree$edge.length[my_tree$edge[,2] > Ntip(my_tree)]
  
	#convert distances to dates
	for (i in 1:nrow(distance_root)) {
		if(i<=root_node){
			distance_root$new_dist[i]=0
		}else {
			distance_root$new_dist[i]=get_branch_length[i-root_node]/2
		}
    distance_root$new_date[i]=(decimal2Date(tmrca+distance_root[i,1])+distance_root$new_dist[i])
	}
	distance_root$new_date = as.Date(distance_root$new_date)
	colnames(distance_root)[5]<- "dates"
  
	#add infos such as localization and clade in tree_table
	tree_table %>% left_join(select(metadata_sub, c('gisaid_id', 'geo_loc', 'Clade')), by = c('label' = 'gisaid_id')) %>%
    identity() -> tree_table
  
	#retain only important data
	tree_table$dates = distance_root$dates
	tree_table = tree_table[,c(1,2,3,4,5,6,9)]
	names(tree_table)[names(tree_table) == "geo_loc.x"] <- "geo_loc"
	names(tree_table)[names(tree_table) == "Clade.x"] <- "clade"

	#Solving the branch with length 0
	zero_offset <- 1e-2/365
	is_zero_branch <- tree_table$branch.length == 0
	tree_table$branch.length[is_zero_branch] <- zero_offset
	summary(tree_table$branch.length)
  
	#ancestrat state reconstruction, here using the SYM model
	ancest <- ace(x = tree_table$geo_loc[1:Ntip(my_tree)], phy = as.phylo(tree_table), type = 'discrete', model = 'SYM')
  
	#get likelihood values and detect the max
	lik_SYM = ancest$lik.anc
	df_rates_SYM = Re(lik_SYM)
	df_rates_SYM = as.data.frame(df_rates_SYM)
  
	df_rates_SYM$max = apply(df_rates_SYM, 1, max)
	origin_SYM = colnames(df_rates_SYM)[apply(df_rates_SYM,1,which.max)]
  
	#update the nodes
	tree_table[(Ntip(my_tree)+1):(Ntip(my_tree) + Nnode(my_tree)),'geo_loc'] <- origin_SYM
  
	#convert to dataframe
	df_tree = as.data.frame(tree_table)
  
	#table initialization
	df_inter_SYM_all = data.frame(In=character(), Out=character(), Date=as.Date(character()), Replicat=integer())
	df_intra_SYM_all = data.frame(In=character(), Out=character(), Date=as.Date(character()), Replicat=integer())
  
	#calculation of transmissions
	for (i in root_node:nrow(tree_table)) {
		child1 = child(tree_table, i)[[2]][1] #child node1
		child2 = child(tree_table, i)[[2]][2] #child node2
    
		#if child nodes are identical and identical to parent = intra, else = inter
		if(df_tree[child1,]$geo_loc==df_tree[child2,]$geo_loc){
			df_intra_SYM_all[nrow(df_intra_SYM_all) + 1,] = list(df_tree[child1,]$geo_loc, df_tree[child2,]$geo_loc, df_tree[i,]$dates, replicat)
		}
		else{
			if(df_tree[i,]$geo_loc!=df_tree[child1,]$geo_loc){
				df_inter_SYM_all[nrow(df_inter_SYM_all) + 1,] = list(df_tree[i,]$geo_loc, df_tree[child1,]$geo_loc, df_tree[i,]$dates, replicat)
			}
			else{
				df_inter_SYM_all[nrow(df_inter_SYM_all) + 1,] = list(df_tree[i,]$geo_loc, df_tree[child2,]$geo_loc, df_tree[i,]$dates, replicat)
			}
		}
	}
  
	#Write table for each replicate inter/intra
	write.table(x = df_inter_SYM_all, file = paste(replicat, "_inter_NEW.sym", sep=""), sep="\t")
	write.table(x = df_intra_SYM_all, file = paste(replicat, "_intra_NEW.sym", sep=""), sep="\t")
	
}           


 
 
 