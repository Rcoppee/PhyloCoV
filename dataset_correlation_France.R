# clear workspace
rm(list = ls())

# set working directory
setwd("C:/Users/romain.coppee/Documents/PhyloCoV/France/")

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

#for wave 1

n_replicates <- 100

#initialize a vector of samples
list_sample = c()

#import metadata
metadata <- readxl::read_xlsx("metadata_file4.xlsx", sheet = "Feuil1") %>%
  mutate(geo_loc=region)

#vector of the French regions investigated
list_region = c("ARA","IDF","PACA")

for(replicat in 1:n_replicates){
	print(replicat)
	#import the tree
	my_tree = read.nexus(file = gsub(" ", "", paste(replicat,"_tree_dated.nexus")))
	#binarize the tree
	my_tree = multi2di(my_tree)
	#get tips and update get corresponding metadata
	my_tips = my_tree$tip.label
  
	metadata2 = metadata[metadata$gisaid_id %in% my_tips,]
  
	#count total of sequences per region
	meta2 = metadata2 %>%
		group_by(region) %>%
		summarize(total=n())
  
	#create the dataframe if not exists
	if(exists("final_data")){
		final_data = rbind(final_data, meta2)
	}else{
		final_data = data.frame(region=character(), total=integer())
		final_data = meta2
	}
}

#calculate the mean of sequences per region across the replicates
final_data2 = final_data %>%
	group_by(region) %>%
	summarize(mean=mean(total))

library(ggpubr)

#number of deaths for ARA, IDF and PACA by the end of the first wave
final_data2$epidemio_total = c(1804,7771,0,980)

#plot our results with epidemiological data
ggplot(final_data2, aes(x = epidemio_total, y = mean))+
	geom_point()+
	geom_smooth(method="lm", formula = y~x)+
	stat_cor(label.y.npc = 1, digits=10)+
	theme(axis.line.x = element_line(size = 1, colour = "black"),
        axis.line.y = element_line(size = 1, colour = "black"),
        axis.line = element_line(size=1, colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        text=element_text(size = 16, family="Arial"),
        axis.text.x=element_text(colour="black", size = 16),
        axis.text.y=element_text(colour="black", size = 16),
        legend.position = "none")

#import epidemiological for weekly analyses
epidemio = read.table(file = "epidemio_france.txt", sep="\t", header=T, quote="\"")

#we used the first replicate as an example
replicat=1
#import the tree
my_tree = read.nexus(file = gsub(" ", "", paste(replicat,"_tree_dated.nexus")))
#binarize the tree
my_tree = multi2di(my_tree)
#get tips and retain only metadata is tips
my_tips = my_tree$tip.label
metadata2 = metadata[metadata$gisaid_id %in% my_tips,]
#convert date as week
metadata2$week = week(metadata2$`Collection date`)

#number of sequences per region and week
metadata3 = metadata2 %>%
	group_by(region, week) %>%
	summarize(total=n())

#combine epidemio and our data
epidemio$nb_seq = NA
for(i in 1:nrow(epidemio)){
	for(j in 1:nrow(metadata3)){
		if(epidemio$week[i]==metadata3$week[j]){
			if(epidemio$region[i]==metadata3$region[j]){
				epidemio$nb_seq[i] = metadata3$total[j]
			}
		}
	}
}

#restrict epidemiological data to the wave 1
epidemio = epidemio[epidemio$week>=5,]
epidemio = epidemio[epidemio$week<=30,]

#remove empty data
for(i in 1:nrow(epidemio)){
	if(is.na(epidemio$nb_seq[i])){
		epidemio$nb_seq[i]=0
	}
}

library(data.table)
#calculate the sum per region of the sequences (cs) and the weekly deaths
setDT(epidemio)[, cs :=cumsum(nb_seq), region]
setDT(epidemio)[, nb_death_week :=cumsum(nb_deaths), region]

#plot the number of deaths per region over time
deaths_time = ggplot(epidemio,aes(x=week,y=nb_death_week, color=region)) +
	geom_point()+
	geom_line()+
	scale_color_manual(values=c("#FF4D31", "#2EC245", "#1790F5"))+
	theme(axis.line.x = element_line(size = 1, colour = "black"),
        axis.line.y = element_line(size = 1, colour = "black"),
        axis.line = element_line(size=1, colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        text=element_text(size = 16, family="Arial"),
        axis.text.x=element_text(colour="black", size = 16),
        axis.text.y=element_text(colour="black", size = 16),
        legend.position = "none")

#plot the number of sequences investigated per region over time
seq_time = ggplot(epidemio,aes(x=week,y=cs, color=region)) +
	geom_point()+
	geom_line()+
	scale_color_manual(values=c("#FF4D31", "#2EC245", "#1790F5"))+
	theme(axis.line.x = element_line(size = 1, colour = "black"),
        axis.line.y = element_line(size = 1, colour = "black"),
        axis.line = element_line(size=1, colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        text=element_text(size = 16, family="Arial"),
        axis.text.x=element_text(colour="black", size = 16),
        axis.text.y=element_text(colour="black", size = 16),
        legend.position = "none")

#view the plots
plot_grid(deaths_time, seq_time, nrow = 1, ncol = 2)


#for wave 2

n_replicates <- 100

#initialize a vector of samples
list_sample = c()

#import metadata
metadata <- readxl::read_xlsx("metadata_file4.xlsx", sheet = "Feuil1") %>%
  mutate(geo_loc=region)

#we focus only on the second wave
metadata = metadata[metadata$`Collection date`>="2020-07-20",]

#vector of the French regions investigated
list_region = c("ARA","IDF","PACA")

for(replicat in 1:n_replicates){
	print(replicat)
	#import the tree
	my_tree = read.nexus(file = gsub(" ", "", paste(replicat,"_tree_dated.nexus")))
	#binarize the tree
	my_tree = multi2di(my_tree)
	#get tips and update get corresponding metadata
	my_tips = my_tree$tip.label
  
	metadata2 = metadata[metadata$gisaid_id %in% my_tips,]
  
	#count total of sequences per region
	meta2 = metadata2 %>%
		group_by(region) %>%
		summarize(total=n())
  
	#create the dataframe if not exists
	if(exists("final_data")){
		final_data = rbind(final_data, meta2)
	}else{
		final_data = data.frame(region=character(), total=integer())
		final_data = meta2
	}
}

#calculate the mean of sequences per region across the replicates
final_data2 = final_data %>%
	group_by(region) %>%
	summarize(mean=mean(total))

library(ggpubr)

#number of deaths for ARA, IDF and PACA by the end of the first wave
final_data2$epidemio_total = c(5461,4916,2844)

#plot our results with epidemiological data
ggplot(final_data2, aes(x = epidemio_total, y = mean))+
	geom_point()+
	geom_smooth(method="lm", formula = y~x)+
	stat_cor(label.y.npc = 1, digits=10)+
	theme(axis.line.x = element_line(size = 1, colour = "black"),
        axis.line.y = element_line(size = 1, colour = "black"),
        axis.line = element_line(size=1, colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        text=element_text(size = 16, family="Arial"),
        axis.text.x=element_text(colour="black", size = 16),
        axis.text.y=element_text(colour="black", size = 16),
        legend.position = "none")

#import epidemiological for weekly analyses
epidemio = read.table(file = "epidemio_france.txt", sep="\t", header=T, quote="\"")

#we used the first replicate as an example
replicat=1
#import the tree
my_tree = read.nexus(file = gsub(" ", "", paste(replicat,"_tree_dated.nexus")))
#binarize the tree
my_tree = multi2di(my_tree)
#get tips and retain only metadata is tips
my_tips = my_tree$tip.label
metadata2 = metadata[metadata$gisaid_id %in% my_tips,]
#convert date as week
metadata2$week = week(metadata2$`Collection date`)

#number of sequences per region and week
metadata3 = metadata2 %>%
	group_by(region, week) %>%
	summarize(total=n())

#combine epidemio and our data
epidemio$nb_seq = NA
for(i in 1:nrow(epidemio)){
	for(j in 1:nrow(metadata3)){
		if(epidemio$week[i]==metadata3$week[j]){
			if(epidemio$region[i]==metadata3$region[j]){
				epidemio$nb_seq[i] = metadata3$total[j]
			}
		}
	}
}

#restrict epidemiological data to the wave 1
epidemio = epidemio[epidemio$week>=30,]
epidemio = epidemio[epidemio$week<=52,]

#remove empty data
for(i in 1:nrow(epidemio)){
	if(is.na(epidemio$nb_seq[i])){
		epidemio$nb_seq[i]=0
	}
}

library(data.table)
#calculate the sum per region of the sequences (cs) and the weekly deaths
setDT(epidemio)[, cs :=cumsum(nb_seq), region]
setDT(epidemio)[, nb_death_week :=cumsum(nb_deaths), region]

#plot the number of deaths per region over time
deaths_time = ggplot(epidemio,aes(x=week,y=nb_death_week, color=region)) +
	geom_point()+
	geom_line()+
	scale_color_manual(values=c("#FF4D31", "#2EC245", "#1790F5"))+
	theme(axis.line.x = element_line(size = 1, colour = "black"),
        axis.line.y = element_line(size = 1, colour = "black"),
        axis.line = element_line(size=1, colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        text=element_text(size = 16, family="Arial"),
        axis.text.x=element_text(colour="black", size = 16),
        axis.text.y=element_text(colour="black", size = 16),
        legend.position = "none")

#plot the number of sequences investigated per region over time
seq_time = ggplot(epidemio,aes(x=week,y=cs, color=region)) +
	geom_point()+
	geom_line()+
	scale_color_manual(values=c("#FF4D31", "#2EC245", "#1790F5"))+
	theme(axis.line.x = element_line(size = 1, colour = "black"),
        axis.line.y = element_line(size = 1, colour = "black"),
        axis.line = element_line(size=1, colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        text=element_text(size = 16, family="Arial"),
        axis.text.x=element_text(colour="black", size = 16),
        axis.text.y=element_text(colour="black", size = 16),
        legend.position = "none")

#view the plots
plot_grid(deaths_time, seq_time, nrow = 1, ncol = 2)