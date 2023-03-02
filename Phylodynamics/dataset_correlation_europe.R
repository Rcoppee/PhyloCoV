#clear workspace
rm(list = ls())

#set working directory
setwd("G:/phylocov/DEFINITIVE/Europe_W1")

#import all packages
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
  mutate(geo_loc=country)

#vector of the territories investigated
list_region=c("Belgium","France","Germany","Italy","Netherlands",
              "Poland","Romania","Russia","Spain","Sweden","United Kingdom")

for(replicat in 1:n_replicates){
  print(replicat)
  #import the tree
  my_tree = read.nexus(file = gsub(" ", "", paste(replicat, ".nexus")))
  #binarize the tree
  my_tree = multi2di(my_tree)
  #get tips and update get corresponding metadata
  my_tips = my_tree$tip.label
  
  metadata2 = metadata[metadata$gisaid_id %in% my_tips,]
  
  #count total of sequences per territory
  meta2 = metadata2 %>%
    group_by(geo_loc) %>%
    summarize(total=n())
  
  #create the dataframe if not exists
  if(exists("final_data")){
    final_data = rbind(final_data, meta2)
  }else{
    final_data = data.frame(geo_loc=character(), total=integer())
    final_data = meta2
  }
}

#calculate the mean of sequences per region across the replicates
final_data2 = final_data %>%
  group_by(geo_loc) %>%
  summarize(mean=mean(total))

library(ggpubr)

final_data2 = final_data2[final_data2$geo_loc!="NA",]
final_data2 = final_data2[final_data2$geo_loc!="China",]
#number of deaths for the different countries (alphabetical order)
final_data2$epidemio_total = c(9805,30182,9054,34938,6151,1562,1847,
                               11000,28403,5526,40866)

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
epidemio = read.table(file = "death_epidemio_data.txt", sep="\t", header=T, quote="\"")
epidemio = epidemio[epidemio$date<="2020-12-31",]
epidemio[is.na(epidemio)] <- 0

#calculate the number of deaths per day
epidemio_fix = epidemio %>%
  group_by(location, date) %>%
  summarise(nb_deaths=sum(total_deaths))

#remove possible NA data
epidemio_fix = epidemio_fix[epidemio_fix$location!="",]

#convert date to week
epidemio_fix$week = week(epidemio_fix$date)

#calculate the number of deaths per week
epidemio_fix2 = epidemio_fix %>%
  group_by(location, week) %>%
  summarize(nb_death_week=max(nb_deaths))

#we used the first replicate as an example
replicat=1
#import the tree
my_tree = read.nexus(file = gsub(" ", "", paste(replicat, ".nexus")))
#binarize the tree
my_tree = multi2di(my_tree)
#get tips and retain only metadata is tips
my_tips = my_tree$tip.label
metadata2 = metadata[metadata$gisaid_id %in% my_tips,]
#convert date as week
metadata2$week = week(metadata2$`Collection date`)

#number of sequences per region and week
metadata3 = metadata2 %>%
  group_by(country, week) %>%
  summarize(total=n())

#update United Kingdom name
for(i in 1:nrow(metadata3)){
  if(metadata3$country[i]=="United.Kingdom"){
    metadata3$country[i]="United Kingdom"
  }
}

#combine epidemio and our data
epidemio_fix2$nb_seq = NA
for(i in 1:nrow(epidemio_fix2)){
  for(j in 1:nrow(metadata3)){
    if(epidemio_fix2$week[i]==metadata3$week[j]){
      if(epidemio_fix2$location[i]==metadata3$country[j]){
        epidemio_fix2$nb_seq[i] = metadata3$total[j]
      }
    }
  }
}

#restrict epidemiological data to the wave 1
epidemio_fix2 = epidemio_fix2[epidemio_fix2$week>=5,]
epidemio_fix2 = epidemio_fix2[epidemio_fix2$week<=30,]

#remove empty data
for(i in 1:nrow(epidemio_fix2)){
  if(is.na(epidemio_fix2$nb_seq[i])){
    epidemio_fix2$nb_seq[i]=0
  }
}

epidemio_fix2 = filter(epidemio_fix2, location %in% list_region)

library(data.table)
#calculate the sum per region of the sequences (cs) and the weekly deaths
setDT(epidemio_fix2)[, cs :=cumsum(nb_seq), location]

#plot the number of deaths per region over time
deaths_time = ggplot(epidemio_fix2,aes(x=week,y=nb_death_week, color=location)) +
  geom_point()+
  geom_line()+
  scale_color_manual(values=c("#FDD93B", "#57C540", "#FD832F", "#558EC8","#872BBA", "#813008", "#C485E6",
                              "#16999F", "#10791E", "#867017", "#ED313C"))+
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
seq_time = ggplot(epidemio_fix2,aes(x=week,y=cs, color=location)) +
  geom_point()+
  geom_line()+
  scale_color_manual(values=c("#FDD93B", "#57C540", "#FD832F", "#558EC8","#872BBA", "#813008", "#C485E6",
                              "#16999F", "#10791E", "#867017", "#ED313C"))+
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






#clear workspace
rm(list = ls())

#set working directory
setwd("G:/phylocov/DEFINITIVE/Europe_W2")

#import all packages
list.of.packages <- c("adephylo","tidyverse", "lubridate", "glue",
                      "ape", "phytools", "tidytree", "ggplot2",
                      "treeio", "data.table", "dplyr",
                      "parallel", "foreach", "doParallel", "sqldf",
                      "ggtree", "zoo")

new.packages <- list.of.packages[!(list.of.packages %in%
                                     installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages)

lapply(list.of.packages, require, character.only = TRUE)

#for wave 2

n_replicates <- 100

#initialize a vector of samples
list_sample = c()

#import metadata
metadata <- readxl::read_xlsx("metadata_file4.xlsx", sheet = "Feuil1") %>%
  mutate(geo_loc=country)

#vector of the territories investigated
list_region=c("Belgium","France","Germany","Italy","Netherlands",
              "Poland","Romania","Russia","Spain","Sweden","United Kingdom")

for(replicat in 1:n_replicates){
  print(replicat)
  #import the tree
  my_tree = read.nexus(file = gsub(" ", "", paste(replicat, ".nexus")))
  #binarize the tree
  my_tree = multi2di(my_tree)
  #get tips and update get corresponding metadata
  my_tips = my_tree$tip.label
  
  metadata2 = metadata[metadata$gisaid_id %in% my_tips,]
  
  #count total of sequences per territory
  meta2 = metadata2 %>%
    group_by(geo_loc) %>%
    summarize(total=n())
  
  #create the dataframe if not exists
  if(exists("final_data")){
    final_data = rbind(final_data, meta2)
  }else{
    final_data = data.frame(geo_loc=character(), total=integer())
    final_data = meta2
  }
}

#calculate the mean of sequences per region across the replicates
final_data2 = final_data %>%
  group_by(geo_loc) %>%
  summarize(mean=mean(total))

library(ggpubr)

final_data2 = final_data2[final_data2$geo_loc!="NA",]
final_data2 = final_data2[final_data2$geo_loc!="China",]

#number of deaths for the different countries (alphabetical order)
final_data2$epidemio_total = c(9726,34462,24017,38221,5308,26992,
                               13920,45271,22434,3201,32703)

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
epidemio = read.table(file = "death_epidemio_data.txt", sep="\t", header=T, quote="\"")
epidemio = epidemio[epidemio$date<="2020-12-31",]
epidemio[is.na(epidemio)] <- 0

#calculate the number of deaths per day
epidemio_fix = epidemio %>%
  group_by(location, date) %>%
  summarise(nb_deaths=sum(total_deaths))

#remove possible NA data
epidemio_fix = epidemio_fix[epidemio_fix$location!="",]

#convert date to week
epidemio_fix$week = week(epidemio_fix$date)

#calculate the number of deaths per week
epidemio_fix2 = epidemio_fix %>%
  group_by(location, week) %>%
  summarize(nb_death_week=max(nb_deaths))

#we used the first replicate as an example
replicat=1
#import the tree
my_tree = read.nexus(file = gsub(" ", "", paste(replicat, ".nexus")))
#binarize the tree
my_tree = multi2di(my_tree)
#get tips and retain only metadata is tips
my_tips = my_tree$tip.label
metadata2 = metadata[metadata$gisaid_id %in% my_tips,]
#convert date as week
metadata2$week = week(metadata2$`Collection date`)

#number of sequences per region and week
metadata3 = metadata2 %>%
  group_by(country, week) %>%
  summarize(total=n())

#update United Kingdom name
for(i in 1:nrow(metadata3)){
  if(metadata3$country[i]=="United.Kingdom"){
    metadata3$country[i]="United Kingdom"
  }
}

#combine epidemio and our data
epidemio_fix2$nb_seq = NA
for(i in 1:nrow(epidemio_fix2)){
  for(j in 1:nrow(metadata3)){
    if(epidemio_fix2$week[i]==metadata3$week[j]){
      if(epidemio_fix2$location[i]==metadata3$country[j]){
        epidemio_fix2$nb_seq[i] = metadata3$total[j]
      }
    }
  }
}

#restrict epidemiological data to the wave 2
epidemio_fix2 = epidemio_fix2[epidemio_fix2$week>=30,]

#remove empty data
for(i in 1:nrow(epidemio_fix2)){
  if(is.na(epidemio_fix2$nb_seq[i])){
    epidemio_fix2$nb_seq[i]=0
  }
}

epidemio_fix2 = filter(epidemio_fix2, location %in% list_region)


#we restart the number of deaths at 0 from week 30
epidemio_fix2$maj_death=0
for(i in 2:nrow(epidemio_fix2)){
  if(epidemio_fix2$location[i]==epidemio_fix2$location[i-1]){
    epidemio_fix2$maj_death[i]=epidemio_fix2$nb_death_week[i]-epidemio_fix2$nb_death_week[i-1]
  }else{
    epidemio_fix2$maj_death[i]=0
  }
}

library(data.table)
#calculate the sum per region of the sequences (cs) and the weekly deaths
setDT(epidemio_fix2)[, cum_death :=cumsum(maj_death), location]
setDT(epidemio_fix2)[, cs :=cumsum(nb_seq), location]


#plot the number of deaths per region over time
deaths_time = ggplot(epidemio_fix2,aes(x=week,y=cum_death, color=location)) +
  geom_point()+
  geom_line()+
  scale_color_manual(values=c("#FDD93B", "#57C540", "#FD832F", "#558EC8","#872BBA", "#813008", "#C485E6",
                              "#16999F", "#10791E", "#867017", "#ED313C"))+
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
seq_time = ggplot(epidemio_fix2,aes(x=week,y=cs, color=location)) +
  geom_point()+
  geom_line()+
  scale_color_manual(values=c("#FDD93B", "#57C540", "#FD832F", "#558EC8","#872BBA", "#813008", "#C485E6",
                              "#16999F", "#10791E", "#867017", "#ED313C"))+
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
