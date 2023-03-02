#clear workspace
rm(list = ls())

# set working directory
setwd("F:/phylocov/DEFINITIVE/World_W1")

# import all packages
library("tidyverse")
library("dplyr")
library("ape")
library("ggtree")
library("tidytree")
library("lubridate")
library("ggplot2")
library("cowplot")
library("circlize")

##############################
##########FUNCTIONS###########
##############################

#function to get all transition events based on location and variants
get_trans <- function(phylo, tip_loc, ace, tMRCA, ace2, tip_clade) {
  tip_loc <- as.character(tip_loc)
  df_phylo <- as_tibble(phylo)
  
  #determine the location and lineages of the nodes
  node_loc <- apply(ace$lik.anc, MARGIN = 1, function(x) names(x)[x == max(x)])
  node_clade <- apply(ace2$lik.anc, MARGIN = 1, function(x) names(x)[x == max(x)])
  
  #update the tibble with the results
  df_phylo$node_loc <- 
    c(tip_loc, node_loc)
  
  df_phylo$node_clade <- 
    c(tip_clade, node_clade)
  
  df_phylo$date <- 
    dist.nodes(phylo)[Ntip(phylo) + 1, ]
  
  df_phylo$parent_loc <-
    df_phylo$node_loc[df_phylo$parent]
  
  df_phylo$parent_date <-
    df_phylo$date[df_phylo$parent]
  
  #determine inter or intra event
  df_phylo$trans <- 
    if_else(df_phylo$node_loc == df_phylo$parent_loc, 'INTRA', 'INTER')
  
  #define the middle of the branch as the date of the event
  df_phylo$trans_date <-
    df_phylo$parent_date + (df_phylo$branch.length / 2)
  
  df_phylo$parent_clade <-
    df_phylo$node_clade[df_phylo$parent]
  
  #return the table of results
  return(df_phylo)
}

#function to determine items not present in a vector
'%!in%'<- function(x,y)!('%in%'(x,y))


##############################
##CALCULATE TRANSITION EVENTS#
##############################

#import metadata using readxl
metadata <- readxl::read_xlsx("metadata_file4.xlsx", sheet = "Feuil1") %>%
  mutate(geo_loc = continent)

#fix the number of replicates
n_replicates <- 100

#read a list of tmrca (one tmrca per replicate)
list_tmrca <- read.table("tmrca.txt", header = T, sep="\t")

for(replicat in 1:n_replicates){
  
  print(replicat)
  #import tree
  my_tree = read.nexus(file = gsub(" ", "", paste(replicat, ".nexus")))
  #resolve all nodes
  my_tree = multi2di(my_tree)
  
  my_tips = my_tree$tip.label
  
  #keep only data from the tree
  metadata_sub = metadata[metadata$gisaid_id %in% my_tips,]
  
  #add France in the continent category
  for(i in 1:nrow(metadata_sub)){
    if(metadata_sub$country[i]=="France"){
      metadata_sub$geo_loc[i]="France"
    }
  }
  
  #get the tmrca of the current replicate
  tMRCA = Date2decimal(as.Date(list_tmrca$tmrca_date[replicat]))

  #generate a tibble of metadata
  tree_table <- as_tibble(my_tree) %>%
    left_join(select(metadata_sub, c('gisaid_id', 'geo_loc', 'Pango_lineage')), by = c('label' = 'gisaid_id'))
  
  #produce a treedata of metadata
  as.treedata(tree_table) -> tree_data
  
  #calculate distances from the root
  root_node <- max(tree_table$node[tree_table$node <= Ntip(tree_data)])+1
  distances = dist.nodes(as.phylo(tree_table))
  distance_root = distances[,root_node]
  distance_root = as.data.frame(distance_root)
  
  #convert distances to dates
  for (i in 1:nrow(distance_root)) {
    distance_root[i,2] = (decimal2Date(tMRCA+distance_root[i,1]))
  }
  colnames(distance_root)[2] <- "dates"
  tree_table$new_date = distance_root$dates
  
  #get the location and pango of all tips
  tip_loc = tree_table$geo_loc[1:Ntip(my_tree)]
  tip_clade = tree_table$Pango_lineage[1:Ntip(my_tree)]
  
  #Solving the branch with length 0
  zero_offset <- 1e-2/365
  is_zero_branch <- tree_table$branch.length == 0
  tree_table$branch.length[is_zero_branch] <- zero_offset
  summary(tree_table$branch.length)
  
  #ancestrat state reconstruction, here using the ER model for location and for lineage
  ace <- ace(x = tree_table$geo_loc[1:Ntip(my_tree)], phy = as.phylo(tree_table), type = 'discrete', model = 'ER')
  ace2 <- ace(x = tree_table$Pango_lineage[1:Ntip(my_tree)], phy = as.phylo(tree_table), type = 'discrete', model = 'ER')
  
  #get all the transition events
  df_transmission = get_trans(my_tree, tip_loc, ace, tMRCA, ace2, tip_clade)
  
  #get the dates from the tree table
  df_transmission$old_dates = tree_table$new_date
  df_transmission$trans_date[Ntip(my_tree)+1] = 0
  
  #update the dates using the length [root to current node] and the tmrca
  df_transmission$new_date=as.Date("0000-01-01")
  for (i in 1:nrow(df_transmission)) {
    df_transmission$new_date[i] = (decimal2Date(tMRCA+df_transmission$trans_date[i]))
  }
  
  #if it is the first replicate, generate a table, else, concatenate the data
  if(replicat==1){
    df_transmission$replicat = replicat
    total_transmission = df_transmission
  }else{
    df_transmission$replicat = replicat
    total_transmission = rbind(total_transmission, df_transmission)
  }
  
}           


##############################
####VARIATION IN REPLICATES###
##############################

#Get all inter- and intra-region events in separate variables
df_intra_SYM = total_transmission[total_transmission$trans=="INTRA",]
df_inter_SYM = total_transmission[total_transmission$trans=="INTER",]

#list of the regions included in the study
list_continent = c("Asia", "Europe", "France", "North.America", "South.America", "Africa", "Oceania")

#We first explore variation in replicates
df_inter_count_SYM = data.frame(Country=character(), In=integer(), Out=integer(), Replicat=integer())

#we count the total on introduction and exportation events per territory
for(i in 1:n_replicates){
  for(j in list_continent){
    sub_df_inter_SYM_IN = df_inter_SYM[df_inter_SYM$parent_loc==j,]
    sub_df_inter_SYM_IN = sub_df_inter_SYM_IN[sub_df_inter_SYM_IN$replicat==i,]
    sub_df_inter_SYM_OUT = df_inter_SYM[df_inter_SYM$node_loc==j,]
    sub_df_inter_SYM_OUT = sub_df_inter_SYM_OUT[sub_df_inter_SYM_OUT$replicat==i,]   
    df_inter_count_SYM[nrow(df_inter_count_SYM) + 1,] = list(j, dim(sub_df_inter_SYM_IN)[1], dim(sub_df_inter_SYM_OUT)[1], i)
  }
}

#plot introduction vs exportation transmissions across replicates
ggplot(df_inter_count_SYM)+
  geom_point(aes(x=In, y=Out, fill=Country), shape=21, size=3, color="black", alpha=1)+
  scale_fill_manual(values=c("#EF3437", "#5F4197", "#538DCA", "#54C538", "#1A7E41", "#FF8421", "#FEDA27"))+
  xlim(0,130)+
  ylim(0,130)+
  xlab("exportation events")+
  ylab("introduction events")+
  geom_abline(intercept = 0, slope = 1, size = 1)+
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

#calculate the mean across the replicates, then plot it
inter_count_mean = df_inter_count_SYM %>%
  group_by(Country) %>%
  summarize(In_m=mean(In),Out_m=mean(Out))

ggplot(inter_count_mean) +
  geom_point(aes(x=In_m, y=Out_m, fill=Country), shape=21, size=10, color="black", alpha=1)+
  scale_fill_manual(values=c("#EF3437", "#5F4197", "#538DCA", "#54C538", "#1A7E41", "#FF8421", "#FEDA27"))+
  xlim(0,130) +
  ylim(0,130) +
  xlab("exportation events")+
  ylab("introduction events")+
  geom_abline(intercept = 0, slope = 1, size = 1)+
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

#initialization of a dataframe
over_replicates = df_inter_count_SYM[order(df_inter_count_SYM$Country, df_inter_count_SYM$Replicat), ]
current_replicat=1

#calculate mean and sd from all the previous replicates
#introduction analysis
for(i in 1:nrow(over_replicates)){
  if(over_replicates$Replicat[i]==1){
    first_r = i
    over_replicates$in_mean[i] = over_replicates$In[i]
    over_replicates$in_mean_sd[i] = 0
    current_replicat = current_replicat+1
  }else{
    over_replicates$in_mean[i] = mean(over_replicates$In[first_r:current_replicat])
    over_replicates$in_mean_sd[i] = sd(over_replicates$in_mean[first_r:current_replicat])
    current_replicat = current_replicat+1
  }
}

current_replicat=1
#exportation analysis
for(i in 1:nrow(over_replicates)){
  if(over_replicates$Replicat[i]==1){
    first_r = i
    over_replicates$out_mean[i] = over_replicates$Out[i]
    over_replicates$out_mean_sd[i] = 0
    current_replicat = current_replicat+1 
  }else{
    over_replicates$out_mean[i] = mean(over_replicates$Out[first_r:current_replicat])
    over_replicates$out_mean_sd[i] = sd(over_replicates$out_mean[first_r:current_replicat])
    current_replicat = current_replicat+1
  } 
}

#plot the tendency of introduction events over replicates
p1=ggplot(over_replicates, aes(x=Replicat, y=in_mean, color=Country)) +
  geom_point()+
  geom_errorbar(aes(ymin=in_mean-in_mean_sd, ymax=in_mean+in_mean_sd, alpha=.5), width=.1) +
  scale_color_manual(values=c("#EF3437", "#5F4197", "#538DCA", "#54C538", "#1A7E41", "#FF8421", "#FEDA27"))+
  geom_line()+
  ylim(0,150) +
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

#plot the tendency of exportation events over replicates
p2=ggplot(over_replicates, aes(x=Replicat, y=out_mean, color=Country)) +
  geom_point()+
  geom_errorbar(aes(ymin=out_mean-out_mean_sd, ymax=out_mean+out_mean_sd, alpha=.5), width=.1) +
  scale_color_manual(values=c("#EF3437", "#5F4197", "#538DCA", "#54C538", "#1A7E41", "#FF8421", "#FEDA27"))+
  geom_line()+
  ylim(0,150) +
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

#viewing the two plots
plot_grid(p1, p2, nrow = 1, ncol = 2)


##############################
#####CIRCOS VISUALISATION#####
##############################

#calculate frequency per replicat, in and out (then for all transmissions)
flux_m_SYM = rename(count(df_inter_SYM, parent_loc, node_loc), Freq = n)

#create a dataframe for circos plot
circos_SYM = flux_m_SYM %>% group_by(parent_loc,node_loc) %>% summarize(B = sum(Freq))
circos_SYM = as.data.frame(circos_SYM)

#remove any incomplete rows
circos_SYM = circos_SYM[complete.cases(circos_SYM), ]

#remove intra_transmission if present
for (row in 1:nrow(circos_SYM)){
  if (circos_SYM$parent_loc[row]==circos_SYM$node_loc[row]){
    circos_SYM = circos_SYM[-row,]
  }
}

#calculate the percentage of each flow (for map production)
#circos_SYM$perc=circos_SYM$B/sum(circos_SYM$B)*100
#write.table(circos_SYM, "for_map.txt", sep="\t")

#if we want to include intra-territory transmission events, use also these lines
# df_SYM = rbind(df_inter_SYM, df_intra_SYM)
# flux_m_SYM_all = rename(count(df_SYM, parent_loc, node_loc, replicat), Freq = n)
# circos_SYM = flux_m_SYM_all %>% group_by(parent_loc, node_loc) %>% summarize(B = sum(Freq))
# circos_SYM = as.data.frame(circos_SYM)
# circos_SYM = circos_SYM[complete.cases(circos_SYM), ]

#parameters of circos plot
circos.clear()
circos.par(start.degree = 90, gap.degree = 1, track.margin = c(-0.1, 0.1), points.overflow.warning = FALSE)
par(mar = rep(0, 4))
# color palette : africa, asia,europe,france, northamerica, southamerica, oceania
mycolor <- c("#EF3437", "#5F4197", "#538DCA", "#54C538", "#1A7E41", "#FEDA27", "#FF8421")

#Base plot
chordDiagram(
  x = circos_SYM,
  grid.col = mycolor,
  transparency = 0.25,
  directional = 1,
  direction.type = c("arrows", "diffHeight"), 
  diffHeight  = -0.04,
  annotationTrack = "grid", 
  annotationTrackHeight = c(0.1, 1),
  link.arr.type = "big.arrow", 
  link.sort = TRUE, 
  link.largest.ontop = TRUE,
  grid.border	= "black"
)

#Add text and axis
circos.trackPlotRegion(
  track.index = 1, 
  bg.border = NA, 
  panel.fun = function(x, y) {
    xlim = get.cell.meta.data("xlim")
    sector.index = get.cell.meta.data("sector.index")
    #Add names to the sector. 
    circos.text(
      x = mean(xlim), 
      y = 1.6,
      labels = sector.index, 
      facing = "bending.inside", 
      adj = par("adj"),
      font= 2,
      cex = 1)
  }
)


##############################
##TRANSITION EVENTS OVER TIME#
##############################

#introduction analysis
current_replicat = 1

#pipeline to calculate the mean over the replicates (per week)
for(continent in list_continent){
  print(continent)
  
  #define the transition for the current region, and calculate the week and year for each event
  mean_df_inter_SYM = df_inter_SYM[df_inter_SYM$node_loc==continent,]
  mean_df_inter_SYM$week = strftime(mean_df_inter_SYM$new_date, format = "%V")
  mean_df_inter_SYM$week = as.integer(mean_df_inter_SYM$week)
  mean_df_inter_SYM$year = NA
  
  #for wave_1, week>=40 suggests the year 2019, since the max date for 2020 was week 30
  for(list_week in 1:nrow(mean_df_inter_SYM)){
    if(mean_df_inter_SYM$week[list_week]>=40){
      mean_df_inter_SYM$year[list_week]=2019
    }else{
      mean_df_inter_SYM$year[list_week]=2020
    }
  }
  
  #regroup the data per week, year and replicat
  mean_df_inter_SYM2 <- mean_df_inter_SYM %>% 
    group_by(week, year, replicat) %>%
    summarise(total=n())
  
  #we complete every lacking data for each week, region and replicate
  for(a_week in 1:30){
    sub_table = mean_df_inter_SYM2[mean_df_inter_SYM2$week==a_week,]
    if(nrow(sub_table)==0){
      for(j in 1:n_replicates){
        mean_df_inter_SYM2[nrow(mean_df_inter_SYM2) + 1,] = list(a_week, 2020, j, 0)
      }
    }else{
      for(i in 1:n_replicates){
        if(i %!in% sub_table$replicat){
          mean_df_inter_SYM2[nrow(mean_df_inter_SYM2) + 1,] = list(a_week, 2020, i, 0)
        }
      }
    }
  }
  
  #we update the frequency for each week and year
  mean_df_inter_SYM3 <- mean_df_inter_SYM2 %>%
    group_by(week, year) %>%
    summarise(total2=(sum(total)/n_replicates))
  
  mean_df_inter_SYM3 <- mean_df_inter_SYM3[order(mean_df_inter_SYM3$year, mean_df_inter_SYM3$week),]
  
  mean_df_inter_SYM3$region = continent
  mean_df_inter_SYM3$cum_freq <- cumsum(mean_df_inter_SYM3$total2)
  
  mean_df_inter_SYM3$week_2 <- as.Date("0000-01-01")
  
  #we convert the week as a date
  for(week2 in 1:nrow(mean_df_inter_SYM3)){
    if(mean_df_inter_SYM3$year[week2]==2019){
      mean_df_inter_SYM3$week_2[week2] = as.Date(lubridate::ymd("2019-01-01") + lubridate::weeks(mean_df_inter_SYM3$week[week2] - 1 ))
    }else{
      mean_df_inter_SYM3$week_2[week2] = as.Date(lubridate::ymd("2020-01-01") + lubridate::weeks(mean_df_inter_SYM3$week[week2] - 1 ))
    }
  }
  
  #concatenation of the results
  if(current_replicat == 1){
    mean_transmission = mean_df_inter_SYM3
  }else{
    mean_transmission = rbind(mean_transmission, mean_df_inter_SYM3)
  }
  
  current_replicat = current_replicat + 1
}

#we order the data for the plot
mean_transmission$region %<>% factor(levels=c("North.America", "South.America", "Europe",
                                              "Asia" , "France", "Africa", "Oceania"))

#stacked area chart
ggplot(mean_transmission, aes(x=week_2, y=total2)) + 
  geom_area(aes(fill=region), size=.5, colour="black", position="stack")+
  scale_fill_manual(values = c("#1A7E41", "#FEDA27","#538DCA", "#5F4197", "#54C538","#EF3437","#FF8421"))+
  theme(axis.line.x = element_line(size = 1, colour = "black"),
        axis.line.y = element_line(size = 1, colour = "black"),
        axis.line = element_line(size=1, colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        text=element_text(size = 16, family="Arial"),
        axis.text.x=element_text(colour="black", size = 12, angle = 45, vjust=0.6),
        axis.text.y=element_text(colour="black", size = 12),
        legend.position = "none")+
  ylim(0,30)+
  scale_x_date(date_breaks = "1 month", limits = as.Date(c('2019-10-01','2020-07-01')))+
  ylab("Number of introductions")+
  xlab("Date")

#exportation analysis
current_replicat = 1

#pipeline to calculate the mean over the replicates (per week)
for(continent in list_continent){
  print(continent)
  
  #define the transition for the current region, and calculate the week and year for each event
  mean_df_inter_SYM = df_inter_SYM[df_inter_SYM$parent_loc==continent,]
  mean_df_inter_SYM$week = strftime(mean_df_inter_SYM$new_date, format = "%V")
  mean_df_inter_SYM$week = as.integer(mean_df_inter_SYM$week)
  mean_df_inter_SYM$year = NA
  
  #for wave_1, week>=40 suggests the year 2019, since the max date for 2020 was week 30
  for(list_week in 1:nrow(mean_df_inter_SYM)){
    if(mean_df_inter_SYM$week[list_week]>=40){
      mean_df_inter_SYM$year[list_week]=2019
    }else{
      mean_df_inter_SYM$year[list_week]=2020
    }
  }
  
  #regroup the data per week, year and replicat
  mean_df_inter_SYM2 <- mean_df_inter_SYM %>% 
    group_by(week, year, replicat) %>%
    summarise(total=n())
  
  #we complete every lacking data for each week, region and replicate
  for(a_week in 1:30){
    sub_table = mean_df_inter_SYM2[mean_df_inter_SYM2$week==a_week,]
    if(nrow(sub_table)==0){
      for(j in 1:n_replicates){
        mean_df_inter_SYM2[nrow(mean_df_inter_SYM2) + 1,] = list(a_week, 2020, j, 0)
      }
    }else{
      for(i in 1:n_replicates){
        if(i %!in% sub_table$replicat){
          mean_df_inter_SYM2[nrow(mean_df_inter_SYM2) + 1,] = list(a_week, 2020, i, 0)
        }
      }
    }
  }
  
  #we update the frequency for each week and year
  mean_df_inter_SYM3 <- mean_df_inter_SYM2 %>%
    group_by(week, year) %>%
    summarise(total2=(sum(total)/n_replicates))
  
  mean_df_inter_SYM3 <- mean_df_inter_SYM3[order(mean_df_inter_SYM3$year, mean_df_inter_SYM3$week),]
  
  mean_df_inter_SYM3$region = continent
  mean_df_inter_SYM3$cum_freq <- cumsum(mean_df_inter_SYM3$total2)
  
  mean_df_inter_SYM3$week_2 <- as.Date("0000-01-01")
  
  #we convert the week as a date
  for(week2 in 1:nrow(mean_df_inter_SYM3)){
    if(mean_df_inter_SYM3$year[week2]==2019){
      mean_df_inter_SYM3$week_2[week2] = as.Date(lubridate::ymd("2019-01-01") + lubridate::weeks(mean_df_inter_SYM3$week[week2] - 1 ))
    }else{
      mean_df_inter_SYM3$week_2[week2] = as.Date(lubridate::ymd("2020-01-01") + lubridate::weeks(mean_df_inter_SYM3$week[week2] - 1 ))
    }
  }
  
  #concatenation of the results
  if(current_replicat == 1){
    mean_transmission = mean_df_inter_SYM3
  }else{
    mean_transmission = rbind(mean_transmission, mean_df_inter_SYM3)
  }
  
  current_replicat = current_replicat + 1
}

#we order the data for the plot
mean_transmission$region %<>% factor(levels=c("Europe", "North.America","Asia",
                                              "South.America", "France","Africa"))

# stacked area chart
ggplot(mean_transmission, aes(x=week_2, y=total2)) + 
  geom_area(aes(fill=region), size=.5, colour="black", position="stack")+
  #scale_fill_manual(values = c("#EF3437", "#54C538", "#FEDA27", "#5F4197", "#1A7E41", "#538DCA"))+
  scale_fill_manual(values = c("#538DCA", "#1A7E41", "#5F4197", "#FEDA27", "#54C538", "#EF3437"))+
  theme(axis.line.x = element_line(size = 1, colour = "black"),
        axis.line.y = element_line(size = 1, colour = "black"),
        axis.line = element_line(size=1, colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        text=element_text(size = 16, family="Arial"),
        axis.text.x=element_text(colour="black", size = 12, angle = 45, vjust=0.6),
        axis.text.y=element_text(colour="black", size = 12),
        legend.position = "none")+
  ylim(0,30)+
  scale_x_date(date_breaks = "1 month", limits = as.Date(c('2019-10-01','2020-07-01')))+
  ylab("Number of exportations")+
  xlab("Date")


##############################
#########FRANCE STUDY#########
##############################

#count the number of identical events for a same replicat (without dating)
df_inter_SYM$number = 1
df_inter_SYM2 = df_inter_SYM %>%
  group_by(parent_loc, node_loc, replicat) %>%
  summarise(a_sum=sum(number), a_mean=(mean(number)))

#convert to dataframe
df_inter_SYM2 = as.data.frame(df_inter_SYM2)

#inialisation dataframe
completion = data.frame(parent_loc=integer(), node_loc=character(), replicat=integer(), a_sum=integer(), a_mean=integer())

#update the transmission table when no transmission between two territories is observed
for (payss in list_continent) {
  print(payss)
  temp_SYM = df_inter_SYM2[df_inter_SYM2$parent_loc == payss,]
  for (payss_out in list_continent){
    if(payss==payss_out){
      next
    } else {
      for(i in 1:n_replicates){
        if (i %!in% temp_SYM[temp_SYM[,2] == payss_out,]$replicat){
          completion[nrow(completion) + 1,] = list(payss, payss_out, i,0,0)
        }
      }
    }
  }
}
df_inter_SYM2 = rbind (df_inter_SYM2, completion)

#visualizing introduction events to France
France_SYM = df_inter_SYM2[df_inter_SYM2$node_loc == "France",]
FR_SYM_intro = ggplot(France_SYM, aes(x=parent_loc, y=a_sum, fill=parent_loc)) + 
  geom_boxplot(outlier.shape = NA) +
  scale_fill_manual(values=c("#EF3437", "#5F4197", "#538DCA", "#1A7E41", "#FF8421", "#FEDA27"))+
  geom_line()+
  theme(axis.line.x = element_line(size = 1, colour = "black"),
        axis.line.y = element_line(size = 1, colour = "black"),
        axis.line = element_line(size=1, colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        text=element_text(size = 16, family="Arial"),
        axis.text.x=element_text(colour="black", size = 16, angle = 45),
        axis.text.y=element_text(colour="black", size = 16),
        legend.position = "none")+
  ylab("Number of introductions to France")+
  xlab("Country")

#view plot
FR_SYM_intro

#visualizing exportation events from France
France_SYM = df_inter_SYM2[df_inter_SYM2$parent_loc == "France",]
FR_SYM_exp = ggplot(France_SYM, aes(x=node_loc, y=a_sum, fill=node_loc)) + 
  geom_boxplot(outlier.shape = NA) +
  scale_fill_manual(values=c("#EF3437", "#5F4197", "#538DCA", "#1A7E41", "#FF8421", "#FEDA27"))+
  geom_line()+
  theme(axis.line.x = element_line(size = 1, colour = "black"),
        axis.line.y = element_line(size = 1, colour = "black"),
        axis.line = element_line(size=1, colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        text=element_text(size = 16, family="Arial"),
        axis.text.x=element_text(colour="black", size = 16, angle = 45),
        axis.text.y=element_text(colour="black", size = 16),
        legend.position = "none")+
  ylab("Number of exportations from France")+
  xlab("Continent")

#view plot
FR_SYM_exp

#introduction events over time
df_inter_SYM_fr = df_inter_SYM[df_inter_SYM$node_loc=="France",]
mean_df_inter_SYM = df_inter_SYM_fr

#calculate the week and year for each event
mean_df_inter_SYM$week = strftime(mean_df_inter_SYM$new_date, format = "%V")
mean_df_inter_SYM$week = as.integer(mean_df_inter_SYM$week)
mean_df_inter_SYM$year = NA

#for wave_1, week>=40 suggests the year 2019, since the max date for 2020 was week 30
for(list_week in 1:nrow(mean_df_inter_SYM)){
  if(mean_df_inter_SYM$week[list_week]>=40){
    mean_df_inter_SYM$year[list_week]=2019
  }else{
    mean_df_inter_SYM$year[list_week]=2020
  }
}

#regroup the data per week, year and replicat
mean_df_inter_SYM2 <- mean_df_inter_SYM %>% 
  group_by(week, year, replicat, parent_loc) %>%
  summarise(total=n())

mean_df_inter_SYM2$year = as.integer(mean_df_inter_SYM2$year)

#location exluding France
locations = c("North.America", "South.America", "Asia", "Europe", "Africa", "Oceania")
mean_df_inter_SYM2 = as.data.frame(mean_df_inter_SYM2)

#we complete every lacking data for each week, region and replicate
for(a_location in locations){
  print(a_location)
  for(a_week in 1:30){
    sub_table = mean_df_inter_SYM2[mean_df_inter_SYM2$week==a_week,]
    sub_table = sub_table[sub_table$parent_loc==a_location,]
    if(nrow(sub_table)==0){
      for(j in 1:n_replicates){
        mean_df_inter_SYM2[nrow(mean_df_inter_SYM2) + 1,] = list(a_week, 2020, j, a_location, 0)
      }
    }else{
      for(i in 1:n_replicates){
        if(i %!in% sub_table$replicat){
          mean_df_inter_SYM2[nrow(mean_df_inter_SYM2) + 1,] = list(a_week, 2020, i, a_location, 0)
        }
      }
    }
  }  
}

#we update the frequency for each week and year
mean_df_inter_SYM3 <- mean_df_inter_SYM2 %>%
  group_by(week, year, parent_loc) %>%
  summarise(total2=(sum(total)/n_replicates))

mean_df_inter_SYM3 <- mean_df_inter_SYM3[order(mean_df_inter_SYM3$year, mean_df_inter_SYM3$week),]

mean_df_inter_SYM3$week_2 <- as.Date("0000-01-01")

#we convert the week as a date
for(week2 in 1:nrow(mean_df_inter_SYM3)){
  if(mean_df_inter_SYM3$year[week2]==2019){
    mean_df_inter_SYM3$week_2[week2] = as.Date(lubridate::ymd("2019-01-01") + lubridate::weeks(mean_df_inter_SYM3$week[week2] - 1 ))
  }else{
    mean_df_inter_SYM3$week_2[week2] = as.Date(lubridate::ymd("2020-01-01") + lubridate::weeks(mean_df_inter_SYM3$week[week2] - 1 ))
  }
}

#we order the data for the plot
mean_df_inter_SYM3$parent_loc %<>% factor(levels=c("North.America","Europe", "Asia",  "South.America", "Africa", "Oceania"))

#stacked area chart
ggplot(mean_df_inter_SYM3, aes(x=week_2, y=total2)) + 
  geom_area(aes(fill=parent_loc), size=.5, colour="black", position="stack")+
  #scale_fill_manual(values = c("#FF8421", "#EF3437", "#FEDA27", "#5F4197", "#538DCA", "#1A7E41"))+
  scale_fill_manual(values = c("#1A7E41", "#538DCA", "#5F4197", "#FEDA27", "#EF3437", "#FF8421"))+
  theme(axis.line.x = element_line(size = 1, colour = "black"),
        axis.line.y = element_line(size = 1, colour = "black"),
        axis.line = element_line(size=1, colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        text=element_text(size = 16, family="Arial"),
        axis.text.x=element_text(colour="black", size = 12, angle = 45, vjust=0.6),
        axis.text.y=element_text(colour="black", size = 12),
        legend.position = "none")+
  ylim(0,5)+
  scale_x_date(date_breaks = "1 month", limits = as.Date(c('2020-01-01','2020-07-31')))+
  ylab("Number of introductions")+
  xlab("Date")

#exportation events over time
df_inter_SYM_fr = df_inter_SYM[df_inter_SYM$parent_loc=="France",]
mean_df_inter_SYM = df_inter_SYM_fr

#calculate the week and year for each event
mean_df_inter_SYM$week = strftime(mean_df_inter_SYM$new_date, format = "%V")
mean_df_inter_SYM$week = as.integer(mean_df_inter_SYM$week)
mean_df_inter_SYM$year = NA

#for wave_1, week>=40 suggests the year 2019, since the max date for 2020 was week 30
for(list_week in 1:nrow(mean_df_inter_SYM)){
  if(mean_df_inter_SYM$week[list_week]>=40){
    mean_df_inter_SYM$year[list_week]=2019
  }else{
    mean_df_inter_SYM$year[list_week]=2020
  }
}

#regroup the data per week, year and replicat
mean_df_inter_SYM2 <- mean_df_inter_SYM %>% 
  group_by(week, year, replicat, node_loc) %>%
  summarise(total=n())

#we complete every lacking data for each week, region and replicate
for(a_location in locations){
  print(a_location)
  for(a_week in 1:30){
    sub_table = mean_df_inter_SYM2[mean_df_inter_SYM2$week==a_week,]
    sub_table = sub_table[sub_table$node_loc==a_location,]
    if(nrow(sub_table)==0){
      for(j in 1:n_replicates){
        mean_df_inter_SYM2[nrow(mean_df_inter_SYM2) + 1,] = list(a_week, 2020, j, a_location, 0)
      }
    }else{
      for(i in 1:n_replicates){
        if(i %!in% sub_table$replicat){
          mean_df_inter_SYM2[nrow(mean_df_inter_SYM2) + 1,] = list(a_week, 2020, i, a_location, 0)
        }
      }
    }
  }  
}

#we update the frequency for each week and year
mean_df_inter_SYM3 <- mean_df_inter_SYM2 %>%
  group_by(week, year, node_loc) %>%
  summarise(total2=(sum(total)/n_replicates))

mean_df_inter_SYM3 <- mean_df_inter_SYM3[order(mean_df_inter_SYM3$year, mean_df_inter_SYM3$week),]

mean_df_inter_SYM3$week_2 <- as.Date("0000-01-01")

#we convert the week as a date
for(week2 in 1:nrow(mean_df_inter_SYM3)){
  if(mean_df_inter_SYM3$year[week2]==2019){
    mean_df_inter_SYM3$week_2[week2] = as.Date(lubridate::ymd("2019-01-01") + lubridate::weeks(mean_df_inter_SYM3$week[week2] - 1 ))
  }else{
    mean_df_inter_SYM3$week_2[week2] = as.Date(lubridate::ymd("2020-01-01") + lubridate::weeks(mean_df_inter_SYM3$week[week2] - 1 ))
  }
}

#we order the data for the plot
mean_df_inter_SYM3$node_loc %<>% factor(levels=c("Europe","North.America", "Asia",  "South.America", "Africa", "Oceania"))

#stacked area chart
ggplot(mean_df_inter_SYM3, aes(x=week_2, y=total2)) + 
  geom_area(aes(fill=node_loc), size=.5, colour="black", position="stack")+
  scale_fill_manual(values = c("#538DCA", "#1A7E41", "#5F4197", "#FEDA27", "#EF3437", "#FF8421"))+
  theme(axis.line.x = element_line(size = 1, colour = "black"),
        axis.line.y = element_line(size = 1, colour = "black"),
        axis.line = element_line(size=1, colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        text=element_text(size = 16, family="Arial"),
        axis.text.x=element_text(colour="black", size = 12, angle = 45, vjust=0.6),
        axis.text.y=element_text(colour="black", size = 12),
        legend.position = "none")+
  ylim(0,2)+
  scale_x_date(date_breaks = "1 month", limits = as.Date(c('2020-01-01','2020-07-31')))+
  ylab("Number of exportations")+
  xlab("Date")


##############################
#########LINEAGE STUDY########
##############################

#get all introduction and exportation events for France
introduced_france = df_inter_SYM[df_inter_SYM$node_loc=="France",]
france_export = df_inter_SYM[df_inter_SYM$parent_loc=="France",]

#group the events for each lineage
introduced_france2 = introduced_france %>%
  group_by(parent_clade) %>%
  summarise(total=n())

france_export2 = france_export %>%
  group_by(parent_clade) %>%
  summarise(total=n())

#introduction analysis

introduced_france2$type="all"
#calculate the frequency
introduced_france2$freq = introduced_france2$total/sum(introduced_france2$total)*100
#we only label lineages with a frequency >= 3%
for(i in 1:nrow(introduced_france2)){
  if(introduced_france2$freq[i]<3){
    introduced_france2$final_clade[i]="other"
  }else{
    introduced_france2$final_clade[i] = introduced_france2$parent_clade[i]
  }
}

#regroup "other" variants together
introduced_france3 = introduced_france2 %>%
  group_by(type, final_clade) %>%
  summarise(total2=sum(total))

#recalculate the frequency of variants
introduced_france3$freq = introduced_france3$total2/sum(introduced_france3$total2)*100
#order the variants according to freq
introduced_france3 <-introduced_france3[order(introduced_france3$freq, decreasing = F),]
level_order <- c(introduced_france3$final_clade)
introduced_france3$final_clade <- factor(introduced_france3$final_clade, levels=level_order)

#show the variants introduced to France
introduction = ggplot(data = introduced_france3, aes(x=type, y=freq, fill=final_clade))+
  geom_bar(stat="identity")+
  geom_label(aes(label = final_clade),
             position = position_stack(vjust = 0.5),
             size = 3,
             colour = 'white')+
  scale_fill_manual(values = c("#ff7c43", "#d45087", "#ffa600",  "#2f4b7c", 
                               "#003f5c", "#f95d6a", "#665191", "#a05195"))+
  theme(axis.line.x = element_line(size = 1, colour = "black"),
        axis.line.y = element_line(size = 1, colour = "black"),
        axis.line = element_line(size=1, colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        text=element_text(size = 12, family="Arial"),
        axis.text.x=element_text(colour="black", size = 12, angle=45),
        axis.text.y=element_text(colour="black", size = 12),
        legend.position = "none")


#exportation analysis

france_export2$type="all"
#calculate the frequency
france_export2$freq = france_export2$total/sum(france_export2$total)*100
#we only label lineages with a frequency >= 3%
for(i in 1:nrow(france_export2)){
  if(france_export2$freq[i]<3){
    france_export2$final_clade[i]="other"
  }else{
    france_export2$final_clade[i] = france_export2$parent_clade[i]
  }
}

#regroup "other" variants together
france_export3 = france_export2 %>%
  group_by(type, final_clade) %>%
  summarise(total2=sum(total))

#recalculate the frequency of variants
france_export3$freq = france_export3$total2/sum(france_export3$total2)*100
#order according to freq
france_export3 <-france_export3[order(france_export3$freq, decreasing = F),]
level_order <- c(france_export3$final_clade)
france_export3$final_clade <- factor(france_export3$final_clade, levels=level_order)

#show the variants exported from France
exportation = ggplot(data = france_export3, aes(x=type, y=freq, fill=final_clade))+
  geom_bar(stat="identity")+
  scale_fill_manual(values = c("#003f5c", "#a05195", "#f95d6a", "#ffa600"))+
  geom_label(aes(label = final_clade),
             position = position_stack(vjust = 0.5),
             size = 3,
             colour = 'white')+
  theme(axis.line.x = element_line(size = 1, colour = "black"),
        axis.line.y = element_line(size = 1, colour = "black"),
        axis.line = element_line(size=1, colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        text=element_text(size = 12, family="Arial"),
        axis.text.x=element_text(colour="black", size = 12, angle=45),
        axis.text.y=element_text(colour="black", size = 12),
        legend.position = "none")

#show the plots
plot_grid(exportation, introduction, nrow = 1, ncol = 2)
