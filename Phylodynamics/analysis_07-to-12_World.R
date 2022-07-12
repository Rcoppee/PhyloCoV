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

#import metadata
metadata <- readxl::read_xlsx("metadata_file4.xlsx", sheet = "Feuil1") %>%
	mutate(geo_loc = continent_wave)

#import all files (n = replicates) with inter transmissions
temp = list.files(pattern="*inter_NEW.sym")
myfiles = lapply(temp, read.delim)
df_inter_SYM = rbindlist(myfiles, fill=FALSE, idcol=NULL)
df_inter_SYM$Date = as.Date(df_inter_SYM$Date)
df_inter_SYM = df_inter_SYM[df_inter_SYM$Date>="2020-07-21",]

#import all files (n = replicates) with intra transmissions
temp = list.files(pattern="*intra_NEW.sym")
myfiles = lapply(temp, read.delim)
df_intra_SYM = rbindlist(myfiles, fill=FALSE, idcol=NULL)
df_intra_SYM$Date = as.Date(df_intra_SYM$Date)
df_intra_SYM = df_intra_SYM[df_intra_SYM$Date>="2020-07-21",]

#merge all transmissions
df_SYM = rbind(df_inter_SYM, df_intra_SYM)
df_inter_SYM$number = 1

#number of replicates in the study
n_replicat = 100

#list of the regions included in the study
list_continent = c("Asia", "Europe", "France", "North.America", "South.America", "Africa", "Oceania")
   
#We first explore variation in replicates
df_inter_count_SYM = data.frame(Country=character(), In=integer(), Out=integer(), Replicat=integer())

#we count the total on introduction and exportation events per territory
for(i in 1:n_replicat){
	for(j in list_continent){
		sub_df_inter_SYM_IN = df_inter_SYM[df_inter_SYM$In==j,]
		sub_df_inter_SYM_IN = sub_df_inter_SYM_IN[sub_df_inter_SYM_IN$Replicat==i,]
		sub_df_inter_SYM_OUT = df_inter_SYM[df_inter_SYM$Out==j,]
		sub_df_inter_SYM_OUT = sub_df_inter_SYM_OUT[sub_df_inter_SYM_OUT$Replicat==i,]   
		df_inter_count_SYM[nrow(df_inter_count_SYM) + 1,] = list(j, dim(sub_df_inter_SYM_IN)[1], dim(sub_df_inter_SYM_OUT)[1], i)
	}
}

#plot introduction vs exportation transmissions
ggplot(df_inter_count_SYM)+
	geom_point(aes(x=In, y=Out, fill=Country), shape=21, size=3, color="black", alpha=1)+
	scale_fill_manual(values=c("#EF3437", "#5F4197", "#538DCA", "#54C538", "#1A7E41", "#FF8421", "#FEDA27"))+
	xlim(0,50)+
	ylim(0,50)+
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

#calculate the mean across the replicates, and plot it
inter_count_mean = df_inter_count_SYM %>%
  group_by(Country) %>%
  summarize(In_m=mean(In),Out_m=mean(Out))

ggplot(inter_count_mean) +
  geom_point(aes(x=In_m, y=Out_m, fill=Country), shape=21, size=10, color="black", alpha=1)+
  scale_fill_manual(values=c("#EF3437", "#5F4197", "#538DCA", "#54C538", "#1A7E41", "#FF8421", "#FEDA27"))+
  xlim(0,50) +
  ylim(0,50) +
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
	ylim(0,50) +
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
	ylim(0,50) +
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

#count the number of identical events for a same replicat (without dating)
df_inter_SYM2 = df_inter_SYM %>%
	group_by(In, Out, Replicat) %>%
	summarise(a_sum=sum(number), a_mean=(mean(number)))

#convert to dataframe
df_inter_SYM2 = as.data.frame(df_inter_SYM2)

#inialisation dataframe
completion = data.frame(In=integer(), Out=character(), Replicat=integer(), a_sum=integer(), a_mean=integer())
'%!in%'<- function(x,y)!('%in%'(x,y))

#update the transmission table when no transmission between two territories is observed
for (payss in list_continent) {
	print(payss)
	temp_SYM = df_inter_SYM2[df_inter_SYM2$In == payss,]
	for (payss_out in list_continent){
		if(payss==payss_out){
			next
		} else {
		for(i in 1:n_replicat){
			if (i %!in% temp_SYM[temp_SYM[,2] == payss_out,]$Replicat){
				completion[nrow(completion) + 1,] = list(payss, payss_out, i,0,0)
			}
		}
		}
	}
}
df_inter_SYM2 = rbind (df_inter_SYM2, completion)

#visualizing exportation events from France
France_SYM = df_inter_SYM2[df_inter_SYM2$In == "France",]
FR_SYM_exp = ggplot(France_SYM, aes(x=Out, y=a_sum, fill=Out)) + 
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
        legend.position = "none")

#view plot
FR_SYM_exp

#visualizing introduction events to France
France_SYM = df_inter_SYM2[df_inter_SYM2$Out == "France",]
FR_SYM_intro = ggplot(France_SYM, aes(x=In, y=a_sum, fill=In)) + 
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
        legend.position = "none")

#view plot
FR_SYM_intro

#calculate frequency per replicat, in and out (then for all transmissions)
flux_m_SYM = rename(count(df_inter_SYM, In, Out, Replicat), Freq = n)
flux_m_SYM_all = rename(count(df_SYM, In, Out, Replicat), Freq = n)

#create a dataframe for circos plot
circos_SYM = flux_m_SYM %>% group_by(In, Out) %>% summarize(B = sum(Freq))
circos_SYM = as.data.frame(circos_SYM)

#remove any incomplete rows
circos_SYM = circos_SYM[complete.cases(circos_SYM), ]

#remove intra_transmission if present
for (row in 1:nrow(circos_SYM)){
	if (circos_SYM$In[row]==circos_SYM$Out[row]){
		circos_SYM = circos_SYM[-row,]
	}
}

#calculate the percentage of each flow. To not use here
#circos_SYM$perc=circos_SYM$B/sum(circos_SYM$B)*100
#write.table(circos_SYM, "for_map.txt", sep="\t")

#if incluse of intra-territory transmissions, update with these lines
#circos_SYM = flux_m_SYM_all %>% group_by(In, Out) %>% summarize(B = sum(Freq))
#circos_SYM = as.data.frame(circos_SYM)
#circos_SYM = circos_SYM[complete.cases(circos_SYM), ]

#parameters of circos plot
circos.clear()
circos.par(start.degree = 90, gap.degree = 1, track.margin = c(-0.1, 0.1), points.overflow.warning = FALSE)
par(mar = rep(0, 4))
# color palette
mycolor <- c("#EF3437", "#5F4197", "#538DCA", "#54C538", "#1A7E41", "#FF8421", "#FEDA27")

#Base plot
chordDiagram(
	x = circos_SYM,
	grid.col = mycolor,
	transparency = 0.25,
	directional = 1,
	direction.type = c("arrows", "diffHeight"), 
	diffHeight  = -0.04,
	annotationTrack = "grid", 
	annotationTrackHeight = c(0.1, 1), #HAUTEUR DES ABSCISSES 0.01
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


#plot cumulative exportations
P_exportation_SYM = ggplot(df_inter_SYM,aes(x=Date,color=In)) +
	stat_bin(data=subset(df_inter_SYM,In=="France"),aes(y=cumsum(..count..)),geom="step", size=1,bins = 100)+
	stat_bin(data=subset(df_inter_SYM,In=="Europe"),aes(y=cumsum(..count..)),geom="step", size=1,bins = 100)+
	stat_bin(data=subset(df_inter_SYM,In=="North.America"),aes(y=cumsum(..count..)),geom="step", size=1,bins = 100)+
	stat_bin(data=subset(df_inter_SYM,In=="South.America"),aes(y=cumsum(..count..)),geom="step", size=1,bins = 100)+
	stat_bin(data=subset(df_inter_SYM,In=="Africa"),aes(y=cumsum(..count..)),geom="step", size=1,bins = 100)+
	stat_bin(data=subset(df_inter_SYM,In=="Asia"),aes(y=cumsum(..count..)),geom="step", size=1,bins = 100)+
	ylim(0,4000)+
	scale_color_manual(values = c("#EF3437", "#5F4197", "#538DCA", "#54C538", "#1A7E41", "#FEDA27", "#FF8421"))+
	scale_x_date(date_breaks = "1 month", limits = as.Date(c('2020-07-21','2020-12-31')))+
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
        legend.position = "none")

#view plot
P_exportation_SYM

#plotting the dynamics of exportation events
#P_exportation_SYM = ggplot(df_inter_SYM,aes(x=Date,color=In)) +
#	stat_bin(data=subset(df_inter_SYM,In=="France"),geom="step", size=1,bins = 100,alpha = 1)+
#	stat_bin(data=subset(df_inter_SYM,In=="Europe"),geom="step", size=1,bins = 100,alpha = 1)+
#	stat_bin(data=subset(df_inter_SYM,In=="North.America"),geom="step", size=1,bins = 100,alpha = 1)+
#	stat_bin(data=subset(df_inter_SYM,In=="South.America"),geom="step", size=1,bins = 100,alpha = 1)+
#	stat_bin(data=subset(df_inter_SYM,In=="Africa"),geom="step", size=1,bins = 100,alpha = 1)+
#	stat_bin(data=subset(df_inter_SYM,In=="Asia"),geom="step", size=1,bins = 100,alpha = 1)+
#	ylim(0,150)+
#	scale_x_date(date_breaks = "1 month", limits = as.Date(c('2020-07-21','2020-12-31')))+
#	scale_color_manual(values = c("#EF3437", "#5F4197", "#538DCA", "#54C538", "#1A7E41", "#FEDA27", "#FF8421"))+
#	theme(axis.line.x = element_line(size = 1, colour = "black"),
#		axis.line.y = element_line(size = 1, colour = "black"),
#       axis.line = element_line(size=1, colour = "black"),
#       panel.grid.major = element_blank(),
#       panel.grid.minor = element_blank(),
#       panel.border = element_blank(),
#       panel.background = element_blank(),
#       text=element_text(size = 16, family="Arial"),
#       axis.text.x=element_text(colour="black", size = 12, angle = 45, vjust=0.6),
#       axis.text.y=element_text(colour="black", size = 12),
#       legend.position = "none")

#plot cumulative introductions
P_importation_SYM = ggplot(df_inter_SYM,aes(x=Date,color=Out)) +
	stat_bin(data=subset(df_inter_SYM ,Out=="France"),aes(y=cumsum(..count..)),geom="step", size=1,bins = 100)+
	stat_bin(data=subset(df_inter_SYM ,Out=="Europe"),aes(y=cumsum(..count..)),geom="step", size=1,bins = 100)+
	stat_bin(data=subset(df_inter_SYM ,Out=="North.America"),aes(y=cumsum(..count..)),geom="step", size=1,bins = 100)+
	stat_bin(data=subset(df_inter_SYM ,Out=="South.America"),aes(y=cumsum(..count..)),geom="step", size=1,bins = 100)+
	stat_bin(data=subset(df_inter_SYM ,Out=="Africa"),aes(y=cumsum(..count..)),geom="step", size=1,bins = 100)+
	stat_bin(data=subset(df_inter_SYM ,Out=="Asia"),aes(y=cumsum(..count..)),geom="step", size=1,bins = 100)+
	ylim(0,4000)+
	scale_color_manual(values = c("#EF3437", "#5F4197", "#538DCA", "#54C538", "#1A7E41", "#FEDA27", "#FF8421"))+
	scale_x_date(date_breaks = "1 month", limits = as.Date(c('2020-07-21','2020-12-31')))+
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
        legend.position = "none")

#view plot
P_importation_SYM

#plotting the dynamics of introduction events
#P_importation_SYM = ggplot(df_inter_SYM,aes(x=Date,color=Out)) +
#	stat_bin(data=subset(df_inter_SYM ,Out=="France"),geom="step", size=1,bins = 100,alpha = 1)+
#   stat_bin(data=subset(df_inter_SYM ,Out=="Europe"),geom="step", size=1,bins = 100,alpha = 1)+
#   stat_bin(data=subset(df_inter_SYM ,Out=="North.America"),geom="step", size=1,bins = 100,alpha = 1)+
#   stat_bin(data=subset(df_inter_SYM ,Out=="South.America"),geom="step", size=1,bins = 100,alpha = 1)+
#   stat_bin(data=subset(df_inter_SYM ,Out=="Africa"),geom="step", size=1,bins = 100,alpha = 1)+
#   stat_bin(data=subset(df_inter_SYM ,Out=="Asia"),geom="step", size=1,bins = 100,alpha = 1)+
#   ylim(0,100)+
#   scale_x_date(date_breaks = "1 month", limits = as.Date(c('2020-07-21','2020-12-31')))+
#	scale_color_manual(values = c("#EF3437", "#5F4197", "#538DCA", "#54C538", "#1A7E41", "#FEDA27", "#FF8421"))+
#	theme(axis.line.x = element_line(size = 1, colour = "black"),
#		axis.line.y = element_line(size = 1, colour = "black"),
#       axis.line = element_line(size=1, colour = "black"),
#       panel.grid.major = element_blank(),
#       panel.grid.minor = element_blank(),
#       panel.border = element_blank(),
#       panel.background = element_blank(),
#       text=element_text(size = 16, family="Arial"),
#       axis.text.x=element_text(colour="black", size = 12, angle = 45, vjust=0.6),
#       axis.text.y=element_text(colour="black", size = 12),
#       legend.position = "none")

#plot cumulative intra-territory transmissions
P_intra_SYM = ggplot(df_intra_SYM,aes(x=Date,color=In)) +
	stat_bin(data=subset(df_intra_SYM ,In=="France"),aes(y=cumsum(..count..)),geom="step", size=1,bins = 100)+
	stat_bin(data=subset(df_intra_SYM ,In=="Europe"),aes(y=cumsum(..count..)),geom="step", size=1,bins = 100)+
	stat_bin(data=subset(df_intra_SYM ,In=="North.America"),aes(y=cumsum(..count..)),geom="step", size=1,bins = 100)+
	stat_bin(data=subset(df_intra_SYM ,In=="South.America"),aes(y=cumsum(..count..)),geom="step", size=1,bins = 100)+
	stat_bin(data=subset(df_intra_SYM ,In=="Africa"),aes(y=cumsum(..count..)),geom="step", size=1,bins = 100)+
	stat_bin(data=subset(df_intra_SYM ,In=="Asia"),aes(y=cumsum(..count..)),geom="step", size=1,bins = 100)+
	ylim(0,25000)+
	scale_color_manual(values = c("#EF3437", "#5F4197", "#538DCA", "#54C538", "#1A7E41", "#FEDA27", "#FF8421"))+
	scale_x_date(date_breaks = "1 month", limits = as.Date(c('2020-07-21','2020-12-31')))+
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
		legend.position = "none")
 
#view plot
P_intra_SYM
 
#plotting the dynamics of intra-territory transmission events
#P_intra_SYM = ggplot(df_intra_SYM,aes(x=Date,color=In)) +
#	stat_bin(data=subset(df_intra_SYM ,In=="France"),geom="step", size=1,bins = 100)+
#	stat_bin(data=subset(df_intra_SYM ,In=="Europe"),geom="step", size=1,bins = 100)+
#	stat_bin(data=subset(df_intra_SYM ,In=="North.America"),geom="step", size=1,bins = 100)+
#	stat_bin(data=subset(df_intra_SYM ,In=="South.America"),geom="step", size=1,bins = 100)+
#	stat_bin(data=subset(df_intra_SYM ,In=="Africa"),geom="step", size=1,bins = 100)+
#	stat_bin(data=subset(df_intra_SYM ,In=="Asia"),geom="step", size=1,bins = 100)+
#	ylim(0,1500)+
#	scale_x_date(date_breaks = "1 month", limits = as.Date(c('2020-01-01','2020-07-20')))+
#	scale_color_manual(values = c("#EF3437", "#5F4197", "#538DCA", "#54C538", "#1A7E41", "#FEDA27", "#FF8421"))+
#	theme(axis.line.x = element_line(size = 1, colour = "black"),
#		axis.line.y = element_line(size = 1, colour = "black"),
#       axis.line = element_line(size=1, colour = "black"),
#       panel.grid.major = element_blank(),
#       panel.grid.minor = element_blank(),
#       panel.border = element_blank(),
#       panel.background = element_blank(),
#       text=element_text(size = 16, family="Arial"),
#       axis.text.x=element_text(colour="black", size = 12, angle = 45, vjust=0.6),
#       axis.text.y=element_text(colour="black", size = 12),
#       legend.position = "none")

#to smooth, update the table using this option 
require(zoo)
myData_smooth <- df_intra_SYM %>%
	mutate(mean7_num_cases = rollmean(Date, k=7, fill=NA, align="right"))

myData_smooth = myData_smooth[complete.cases(myData_smooth), ]

#focusing on France transmissions
#note that we can smooth again if needed using the same strategy as previously
df_from_FR_to = df_inter_SYM[df_inter_SYM$In=="France",]
df_to_FR = df_inter_SYM[df_inter_SYM$Out=="France",]

#plot exportation events from France
plot_from_France = ggplot(df_from_FR_to,aes(x=Date,color=Out)) +
	stat_bin(data=subset(df_from_FR_to ,Out=="Europe"),geom="step", size=1,bins = 100,)+
	stat_bin(data=subset(df_from_FR_to ,Out=="North.America"),geom="step", size=1,bins = 100)+
	stat_bin(data=subset(df_from_FR_to ,Out=="South.America"),geom="step", size=1,bins = 100)+
	stat_bin(data=subset(df_from_FR_to ,Out=="Africa"),geom="step", size=1,bins = 100)+
	stat_bin(data=subset(df_from_FR_to ,Out=="Asia"),geom="step", size=1,bins = 100)+
	ylim(0,60)+
	scale_color_manual(values = c("#EF3437", "#5F4197", "#538DCA", "#1A7E41", "#FEDA27", "#FF8421"))+
	scale_x_date(date_breaks = "1 month", limits = as.Date(c('2020-07-21','2020-12-31')))+
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
		labs(x="Date", y="Exportation across replicates")

#plot introduction events to France
plot_to_France = ggplot(df_to_FR,aes(x=Date,color=In)) +
	stat_bin(data=subset(df_to_FR ,In=="Europe"),geom="step", size=1,bins = 100)+
	stat_bin(data=subset(df_to_FR ,In=="North.America"),geom="step", size=1,bins = 100)+
	stat_bin(data=subset(df_to_FR ,In=="South.America"),geom="step", size=1,bins = 100)+
	stat_bin(data=subset(df_to_FR ,In=="Africa"),geom="step", size=1,bins = 100)+
	stat_bin(data=subset(df_to_FR ,In=="Asia"),geom="step", size=1,bins = 100)+
	ylim(0,50)+
	scale_color_manual(values = c("#EF3437", "#5F4197", "#538DCA", "#1A7E41", "#FEDA27", "#FF8421"))+
	scale_x_date(date_breaks = "1 month", limits = as.Date(c('2020-07-21','2020-12-31')))+
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
	labs(x="Date", y="Importation across replicates")

#plotting the results
plot_grid(plot_from_France, plot_to_France, ncol=2, nrow=1)