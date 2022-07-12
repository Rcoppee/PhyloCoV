#import libraries
library(ggplot2)
library(ggpubr)

# set working directory
setwd("C:/Users/romain.coppee/Documents/PhyloCoV/France/")

#define the max week of the first wave
max_w1=31

#import all files (n = replicates) with intra transmissions
temp = list.files(pattern="*intra_NEW.sym")
myfiles = lapply(temp, read.delim)
df_intra_SYM = rbindlist(myfiles, fill=FALSE, idcol=NULL)
df_intra_SYM$Date = as.Date(df_intra_SYM$Date)

#count the total of transmissions per region
df_intra_SYM$number = 1
nb_transmissions = df_intra_SYM %>%
	group_by(In, Out, Date) %>%
	summarise(a_sum=n())

#wave1
nb_transmissions = nb_transmissions[nb_transmissions$Date<="2020-07-01",]
nb_transmissions = nb_transmissions[nb_transmissions$Date>="2020-03-01",]

#remove possible incomplete data
nb_transmissions = nb_transmissions[complete.cases(nb_transmissions),]
#convert date to week
nb_transmissions$week = week(nb_transmissions$Date)

nb_transmissions_weekly = nb_transmissions %>%
	group_by(In, Out, week) %>%
	summarize(nb_deaths=sum(a_sum))

nb_transmissions_weekly$total = 0
for(i in 1:nrow(nb_transmissions_weekly)){
	if(i == 1 || nb_transmissions_weekly$In[i] != nb_transmissions_weekly$In[i-1]){
		nb_transmissions_weekly$total[i] = nb_transmissions_weekly$nb_deaths[i]
		next
	}
	if(nb_transmissions_weekly$In[i] == nb_transmissions_weekly$In[i-1]){
		nb_transmissions_weekly$total[i] = nb_transmissions_weekly$total[i-1] + nb_transmissions_weekly$nb_deaths[i]
	}
}

#importation of epidemiological data
epidemio = read.table(file = "epidemio_france.txt", sep="\t", header=T, quote="\"")
epidemio = epidemio[epidemio$week<="30",]
epidemio[is.na(epidemio)] <- 0

epidemio$total = 0
#calculate the cumul per region
for(i in 1:nrow(epidemio)){
	if(i == 1 || epidemio$region[i] != epidemio$region[i-1]){
		epidemio$total[i] = epidemio$nb_deaths[i]
		next
	}
	if(epidemio$region[i] == epidemio$region[i-1]){
		epidemio$total[i] = epidemio$total[i-1] + epidemio$nb_deaths[i]
	}
}

#obtain subdata of regions
idf = epidemio[epidemio$region=="IDF",]
ara = epidemio[epidemio$region=="ARA",]
paca = epidemio[epidemio$region=="PACA",]

#restrict to wave 1
ara_w1 = ara[ara$week<max_w1,]
idf_w1 = idf[idf$week<max_w1,]
paca_w1 = paca[paca$week<max_w1,]

#combine observations with epidemio
my_ara_w1 = nb_transmissions_weekly[nb_transmissions_weekly$In=="ARA",]
combine_ara_w1 = merge(ara_w1,my_ara_w1, by="week") 

my_idf_w1 = nb_transmissions_weekly[nb_transmissions_weekly$In=="IDF",]
combine_idf_w1 = merge(idf_w1,my_idf_w1, by="week") 

my_paca_w1 = nb_transmissions_weekly[nb_transmissions_weekly$In=="PACA",]
combine_paca_w1 = merge(paca_w1,my_paca_w1, by="week") 

#plot and regression
plot_ara = ggplot(combine_ara_w1, aes(x=total.x, y=total.y)) + geom_point(colour="#FF4D31") +
	geom_smooth(method=lm, se=FALSE, color="black")+
	stat_cor(label.y.npc = 1, digits=10)+
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

plot_idf = ggplot(combine_idf_w1, aes(x=total.x, y=total.y)) + geom_point(colour="#2EC245") +
	geom_smooth(method=lm, se=FALSE, color="black")+
	stat_cor(label.y.npc = 1, digits=10)+
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

plot_paca = ggplot(combine_paca_w1, aes(x=total.x, y=total.y)) + geom_point(colour="#1790F5") +
	geom_smooth(method=lm, se=FALSE, color="black")+
	stat_cor(label.y.npc = 1, digits=10)+
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

#plot results
plot_grid(plot_ara, plot_idf, plot_paca, nrow = 1, ncol = 3)



#wave2
#define the max week of the first wave
max_w1=31

#import all files (n = replicates) with intra transmissions
temp = list.files(pattern="*intra_NEW.sym")
myfiles = lapply(temp, read.delim)
df_intra_SYM = rbindlist(myfiles, fill=FALSE, idcol=NULL)
df_intra_SYM$Date = as.Date(df_intra_SYM$Date)


#count the total of transmissions per region
df_intra_SYM$number = 1
nb_transmissions = df_intra_SYM %>%
	group_by(In, Out, Date) %>%
	summarise(a_sum=n())

#wave1
nb_transmissions = nb_transmissions[nb_transmissions$Date<="2020-12-31",]
nb_transmissions = nb_transmissions[nb_transmissions$Date>="2020-07-21",]

#remove possible incomplete data
nb_transmissions = nb_transmissions[complete.cases(nb_transmissions),]
#convert date to week
nb_transmissions$week = week(nb_transmissions$Date)

nb_transmissions_weekly = nb_transmissions %>%
	group_by(In, Out, week) %>%
	summarize(nb_deaths=sum(a_sum))

nb_transmissions_weekly$total = 0
for(i in 1:nrow(nb_transmissions_weekly)){
	if(i == 1 || nb_transmissions_weekly$In[i] != nb_transmissions_weekly$In[i-1]){
		nb_transmissions_weekly$total[i] = nb_transmissions_weekly$nb_deaths[i]
		next
	}
	if(nb_transmissions_weekly$In[i] == nb_transmissions_weekly$In[i-1]){
		nb_transmissions_weekly$total[i] = nb_transmissions_weekly$total[i-1] + nb_transmissions_weekly$nb_deaths[i]
	}
}

#importation of epidemiological data
epidemio = read.table(file = "epidemio_france.txt", sep="\t", header=T, quote="\"")
epidemio = epidemio[epidemio$week>="31",]
epidemio[is.na(epidemio)] <- 0

epidemio$total = 0
#calculate the cumul per region
for(i in 1:nrow(epidemio)){
	if(i == 1 || epidemio$region[i] != epidemio$region[i-1]){
		epidemio$total[i] = epidemio$nb_deaths[i]
		next
	}
	if(epidemio$region[i] == epidemio$region[i-1]){
		epidemio$total[i] = epidemio$total[i-1] + epidemio$nb_deaths[i]
	}
}

#obtain subdata of regions
idf = epidemio[epidemio$region=="IDF",]
ara = epidemio[epidemio$region=="ARA",]
paca = epidemio[epidemio$region=="PACA",]

#restrict to wave 2
ara_w2 = ara[ara$week>=max_w1,]
idf_w2 = idf[idf$week>=max_w1,]
paca_w2 = paca[paca$week>=max_w1,]

#combine observations with epidemio
my_ara_w2 = nb_transmissions_weekly[nb_transmissions_weekly$In=="ARA",]
combine_ara_w2 = merge(ara_w2,my_ara_w2, by="week") 

my_idf_w2 = nb_transmissions_weekly[nb_transmissions_weekly$In=="IDF",]
combine_idf_w2 = merge(idf_w2,my_idf_w2, by="week") 

my_paca_w2 = nb_transmissions_weekly[nb_transmissions_weekly$In=="PACA",]
combine_paca_w2 = merge(paca_w2,my_paca_w2, by="week") 

#plot and regression
plot_ara = ggplot(combine_ara_w2, aes(x=total.x, y=total.y)) + geom_point(colour="#FF4D31") +
	geom_smooth(method=lm, se=FALSE, color="black")+
	stat_cor(label.y.npc = 1, digits=10)+
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

plot_idf = ggplot(combine_idf_w2, aes(x=total.x, y=total.y)) + geom_point(colour="#2EC245") +
	geom_smooth(method=lm, se=FALSE, color="black")+
	stat_cor(label.y.npc = 1, digits=10)+
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

plot_paca = ggplot(combine_paca_w2, aes(x=total.x, y=total.y)) + geom_point(colour="#1790F5") +
	geom_smooth(method=lm, se=FALSE, color="black")+
	stat_cor(label.y.npc = 1, digits=10)+
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

#plot results
plot_grid(plot_ara, plot_idf, plot_paca, nrow = 1, ncol = 3)