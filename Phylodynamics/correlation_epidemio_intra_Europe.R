#import libraries
library(ggplot2)
library(ggpubr)
library(cowplot)

# set working directory
setwd("C:/Users/romain.coppee/Documents/PhyloCoV/Europe/")

#define the max of the first wave
max_w1="2020-07-01"

#import all files (n = replicates) with intra transmissions
temp = list.files(pattern="*intra_NEW.sym")
myfiles = lapply(temp, read.delim)
df_intra_SYM = rbindlist(myfiles, fill=FALSE, idcol=NULL)
df_intra_SYM$Date = as.Date(df_intra_SYM$Date)

#count the total of transmissions per region
df_intra_SYM$number = 1
nb_transmissions = df_intra_SYM %>%
	group_by(In, Out, Date) %>%
	summarise(a_sum=sum(number))

#wave1
#we fix the start to april to be sure that all region have data
nb_transmissions = nb_transmissions[nb_transmissions$Date<="2020-07-01",]
nb_transmissions = nb_transmissions[nb_transmissions$Date>="2020-04-05",]

#wave2. Remove the two previous lines and use the two followings:
#nb_transmissions = nb_transmissions[nb_transmissions$Date<="2020-12-31",]
#nb_transmissions = nb_transmissions[nb_transmissions$Date>="2020-07-01",]

nb_transmissions = nb_transmissions[!is.na(nb_transmissions$In),]

nb_transmissions$total = 0
#we calculate the cumulative of transmissions
for(i in 1:nrow(nb_transmissions)){
	if(i == 1 || nb_transmissions$In[i] != nb_transmissions$In[i-1]){
		nb_transmissions$total[i] = nb_transmissions$a_sum[i]
		next
	}
	if(nb_transmissions$In[i] == nb_transmissions$In[i-1]){
		nb_transmissions$total[i] = nb_transmissions$total[i-1] + nb_transmissions$a_sum[i]
	}
}

names(nb_transmissions)[names(nb_transmissions) == "Date"] <- "date"

#importation of epidemiological data
epidemio = read.table(file = "death_epidemio_data.txt", sep="\t", header=T, quote="\"")
epidemio = epidemio[epidemio$date<="2020-12-31",]
epidemio[is.na(epidemio)] <- 0

#obtain subdata of continents
belgium = epidemio[epidemio$location=="Belgium",]
france = epidemio[epidemio$location=="France",]
germany = epidemio[epidemio$location=="Germany",]
italy = epidemio[epidemio$location=="Italy",]
netherlands = epidemio[epidemio$location=="Netherlands",]
russia = epidemio[epidemio$location=="Russia",]
spain = epidemio[epidemio$location=="Spain",]
uk = epidemio[epidemio$location=="United Kingdom",]

#cutting into waves
#we calculate the number of observation for each date
#either wave 1 (w1) or wave 2 (w2) will be generated
belgium_w1 = belgium[belgium$date<max_w1,]
belgium_w2 = belgium[belgium$date>=max_w1,]
belgium_w1_fix = belgium_w1 %>%
	group_by(location, date) %>%
	summarise(nb_deaths=sum(total_deaths))
belgium_w2_fix = belgium_w2 %>%
	group_by(location, date) %>%
	summarise(nb_deaths=sum(total_deaths))

france_w1 = france[france$date<max_w1,]
france_w2 = france[france$date>=max_w1,]
france_w1_fix = france_w1 %>%
	group_by(location, date) %>%
	summarise(nb_deaths=sum(total_deaths))
france_w2_fix = france_w2 %>%
	group_by(location, date) %>%
	summarise(nb_deaths=sum(total_deaths))

germany_w1 = germany[germany$date<max_w1,]
germany_w2 = germany[germany$date>=max_w1,]
germany_w1_fix = germany_w1 %>%
	group_by(location, date) %>%
	summarise(nb_deaths=sum(total_deaths))
germany_w2_fix = germany_w2 %>%
	group_by(location, date) %>%
	summarise(nb_deaths=sum(total_deaths))

netherlands_w1 = netherlands[netherlands$date<max_w1,]
netherlands_w2 = netherlands[netherlands$date>=max_w1,]
netherlands_w1_fix = netherlands_w1 %>%
	group_by(location, date) %>%
	summarise(nb_deaths=sum(total_deaths))
netherlands_w2_fix = netherlands_w2 %>%
	group_by(location, date) %>%
	summarise(nb_deaths=sum(total_deaths))

italy_w1 = italy[italy$date<max_w1,]
italy_w2 = italy[italy$date>=max_w1,]
italy_w1_fix = italy_w1 %>%
	group_by(location, date) %>%
	summarise(nb_deaths=sum(total_deaths))
italy_w2_fix = italy_w2 %>%
	group_by(location, date) %>%
	summarise(nb_deaths=sum(total_deaths))

russia_w1 = russia[russia$date<max_w1,]
russia_w2 = russia[russia$date>=max_w1,]
russia_w1_fix = russia_w1 %>%
	group_by(location, date) %>%
	summarise(nb_deaths=sum(total_deaths))
russia_w2_fix = russia_w2 %>%
	group_by(location, date) %>%
	summarise(nb_deaths=sum(total_deaths))

spain_w1 = spain[spain$date<max_w1,]
spain_w2 = spain[spain$date>=max_w1,]
spain_w1_fix = spain_w1 %>%
	group_by(location, date) %>%
	summarise(nb_deaths=sum(total_deaths))
spain_w2_fix = spain_w2 %>%
	group_by(location, date) %>%
	summarise(nb_deaths=sum(total_deaths))

uk_w1 = uk[uk$date<max_w1,]
uk_w2 = uk[uk$date>=max_w1,]
uk_w1_fix = uk_w1 %>%
	group_by(location, date) %>%
	summarise(nb_deaths=sum(total_deaths))
uk_w2_fix = uk_w2 %>%
	group_by(location, date) %>%
	summarise(nb_deaths=sum(total_deaths))
	
#wave 1
#combine observations with epidemio
my_belgium_w1 = nb_transmissions[nb_transmissions$In=="Belgium",]
my_belgium_w1$date = as.character(my_belgium_w1$date)
combine_belgium_w1 = merge(belgium_w1_fix,my_belgium_w1, by="date") 

my_france_w1 = nb_transmissions[nb_transmissions$In=="France",]
my_france_w1$date = as.character(my_france_w1$date)
combine_france_w1 = merge(france_w1_fix,my_france_w1, by="date") 

my_germany_w1 = nb_transmissions[nb_transmissions$In=="Germany",]
my_germany_w1$date = as.character(my_germany_w1$date)
combine_germany_w1 = merge(germany_w1_fix,my_germany_w1, by="date") 

my_italy_w1 = nb_transmissions[nb_transmissions$In=="Italy",]
my_italy_w1$date = as.character(my_italy_w1$date)
combine_italy_w1 = merge(italy_w1_fix,my_italy_w1, by="date") 

my_netherlands_w1 = nb_transmissions[nb_transmissions$In=="Netherlands",]
my_netherlands_w1$date = as.character(my_netherlands_w1$date)
combine_netherlands_w1 = merge(netherlands_w1_fix,my_netherlands_w1, by="date") 

my_russia_w1 = nb_transmissions[nb_transmissions$In=="Russia",]
my_russia_w1$date = as.character(my_russia_w1$date)
combine_russia_w1 = merge(russia_w1_fix,my_russia_w1, by="date") 

my_spain_w1 = nb_transmissions[nb_transmissions$In=="Spain",]
my_spain_w1$date = as.character(my_spain_w1$date)
combine_spain_w1 = merge(spain_w1_fix,my_spain_w1, by="date") 

my_uk_w1 = nb_transmissions[nb_transmissions$In=="United.Kingdom",]
my_uk_w1$date = as.character(my_uk_w1$date)
combine_uk_w1 = merge(uk_w1_fix,my_uk_w1, by="date") 

#plot and regression
plot_belgium = ggplot(combine_belgium_w1, aes(x=nb_deaths, y=total)) +
	geom_point(colour="#FDD93B") +
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

plot_france = ggplot(combine_france_w1, aes(x=nb_deaths, y=total)) +
	geom_point(colour="#54C538") +
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

plot_germany = ggplot(combine_germany_w1, aes(x=nb_deaths, y=total)) +
	geom_point(colour="#FD832F") +
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

plot_italy = ggplot(combine_italy_w1, aes(x=nb_deaths, y=total)) +
	geom_point(colour="#558EC8") +
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

plot_netherlands = ggplot(combine_netherlands_w1, aes(x=nb_deaths, y=total)) +
	geom_point(colour="#872BBA") +
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

plot_russia = ggplot(combine_russia_w1, aes(x=nb_deaths, y=total)) +
	geom_point(colour="#16999F") +
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

plot_spain = ggplot(combine_spain_w1, aes(x=nb_deaths, y=total)) +
	geom_point(colour="#10791E") +
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

plot_uk = ggplot(combine_uk_w1, aes(x=nb_deaths, y=total)) +
	geom_point(colour="#ED313C") +
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
plot_grid(plot_belgium, plot_france, plot_germany, plot_italy,
          plot_netherlands, plot_russia, plot_spain, plot_uk, nrow = 2, ncol = 4)
  

#wave 2
#combine observations with epidemio
my_belgium_w2 = nb_transmissions[nb_transmissions$In=="Belgium",]
my_belgium_w2$date = as.character(my_belgium_w2$date)
combine_belgium_w2 = merge(belgium_w2_fix,my_belgium_w2, by="date") 

my_france_w2 = nb_transmissions[nb_transmissions$In=="France",]
my_france_w2$date = as.character(my_france_w2$date)
combine_france_w2 = merge(france_w2_fix,my_france_w2, by="date") 

my_germany_w2 = nb_transmissions[nb_transmissions$In=="Germany",]
my_germany_w2$date = as.character(my_germany_w2$date)
combine_germany_w2 = merge(germany_w2_fix,my_germany_w2, by="date") 

my_italy_w2 = nb_transmissions[nb_transmissions$In=="Italy",]
my_italy_w2$date = as.character(my_italy_w2$date)
combine_italy_w2 = merge(italy_w2_fix,my_italy_w2, by="date") 

my_netherlands_w2 = nb_transmissions[nb_transmissions$In=="Netherlands",]
my_netherlands_w2$date = as.character(my_netherlands_w2$date)
combine_netherlands_w2 = merge(netherlands_w2_fix,my_netherlands_w2, by="date") 

my_russia_w2 = nb_transmissions[nb_transmissions$In=="Russia",]
my_russia_w2$date = as.character(my_russia_w2$date)
combine_russia_w2 = merge(russia_w2_fix,my_russia_w2, by="date") 

my_spain_w2 = nb_transmissions[nb_transmissions$In=="Spain",]
my_spain_w2$date = as.character(my_spain_w2$date)
combine_spain_w2 = merge(spain_w2_fix,my_spain_w2, by="date") 

my_uk_w2 = nb_transmissions[nb_transmissions$In=="United.Kingdom",]
my_uk_w2$date = as.character(my_uk_w2$date)
combine_uk_w2 = merge(uk_w2_fix,my_uk_w2, by="date") 


#plot and regression
plot_belgium = ggplot(combine_belgium_w2, aes(x=nb_deaths, y=total)) +
	geom_point(colour="#FDD93B") +
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

plot_france = ggplot(combine_france_w2, aes(x=nb_deaths, y=total)) +
	geom_point(colour="#54C538") +
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

plot_germany = ggplot(combine_germany_w2, aes(x=nb_deaths, y=total)) +
	geom_point(colour="#FD832F") +
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

plot_italy = ggplot(combine_italy_w2, aes(x=nb_deaths, y=total)) +
	geom_point(colour="#558EC8") +
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

plot_netherlands = ggplot(combine_netherlands_w2, aes(x=nb_deaths, y=total)) +
	geom_point(colour="#872BBA") +
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

plot_russia = ggplot(combine_russia_w2, aes(x=nb_deaths, y=total)) +
	geom_point(colour="#16999F") +
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

plot_spain = ggplot(combine_spain_w2, aes(x=nb_deaths, y=total)) +
	geom_point(colour="#10791E") +
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

plot_uk = ggplot(combine_uk_w2, aes(x=nb_deaths, y=total)) +
	geom_point(colour="#ED313C") +
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
plot_grid(plot_belgium, plot_france, plot_germany, plot_italy,
          plot_netherlands, plot_russia, plot_spain, plot_uk, nrow = 2, ncol = 4)