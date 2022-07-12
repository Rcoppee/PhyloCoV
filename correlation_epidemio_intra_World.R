#import libraries
library(ggplot2)
library(ggpubr)

# set working directory
setwd("C:/Users/romain.coppee/Documents/PhyloCoV/World/")

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
nb_transmissions = nb_transmissions[nb_transmissions$Date>="2020-04-01",]

#wave2. Remove the two previous lines and use the two followings:
#nb_transmissions = nb_transmissions[nb_transmissions$Date<="2020-12-31",]
#nb_transmissions = nb_transmissions[nb_transmissions$Date>="2020-07-01",]

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
asia = epidemio[epidemio$continent=="Asia",]
sa = epidemio[epidemio$continent=="South America",]
na = epidemio[epidemio$continent=="North America",]
europe = epidemio[epidemio$continent=="Europe",]
europe = europe[europe$location!="France",]
france = epidemio[epidemio$location=="France",]
africa = epidemio[epidemio$continent=="Africa",]

#cutting into waves
#we calculate the number of observation for each date
#either wave 1 (w1) or wave 2 (w2) will be generated
asia_w1 = asia[asia$date<max_w1,]
asia_w2 = asia[asia$date>=max_w1,]
asia_w1_fix = asia_w1 %>%
	group_by(continent, date) %>%
	summarise(nb_deaths=sum(total_deaths))
asia_w2_fix = asia_w2 %>%
	group_by(continent, date) %>%
	summarise(nb_deaths=sum(total_deaths))

sa_w1 = sa[sa$date<max_w1,]
sa_w2 = sa[sa$date>=max_w1,]
sa_w1_fix = sa_w1 %>%
	group_by(continent, date) %>%
	summarise(nb_deaths=sum(total_deaths))
sa_w2_fix = sa_w2 %>%
	group_by(continent, date) %>%
	summarise(nb_deaths=sum(total_deaths))

na_w1 = na[na$date<max_w1,]
na_w2 = na[na$date>=max_w1,]
na_w1_fix = na_w1 %>%
	group_by(continent, date) %>%
	summarise(nb_deaths=sum(total_deaths))
na_w2_fix = na_w2 %>%
	group_by(continent, date) %>%
	summarise(nb_deaths=sum(total_deaths))

europe_w1 = europe[europe$date<max_w1,]
europe_w2 = europe[europe$date>=max_w1,]
europe_w1_fix = europe_w1 %>%
	group_by(continent, date) %>%
	summarise(nb_deaths=sum(total_deaths))
europe_w2_fix = europe_w2 %>%
	group_by(continent, date) %>%
	summarise(nb_deaths=sum(total_deaths))

france_w1 = france[france$date<max_w1,]
france_w2 = france[france$date>=max_w1,]
france_w1_fix = france_w1 %>%
	group_by(continent, date) %>%
	summarise(nb_deaths=sum(total_deaths))
france_w2_fix = france_w2 %>%
	group_by(continent, date) %>%
	summarise(nb_deaths=sum(total_deaths))

africa_w1 = africa[africa$date<max_w1,]
africa_w2 = africa[africa$date>=max_w1,]
africa_w1_fix = africa_w1 %>%
	group_by(continent, date) %>%
	summarise(nb_deaths=sum(total_deaths))
africa_w2_fix = africa_w2 %>%
	group_by(continent, date) %>%
	summarise(nb_deaths=sum(total_deaths))

#wave 1
#combine observations with epidemio
my_asia_w1 = nb_transmissions[nb_transmissions$In=="Asia",]
my_asia_w1$date = as.character(my_asia_w1$date)
combine_asia_w1 = merge(asia_w1_fix,my_asia_w1, by="date") 

my_sa_w1 = nb_transmissions[nb_transmissions$In=="South.America",]
my_sa_w1$date = as.character(my_sa_w1$date)
combine_sa_w1 = merge(sa_w1_fix,my_sa_w1, by="date") 

my_na_w1 = nb_transmissions[nb_transmissions$In=="North.America",]
my_na_w1$date = as.character(my_na_w1$date)
combine_na_w1 = merge(na_w1_fix,my_na_w1, by="date") 

my_europe_w1 = nb_transmissions[nb_transmissions$In=="Europe",]
my_europe_w1$date = as.character(my_europe_w1$date)
combine_europe_w1 = merge(europe_w1_fix,my_europe_w1, by="date") 

my_france_w1 = nb_transmissions[nb_transmissions$In=="France",]
my_france_w1$date = as.character(my_france_w1$date)
combine_france_w1 = merge(france_w1_fix,my_france_w1, by="date") 

my_africa_w1 = nb_transmissions[nb_transmissions$In=="Africa",]
my_africa_w1$date = as.character(my_africa_w1$date)
combine_africa_w1 = merge(africa_w1_fix,my_africa_w1, by="date") 

#plot and regression
plot_asia = ggplot(combine_asia_w1, aes(x=nb_deaths, y=total)) + geom_point(colour="#5F4197") +
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

plot_sa = ggplot(combine_sa_w1, aes(x=nb_deaths, y=total)) + geom_point(colour="#FEDA27") +
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

plot_na = ggplot(combine_na_w1, aes(x=nb_deaths, y=total)) + geom_point(colour="#1A7E41") +
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

plot_europe = ggplot(combine_europe_w1, aes(x=nb_deaths, y=total)) + geom_point(colour="#538DCA") +
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

plot_france = ggplot(combine_france_w1, aes(x=nb_deaths, y=total)) + geom_point(colour="#54C538") +
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

plot_africa = ggplot(combine_africa_w1, aes(x=nb_deaths, y=total)) + geom_point(colour="#EF3437") +
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
plot_grid(plot_asia, plot_na, plot_sa, plot_africa, plot_europe, plot_france, nrow = 2, ncol = 3)
  

#wave 2
#combine observations with epidemio
my_asia_w2 = nb_transmissions[nb_transmissions$In=="Asia",]
my_asia_w2$date = as.character(my_asia_w2$date)
combine_asia_w2 = merge(asia_w2_fix,my_asia_w2, by="date") 

my_sa_w2 = nb_transmissions[nb_transmissions$In=="South.America",]
my_sa_w2$date = as.character(my_sa_w2$date)
combine_sa_w2 = merge(sa_w2_fix,my_sa_w2, by="date") 

my_na_w2 = nb_transmissions[nb_transmissions$In=="North.America",]
my_na_w2$date = as.character(my_na_w2$date)
combine_na_w2 = merge(na_w2_fix,my_na_w2, by="date") 

my_europe_w2 = nb_transmissions[nb_transmissions$In=="Europe",]
my_europe_w2$date = as.character(my_europe_w2$date)
combine_europe_w2 = merge(europe_w2_fix,my_europe_w2, by="date") 

my_france_w2 = nb_transmissions[nb_transmissions$In=="France",]
my_france_w2$date = as.character(my_france_w2$date)
combine_france_w2 = merge(france_w2_fix,my_france_w2, by="date") 

my_africa_w2 = nb_transmissions[nb_transmissions$In=="Africa",]
my_africa_w2$date = as.character(my_africa_w2$date)
combine_africa_w2 = merge(africa_w2_fix,my_africa_w2, by="date") 

#plot and regression
plot_asia = ggplot(combine_asia_w2, aes(x=nb_deaths, y=total)) + geom_point(colour="#5F4197") +
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

plot_sa = ggplot(combine_sa_w2, aes(x=nb_deaths, y=total)) + geom_point(colour="#FEDA27") +
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

plot_na = ggplot(combine_na_w2, aes(x=nb_deaths, y=total)) + geom_point(colour="#1A7E41") +
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

plot_europe = ggplot(combine_europe_w2, aes(x=nb_deaths, y=total)) + geom_point(colour="#538DCA") +
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

plot_france = ggplot(combine_france_w2, aes(x=nb_deaths, y=total)) + geom_point(colour="#54C538") +
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

plot_africa = ggplot(combine_africa_w2, aes(x=nb_deaths, y=total)) + geom_point(colour="#EF3437") +
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
plot_grid(plot_asia, plot_na, plot_sa, plot_africa, plot_europe, plot_france, nrow = 2, ncol = 3)