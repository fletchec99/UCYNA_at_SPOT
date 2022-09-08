#Look at what environmental factors are different on the dates UCYN-A is present vs. absent
#Ver. 12.07.2021

setwd("/Users/cfletcherh1/Desktop/DataAnalysis/SPOT_16S_18S_eLSA/")

#1. Load data of environmental variables when UCYN-A and host are present/ absent-----
#Read in taxonomy to set ix's 
proks_tax=read.table("ModifiedFiles/Proks_tax_classified_23072020_SILVA132.tsv", header=T, stringsAsFactors = F, sep="\t")
dim(proks_tax)
UCYNA_ix=proks_tax$Feature.ID[grep("UCYN", proks_tax$Taxon)]

euks_tax=read.table("ModifiedFiles/Euks_tax_classified_19052020_SILVA132.tsv", header=T, stringsAsFactors = F, sep="\t")
dim(euks_tax) #30537 3
head(euks_tax) 
Brad_ix=euks_tax$Feature.ID[grep("Braarudo", euks_tax$Taxon)]

#Read in ASV data (only used to figure out when UCYN-A present/ absent)
#Use interpolated data from 5m
proks_data=read.table("int_CLR_data/proks_AE_5m_int_10.08.2020.tsv", stringsAsFactors = F, header=T, row.names=1)
dim(proks_data) #78953 125
head(names(proks_data))

#Euks data
euks_data=read.table("int_CLR_data/euks_AE_5m_int_10.09.2020.tsv", stringsAsFactors = F, header=T, row.names=1)
dim(euks_data) #28402 125
head(names(euks_data))
head(rownames(euks_data))

which((names(euks_data)==names(proks_data))==FALSE) #None! All the names match, huzzah!

#Format environmental data so that we have upwelling/ PO + nutrient data available for 5m depth in one object
#nutrient data first
nutrient_data=read.table("ModifiedFiles/SPOT_nutrient_CTD_prod_Hammond_data_04.20.2021.tsv", header=T, stringsAsFactors = F, sep="\t")
dim(nutrient_data) #860 18
names(nutrient_data)

#Subset to pull out 5m data only, and only dates for which we have AE filters
length(which(nutrient_data$DepthBin=="5m" & nutrient_data$Year>=2008 & nutrient_data$Year<=2018)) #125
names(proks_data[c(1, ncol(proks_data))]) #Btwn 2008/03 to 2018/07

#Is this data interpolated?
yr_mo=paste(nutrient_data$Year[which(nutrient_data$DepthBin=="5m" & nutrient_data$Year>=2008 & nutrient_data$Year<=2018)], nutrient_data$Month[which(nutrient_data$DepthBin=="5m" & nutrient_data$Year>=2008 & nutrient_data$Year<=2018)], sep="_")
head(yr_mo)
tail(yr_mo)
length(yr_mo) #125
grep("2018.5", yr_mo) #118

full_list=paste(c(rep("2008", times=length(c(3:12))), rep(c(2009:2018), times=12)), c(3:12, rep(1:12, times=length(2009:2018))), sep="_") #130 elements long
yr_mo==full_list
length(full_list) #130 #And there are 125 AE samples in asv data #Last date is May 2018 
tail(yr_mo) #Goes till Dec 2018 - 7 elements too many in nutrient data, will need to remove
grep("2018.5", yr_mo) #118

#Which dates are missing?
full_list %in% yr_mo
which(full_list %in% yr_mo==FALSE)
full_list[which(full_list %in% yr_mo==FALSE)] #DAMN, there are the missing dates!

#Oh right, need to actually subset to just take 5m data
nutrient_subs=nutrient_data[which(nutrient_data$DepthBin=="5m" & nutrient_data$Year>=2008 & nutrient_data$Year<=2018),]
dim(nutrient_subs) #125 18
nutrient_subs$Year_Month=paste(nutrient_subs$Year, nutrient_subs$Month, sep="_")
nutrient_subs$Year_Month #Ends on Dec 2018
grep("2018_7", nutrient_subs$Year_Month) #120
nutrient_subs=nutrient_subs[c(1:grep("2018_7", nutrient_subs$Year_Month)),]
dim(nutrient_subs) #120 19
head(nutrient_subs$Year_Month)
nutrient_subs=nutrient_subs[-1,]
dim(nutrient_subs) #119 19

#PO data
PO_data=read.table("ModifiedFiles/Upwelling_Chl_NO2_Data_CA_Coast_2000_2018.tsv", sep="\t", stringsAsFactors = F, header=T)
PO_data$Year_Month[which(PO_data$Year>=2008 & PO_data$Year<=2018)]
dim(PO_data) #245 19
#Are there any dates missing?
full_list=paste(c(rep("2008", times=length(c(3:12))), rep(c(2009:2018), times=12)), c(3:12, rep(1:12, times=length(2009:2018))), sep="_")
length(full_list) #130
length(PO_data$Year_Month[which(PO_data$Year>=2008 & PO_data$Year<=2018)]) #132
full_list %in% PO_data$Year_Month[which(PO_data$Year>=2008 & PO_data$Year<=2018)] #ALL TRUE! :) no dates missing!

#Paste into PO_data, which is complete/does not need to be linearly interpolated (see below)
names(PO_data)
nutrient_data$Year_Month=paste(nutrient_data$Year, nutrient_data$Month, sep="_")
dim(nutrient_data) #860 19
head(nutrient_data$Year_Month)
length(unique(nutrient_data$Year_Month)) #189

#Subset to only pull out dates I want
PO_subs=PO_data[which(PO_data$Year>=2008 & PO_data$Year<=2018),]
dim(PO_subs) #132 19
head(PO_subs$Year_Month) #Starts on 1/2008
head(nutrient_subs$Year_Month) #Starts on 3/2008, subset further, above

tail(nutrient_subs$Year_Month) #07/2018
tail(PO_subs$Year_Month) #12/2018, subset further

PO_subs=PO_subs[grep("2008_4", PO_subs$Year_Month):grep("2018_7",PO_subs$Year_Month),]
dim(PO_subs) #124 19
dim(nutrient_subs) #119 19

library("dplyr")
env_data=full_join(PO_subs, nutrient_subs, by="Year_Month")
dim(env_data) #124 37
names(env_data)
head(env_data$Year_Month) #2008_4
tail(env_data$Year_Month)
#env_data has two too many dates 

asv_yr_mo=c()
for(i in c(1:ncol(proks_data))){
  a=strsplit(names(proks_data)[i], fixed=T, split=".")
  yr=a[[1]][2]
  if(a[[1]][3]<10){
    b=strsplit(a[[1]][3], split="", fixed=T)
    mo=b[[1]][2]
  }else{
    mo=a[[1]][3]
  }
  d=paste(yr, mo, sep="_")
  asv_yr_mo=c(asv_yr_mo, d)
}
length(asv_yr_mo) #125
length(env_data$Year_Month) #124
head(asv_yr_mo) #2008_3
tail(asv_yr_mo) 
tail(env_data$Year_Month)
length(env_data$Year_Month) #124

asv_yr_mo %in% env_data$Year_Month #First one is false
tail(env_data$Year_Month)

#Subset to remove those last two months from env data, then also remove data with too many NAs
env_subs=env_data[-c(grep("2018_6|2018_7", env_data$Year_Month)),]
tail(env_subs$Year_Month)

#First get rid of obvious data that we do not need
names(env_subs)
unique(env_subs$Depth_m)
unique(env_subs$DepthBin)
env_subs=env_subs[,-c(grep("Year.x|Month.x|Year.y|Month.y|Depth|Day", names(env_subs)))]
names(env_subs)
dim(env_subs) #122 30

#Figure out which variables are missing too much data to interpolate
for(i in c(1:ncol(env_subs))){
  print(names(env_subs)[i])
  #print(summary(as.numeric(env_subs[,i])))
  print(length(which(is.na(env_subs[,i]))))
  if(length(which(is.na(env_subs[,i])))<40){
    print(which(is.na(env_subs[,i])))
  }
}
#Get rid of: NO2, NO3, Oxygen, Salinity, Temperature, Si, pH, Fluorescence
#Interpolate: NO2+NO3, PO4, Bacterial production Leu, Thy, MODIS SST, Model NO2, and Model SST

which(is.na(as.numeric(env_subs$Model_NO2))) #1  6  7  8  9 10 12 19 20 49 53
which(is.na(as.numeric(env_subs$Model_SST))) #Same #Interpolate anyway

env_subs=env_subs[,-c(grep("Oxy|Sal|Temperature|Si|pH|Fluor", names(env_subs)))]
#Now separately take out nitrite and nitrate
grep("NO2_uM", names(env_subs))
env_subs=env_subs[,-grep("NO2_uM", names(env_subs))]
names(env_subs)[grep("NO3_uM", names(env_subs))][1] #Remove the first one
grep("NO3_uM", names(env_subs))[1]
env_subs=env_subs[,-grep("NO3_uM", names(env_subs))[1]]
names(env_subs)
dim(env_subs) #122 22

#Linearly interpolate missing nutrient data
library(zoo)
env_subs$PO4=na.approx(env_subs$PO4_uM)
env_subs$NO2_NO3=na.approx(env_subs$NO2.NO3_uM)
env_subs$BactProd_Leu=na.approx(env_subs$BactProd_Leu)
env_subs$BactProd_Thy=na.approx(env_subs$BactProd_Thy)
env_subs$MODIS_SST=na.approx(env_subs$MODIS_SST)

env_subs$Model_NO2[grep("na", env_subs$Model_NO2)]=NA #Shoot, the first one is missing, na.approx won't work
env_subs$Model_NO2[2:nrow(env_subs)]=na.approx(as.numeric(env_subs$Model_NO2[2:nrow(env_subs)]))

env_subs$Model_SST[grep("na", env_subs$Model_SST)]=NA
env_subs$Model_SST[2:nrow(env_subs)]=na.approx(as.numeric(env_subs$Model_SST[2:nrow(env_subs)]))

#Check if any columns have NAs in them
for(i in c(1:ncol(env_subs))){
  print(colnames(env_subs)[i])
  print(length(which(is.na(env_subs[,i]))))
} #Model SST and Model Chl both missing the first date #Sample ID also missing a few
dim(env_subs)

#Create vectors in env data indicating whether UCYN-A1 and UCYN-A2 present or absent
#Also rows for abundance (absent, low, high abundance) #Low vs. high= is it above or below mean abundance
A1_pres=rep("Absent", times=nrow(env_subs))
length(A1_pres) #122 
env_subs$Year_Month[nrow(env_subs)] #Last one is 5/2018
which(as.numeric(proks_data[grep(UCYNA_ix[1], rownames(proks_data)),2:grep("2018.05", colnames(proks_data))])>0.0001) #52 dates
which(as.numeric(proks_data[grep(UCYNA_ix[1], rownames(proks_data)),2:grep("2018.05", colnames(proks_data))])>0) #51 dates
A1_pres[which(as.numeric(proks_data[grep(UCYNA_ix[1], rownames(proks_data)),2:grep("2018.05", colnames(proks_data))])>0.0001)]="Present"

A2_pres=c(rep("Absent", times=nrow(env_subs)))
length(A2_pres) #122
which(as.numeric(proks_data[grep(UCYNA_ix[5], rownames(proks_data)),2:grep("2018.05", colnames(proks_data))])>0.0001)
A2_pres[which(as.numeric(proks_data[grep(UCYNA_ix[5], rownames(proks_data)),2:grep("2018.05", colnames(proks_data))])>0.0001)]="Present"

#Now relative abundance
A1_abund=rep("a_Absent", times=nrow(env_subs))
A1_abund[which(as.numeric(proks_data[grep(UCYNA_ix[1], rownames(proks_data)),2:grep("2018.05", colnames(proks_data))])>0.0001)]="b_Low"
summary(as.numeric(proks_data[grep(UCYNA_ix[1], rownames(proks_data)),2:grep("2018.05", colnames(proks_data))]))
length(which(as.numeric(proks_data[grep(UCYNA_ix[1], rownames(proks_data)),2:grep("2018.05", colnames(proks_data))])>mean(as.numeric(proks_data[grep(UCYNA_ix[1], rownames(proks_data)),2:grep("2018.05", colnames(proks_data))])))) #26 #The mean relative abundance is 0.0009882487
A1_abund[which(as.numeric(proks_data[grep(UCYNA_ix[1], rownames(proks_data)),2:grep("2018.05", colnames(proks_data))])>mean(as.numeric(proks_data[grep(UCYNA_ix[1], rownames(proks_data)),2:grep("2018.05", colnames(proks_data))])))]="c_High"

A2_abund=rep("a_Absent", times=nrow(env_subs))
A2_abund[which(as.numeric(proks_data[grep(UCYNA_ix[5], rownames(proks_data)),2:grep("2018.05", colnames(proks_data))])>0.0001)]="b_Low"
A2_abund[which(as.numeric(proks_data[grep(UCYNA_ix[5],rownames(proks_data)),2:grep("2018.05",colnames(proks_data))])>mean(as.numeric(proks_data[grep(UCYNA_ix[5], rownames(proks_data)),2:grep("2018.05", colnames(proks_data))])))]="c_High" #0.000264276=mean of A2 relabund

#Look at presence/ absence of Brad
Brad7_pres=rep("Absent", times=length(2:grep("2018.05", colnames(euks_data))))
Brad7_pres[which(as.numeric(euks_data[grep(Brad_ix[7], rownames(euks_data)),2:grep("2018.05",colnames(euks_data))])>0.0001)]="Present"

Brad3_pres=rep("Absent", times=length(2:grep("2018.05", colnames(euks_data))))
Brad3_pres[which(as.numeric(euks_data[grep(Brad_ix[4], rownames(euks_data)),2:grep("2018.05",colnames(euks_data))])>0.0001)]="Present"

#Also add in a column for Lepidodinium
#First, find the hash 
eLSA_data=read.table("eLSA_output/piecewise_eLSA/CLR_5m/parsed_eLSA_output_CLR_5m_w_Hammond_abio_02.27.2021.tsv", sep="\t", header=T, stringsAsFactors = F)
dim(eLSA_data) #5240 33
names(eLSA_data)
unique(eLSA_data$Y[grep("Lepido", eLSA_data$Y_tax)])
unique(eLSA_data$X[grep("Lepido", eLSA_data$X_tax)]) #FOUND IT! :) 
lepido=unique(eLSA_data$X[grep("Lepido", eLSA_data$X_tax)])

grep(lepido, rownames(euks_data))
Lepido_pres=rep("Absent", times=length(2:grep("2018.05", colnames(euks_data))))
Lepido_pres[which(as.numeric(euks_data[grep(lepido, rownames(euks_data)), 2:grep("2018.05", colnames(euks_data))])>0.0001)]="Present"
length(which(Lepido_pres=="Present"))

#Add them in as columns
env_subs$A1_pres=A1_pres
env_subs$A2_pres=A2_pres
env_subs$A1_abund=A1_abund
env_subs$A2_abund=A2_abund
env_subs$Brad7_pres=Brad7_pres
env_subs$Brad3_pres=Brad3_pres
env_subs$Lepido_pres=Lepido_pres

names(env_subs)

#WRITE THIS TABLE OUT
write.table(file="int_CLR_data/int_env_data_UCYNA_Brad_Lepido_pres_abs_cutoff_0.0001_06.27.2022.tsv", quote=F, row.names=F, x=env_subs, sep="\t")

#2. Plot mean environmental variable on dates UCYN-A present/ absent------
#install.packages("tidyverse")
library(tidyverse)
#install.packages("ggplot2")
library(ggplot2)
#install.packages("dplyr")
library(dplyr)

#Chemicals/ nutrients: PO4, NO2+NO3-----
#UCYN-A1 N species
N_table = env_subs %>% group_by(A1_pres) %>% summarize(avgN=mean(NO2_NO3),
                                                       SD=sd(NO2_NO3), n=n(), 
                                                       SE=sd(NO2_NO3/sqrt(n())))
#Try plotting just the mean, using a circle
N_plot=ggplot(N_table, aes(x=A1_pres, y=avgN)) + 
  geom_point(size=15, color="#0099f9") +
  geom_errorbar(aes(ymin=avgN-SE, ymax=avgN+SE), width=0.2) + 
  labs(title="UCYN-A1 vs NO2+NO3", x="", y="[NO2+NO3] (uM)", size=18) +
  theme_classic() +
  theme(axis.title.y=element_text(size=20, face="bold", color="black"), axis.text=element_text(size=20, face="bold", color="black"))
N_plot #YEAH!! No difference, though, lol 

ggsave("/Users/cfletcherh1/Desktop/DataAnalysis/SPOT_16S_18S_eLSA/Plots/04.20.2022_abio_brad_pres_abs/N_1_plot.pdf", height = 4, width = 3.5, units = "in", dpi=300)

#UCYN-A2 N species
N_table = env_subs %>% group_by(A2_pres) %>% summarize(avgN=mean(NO2_NO3),
                                                       SD=sd(NO2_NO3), n=n(), 
                                                       SE=sd(NO2_NO3/sqrt(n())))
#Just the mean
N_plot=ggplot(N_table, aes(x=A2_pres, y=avgN)) + 
  geom_point(size=15, color="#6699CC") +
  geom_errorbar(aes(ymin=avgN-SE, ymax=avgN+SE), width=0.2) + 
  labs(title="UCYN-A2 vs NO2+NO3", x="", y="[NO2+NO3] (uM)", size=18) +
  #theme_classic() + 
  #theme(axis.title.y=element_text(size=20, face="bold", color="black"), axis.text=element_text(size=20, face="bold", color="black")) 
N_plot + theme_classic() + theme(axis.title.y=element_text(size=15), axis.text=element_text(size=15, face="bold"))  #No statistical difference 

ggsave("/Users/cfletcherh1/Desktop/DataAnalysis/SPOT_16S_18S_eLSA/Plots/04.20.2022_abio_brad_pres_abs/N_2_plot.pdf", height = 4, width = 3.50, units = "in", dpi=300)

#UCYN-A1 PO4
P_table = env_subs %>% group_by(A1_pres) %>% summarize(avgP=mean(PO4),
                                                       SD=sd(PO4), n=n(), 
                                                       SE=sd(PO4/sqrt(n())))
#Just the mean #A1 PO4
P_plot=ggplot(P_table, aes(x=A1_pres, y=avgP)) + 
  geom_point(size=15, color="#990099") +
  geom_errorbar(aes(ymin=avgP-SE, ymax=avgP+SE), width=0.2) + 
  labs(title="UCYN-A1 vs PO4", x="", y="[PO4] (uM)", size=18) +
  theme_classic() + 
  theme(axis.title.y=element_text(size=20, face="bold", color="black"), axis.text=element_text(size=20, face="bold", color="black"))  
P_plot  #No statistical difference

ggsave("/Users/cfletcherh1/Desktop/DataAnalysis/SPOT_16S_18S_eLSA/Plots/04.20.2022_abio_brad_pres_abs/P_1_plot.pdf", height = 4, width = 3.50, units = "in", dpi=300)

#UCYN-A2 P
P_table = env_subs %>% group_by(A2_pres) %>% summarize(avgP=mean(PO4),
                                                       SD=sd(PO4), n=n(), 
                                                       SE=sd(PO4/sqrt(n())))
#Just the mean #A2 PO4
P_plot=ggplot(P_table, aes(x=A2_pres, y=avgP)) + 
  geom_point(size=15, color="#CC99CC") +
  geom_errorbar(aes(ymin=avgP-SE, ymax=avgP+SE), width=0.2) + 
  labs(title="UCYN-A2 vs PO4", x="", y="[PO4] (uM)", size=18) +
  theme_classic() + 
  theme(axis.title.y=element_text(size=20, face="bold", color="black"), axis.text=element_text(size=20, face="bold", color="black"))
P_plot

ggsave("/Users/cfletcherh1/Desktop/DataAnalysis/SPOT_16S_18S_eLSA/Plots/04.20.2022_abio_brad_pres_abs/P_2_plot.pdf", height = 4, width = 3.50, units = "in", dpi=300)

#Upwelling: Indices, MEI, Chl, SST, upwelling-derived NO2, bacterial production-----
#A1 vs CUTI ix
CUTI_table = env_subs %>% group_by(A1_pres) %>% summarize(avgCUTI=mean(CUTI_ix),
                                                          SD=sd(CUTI_ix), n=n(), 
                                                          SE=sd(CUTI_ix/sqrt(n())))

#Just the mean #A1
CUTI_plot=ggplot(CUTI_table, aes(x=A1_pres, y=avgCUTI)) + 
  geom_point(size=15, color="#CC6600") +
  geom_errorbar(aes(ymin=avgCUTI-SE, ymax=avgCUTI+SE), width=0.2) + 
  labs(title="UCYN-A1 vs CUTI", x="", y="CUTI", size=18) +
  theme_classic() +
  theme(axis.title.y=element_text(size=20, face="bold", color="black"), axis.text=element_text(size=20, face="bold", color="black"))
CUTI_plot

ggsave("/Users/cfletcherh1/Desktop/DataAnalysis/SPOT_16S_18S_eLSA/Plots/04.20.2022_abio_brad_pres_abs/CUTI_1_plot.pdf", height = 4, width = 3.50, units = "in", dpi=300)


#A2 vs CUTI ix
CUTI_table = env_subs %>% group_by(A2_pres) %>% summarize(avgCUTI=mean(CUTI_ix),
                                                          SD=sd(CUTI_ix), n=n(), 
                                                          SE=sd(CUTI_ix/sqrt(n())))

#A2 vs CUTI avg
CUTI_plot=ggplot(CUTI_table, aes(x=A2_pres, y=avgCUTI)) + 
  geom_point(size=15, color="#FFCC99") +
  geom_errorbar(aes(ymin=avgCUTI-SE, ymax=avgCUTI+SE), width=0.2) + 
  labs(title="UCYN-A2 vs CUTI", x="", y="CUTI", size=18) +
  theme_classic() + 
  theme(axis.title.y=element_text(size=20, face="bold", color="black"), axis.text=element_text(size=20, face="bold", color="black"))
CUTI_plot

ggsave("/Users/cfletcherh1/Desktop/DataAnalysis/SPOT_16S_18S_eLSA/Plots/04.20.2022_abio_brad_pres_abs/CUTI_2_plot.pdf", height = 4, width = 3.50, units = "in", dpi=300)

#A1 vs BEUTI index
BEUTI_table = env_subs %>% group_by(A1_pres) %>% summarize(avgBEUTI=mean(BEUTI_ix),
                                                           SD=sd(BEUTI_ix), n=n(), 
                                                           SE=sd(BEUTI_ix/sqrt(n())))
#Plot just the mean
BEUTI_plot=ggplot(BEUTI_table, aes(x=A1_pres, y=avgBEUTI)) + 
  geom_point(size=15, color="#FF00FF") +
  geom_errorbar(aes(ymin=avgBEUTI-SE, ymax=avgBEUTI+SE), width=0.2) + 
  labs(title="UCYN-A1 vs BEUTI", x="", y="BEUTI", size=18) + 
  theme_classic() + 
  theme(axis.title.y=element_text(size=20, face="bold", color="black"), axis.text=element_text(size=20, face="bold", color="black"))
BEUTI_plot

ggsave("/Users/cfletcherh1/Desktop/DataAnalysis/SPOT_16S_18S_eLSA/Plots/04.20.2022_abio_brad_pres_abs/BEUTI_1_plot.pdf", height = 4, width = 3.50, units = "in", dpi=300)

#A2 vs BEUTI index
BEUTI_table = env_subs %>% group_by(A2_pres) %>% summarize(avgBEUTI=mean(BEUTI_ix),
                                                           SD=sd(BEUTI_ix), n=n(), 
                                                           SE=sd(BEUTI_ix/sqrt(n())))

#Plot just the mean #A2
BEUTI_plot=ggplot(BEUTI_table, aes(x=A2_pres, y=avgBEUTI)) + 
  geom_point(size=15, color="#FF99FF") +
  geom_errorbar(aes(ymin=avgBEUTI-SE, ymax=avgBEUTI+SE), width=0.2) + 
  labs(title="UCYN-A2 vs BEUTI", x="", y="BEUTI", size=18) +
  theme_classic() + 
  theme(axis.title.y=element_text(size=20, face="bold", color="black"), axis.text=element_text(size=20, face="bold", color="black"))
BEUTI_plot 

ggsave("/Users/cfletcherh1/Desktop/DataAnalysis/SPOT_16S_18S_eLSA/Plots/04.20.2022_abio_brad_pres_abs/BEUTI_2_plot.pdf", height = 4, width = 3.50, units = "in", dpi=300)

#A1 vs MEI
MEI_table = env_subs %>% group_by(A1_pres) %>% summarize(avgMEI=mean(MEI_ix),
                                                         SD=sd(MEI_ix), n=n(), 
                                                         SE=sd(MEI_ix/sqrt(n())))
#Plot just the mean
MEI_plot=ggplot(MEI_table, aes(x=A1_pres, y=avgMEI)) + 
  geom_point(size=8, color="#CC0033") +
  geom_errorbar(aes(ymin=avgMEI-SE, ymax=avgMEI+SE), width=0.2, size=1.2) + labs(title="UCYN-A1 vs MEI", x="", y="MEI") 
MEI_plot + theme_classic() + theme(axis.title.y=element_text(size=9), axis.text=element_text(size=7, face="bold", color="black"))

ggsave("/Users/colettef/Desktop/DataAnalysis/SPOT_16S_18S_eLSA/Plots/04.20.2022_abio_brad_pres_abs/formatted/MEI_1_plot.pdf", height = 2, width = 1.5, units = "in", dpi=300)
#ggsave("/Users/cfletcherh1/Desktop/DataAnalysis/SPOT_16S_18S_eLSA/Plots/04.20.2022_abio_brad_pres_abs/MEI_1_plot.pdf", height = 4, width = 3.5, units = "in", dpi=300)

#Now for UCYN-A2 vs MEI
MEI_table = env_subs %>% group_by(A2_pres) %>% summarize(avgMEI=mean(MEI_ix),
                                                         SD=sd(MEI_ix), n=n(), 
                                                         SE=sd(MEI_ix/sqrt(n())))

#Plot just the mean #A2
MEI_plot=ggplot(MEI_table, aes(x=A2_pres, y=avgMEI)) + 
  geom_point(size=8, color="#FF9999") +
  geom_errorbar(aes(ymin=avgMEI-SE, ymax=avgMEI+SE), width=0.2, size=1.2) + labs(title="UCYN-A2 vs MEI", x="", y="MEI") 
MEI_plot + theme_classic() + theme(axis.title.y=element_text(size=9), axis.text=element_text(size=7, face="bold", color="black"))

ggsave("/Users/colettef/Desktop/DataAnalysis/SPOT_16S_18S_eLSA/Plots/04.20.2022_abio_brad_pres_abs/formatted/MEI_2_plot.pdf", height = 2, width = 1.5, units = "in", dpi=300)
#ggsave("/Users/cfletcherh1/Desktop/DataAnalysis/SPOT_16S_18S_eLSA/Plots/04.20.2022_abio_brad_pres_abs/MEI_2_plot.pdf", height = 4, width = 3.5, units = "in", dpi=300)

#MODIS SST #A1
SST_table= env_subs %>% group_by(A1_pres) %>% summarize(avg_SST=mean(MODIS_SST),
                                                        SD=sd(MODIS_SST), n=n(),
                                                        SE=sd(MODIS_SST)/sqrt(n()))
SST_plot=ggplot(SST_table, aes(x=A1_pres, y=avg_SST))+
  geom_point(size=8, color="#0000FF") + 
  geom_errorbar(aes(ymin=avg_SST-SE, ymax=avg_SST+SE), width=0.2, size=1.2) + labs(title="UCYN-A1 vs MODIS SST", x="", y="SST (ºC)")
SST_plot + theme_classic() + theme(axis.title.y=element_text(size=9), axis.text=element_text(size=7, face="bold", color="black"))

ggsave("/Users/colettef/Desktop/DataAnalysis/SPOT_16S_18S_eLSA/Plots/04.20.2022_abio_brad_pres_abs/formatted/MODIS_SST_1_plot.pdf", height=2, width=1.5, units="in", dpi=300)
#ggsave("/Users/cfletcherh1/Desktop/DataAnalysis/SPOT_16S_18S_eLSA/Plots/04.20.2022_abio_brad_pres_abs/MODIS_SST_1_plot.pdf", height = 4, width = 3.5, units = "in", dpi=300)

#A2
SST_table= env_subs %>% group_by(A2_pres) %>% summarize(avg_SST=mean(MODIS_SST),
                                                        SD=sd(MODIS_SST), n=n(),
                                                        SE=sd(MODIS_SST)/sqrt(n()))
SST_plot=ggplot(SST_table, aes(x=A2_pres, y=avg_SST))+
  geom_point(size=8, color="#66CCFF") + 
  geom_errorbar(aes(ymin=avg_SST-SE, ymax=avg_SST+SE), width=0.2, size=1.2) + labs(title="UCYN-A2 vs MODIS SST", x="", y="SST (ºC)")
SST_plot + theme_classic() + theme(axis.title.y=element_text(size=9), axis.text=element_text(size=7, face="bold", color="black"))

ggsave("/Users/colettef/Desktop/DataAnalysis/SPOT_16S_18S_eLSA/Plots/04.20.2022_abio_brad_pres_abs/formatted/MODIS_SST_2_plot.pdf", height=2, width=1.5, units="in", dpi=300)
#ggsave("/Users/cfletcherh1/Desktop/DataAnalysis/SPOT_16S_18S_eLSA/Plots/04.20.2022_abio_brad_pres_abs/MODIS_SST_2_plot.pdf", height = 4, width = 3.5, units = "in", dpi=300)

#MODIS Chl #A1
Chl_table=env_subs %>% group_by(A1_pres) %>% summarize(avgChl=mean(as.numeric(MODIS_Chl)),
                                                       SD=sd(as.numeric(MODIS_Chl)), n=n(),
                                                       SE=sd(as.numeric(MODIS_Chl))/sqrt(n()))

Chl_plot=ggplot(Chl_table, aes(x=A1_pres, y=avgChl)) +
  geom_point(size=8, color="#009E73") + 
  geom_errorbar(aes(ymin=avgChl-SE, ymax=avgChl+SE), width=0.2, size=1.2) + labs(title="UCYN-A1 vs MODIS Chl", x="", y="[Chl] (mg*m-3)")
Chl_plot + theme_classic() + theme(axis.title.y=element_text(size=9), axis.text=element_text(size=7, face="bold", color="black"))

ggsave("/Users/colettef/Desktop/DataAnalysis/SPOT_16S_18S_eLSA/Plots/04.20.2022_abio_brad_pres_abs/formatted/MODIS_Chl_1_plot.pdf", height=2, width=1.5, units="in", dpi=300)
#ggsave("/Users/cfletcherh1/Desktop/DataAnalysis/SPOT_16S_18S_eLSA/Plots/04.20.2022_abio_brad_pres_abs/MODIS_Chl_1_plot.pdf", height = 4, width = 3.5, units = "in", dpi=300)

#A2
Chl_table=env_subs %>% group_by(A2_pres) %>% summarize(avgChl=mean(as.numeric(MODIS_Chl)),
                                                       SD=sd(as.numeric(MODIS_Chl)), n=n(),
                                                       SE=sd(as.numeric(MODIS_Chl))/sqrt(n()))

Chl_plot=ggplot(Chl_table, aes(x=A2_pres, y=avgChl))+
  geom_point(size=8, color="#66CC99") + 
  geom_errorbar(aes(ymin=avgChl-SE, ymax=avgChl+SE), width=0.2, size=1.2) + labs(title="UCYN-A2 vs MODIS Chl", x="", y="[Chl] (mg*m-3)")
Chl_plot + theme_classic() + theme(axis.title.y=element_text(size=9), axis.text=element_text(size=7, face="bold", color="black"))

ggsave("/Users/colettef/Desktop/DataAnalysis/SPOT_16S_18S_eLSA/Plots/04.20.2022_abio_brad_pres_abs/formatted/MODIS_Chl_2_plot.pdf", height=2, width=1.5, units="in", dpi=300)
#ggsave("/Users/cfletcherh1/Desktop/DataAnalysis/SPOT_16S_18S_eLSA/Plots/04.20.2022_abio_brad_pres_abs/MODIS_Chl_2_plot.pdf", height = 4, width = 3.5, units = "in", dpi=300)

#Bacterial production leucine 
leu_table=env_subs %>% group_by(A1_pres) %>% summarize(avg_leu=mean(BactProd_Leu), 
                                                       SD=sd(BactProd_Leu), n=n(),
                                                       SE=sd(BactProd_Leu)/sqrt(n()))

Leu_plot=ggplot(leu_table, aes(x=A1_pres, y=avg_leu)) + 
  geom_point(size=8, color="#999999") + 
  geom_errorbar(aes(ymin=avg_leu-SE, ymax=avg_leu+SE), width=0.2, size=1.2) + labs(title="UCYN-A1 vs Bacterial Production - Leucine", x="", y="Bact Prod. (cells/mL/day)")
Leu_plot + theme_classic() + theme(axis.title.y=element_text(size=9), axis.text=element_text(size=7, face="bold"))

ggsave("/Users/colettef/Desktop/DataAnalysis/SPOT_16S_18S_eLSA/Plots/04.20.2022_abio_brad_pres_abs/formatted/BactProd_Leu_1_plot.pdf", height=2, width=1.5, units="in", dpi=300)
ggsave("/Users/cfletcherh1/Desktop/DataAnalysis/SPOT_16S_18S_eLSA/Plots/04.20.2022_abio_brad_pres_abs/BactProd_Leu_1_plot.pdf", height = 4, width = 3.5, units = "in", dpi=300)

#A2 
leu_table=env_subs %>% group_by(A2_pres) %>% summarize(avg_leu=mean(BactProd_Leu), 
                                                       SD=sd(BactProd_Leu), n=n(),
                                                       SE=sd(BactProd_Leu)/sqrt(n()))

Leu_plot=ggplot(leu_table, aes(x=A2_pres, y=avg_leu)) + 
  geom_point(size=10, col="#999999") + 
  geom_errorbar(aes(ymin=avg_leu-SE, ymax=avg_leu+SE), width=0.2) + labs(title="UCYN-A2 vs Bacterial Production - Leucine", x="", y="Bact Prod. (cells/mL/day)")
Leu_plot + theme_classic() + theme(axis.title.y=element_text(size=15), axis.text=element_text(size=15, face="bold"))

ggsave("/Users/cfletcherh1/Desktop/DataAnalysis/SPOT_16S_18S_eLSA/Plots/04.20.2022_abio_brad_pres_abs/BactProd_Leu_2_plot.pdf", height = 4, width = 3.5, units = "in", dpi=300)

#Bacterial production thymidine 
thy_table=env_subs %>% group_by(A1_pres) %>% summarize(avg_thy=mean(BactProd_Thy), 
                                                       SD=sd(BactProd_Thy), n=n(),
                                                       SE=sd(BactProd_Thy)/sqrt(n()))

Thy_plot=ggplot(thy_table, aes(x=A1_pres, y=avg_thy)) + 
  geom_point(size=10, color="#999999") + 
  geom_errorbar(aes(ymin=avg_thy-SE, ymax=avg_thy+SE), width=0.2) + labs(title="UCYN-A1 vs Bacterial Production - Thymidine", x="", y="Bact Prod. (cells/mL/day)")
Thy_plot + theme_classic() + theme(axis.title.y=element_text(size=15), axis.text=element_text(size=15, face="bold"))

ggsave("/Users/cfletcherh1/Desktop/DataAnalysis/SPOT_16S_18S_eLSA/Plots/04.20.2022_abio_brad_pres_abs/BactProd_Thy_1_plot.pdf", height = 4, width = 3.5, units = "in", dpi=300)

#A2 
thy_table=env_subs %>% group_by(A2_pres) %>% summarize(avg_thy=mean(BactProd_Thy), 
                                                       SD=sd(BactProd_Thy), n=n(),
                                                       SE=sd(BactProd_Thy)/sqrt(n()))

Thy_plot=ggplot(thy_table, aes(x=A2_pres, y=avg_thy)) + 
  geom_point(size=10, color="#999999") + 
  geom_errorbar(aes(ymin=avg_thy-SE, ymax=avg_thy+SE), width=0.2) + labs(title="UCYN-A2 vs Bacterial Production - Thymidine", x="", y="Bact Prod. (cells/mL/day)")
Thy_plot + theme_classic() + theme(axis.title.y=element_text(size=15), axis.text=element_text(size=15, face="bold"))

ggsave("/Users/cfletcherh1/Desktop/DataAnalysis/SPOT_16S_18S_eLSA/Plots/04.20.2022_abio_brad_pres_abs/BactProd_Thy_2_plot.pdf", height = 4, width = 3.5, units = "in", dpi=300)

#Upwelling derived NO2 #A1
no2_table=env_subs %>% group_by(A1_pres) %>% summarize(avg_no2=mean(as.numeric(env_subs$Model_NO2[2:nrow(env_subs)])), 
                                                       SD=sd(as.numeric(env_subs$Model_NO2[2:nrow(env_subs)]), n=n()), 
                                                       SE=sd(as.numeric(env_subs$Model_NO2[2:nrow(env_subs)])/sqrt(n())))
no2_plot=ggplot(no2_table, aes(x=A1_pres, y=avg_no2)) + 
  geom_point(size=10, color="#999999") + 
  geom_errorbar(aes(ymin=avg_no2-SE, ymax=avg_no2+SE), width=0.2) + labs(title="UCYN-A1 vs Upwelling-derived NO2", x="", y="[NO2] (uM)")
no2_plot + theme_classic() + theme(axis.title.y=element_text(size=15), axis.text=element_text(size=15, face="bold"))

ggsave("/Users/cfletcherh1/Desktop/DataAnalysis/SPOT_16S_18S_eLSA/Plots/04.20.2022_abio_brad_pres_abs/Model_NO2_1_plot.pdf", height = 4, width = 3.5, units = "in", dpi=300)

#A2
no2_table=env_subs %>% group_by(A2_pres) %>% summarize(avg_no2=mean(as.numeric(env_subs$Model_NO2[2:nrow(env_subs)])), 
                                                       SD=sd(as.numeric(env_subs$Model_NO2[2:nrow(env_subs)]), n=n()), 
                                                       SE=sd(as.numeric(env_subs$Model_NO2[2:nrow(env_subs)])/sqrt(n())))
no2_plot=ggplot(no2_table, aes(x=A2_pres, y=avg_no2)) + 
  geom_point(size=10, color="#999999") + 
  geom_errorbar(aes(ymin=avg_no2-SE, ymax=avg_no2+SE), width=0.2) + labs(title="UCYN-A2 vs Upwelling-derived NO2", x="", y="[NO2] (uM)")
no2_plot + theme_classic() + theme(axis.title.y=element_text(size=15), axis.text=element_text(size=15, face="bold"))

ggsave("/Users/cfletcherh1/Desktop/DataAnalysis/SPOT_16S_18S_eLSA/Plots/04.20.2022_abio_brad_pres_abs/Model_NO2_2_plot.pdf", height = 4, width = 3.5, units = "in", dpi=300)

#Model-derived SST 
sst_table=env_subs %>% group_by(A1_pres) %>% summarize(avg_sst=mean(as.numeric(env_subs$Model_SST[2:nrow(env_subs)])),
                                                       SD=sd(as.numeric(env_subs$Model_SST[2:nrow(env_subs)]), n=n()),
                                                       SE=sd(as.numeric(env_subs$Model_SST[2:nrow(env_subs)])/sqrt(n())))
sst_plot=ggplot(sst_table, aes(x=A1_pres, y=avg_sst)) + 
  geom_point(size=10, color="#999999") + 
  geom_errorbar(aes(ymin=avg_sst-SE, ymax=avg_sst+SE), width=0.2) + labs(title="UCYN-A1 vs model-derived SST", x="", y="SST (ºC)")
sst_plot + theme_classic() + theme(axis.title.y=element_text(size=15), axis.text=element_text(size=15, face="bold"))

ggsave("/Users/cfletcherh1/Desktop/DataAnalysis/SPOT_16S_18S_eLSA/Plots/04.20.2022_abio_brad_pres_abs/Model_SST_1_plot.pdf", height = 4, width = 3.5, units = "in", dpi=300)

#A2
sst_table=env_subs %>% group_by(A2_pres) %>% summarize(avg_sst=mean(as.numeric(env_subs$Model_SST[2:nrow(env_subs)])),
                                                       SD=sd(as.numeric(env_subs$Model_SST[2:nrow(env_subs)]), n=n()),
                                                       SE=sd(as.numeric(env_subs$Model_SST[2:nrow(env_subs)])/sqrt(n())))
sst_plot=ggplot(sst_table, aes(x=A2_pres, y=avg_sst)) + 
  geom_point(size=10, color="#999999") + 
  geom_errorbar(aes(ymin=avg_sst-SE, ymax=avg_sst+SE), width=0.2) + labs(title="UCYN-A2 vs model-derived SST", x="", y="SST (ºC)")
sst_plot + theme_classic() + theme(axis.title.y=element_text(size=15), axis.text=element_text(size=15, face="bold"))

ggsave("/Users/cfletcherh1/Desktop/DataAnalysis/SPOT_16S_18S_eLSA/Plots/04.20.2022_abio_brad_pres_abs/Model_SST_2_plot.pdf", height = 4, width = 3.5, units = "in", dpi=300)


#Wind factors: pmsl, ektry, uv_mean, etc.----
#pmsl 
pmsl_table= env_subs %>% group_by(A1_pres) %>% summarize(avgpmsl=mean(as.numeric(pmsl)),
                                                         SD=sd(as.numeric(pmsl)), n=n(),
                                                         SE=sd(as.numeric(pmsl)/sqrt(n())))

#Plot just the mean
pmsl_plot=ggplot(pmsl_table, aes(x=A1_pres, y=avgpmsl)) + 
  geom_point(size=10, color="#999999") + 
  geom_errorbar(aes(ymin=avgpmsl-SE, ymax=avgpmsl+SE), width=0.2) + labs(title="UCYN-A1 vs pmsl", x="", y="pmsl")
pmsl_plot + theme_classic() + theme(axis.title.y=element_text(size=15), axis.text=element_text(size=15, face="bold"))

ggsave("/Users/cfletcherh1/Desktop/DataAnalysis/SPOT_16S_18S_eLSA/Plots/04.20.2022_abio_brad_pres_abs/pmsl_1_plot.pdf", height = 4, width = 3.5, units = "in", dpi=300)

#Now for UCYN-A2
pmsl_table= env_subs %>% group_by(A2_pres) %>% summarize(avgpmsl=mean(as.numeric(pmsl)),
                                                         SD=sd(as.numeric(pmsl)), n=n(),
                                                         SE=sd(as.numeric(pmsl)/sqrt(n())))

pmsl_plot=ggplot(pmsl_table, aes(x=A2_pres, y=avgpmsl)) + 
  geom_point(size=10, color="#999999") + 
  geom_errorbar(aes(ymin=avgpmsl-SE, ymax=avgpmsl+SE), width=0.2) + labs(title="UCYN-A2 vs pmsl", x="", y="pmsl")
pmsl_plot + theme_classic() + theme(axis.title.y=element_text(size=15), axis.text=element_text(size=15, face="bold"))

ggsave("/Users/cfletcherh1/Desktop/DataAnalysis/SPOT_16S_18S_eLSA/Plots/04.20.2022_abio_brad_pres_abs/pmsl_2_plot.pdf", height = 4, width = 3.5, units = "in", dpi=300)

#u_mean
names(env_subs)
u_table= env_subs %>% group_by(A1_pres) %>% summarize(avg_u=mean(as.numeric(u_mean)), 
                                                      SD=sd(as.numeric(u_mean)), n=n(),
                                                      SE=sd(as.numeric(u_mean))/sqrt(n()))
u_plot=ggplot(u_table, aes(x=A1_pres, y=avg_u)) + 
  geom_point(size=10, color="#999999") + 
  geom_errorbar(aes(ymin=avg_u-SE, ymax=avg_u+SE), width=0.2) + labs(title="UCYN-A1 vs u_mean", x="", y="E-W wind (m/s)")
u_plot + theme_classic() + theme(axis.title.y=element_text(size=15), axis.text=element_text(size=15, face="bold"))

ggsave("/Users/cfletcherh1/Desktop/DataAnalysis/SPOT_16S_18S_eLSA/Plots/04.20.2022_abio_brad_pres_abs/u_mean_1_plot.pdf", height=4, width=3.5, units="in", dpi=300)

#HMM
sort(env_subs$A1_pres[which(as.numeric(env_subs$u_mean)<0)]) #When u_mean is negative (blowing E->W, offshore from SPOT), UCYN-A1 mostly present

#A2
u_table= env_subs %>% group_by(A2_pres) %>% summarize(avg_u=mean(as.numeric(u_mean)), 
                                                      SD=sd(as.numeric(u_mean)), n=n(),
                                                      SE=sd(as.numeric(u_mean))/sqrt(n()))
u_plot=ggplot(u_table, aes(x=A2_pres, y=avg_u)) + 
  geom_point(size=10, color="#999999") + 
  geom_errorbar(aes(ymin=avg_u-SE, ymax=avg_u+SE), width=0.2) + labs(title="UCYN-A2 vs u_mean", x="", y="E-W wind (m/s)")
u_plot + theme_classic() + theme(axis.title.y=element_text(size=15), axis.text=element_text(size=15, face="bold"))

ggsave("/Users/cfletcherh1/Desktop/DataAnalysis/SPOT_16S_18S_eLSA/Plots/04.20.2022_abio_brad_pres_abs/u_mean_2_plot.pdf", height=4, width=3.5, units="in", dpi=300)

#V mean #N-S wind component
#A1
v_table= env_subs %>% group_by(A1_pres) %>% summarize(avg_v=mean(as.numeric(v_mean)),
                                                      SD=sd(as.numeric(v_mean)), n=n(),
                                                      SE=SD/sqrt(n()))
v_plot=ggplot(v_table, aes(A1_pres, y=avg_v)) + 
  geom_point(size=10, color="#999999") + 
  geom_errorbar(aes(ymin=avg_v-SE, ymax=avg_v+SE), width=0.2) + labs(title="UCYN-A1 vs v_mean", x="", y="N-S wind (m/s)")
v_plot + theme_classic() + theme(axis.title.y=element_text(size=15), axis.text=element_text(size=15, face="bold"))

ggsave("/Users/cfletcherh1/Desktop/DataAnalysis/SPOT_16S_18S_eLSA/Plots/04.20.2022_abio_brad_pres_abs/v_mean_1_plot.pdf", height=4, width=3.5, units="in", dpi=300) #North= more positive

#A2
v_table= env_subs %>% group_by(A2_pres) %>% summarize(avg_v=mean(as.numeric(v_mean)),
                                                      SD=sd(as.numeric(v_mean)), n=n(),
                                                      SE=SD/sqrt(n()))
v_plot=ggplot(v_table, aes(A2_pres, y=avg_v)) + 
  geom_point(size=10, color="#999999") + 
  geom_errorbar(aes(ymin=avg_v-SE, ymax=avg_v+SE), width=0.2) + labs(title="UCYN-A2 vs v_mean", x="", y="N-S wind (m/s)")
v_plot + theme_classic() + theme(axis.title.y=element_text(size=15), axis.text=element_text(size=15, face="bold"))

ggsave("/Users/cfletcherh1/Desktop/DataAnalysis/SPOT_16S_18S_eLSA/Plots/04.20.2022_abio_brad_pres_abs/v_mean_2_plot.pdf", height=4, width=3.5, units="in", dpi=300)

#UV mean
uv_table= env_subs %>% group_by(A1_pres) %>% summarize(avg_uv=mean(as.numeric(uv_mag_mean)), 
                                                       SD=sd(as.numeric(uv_mag_mean)), n=n(), 
                                                       SE=sd(as.numeric(uv_mag_mean))/sqrt(n()))

uv_plot=ggplot(uv_table, aes(A1_pres, y=avg_uv)) +
  geom_point(size=10, color="#999999") + 
  geom_errorbar(aes(ymin=avg_uv-SE, ymax=avg_uv+SE), width=0.2) + labs(title="UCYN-A1 vs uv_mean", x="", y="uv_mean (m/s")
uv_plot + theme_classic() + theme(axis.title.y=element_text(size=15), axis.text=element_text(size=15, face="bold"))

ggsave("/Users/cfletcherh1/Desktop/DataAnalysis/SPOT_16S_18S_eLSA/Plots/04.20.2022_abio_brad_pres_abs/uv_mean_1_plot.pdf", height=4, width=3.5, units="in", dpi=300)

#A2
uv_table= env_subs %>% group_by(A2_pres) %>% summarize(avg_uv=mean(as.numeric(uv_mag_mean)), 
                                                       SD=sd(as.numeric(uv_mag_mean)), n=n(), 
                                                       SE=sd(as.numeric(uv_mag_mean))/sqrt(n()))

uv_plot=ggplot(uv_table, aes(A2_pres, y=avg_uv)) +
  geom_point(size=10, color="#999999") + 
  geom_errorbar(aes(ymin=avg_uv-SE, ymax=avg_uv+SE), width=0.2) + labs(title="UCYN-A2 vs uv_mean", x="", y="uv_mean (m/s")
uv_plot + theme_classic() + theme(axis.title.y=element_text(size=15), axis.text=element_text(size=15, face="bold"))

ggsave("/Users/cfletcherh1/Desktop/DataAnalysis/SPOT_16S_18S_eLSA/Plots/04.20.2022_abio_brad_pres_abs/uv_mean_2_plot.pdf", height=4, width=3.5, units="in", dpi=300)

#taux #A1
taux_table= env_subs %>% group_by(A1_pres) %>% summarize(avg_taux=mean(as.numeric(taux_mean)),
                                                        SD=sd(as.numeric(taux_mean)), n=n(),
                                                        SE=sd(as.numeric(taux_mean))/sqrt(n()))
taux_plot=ggplot(taux_table, aes(x=A1_pres, y=avg_taux)) + 
  geom_point(size=10, color="#999999") + 
  geom_errorbar(aes(ymin=avg_taux-SE, ymax=avg_taux+SE), width=0.2) + labs(title="UCYN-A1 vs taux_mean", x="", y="E-W wind stress (N/m^2)")
taux_plot + theme_classic() + theme(axis.title.y=element_text(size=15), axis.text=element_text(size=15, face="bold"))

ggsave("/Users/cfletcherh1/Desktop/DataAnalysis/SPOT_16S_18S_eLSA/Plots/04.20.2022_abio_brad_pres_abs/taux_1_plot.pdf", height=4, width=3.5, units="in", dpi=300)

#A2
taux_table= env_subs %>% group_by(A2_pres) %>% summarize(avg_taux=mean(as.numeric(taux_mean)),
                                                         SD=sd(as.numeric(taux_mean)), n=n(),
                                                         SE=sd(as.numeric(taux_mean))/sqrt(n()))
taux_plot=ggplot(taux_table, aes(x=A2_pres, y=avg_taux)) + 
  geom_point(size=10, color="#999999") + 
  geom_errorbar(aes(ymin=avg_taux-SE, ymax=avg_taux+SE), width=0.2) + labs(title="UCYN-A2 vs taux_mean", x="", y="E-W wind stress (N/m^2)")
taux_plot + theme_classic() + theme(axis.title.y=element_text(size=15), axis.text=element_text(size=15, face="bold"))

ggsave("/Users/cfletcherh1/Desktop/DataAnalysis/SPOT_16S_18S_eLSA/Plots/04.20.2022_abio_brad_pres_abs/taux_2_plot.pdf", height=4, width=3.5, units="in", dpi=300)

#tauy #A1
tauy_table= env_subs %>% group_by(A1_pres) %>% summarize(avg_tauy=mean(as.numeric(tauy_mean)),
                                                         SD=sd(as.numeric(tauy_mean)), n=n(),
                                                         SE=sd(as.numeric(tauy_mean))/sqrt(n()))
tauy_plot=ggplot(tauy_table, aes(x=A1_pres, y=avg_tauy)) + 
  geom_point(size=10, color="#999999") + 
  geom_errorbar(aes(ymin=avg_tauy-SE, ymax=avg_tauy+SE), width=0.2) + labs(title="UCYN-A1 vs tauy_mean", x="", y="E-W wind stress (N/m^2)")
tauy_plot + theme_classic() + theme(axis.title.y=element_text(size=15), axis.text=element_text(size=15, face="bold"))

ggsave("/Users/cfletcherh1/Desktop/DataAnalysis/SPOT_16S_18S_eLSA/Plots/04.20.2022_abio_brad_pres_abs/tauy_1_plot.pdf", height=4, width=3.5, units="in", dpi=300)

#A2
tauy_table= env_subs %>% group_by(A2_pres) %>% summarize(avg_tauy=mean(as.numeric(tauy_mean)),
                                                         SD=sd(as.numeric(tauy_mean)), n=n(),
                                                         SE=sd(as.numeric(tauy_mean))/sqrt(n()))
tauy_plot=ggplot(tauy_table, aes(x=A2_pres, y=avg_tauy)) + 
  geom_point(size=10, color="#999999") + 
  geom_errorbar(aes(ymin=avg_tauy-SE, ymax=avg_tauy+SE), width=0.2) + labs(title="UCYN-A2 vs tauy_mean", x="", y="E-W wind stress (N/m^2)")
tauy_plot + theme_classic() + theme(axis.title.y=element_text(size=15), axis.text=element_text(size=15, face="bold"))

ggsave("/Users/cfletcherh1/Desktop/DataAnalysis/SPOT_16S_18S_eLSA/Plots/04.20.2022_abio_brad_pres_abs/tauy_2_plot.pdf", height=4, width=3.5, units="in", dpi=300)

#Ektrx (E-W component of Ekman transport)
ektrx_table=env_subs %>% group_by(A1_pres) %>% summarize(avg_ektrx=mean(as.numeric(ektrx)),
                                                         SD=sd(as.numeric(ektrx)), n=n(),
                                                         SE=sd(as.numeric(ektrx))/sqrt(n()))
ektrx_plot=ggplot(ektrx_table, aes(x=A1_pres, y=avg_ektrx)) + 
  geom_point(size=8, color="#999999") + 
  geom_errorbar(aes(ymin=avg_ektrx-SE, ymax=avg_ektrx+SE), width=0.2, size=1.2) + labs(title="UCYN-A1 vs. ektrx", x="", y="E-W component of Ekman \n transport (kg/m/s)")
ektrx_plot + theme_classic() + theme(axis.title.y=element_text(size=9), axis.text=element_text(size=7, face="bold"))

ggsave("/Users/colettef/Desktop/DataAnalysis/SPOT_16S_18S_eLSA/Plots/04.20.2022_abio_brad_pres_abs/formatted/ektrx_1_plot.pdf", height=2, width=1.5, units="in", dpi=300)
ggsave("/Users/cfletcherh1/Desktop/DataAnalysis/SPOT_16S_18S_eLSA/Plots/04.20.2022_abio_brad_pres_abs/ektrx_1_plot.pdf", height=4, width=3.5, units="in", dpi=300)

#A2
#Ektrx (E-W component of Ekman transport)
ektrx_table=env_subs %>% group_by(A2_pres) %>% summarize(avg_ektrx=mean(as.numeric(ektrx)),
                                                         SD=sd(as.numeric(ektrx)), n=n(),
                                                         SE=sd(as.numeric(ektrx))/sqrt(n()))
ektrx_plot=ggplot(ektrx_table, aes(x=A2_pres, y=avg_ektrx)) + 
  geom_point(size=10, color="#999999") + 
  geom_errorbar(aes(ymin=avg_ektrx-SE, ymax=avg_ektrx+SE), width=0.2) + labs(title="UCYN-A2 vs. ektrx", x="", y="E-W component of Ekman \n transport (kg/m/s)")
ektrx_plot + theme_classic() + theme(axis.title.y=element_text(size=15), axis.text=element_text(size=15, face="bold"))

ggsave("/Users/cfletcherh1/Desktop/DataAnalysis/SPOT_16S_18S_eLSA/Plots/04.20.2022_abio_brad_pres_abs/ektrx_2_plot.pdf", height=4, width=3.5, units="in", dpi=300)

#ektry (E-W component of Ekman transport)
ektry_table=env_subs %>% group_by(A1_pres) %>% summarize(avg_ektry=mean(as.numeric(ektry)),
                                                         SD=sd(as.numeric(ektry)), n=n(),
                                                         SE=sd(as.numeric(ektry))/sqrt(n()))
ektry_plot=ggplot(ektry_table, aes(x=A1_pres, y=avg_ektry)) + 
  geom_point(size=8, color="#999999") + 
  geom_errorbar(aes(ymin=avg_ektry-SE, ymax=avg_ektry+SE), width=0.2, size=1.2) + labs(title="UCYN-A1 vs. ektry", x="", y="E-W component of Ekman \n transport (kg/m/s)")
ektry_plot + theme_classic() + theme(axis.title.y=element_text(size=9), axis.text=element_text(size=7, face="bold"))

ggsave("/Users/colettef/Desktop/DataAnalysis/SPOT_16S_18S_eLSA/Plots/04.20.2022_abio_brad_pres_abs/formatted/ektry_1_plot.pdf", height=2, width=1.5, units="in", dpi=300)
ggsave("/Users/cfletcherh1/Desktop/DataAnalysis/SPOT_16S_18S_eLSA/Plots/04.20.2022_abio_brad_pres_abs/ektry_1_plot.pdf", height=4, width=3.5, units="in", dpi=300)

#A2
#ektry (E-W component of Ekman transport)
ektry_table=env_subs %>% group_by(A2_pres) %>% summarize(avg_ektry=mean(as.numeric(ektry)),
                                                         SD=sd(as.numeric(ektry)), n=n(),
                                                         SE=sd(as.numeric(ektry))/sqrt(n()))
ektry_plot=ggplot(ektry_table, aes(x=A2_pres, y=avg_ektry)) + 
  geom_point(size=10, color="#999999") + 
  geom_errorbar(aes(ymin=avg_ektry-SE, ymax=avg_ektry+SE), width=0.2) + labs(title="UCYN-A2 vs. ektry", x="", y="E-W component of Ekman \n transport (kg/m/s)")
ektry_plot + theme_classic() + theme(axis.title.y=element_text(size=15), axis.text=element_text(size=15, face="bold"))

ggsave("/Users/cfletcherh1/Desktop/DataAnalysis/SPOT_16S_18S_eLSA/Plots/04.20.2022_abio_brad_pres_abs/ektry_2_plot.pdf", height=4, width=3.5, units="in", dpi=300)

#3. Try logistic regression to see if env variables different on dates when UCYN-A present/absent -----
names(env_subs)
#Add in dummy variables for when UCYN-A1 and UCYN-A2 are present vs absent (0/1)
env_subs$A1_dummy=0
head(env_subs$A1_dummy)
env_subs$A1_dummy[which(env_subs$A1_pres=="Present")]=1
mean(env_subs$A1_dummy) #0.426

#A2
env_subs$A2_dummy=0
head(env_subs$A2_dummy)
which(env_subs$A2_pres=="Present")
env_subs$A2_dummy[which(env_subs$A2_pres=="Present")]=1
mean(env_subs$A2_dummy) #0.385

#Add in dummy variables for when Brad ASV7x and Brad ASV3 are present vs. absent (0/1)
env_subs$Brad7_dummy=0
env_subs$Brad7_dummy[which(env_subs$Brad7_pres=="Present")]=1
mean(env_subs$Brad7_dummy) #0.361

#Brad3
env_subs$Brad3_dummy=0
env_subs$Brad3_dummy[which(env_subs$Brad3_pres=="Present")]=1
mean(env_subs$Brad3_dummy) #0.286

#Now for lepidodinium
env_subs$Lepido_dummy=0
env_subs$Lepido_dummy[which(env_subs$Lepido_pres=="Present")]=1
mean(env_subs$Lepido_dummy) #0.336

#Define functions
A1_pval <- function(variable){
  A1_model=glm(formula=A1_dummy~variable, data=env_subs, family=binomial)
  A1_pval=summary(A1_model)$coef[2,4]
  return(A1_pval)
}

A1_coeficient <- function(variable){
  A1_model=glm(formula=A1_dummy~variable, data=env_subs, family=binomial)
  A1_coeficient=summary(A1_model)$coef[2,1]
  return(A1_coeficient)
}

A2_pval <- function(variable){
  A2_model=glm(formula=A2_dummy~variable, data=env_subs, family=binomial)
  A2_pval=summary(A2_model)$coef[2,4]
  return(A2_pval)
}

A2_coeficient <- function(variable){
  A2_model=glm(formula=A2_dummy~variable, data=env_subs, family=binomial)
  A2_coeficient=summary(A2_model)$coef[2,1]
  return(A2_coeficient)
}

Brad7_pval <- function(variable){
  Brad7_model=glm(formula=Brad7_dummy~variable, data=env_subs, family=binomial)
  Brad7_pval=summary(Brad7_model)$coef[2,4]
  return(Brad7_pval)
}

Brad7_coefficient <- function(variable){
  Brad7_model=glm(formula=Brad7_dummy~variable, data=env_subs, family=binomial)
  Brad7_coefficient=summary(Brad7_model)$coef[2,1]
  return(Brad7_coefficient)
}

Brad3_pval <- function(variable){
  Brad3_model=glm(formula=Brad3_dummy~variable, data=env_subs, family=binomial)
  Brad3_pval=summary(Brad3_model)$coef[2,4]
  return(Brad3_pval)
}

Brad3_coefficient <- function(variable){
  Brad3_model=glm(formula=Brad3_dummy~variable, data=env_subs, family=binomial)
  Brad3_coefficient=summary(Brad3_model)$coef[2,1]
  return(Brad3_coefficient)
}

L_pval <- function(variable){
  L_model=glm(formula=Lepido_dummy~variable, data=env_subs, family=binomial)
  L_pval=summary(L_model)$coef[2,4]
  return(L_pval)
}

L_coefficient <- function(variable){
  L_model=glm(formula=Lepido_dummy~variable, data=env_subs, family=binomial)
  L_coef=summary(L_model)$coef[2,1]
  return(L_coef)
}

#Put the important terms in a vector
A1_Pr=c()
A1_coef=c()
A2_Pr=c()
A2_coef=c()

Brad7_Pr=c()
Brad7_coef=c()
Brad3_Pr=c()
Brad3_coef=c()

L_Pr=c()
L_coef=c()

#Cycle through the variables 
names(env_subs)

rownames=c("CUTI_ix", "BEUTI_ix", "MEI_ix",  "pmsl", "u_mean", "v_mean", "UV_mag_mean", "taux_mean", "tauy_mean", "curl", "ektrx", "ektry", "MODIS_SST", "MODIS_Chl", "Model_NO2", "Model_SST", "NO2+NO3", "PO4", "BactProd_Leu", "BactProd_Thy")

#Remember, Model NO2 and SST won't model properly because they are missing the first row
GLM=glm(formula=as.numeric(A1_dummy[2:nrow(env_subs)])~as.numeric(Model_NO2[2:nrow(env_subs)]), data=env_subs[2:nrow(env_subs),], family=binomial)
no2_P=summary(GLM)$coef[2,4]
no2_coef=summary(GLM)$coef[2,1]

GLM=glm(formula=as.numeric(A1_dummy[2:nrow(env_subs)])~as.numeric(Model_SST[2:nrow(env_subs)]), data=env_subs[2:nrow(env_subs),], family=binomial)
sst_P=summary(GLM)$coef[2,4]
sst_coef=summary(GLM)$coef[2,1]

A1_Pr=c(A1_pval(env_subs$CUTI_ix), A1_pval(env_subs$BEUTI_ix), A1_pval(env_subs$MEI_ix), A1_pval(as.numeric(env_subs$pmsl)), A1_pval(as.numeric(env_subs$u_mean)), A1_pval(as.numeric(env_subs$v_mean)), A1_pval(as.numeric(env_subs$uv_mag_mean)), A1_pval(as.numeric(env_subs$taux_mean)), A1_pval(as.numeric(env_subs$tauy_mean)), A1_pval(as.numeric(env_subs$curl)), A1_pval(as.numeric(env_subs$ektrx)), A1_pval(as.numeric(env_subs$ektry)), A1_pval(env_subs$MODIS_SST), A1_pval(as.numeric(env_subs$MODIS_Chl)), no2_P, sst_P, A1_pval(env_subs$NO2_NO3), A1_pval(env_subs$PO4), A1_pval(env_subs$BactProd_Leu), A1_pval(env_subs$BactProd_Thy))

A1_coef=c(A1_coeficient(env_subs$CUTI_ix), A1_coeficient(env_subs$BEUTI_ix), A1_coeficient(env_subs$MEI_ix), A1_coeficient(as.numeric(env_subs$pmsl)), A1_coeficient(as.numeric(env_subs$u_mean)), A1_coeficient(as.numeric(env_subs$v_mean)), A1_coeficient(as.numeric(env_subs$uv_mag_mean)), A1_coeficient(as.numeric(env_subs$taux_mean)), A1_coeficient(as.numeric(env_subs$tauy_mean)), A1_coeficient(as.numeric(env_subs$curl)), A1_coeficient(as.numeric(env_subs$ektrx)), A1_coeficient(as.numeric(env_subs$ektry)), A1_coeficient(env_subs$MODIS_SST), A1_coeficient(as.numeric(env_subs$MODIS_Chl)), no2_coef, sst_coef, A1_coeficient(env_subs$NO2_NO3), A1_coeficient(env_subs$PO4), A1_coeficient(env_subs$BactProd_Leu), A1_coeficient(env_subs$BactProd_Thy))

GLM=glm(formula=as.numeric(A2_dummy[2:nrow(env_subs)])~as.numeric(Model_NO2[2:nrow(env_subs)]), data=env_subs[2:nrow(env_subs),], family=binomial)
no2_P2=summary(GLM)$coef[2,4]
no2_coef2=summary(GLM)$coef[2,1]

GLM=glm(formula=as.numeric(A2_dummy[2:nrow(env_subs)])~as.numeric(Model_SST[2:nrow(env_subs)]), data=env_subs[2:nrow(env_subs),], family=binomial)
sst_P2=summary(GLM)$coef[2,4]
sst_coef2=summary(GLM)$coef[2,1]

A2_Pr=c(A2_pval(env_subs$CUTI_ix), A2_pval(env_subs$BEUTI_ix), A2_pval(env_subs$MEI_ix), A2_pval(as.numeric(env_subs$pmsl)), A2_pval(as.numeric(env_subs$u_mean)), A2_pval(as.numeric(env_subs$v_mean)), A2_pval(as.numeric(env_subs$uv_mag_mean)), A2_pval(as.numeric(env_subs$taux_mean)), A2_pval(as.numeric(env_subs$tauy_mean)), A2_pval(as.numeric(env_subs$curl)), A2_pval(as.numeric(env_subs$ektrx)), A2_pval(as.numeric(env_subs$ektry)), A2_pval(env_subs$MODIS_SST), A2_pval(as.numeric(env_subs$MODIS_Chl)), no2_P2, sst_P2, A2_pval(env_subs$NO2_NO3), A2_pval(env_subs$PO4), A2_pval(env_subs$BactProd_Leu), A2_pval(env_subs$BactProd_Thy))

A2_coef=c(A2_coeficient(env_subs$CUTI_ix), A2_coeficient(env_subs$BEUTI_ix), A2_coeficient(env_subs$MEI_ix), A2_coeficient(as.numeric(env_subs$pmsl)), A2_coeficient(as.numeric(env_subs$u_mean)), A2_coeficient(as.numeric(env_subs$v_mean)), A2_coeficient(as.numeric(env_subs$uv_mag_mean)), A2_coeficient(as.numeric(env_subs$taux_mean)), A2_coeficient(as.numeric(env_subs$tauy_mean)), A2_coeficient(as.numeric(env_subs$curl)), A2_coeficient(as.numeric(env_subs$ektrx)), A2_coeficient(as.numeric(env_subs$ektry)), A2_coeficient(env_subs$MODIS_SST), A2_coeficient(as.numeric(env_subs$MODIS_Chl)), no2_coef2, sst_coef2, A2_coeficient(env_subs$NO2_NO3), A2_coeficient(env_subs$PO4), A2_coeficient(env_subs$BactProd_Leu), A2_coeficient(env_subs$BactProd_Thy))

#Braarudospharea's 
GLM=glm(formula=as.numeric(Brad7_dummy[2:nrow(env_subs)])~as.numeric(Model_NO2[2:nrow(env_subs)]), data=env_subs[2:nrow(env_subs),], family=binomial)
no2_P=summary(GLM)$coef[2,4]
no2_coef=summary(GLM)$coef[2,1]

GLM=glm(formula=as.numeric(Brad7_dummy[2:nrow(env_subs)])~as.numeric(Model_SST[2:nrow(env_subs)]), data=env_subs[2:nrow(env_subs),], family=binomial)
sst_P=summary(GLM)$coef[2,4]
sst_coef=summary(GLM)$coef[2,1]

Brad7_Pr=c(Brad7_pval(env_subs$CUTI_ix), Brad7_pval(env_subs$BEUTI_ix), Brad7_pval(env_subs$MEI_ix), Brad7_pval(as.numeric(env_subs$pmsl)), Brad7_pval(as.numeric(env_subs$u_mean)), Brad7_pval(as.numeric(env_subs$v_mean)), Brad7_pval(as.numeric(env_subs$uv_mag_mean)), Brad7_pval(as.numeric(env_subs$taux_mean)), Brad7_pval(as.numeric(env_subs$tauy_mean)), Brad7_pval(as.numeric(env_subs$curl)), Brad7_pval(as.numeric(env_subs$ektrx)), Brad7_pval(as.numeric(env_subs$ektry)), Brad7_pval(env_subs$MODIS_SST), Brad7_pval(as.numeric(env_subs$MODIS_Chl)), no2_P, sst_P, Brad7_pval(env_subs$NO2_NO3), Brad7_pval(env_subs$PO4), Brad7_pval(env_subs$BactProd_Leu), Brad7_pval(env_subs$BactProd_Thy))

Brad7_coef=c(Brad7_coefficient(env_subs$CUTI_ix), Brad7_coefficient(env_subs$BEUTI_ix), Brad7_coefficient(env_subs$MEI_ix), Brad7_coefficient(as.numeric(env_subs$pmsl)), Brad7_coefficient(as.numeric(env_subs$u_mean)), Brad7_coefficient(as.numeric(env_subs$v_mean)), Brad7_coefficient(as.numeric(env_subs$uv_mag_mean)), Brad7_coefficient(as.numeric(env_subs$taux_mean)), Brad7_coefficient(as.numeric(env_subs$tauy_mean)), Brad7_coefficient(as.numeric(env_subs$curl)), Brad7_coefficient(as.numeric(env_subs$ektrx)), Brad7_coefficient(as.numeric(env_subs$ektry)), Brad7_coefficient(env_subs$MODIS_SST), Brad7_coefficient(as.numeric(env_subs$MODIS_Chl)), no2_coef, sst_coef, Brad7_coefficient(env_subs$NO2_NO3), Brad7_coefficient(env_subs$PO4), Brad7_coefficient(env_subs$BactProd_Leu), Brad7_coefficient(env_subs$BactProd_Thy))

GLM=glm(formula=as.numeric(Brad3_dummy[2:nrow(env_subs)])~as.numeric(Model_NO2[2:nrow(env_subs)]), data=env_subs[2:nrow(env_subs),], family=binomial)
no2_P=summary(GLM)$coef[2,4]
no2_coef=summary(GLM)$coef[2,1]

GLM=glm(formula=as.numeric(Brad3_dummy[2:nrow(env_subs)])~as.numeric(Model_SST[2:nrow(env_subs)]), data=env_subs[2:nrow(env_subs),], family=binomial)
sst_P=summary(GLM)$coef[2,4]
sst_coef=summary(GLM)$coef[2,1]

Brad3_Pr=c(Brad3_pval(env_subs$CUTI_ix), Brad3_pval(env_subs$BEUTI_ix), Brad3_pval(env_subs$MEI_ix), Brad3_pval(as.numeric(env_subs$pmsl)), Brad3_pval(as.numeric(env_subs$u_mean)), Brad3_pval(as.numeric(env_subs$v_mean)), Brad3_pval(as.numeric(env_subs$uv_mag_mean)), Brad3_pval(as.numeric(env_subs$taux_mean)), Brad3_pval(as.numeric(env_subs$tauy_mean)), Brad3_pval(as.numeric(env_subs$curl)), Brad3_pval(as.numeric(env_subs$ektrx)), Brad3_pval(as.numeric(env_subs$ektry)), Brad3_pval(env_subs$MODIS_SST), Brad3_pval(as.numeric(env_subs$MODIS_Chl)), no2_P, sst_P, Brad3_pval(env_subs$NO2_NO3), Brad3_pval(env_subs$PO4), Brad3_pval(env_subs$BactProd_Leu), Brad3_pval(env_subs$BactProd_Thy))

Brad3_coef=c(Brad3_coefficient(env_subs$CUTI_ix), Brad3_coefficient(env_subs$BEUTI_ix), Brad3_coefficient(env_subs$MEI_ix), Brad3_coefficient(as.numeric(env_subs$pmsl)), Brad3_coefficient(as.numeric(env_subs$u_mean)), Brad3_coefficient(as.numeric(env_subs$v_mean)), Brad3_coefficient(as.numeric(env_subs$uv_mag_mean)), Brad3_coefficient(as.numeric(env_subs$taux_mean)), Brad3_coefficient(as.numeric(env_subs$tauy_mean)), Brad3_coefficient(as.numeric(env_subs$curl)), Brad3_coefficient(as.numeric(env_subs$ektrx)), Brad3_coefficient(as.numeric(env_subs$ektry)), Brad3_coefficient(env_subs$MODIS_SST), Brad3_coefficient(as.numeric(env_subs$MODIS_Chl)), no2_coef, sst_coef, Brad3_coefficient(env_subs$NO2_NO3), Brad3_coefficient(env_subs$PO4), Brad3_coefficient(env_subs$BactProd_Leu), Brad3_coefficient(env_subs$BactProd_Thy))

#Lepidodinium
GLM=glm(formula=as.numeric(Lepido_dummy[2:nrow(env_subs)])~as.numeric(Model_NO2[2:nrow(env_subs)]), data=env_subs[2:nrow(env_subs),], family=binomial)
no2_P=summary(GLM)$coef[2,4]
no2_coef=summary(GLM)$coef[2,1]

GLM=glm(formula=as.numeric(Lepido_dummy[2:nrow(env_subs)])~as.numeric(Model_SST[2:nrow(env_subs)]), data=env_subs[2:nrow(env_subs),], family=binomial)
sst_P=summary(GLM)$coef[2,4]
sst_coef=summary(GLM)$coef[2,1]

L_Pr=c(L_pval(env_subs$CUTI_ix),L_pval(env_subs$BEUTI_ix), L_pval(env_subs$MEI_ix), L_pval(as.numeric(env_subs$pmsl)), L_pval(as.numeric(env_subs$u_mean)), L_pval(as.numeric(env_subs$v_mean)), L_pval(as.numeric(env_subs$uv_mag_mean)), L_pval(as.numeric(env_subs$taux_mean)), L_pval(as.numeric(env_subs$tauy_mean)), L_pval(as.numeric(env_subs$curl)), L_pval(as.numeric(env_subs$ektrx)), L_pval(as.numeric(env_subs$ektry)), L_pval(env_subs$MODIS_SST), L_pval(as.numeric(env_subs$MODIS_Chl)), no2_P, sst_P, L_pval(env_subs$NO2_NO3), L_pval(env_subs$PO4), L_pval(env_subs$BactProd_Leu), L_pval(env_subs$BactProd_Thy))

L_coef=c(L_coefficient(env_subs$CUTI_ix), L_coefficient(env_subs$BEUTI_ix), L_coefficient(env_subs$MEI_ix), L_coefficient(as.numeric(env_subs$pmsl)), L_coefficient(as.numeric(env_subs$u_mean)), L_coefficient(as.numeric(env_subs$v_mean)), L_coefficient(as.numeric(env_subs$uv_mag_mean)), L_coefficient(as.numeric(env_subs$taux_mean)), L_coefficient(as.numeric(env_subs$tauy_mean)), L_coefficient(as.numeric(env_subs$curl)), L_coefficient(as.numeric(env_subs$ektrx)), L_coefficient(as.numeric(env_subs$ektry)), L_coefficient(env_subs$MODIS_SST), L_coefficient(as.numeric(env_subs$MODIS_Chl)), no2_coef, sst_coef, L_coefficient(env_subs$NO2_NO3), L_coefficient(env_subs$PO4), L_coefficient(env_subs$BactProd_Leu), L_coefficient(env_subs$BactProd_Thy))

#Put together in a table

length(rownames) #20
length(A1_Pr) #20
length(A1_coef) #20
length(A2_Pr) #20
length(A2_coef) #20
length(L_Pr) #20
length(L_coef) #20

#Put this into a table
p_values_table=as.data.frame(matrix(nrow=length(rownames), ncol=(4*4)))
rownames(p_values_table)=rownames
colnames(p_values_table)=c("A1_Pr", "A1_Pr_adj", "A1_coef", "A1_inc_dec", "A2_Pr", "A2_Pr_adj", "A2_coef", "A2_inc_dec", "Brad7_Pr", "Brad7_Pr_adj", "Brad7_coef", "Brad7_inc_dec", "Brad3_Pr", "Brad3_Pr_adj", "Brad3_coef", "Brad3_inc_dec")
dim(p_values_table)
rownames(p_values_table)
names(p_values_table)

#Also add in bonferroni correction, increase/decrease column
p_values_table$A1_Pr=A1_Pr
p_values_table$A1_Pr_adj=p.adjust(A1_Pr, method="BH")
p_values_table$A1_coef=A1_coef
p_values_table$A1_inc_dec="Increase"
p_values_table$A1_inc_dec[which(p_values_table$A1_coef<0)]="Decrease"

p_values_table$A2_Pr=A2_Pr
p_values_table$A2_Pr_adj=p.adjust(A2_Pr, method="BH")
p_values_table$A2_coef=A2_coef
p_values_table$A2_inc_dec="Increase"
p_values_table$A2_inc_dec[which(p_values_table$A2_coef<0)]="Decrease"

p_values_table$Brad7_Pr=Brad7_Pr
p_values_table$Brad7_Pr_adj=p.adjust(Brad7_Pr, method="BH")
p_values_table$Brad7_coef=Brad7_coef
p_values_table$Brad7_inc_dec="Increase"
p_values_table$Brad7_inc_dec[which(p_values_table$Brad7_coef<0)]="Decrease"

p_values_table$Brad3_Pr=Brad3_Pr
p_values_table$Brad3_Pr_adj=p.adjust(Brad3_Pr, method="BH")
p_values_table$Brad3_coef=Brad3_coef
p_values_table$Brad3_inc_dec="Increase"
p_values_table$Brad3_inc_dec[which(p_values_table$Brad3_coef<0)]="Decrease"

p_values_table$Lepido_Pr=L_Pr
p_values_table$Lepido_Pr_adj=p.adjust(L_Pr, method="BH")
p_values_table$Lepido_coef=L_coef
p_values_table$Lepido_inc_dec="Increase"
p_values_table$Lepido_inc_dec[which(p_values_table$Lepido_coef<0)]="Decrease"

#Write it out!
write.table(p_values_table, file="UCYN-A_Brad_pres_absence_vs_env_variables_Pvals_06.27.2022.tsv", row.names = T, quote=F, sep="\t")

#4. Put these data in a table----
#Read env_subs data back in, this variable is being weird
env_subs=read.table("int_CLR_data/int_env_data_UCYNA_Brad_Lepido_pres_abs_cutoff_0.0001_06.27.2022.tsv",header=T, stringsAsFactors = F, sep="\t")
dim(env_subs) #122 31
names(env_subs)

#Drop model data
env_subs=env_subs[,-grep("Model", names(env_subs))]

#Drop "Sample ID"
env_subs=env_subs[,-grep("SampleID", names(env_subs))]
dim(env_subs) #122 30

for(i in c(2:grep("NO2_NO3", names(env_subs)))){
  print(names(env_subs)[i])
  print(which(is.na(env_subs[,i])))
  print(mean(env_subs[,i]))
}

names(env_subs)[2:grep("NO2_NO3", names(env_subs))]

table_A1_mean= env_subs %>% group_by(A1_pres) %>% summarize_at(vars(names(env_subs)[2:grep("NO2_NO3", names(env_subs))]), mean)
table_A1_sd=env_subs %>% group_by(A1_pres) %>% summarize_at(vars(names(env_subs)[2:grep("NO2_NO3", names(env_subs))]), sd)
length(which(env_subs$A1_pres=="Present")) #50
length(which(env_subs$A1_pres=="Absent")) #72
#A1 is present on 50 dates, and absent on 72 dates --divide SD by sqrt(n) to get SE

table_B7_mean= env_subs %>% group_by(Brad7_pres) %>% summarize_at(vars(names(env_subs)[2:grep("NO2_NO3", names(env_subs))]), mean)
table_B7_sd=env_subs %>% group_by(Brad7_pres) %>% summarize_at(vars(names(env_subs)[2:grep("NO2_NO3", names(env_subs))]), sd)
length(which(env_subs$Brad7_pres=="Present")) #44
length(which(env_subs$Brad7_pres=="Absent")) #78

table_A2_mean= env_subs %>% group_by(A2_pres) %>% summarize_at(vars(names(env_subs)[2:grep("NO2_NO3", names(env_subs))]), mean)
table_A2_sd=env_subs %>% group_by(A2_pres) %>% summarize_at(vars(names(env_subs)[2:grep("NO2_NO3", names(env_subs))]), sd)
length(which(env_subs$A2_pres=="Present")) #47
length(which(env_subs$A2_pres=="Absent")) #75

table_B3_mean= env_subs %>% group_by(Brad3_pres) %>% summarize_at(vars(names(env_subs)[2:grep("NO2_NO3", names(env_subs))]), mean)
table_B3_sd=env_subs %>% group_by(Brad3_pres) %>% summarize_at(vars(names(env_subs)[2:grep("NO2_NO3", names(env_subs))]), sd)
length(which(env_subs$Brad3_pres=="Present")) #35
length(which(env_subs$Brad3_pres=="Absent")) #87

table_L_mean= env_subs %>% group_by(Lepido_pres) %>% summarize_at(vars(names(env_subs)[2:grep("NO2_NO3", names(env_subs))]), mean)
table_L_sd=env_subs %>% group_by(Lepido_pres) %>% summarize_at(vars(names(env_subs[2:grep("NO2_NO3", names(env_subs))])), sd)
length(which(env_subs$Lepido_pres=="Present")) #41
length(which(env_subs$Lepido_pres=="Absent")) #81

TABLE=as.data.frame(cbind(t(table_A1_mean), t(table_A1_sd), t(table_B7_mean), t(table_B7_sd), t(table_A2_mean), t(table_A2_sd), t(table_B3_mean), t(table_B3_sd), t(table_L_mean), t(table_L_sd)))
names(TABLE)=c("A1_Absent_mean", "A1_Present_mean", "A1_Absent_SD", "A1_Present_SD",
               "B7_Absent_mean", "B7_Present_mean", "B7_Absent_SD", "B7_Present_SD",
               "A2_Absent_mean", "A2_Present_mean", "A2_Absent_SD", "A2_Present_SD", 
               "B3_Absent_mean", "B3_Present_mean", "B3_Absent_SD", "B3_Present_SD",
               "L_Absent_mean", "L_Present_mean", "L_Absent_SD", "L_Present_SD")
dim(TABLE) #13 20
head(TABLE)

#Get rid of that first row 
TABLE=TABLE[-1,]
dim(TABLE) #12 20

#Calculate standard error
TABLE$A1_absent_SE=as.numeric(TABLE$A1_Absent_SD)/sqrt(length(which(env_subs$A1_pres=="Absent")))
TABLE$A1_present_SE=as.numeric(TABLE$A1_Present_SD)/sqrt(length(which(env_subs$A1_pres=="Present")))
TABLE$B7_absent_SE=as.numeric(TABLE$B7_Absent_SD)/sqrt(length(which(env_subs$Brad7_pres=="Absent")))
TABLE$B7_present_SE=as.numeric(TABLE$B7_Present_SD)/sqrt(length(which(env_subs$Brad7_pres=="Present")))
TABLE$A2_absent_SE=as.numeric(TABLE$A2_Absent_SD)/sqrt(length(which(env_subs$A2_pres=="Absent")))
TABLE$A2_present_SE=as.numeric(TABLE$A2_Present_SD)/sqrt(length(which(env_subs$A2_pres=="Present")))
TABLE$B3_absent_SE=as.numeric(TABLE$B3_Present_SD)/sqrt(length(which(env_subs$Brad3_pres=="Absent")))
TABLE$B3_present_SE=as.numeric(TABLE$B3_Absent_SD)/sqrt(length(which(env_subs$Brad3_pres=="Present")))
TABLE$L_absent_SE=as.numeric(TABLE$L_Absent_SD)/sqrt(length(which(env_subs$Lepido_pres=="Absent")))
TABLE$L_present_SE=as.numeric(TABLE$L_Present_SD)/sqrt(length(which(env_subs$Lepido_pres=="Present")))

dim(TABLE) #22 30

#order
ORDER=c("A1_Absent_mean", "A1_Present_mean", "A1_Absent_SD", "A1_Present_SD", "A1_absent_SE", "A1_present_SE",
"B7_Absent_mean", "B7_Present_mean", "B7_Absent_SD", "B7_Present_SD", "B7_absent_SE", "B7_present_SE",
"A2_Absent_mean", "A2_Present_mean", "A2_Absent_SD", "A2_Present_SD", "A2_absent_SE", "A2_present_SE", 
"B3_Absent_mean", "B3_Present_mean", "B3_Absent_SD", "B3_Present_SD", "A2_absent_SE", "A2_present_SE", 
"L_Absent_mean", "L_Present_mean", "L_Absent_SD", "L_Present_SD", "L_absent_SE", "L_present_SE")

TABLE2=select(TABLE, ORDER)

write.table(TABLE2, quote=F, file="UCYNA_hosts_Lepido_pres_abs_abiotic_data_mean_SD+SE_07.05.2022.tsv", sep="\t")


#5. Plot upwelling indices based on relative abundance of UCYN-A (absent, low, high)-----
#Turk Kubo et al. observed upwelling stimulated UCYN-A relative abundance in May 2017 at SIO pier
#SST was 16.82, No2+no3 was 0.1uM, PO4 was 0.01 uM, chl was 0.213 #So check these variables and upwelling indices too

#UCYN-A1 N species
N_table = env_subs %>% group_by(A1_abund) %>% summarize(avgN=mean(NO2_NO3),
                                                       SD=sd(NO2_NO3), n=n(), 
                                                       SE=sd(NO2_NO3/sqrt(n())))
#Try plotting just the mean, using a circle
N_plot=ggplot(N_table, aes(x=A1_abund, y=avgN)) + 
  geom_point(size=15, color="#0099f9") +
  geom_errorbar(aes(ymin=avgN-SE, ymax=avgN+SE), width=0.2) + 
  labs(title="UCYN-A1 vs NO2+NO3", x="", y="[NO2+NO3] (uM)", size=18) 
N_plot + theme_classic() + theme(axis.title.y=element_text(size=20, face="bold", color="black"), axis.text=element_text(size=20, face="bold", color="black"))

ggsave("/Users/cfletcherh1/Desktop/DataAnalysis/SPOT_16S_18S_eLSA/Plots/04.20.2022_abio_brad_pres_abs/abund_N_1_plot.pdf", height = 4, width = 3.5, units = "in", dpi=300)

#UCYN-A2 N species
N_table = env_subs %>% group_by(A2_abund) %>% summarize(avgN=mean(NO2_NO3),
                                                       SD=sd(NO2_NO3), n=n(), 
                                                       SE=sd(NO2_NO3/sqrt(n())))
N_table

#Just the mean
N_plot=ggplot(N_table, aes(x=A2_abund, y=avgN)) + 
  geom_point(size=15, color="#6699CC") +
  geom_errorbar(aes(ymin=avgN-SE, ymax=avgN+SE), width=0.2) + 
  labs(title="UCYN-A1 vs NO2+NO3", x="", y="[NO2+NO3] (uM)", size=18) 
N_plot  + theme_classic() + theme(axis.title.y=element_text(size=20, face="bold", color="black"), axis.text=element_text(size=20, face="bold", color="black"))
#No statistical difference 

ggsave("/Users/cfletcherh1/Desktop/DataAnalysis/SPOT_16S_18S_eLSA/Plots/04.20.2022_abio_brad_pres_abs/abund_N_2_plot.pdf", height = 4, width = 3.50, units = "in", dpi=300)

#UCYN-A1 PO4
P_table = env_subs %>% group_by(A1_abund) %>% summarize(avgP=mean(PO4),
                                                       SD=sd(PO4), n=n(), 
                                                       SE=sd(PO4/sqrt(n())))
#Just the mean #A1 PO4
P_plot=ggplot(P_table, aes(x=A1_abund, y=avgP)) + 
  geom_point(size=8, color="#990099") +
  geom_errorbar(aes(ymin=avgP-SE, ymax=avgP+SE), width=0.2, size=1.2) + 
  labs(title="UCYN-A1 vs PO4", x="", y="[PO4] (uM)", size=18) +
  theme_classic() + 
  theme(axis.title.y=element_text(size=9, face="bold", color="black"), axis.text=element_text(size=7, face="bold", color="black"))  
P_plot  

ggsave("/Users/colettef/Desktop/DataAnalysis/SPOT_16S_18S_eLSA/Plots/04.20.2022_abio_brad_pres_abs/formatted/abund_PO4_1_plot.pdf", height = 2, width = 1.8, units = "in", dpi=300)
ggsave("/Users/cfletcherh1/Desktop/DataAnalysis/SPOT_16S_18S_eLSA/Plots/04.20.2022_abio_brad_pres_abs/abund_P_1_plot.pdf", height = 4, width = 3.50, units = "in", dpi=300)

#UCYN-A2 P
P_table = env_subs %>% group_by(A2_abund) %>% summarize(avgP=mean(PO4),
                                                       SD=sd(PO4), n=n(), 
                                                       SE=sd(PO4/sqrt(n())))
#Just the mean #A2 PO4
P_plot=ggplot(P_table, aes(x=A2_abund, y=avgP)) + 
  geom_point(size=8, color="#CC99CC") +
  geom_errorbar(aes(ymin=avgP-SE, ymax=avgP+SE), width=0.2, size=1.2) + 
  labs(title="UCYN-A2 vs PO4", x="", y="[PO4] (uM)", size=18) +
  theme_classic() + 
  theme(axis.title.y=element_text(size=9, face="bold", color="black"), axis.text=element_text(size=7, face="bold", color="black"))
P_plot

ggsave("/Users/colettef/Desktop/DataAnalysis/SPOT_16S_18S_eLSA/Plots/04.20.2022_abio_brad_pres_abs/formatted/abund_PO4_2_plot.pdf", height = 2, width = 1.8, units = "in", dpi=300)
ggsave("/Users/cfletcherh1/Desktop/DataAnalysis/SPOT_16S_18S_eLSA/Plots/04.20.2022_abio_brad_pres_abs/abund_P_2_plot.pdf", height = 4, width = 3.50, units = "in", dpi=300)

#Upwelling indices
CUTI_table = env_subs %>% group_by(A1_abund) %>% summarize(avgCUTI=mean(CUTI_ix),
                                                          SD=sd(CUTI_ix), n=n(), 
                                                          SE=sd(CUTI_ix/sqrt(n())))

#Just the mean #A1
CUTI_plot=ggplot(CUTI_table, aes(x=A1_abund, y=avgCUTI)) + 
  geom_point(size=8, color="#CC6600") +
  geom_errorbar(aes(ymin=avgCUTI-SE, ymax=avgCUTI+SE), width=0.2, size=1.2) + 
  labs(title="UCYN-A1 vs CUTI", x="", y="CUTI", size=18) + theme_classic() + 
  theme(axis.title.y=element_text(size=9, face="bold", color="black"), axis.text=element_text(size=7, face="bold", color="black"))
CUTI_plot 

ggsave("/Users/colettef/Desktop/DataAnalysis/SPOT_16S_18S_eLSA/Plots/04.20.2022_abio_brad_pres_abs/formatted/abund_CUTI_1_plot.pdf", height = 2, width = 1.8, units = "in", dpi=300)
ggsave("/Users/cfletcherh1/Desktop/DataAnalysis/SPOT_16S_18S_eLSA/Plots/04.20.2022_abio_brad_pres_abs/abund_CUTI_1_plot.pdf", height = 4, width = 3.50, units = "in", dpi=300)

#A2 vs CUTI avg
CUTI_table = env_subs %>% group_by(A2_abund) %>% summarize(avgCUTI=mean(CUTI_ix),
                                    SD=sd(CUTI_ix), n=n(), 
                                    SE=sd(CUTI_ix/sqrt(n())))
CUTI_table

CUTI_plot=ggplot(CUTI_table, aes(x=A2_abund, y=avgCUTI)) + 
  geom_point(size=8, color="#FFCC99") +
  geom_errorbar(aes(ymin=avgCUTI-SE, ymax=avgCUTI+SE), width=0.2, size=1.2) + 
  labs(title="UCYN-A2 vs CUTI", x="", y="CUTI", size=18) +
  theme_classic() + 
  theme(axis.title.y=element_text(size=9, face="bold", color="black"), axis.text=element_text(size=7, face="bold", color="black"))
CUTI_plot

ggsave("/Users/colettef/Desktop/DataAnalysis/SPOT_16S_18S_eLSA/Plots/04.20.2022_abio_brad_pres_abs/formatted/abund_CUTI_2_plot.pdf", height = 2, width = 1.8, units = "in", dpi=300)
ggsave("/Users/cfletcherh1/Desktop/DataAnalysis/SPOT_16S_18S_eLSA/Plots/04.20.2022_abio_brad_pres_abs/abund_CUTI_2_plot.pdf", height = 4, width = 3.50, units = "in", dpi=300)

#A1 vs BEUTI index
BEUTI_table = env_subs %>% group_by(A1_abund) %>% summarize(avgBEUTI=mean(BEUTI_ix),
                                                           SD=sd(BEUTI_ix), n=n(), 
                                                           SE=sd(BEUTI_ix/sqrt(n())))
#Plot just the mean
BEUTI_plot=ggplot(BEUTI_table, aes(x=A1_abund, y=avgBEUTI)) + 
  geom_point(size=8, color="#FF00FF") +
  geom_errorbar(aes(ymin=avgBEUTI-SE, ymax=avgBEUTI+SE), width=0.2, size=1.2) + 
  labs(title="UCYN-A1 vs BEUTI", x="", y="BEUTI", size=18) + 
  theme_classic() + 
  theme(axis.title.y=element_text(size=9, face="bold", color="black"), axis.text=element_text(size=7, face="bold", color="black"))
BEUTI_plot

ggsave("/Users/colettef/Desktop/DataAnalysis/SPOT_16S_18S_eLSA/Plots/04.20.2022_abio_brad_pres_abs/formatted/abund_BEUTI_1_plot.pdf", height = 2, width = 1.8, units = "in", dpi=300)
ggsave("/Users/cfletcherh1/Desktop/DataAnalysis/SPOT_16S_18S_eLSA/Plots/04.20.2022_abio_brad_pres_abs/abund_BEUTI_1_plot.pdf", height = 4, width = 3.50, units = "in", dpi=300)

#A2 vs BEUTI index
BEUTI_table = env_subs %>% group_by(A2_abund) %>% summarize(avgBEUTI=mean(BEUTI_ix),
                                                           SD=sd(BEUTI_ix), n=n(), 
                                                           SE=sd(BEUTI_ix/sqrt(n())))

#Plot just the mean #A2
BEUTI_plot=ggplot(BEUTI_table, aes(x=A2_abund, y=avgBEUTI)) + 
  geom_point(size=8, color="#FF99FF") +
  geom_errorbar(aes(ymin=avgBEUTI-SE, ymax=avgBEUTI+SE), width=0.2, size=1.2) + 
  labs(title="UCYN-A2 vs BEUTI", x="", y="BEUTI", size=18) +
  theme_classic() + 
  theme(axis.title.y=element_text(size=9, face="bold", color="black"), axis.text=element_text(size=7, face="bold", color="black"))
BEUTI_plot

ggsave("/Users/colettef/Desktop/DataAnalysis/SPOT_16S_18S_eLSA/Plots/04.20.2022_abio_brad_pres_abs/formatted/abund_BEUTI_2_plot.pdf", height = 2, width = 1.8, units = "in", dpi=300)
ggsave("/Users/cfletcherh1/Desktop/DataAnalysis/SPOT_16S_18S_eLSA/Plots/04.20.2022_abio_brad_pres_abs/abund_BEUTI_2_plot.pdf", height = 4, width = 3.50, units = "in", dpi=300)

#A1 vs MEI
MEI_table = env_subs %>% group_by(A1_abund) %>% summarize(avgMEI=mean(MEI_ix),
                                                         SD=sd(MEI_ix), n=n(), 
                                                         SE=sd(MEI_ix/sqrt(n())))
#Plot just the mean
MEI_plot=ggplot(MEI_table, aes(x=A1_abund, y=avgMEI)) + 
  geom_point(size=8, color="#CC0033") +
  geom_errorbar(aes(ymin=avgMEI-SE, ymax=avgMEI+SE), width=0.2, size=1.2) + 
  labs(title="UCYN-A1 vs MEI", x="", y="MEI") +
  theme_classic() + 
  theme(axis.title.y=element_text(size=9, face="bold", color="black"), axis.text=element_text(size=7, face="bold", color="black"))
MEI_plot

ggsave("/Users/colettef/Desktop/DataAnalysis/SPOT_16S_18S_eLSA/Plots/04.20.2022_abio_brad_pres_abs/formatted/abund_MEI_1_plot.pdf", height = 2, width = 1.8, units = "in", dpi=300)
ggsave("/Users/cfletcherh1/Desktop/DataAnalysis/SPOT_16S_18S_eLSA/Plots/04.20.2022_abio_brad_pres_abs/abund_MEI_1_plot.pdf", height = 4, width = 3.5, units = "in", dpi=300)

#Now for UCYN-A2 vs MEI
MEI_table = env_subs %>% group_by(A2_abund) %>% summarize(avgMEI=mean(MEI_ix),
                                                         SD=sd(MEI_ix), n=n(), 
                                                         SE=sd(MEI_ix/sqrt(n())))

#Plot just the mean #A2
MEI_plot=ggplot(MEI_table, aes(x=A2_abund, y=avgMEI)) + 
  geom_point(size=8, color="#FF9999") +
  geom_errorbar(aes(ymin=avgMEI-SE, ymax=avgMEI+SE), width=0.2, size=1.2) + 
  labs(title="UCYN-A2 vs MEI", x="", y="MEI") +
  theme_classic() + 
  theme(axis.title.y=element_text(size=9, face="bold", color="black"), axis.text=element_text(size=7, face="bold", color="black"))
MEI_plot

ggsave("/Users/colettef/Desktop/DataAnalysis/SPOT_16S_18S_eLSA/Plots/04.20.2022_abio_brad_pres_abs/formatted/abund_MEI_2_plot.pdf", height = 2, width = 1.8, units = "in", dpi=300)
ggsave("/Users/cfletcherh1/Desktop/DataAnalysis/SPOT_16S_18S_eLSA/Plots/04.20.2022_abio_brad_pres_abs/abund_MEI_2_plot.pdf", height = 4, width = 3.5, units = "in", dpi=300)

#MODIS SST #A1
SST_table= env_subs %>% group_by(A1_abund) %>% summarize(avg_SST=mean(MODIS_SST),
                                                        SD=sd(MODIS_SST), n=n(),
                                                        SE=sd(MODIS_SST)/sqrt(n()))
SST_plot=ggplot(SST_table, aes(x=A1_abund, y=avg_SST))+
  geom_point(size=8, color="#999999") + 
  geom_errorbar(aes(ymin=avg_SST-SE, ymax=avg_SST+SE), width=0.2, size=1.2) + 
  labs(title="UCYN-A1 vs MODIS SST", x="", y="SST (ºC)") + 
  theme_classic() + 
  theme(axis.title.y=element_text(size=9, face="bold", color="black"), axis.text=element_text(size=7, face="bold", color="black"))
SST_plot

ggsave("/Users/colettef/Desktop/DataAnalysis/SPOT_16S_18S_eLSA/Plots/04.20.2022_abio_brad_pres_abs/formatted/abund_SST_1_plot.pdf", height = 2, width = 1.8, units = "in", dpi=300)
ggsave("/Users/cfletcherh1/Desktop/DataAnalysis/SPOT_16S_18S_eLSA/Plots/04.20.2022_abio_brad_pres_abs/abund_MODIS_SST_1_plot.pdf", height = 4, width = 3.5, units = "in", dpi=300)

#A2
SST_table= env_subs %>% group_by(A2_abund) %>% summarize(avg_SST=mean(MODIS_SST),
                                                         SD=sd(MODIS_SST), n=n(),
                                                         SE=sd(MODIS_SST)/sqrt(n()))

SST_plot=ggplot(SST_table, aes(x=A2_abund, y=avg_SST))+
  geom_point(size=8, color="#999999") + 
  geom_errorbar(aes(ymin=avg_SST-SE, ymax=avg_SST+SE), width=0.2, size=1.2) + 
  labs(title="UCYN-A2 vs MODIS SST", x="", y="SST (ºC)") + 
  theme_classic() + 
  theme(axis.title.y=element_text(size=9, face="bold", color="black"), axis.text=element_text(size=7, face="bold", color="black"))
SST_plot

ggsave("/Users/colettef/Desktop/DataAnalysis/SPOT_16S_18S_eLSA/Plots/04.20.2022_abio_brad_pres_abs/formatted/abund_SST_2_plot.pdf", height = 2, width = 1.8, units = "in", dpi=300)
ggsave("/Users/cfletcherh1/Desktop/DataAnalysis/SPOT_16S_18S_eLSA/Plots/04.20.2022_abio_brad_pres_abs/abund_MODIS_SST_2_plot.pdf", height = 4, width = 3.5, units = "in", dpi=300)

#A2
SST_table= env_subs %>% group_by(A2_abund) %>% summarize(avg_SST=mean(MODIS_SST),
                                                        SD=sd(MODIS_SST), n=n(),
                                                        SE=sd(MODIS_SST)/sqrt(n()))
SST_plot=ggplot(SST_table, aes(x=A2_abund, y=avg_SST)) +
  geom_point(size=10, color="#999999") + 
  geom_errorbar(aes(ymin=avg_SST-SE, ymax=avg_SST+SE), width=0.2) + labs(title="UCYN-A2 vs MODIS SST", x="", y="SST (ºC)") + 
  theme_classic() + 
  theme(axis.title.y=element_text(size=20, face="bold", color="black"), axis.text=element_text(size=20, face="bold", color="black"))
SST_plot

ggsave("/Users/cfletcherh1/Desktop/DataAnalysis/SPOT_16S_18S_eLSA/Plots/04.20.2022_abio_brad_pres_abs/abund_MODIS_SST_2_plot.pdf", height = 4, width = 3.5, units = "in", dpi=300)

#MODIS Chl #A1
Chl_table=env_subs %>% group_by(A1_abund) %>% summarize(avgChl=mean(as.numeric(MODIS_Chl)),
                                                       SD=sd(as.numeric(MODIS_Chl)), n=n(),
                                                       SE=sd(as.numeric(MODIS_Chl))/sqrt(n()))

Chl_plot=ggplot(Chl_table, aes(x=A1_abund, y=avgChl)) +
  geom_point(size=10, color="#009E73") + 
  geom_errorbar(aes(ymin=avgChl-SE, ymax=avgChl+SE), width=0.2) + labs(title="UCYN-A1 vs MODIS Chl", x="", y="[Chl]") + 
  theme_classic() + 
  theme(axis.title.y=element_text(size=20, face="bold", color="black"), axis.text=element_text(size=20, face="bold", color="black"))
Chl_plot  

ggsave("/Users/cfletcherh1/Desktop/DataAnalysis/SPOT_16S_18S_eLSA/Plots/04.20.2022_abio_brad_pres_abs/abund_MODIS_Chl_1_plot.pdf", height = 4, width = 3.5, units = "in", dpi=300)

#A2
Chl_table=env_subs %>% group_by(A2_abund) %>% summarize(avgChl=mean(as.numeric(MODIS_Chl)),
                                                       SD=sd(as.numeric(MODIS_Chl)), n=n(),
                                                       SE=sd(as.numeric(MODIS_Chl))/sqrt(n()))

Chl_plot=ggplot(Chl_table, aes(x=A2_abund, y=avgChl))+
  geom_point(size=10, color="#66CC99") + 
  geom_errorbar(aes(ymin=avgChl-SE, ymax=avgChl+SE), width=0.2) + labs(title="UCYN-A2 vs MODIS Chl", x="", y="[Chl]") +
  theme_classic() + 
  theme(axis.title.y=element_text(size=20, face="bold", color="black"), axis.text=element_text(size=20, face="bold", color="black"))
Chl_plot

ggsave("/Users/cfletcherh1/Desktop/DataAnalysis/SPOT_16S_18S_eLSA/Plots/04.20.2022_abio_brad_pres_abs/abund_MODIS_Chl_2_plot.pdf", height = 4, width = 3.5, units = "in", dpi=300)

#Bacterial production leucine 
leu_table=env_subs %>% group_by(A1_abund) %>% summarize(avg_leu=mean(BactProd_Leu), 
                                                       SD=sd(BactProd_Leu), n=n(),
                                                       SE=sd(BactProd_Leu)/sqrt(n()))

Leu_plot=ggplot(leu_table, aes(x=A1_abund, y=avg_leu)) + 
  geom_point(size=10, color="#999999") + 
  geom_errorbar(aes(ymin=avg_leu-SE, ymax=avg_leu+SE), width=0.2) + labs(title="UCYN-A1 vs Bacterial Production - Leucine", x="", y="Bact Prod. (cells/mL/day)") + 
  theme_classic() + 
  theme(axis.title.y=element_text(size=20, face="bold", color="black"), axis.text=element_text(size=20, face="bold", color="black"))
Leu_plot 

ggsave("/Users/cfletcherh1/Desktop/DataAnalysis/SPOT_16S_18S_eLSA/Plots/04.20.2022_abio_brad_pres_abs/abund_BactProd_Leu_1_plot.pdf", height = 4, width = 3.5, units = "in", dpi=300)

#A2 
leu_table=env_subs %>% group_by(A2_abund) %>% summarize(avg_leu=mean(BactProd_Leu), 
                                                       SD=sd(BactProd_Leu), n=n(),
                                                       SE=sd(BactProd_Leu)/sqrt(n()))

Leu_plot=ggplot(leu_table, aes(x=A2_abund, y=avg_leu)) + 
  geom_point(size=10, col="#999999") + 
  geom_errorbar(aes(ymin=avg_leu-SE, ymax=avg_leu+SE), width=0.2) + labs(title="UCYN-A2 vs Bacterial Production - Leucine", x="", y="Bact Prod. (cells/mL/day)") + 
  theme_classic() + 
  theme(axis.title.y=element_text(size=20, face="bold", color="black"), axis.text=element_text(size=20, face="bold", color="black"))
Leu_plot

ggsave("/Users/cfletcherh1/Desktop/DataAnalysis/SPOT_16S_18S_eLSA/Plots/04.20.2022_abio_brad_pres_abs/abund_BactProd_Leu_2_plot.pdf", height = 4, width = 3.5, units = "in", dpi=300)

#6. Plot CLR transformed UCYN-A abundances vs. environmental data----
names(env_subs)
head(env_subs$Year_Month)
head(names(proks_CLR))

#Subset proks_CLR data to get UCYN-A data only
include_rows=c()
for(i in UCYNA_ix){
  include_rows=c(include_rows, grep(i, proks_CLR$OTU_ID))
}

CLR_subs=proks_CLR[include_rows,]
CLR_subs$OTU_ID
dim(CLR_subs) #6 126
head(names(CLR_subs))

#format sampling date like year_month 
year_month=c()
for(i in c(2:ncol(CLR_subs))){
  #print(names(CLR_subs)[i])
  a=strsplit(names(CLR_subs)[i], split=".", fixed=T)
  year=a[[1]][2]
  #print(year)
  month=a[[1]][3]
  if(month<10){
    b=strsplit(month, split="", fixed=T)
    month=b[[1]][2]
  }
  #print(month)
  year_month=c(year_month, paste(year, month, sep="_"))
}
head(year_month)
head(env_subs$Year_Month)

#Transform CLR data into data frame
CLR_t=as.data.frame(t(CLR_subs))
dim(CLR_t) #126 6
head(CLR_t)
names(CLR_t)
CLR_t[1,]
names(CLR_t)=CLR_t[1,]
names(CLR_t)
head(rownames(CLR_t))

#Get rid of that first row
CLR_t=CLR_t[-1,]
dim(CLR_t) #125 6
head(CLR_t)

#Add in year_month as a column 
dim(CLR_t) #125 6
length(year_month) #125
CLR_t$Year_Month=year_month
head(CLR_t)

#Figure out which samples are in both CLR_t and env_subs

env_subs$Year_Month %in% CLR_t$Year_Month #All
CLR_t$Year_Month %in% env_subs$Year_Month #First one and last two are missing
length(env_subs$Year_Month) #122 #Ah, there we go

#Add in transformed abundances of UCYN-A ASV1 and ASV5 into env_subs
env_subs$A1_CLR=CLR_t[c(2:123),1]
env_subs$A2_CLR=CLR_t[c(2:123),5]

#Now plot!! 
#Throw this in a .PDF
pdf(file="Plots/CLR_transformed_UCYNA_ASVs_vs_env_parameters_07.14.2022.pdf")
for(i in c(2:4, 8:13)){
  plot(x=env_subs[,i], y=env_subs$A1_CLR, pch=16, xlab=names(env_subs)[i], ylab="UCYN-A1 CLR Abundance", main="UCYN-A1")
  plot(x=env_subs[,i], y=env_subs$A2_CLR, pch=16, xlab=names(env_subs)[i], ylab="UCYN-A2 CLR Abundance", main="UCYN-A2")
}
dev.off()


