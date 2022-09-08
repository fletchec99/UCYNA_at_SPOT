#Plot UCYN-A ASVs vs depth on 7.9.2009, the date it appears at 890m
#11.19.2021

getwd()
setwd("Desktop/DataAnalysis/SPOT_16S_18S_eLSA/")

#1. Read in data
data=read.table("ModifiedFiles/4.SPOT_16S_no_chloro_proportions.tsv", header=T, stringsAsFactors = F)
dim(data) #72284 1206

proks_tax=read.table("ModifiedFiles/Proks_tax_classified_23072020_SILVA132.tsv", stringsAsFactors = T, header=T, sep="\t")
dim(proks_tax) #78953 3

#2. Set up appropriate variables
UCYNA_ix=proks_tax$Feature.ID[grep("UCYN", proks_tax$Taxon)]

#Graph of 07.2009, when UCYN-A present at depths------

#Subset data, format it 
head(colnames(data))
colnames(data)[grep("2009.07", colnames(data))]

include_rows=c()
for(i in UCYNA_ix){
  include_rows=c(include_rows, grep(i, data$OTU_ID))
}
include_rows

data_subs=data[include_rows, grep("2009.07", colnames(data))]
dim(data_subs) #6 10, looks good
data_subs=as.data.frame(t(data_subs))
colnames(data_subs)=paste("UCYNA", "ASV", c(1:6), sep="_")

#Figure out how to get depth 
data_subs$ID=row.names(data_subs)
data_subs$depth_category="deep"

for(i in c(1:nrow(data_subs))){
  a=data_subs$ID[i]
  data_subs$depth_category[i]=unlist(strsplit(a, split="\\."))[5]
}
data_subs$depth_category

#Wait, check DCM depth
abio_data=read.table("ModifiedFiles/SPOT_nutrient_CTD_prod_Hammond_data_04.20.2021.tsv", header=T, stringsAsFactors = F, sep="\t")
dim(abio_data) #860 18
names(abio_data)
which(abio_data$Year==2009 & abio_data$Month==07)
abio_data[which(abio_data$Year==2009 & abio_data$Month==07), c(1:5)]
#DCM was 22m

#Write in depth 
data_subs$depth=5
data_subs$depth[which(data_subs$depth_category=="5m")]
data_subs$depth[which(data_subs$depth_category=="DCM")]=22
data_subs$depth[which(data_subs$depth_category=="150m")]=150
data_subs$depth[which(data_subs$depth_category=="500m")]=500
data_subs$depth[which(data_subs$depth_category=="890m")]=890

#Write in filter type
data_subs$Filter="AE"
for(i in c(1:nrow(data_subs))){
  a=data_subs$ID[i]
  data_subs$Filter[i]=unlist(strsplit(a, split="\\."))[6]
  #print(b)
}

#Actually, just remove all Durapore samples bc UCYN-A abundances are zero anyway
data_AE=data_subs[which(data_subs$Filter=="AE"),]
dim(data_AE) #5 10
dim(data_subs) #10 10 

#Now PLOT! 
library(ggplot2)

ggplot(data_AE, aes(x=depth, y=UCYNA_ASV_1))+
  geom_line(color="red", size=1, linetype=2) + 
  geom_point() + 
  scale_x_reverse() + 
  scale_y_continuous(position="right") + 
  coord_flip() 

ggplot(data_AE, aes(x=depth)) +
  geom_line(aes(y=UCYNA_ASV_1), color="red", size=1, linetype=1) + 
  geom_line(aes(y=UCYNA_ASV_2), color="orange", size=1, linetype=2) + 
  
  geom_point(aes(y=UCYNA_ASV_1), color="red", size=3, shape=15) + 
  geom_point(aes(y=UCYNA_ASV_2), color="orange", size=3, shape=16) + 
  
  scale_color_manual(name="Line Color", values=c(ASV1="red", ASV2="orange")) + 
  
  scale_x_reverse("Depth", breaks=c(5, 22, 150, 500, 890)) + 
  scale_y_continuous("Relative abundance out of 16S sequences in 1-80um size fraction", position="right") + 
  coord_flip() + 
  theme(legend.position=c(5, 0.005), legend.justification = c("right", "top")) + 
  theme_bw()

#Spruced up plot! :) 
ggplot(data_AE, aes(x=depth)) +
  geom_line(aes(y=UCYNA_ASV_1, color="ASV1"), size=1, linetype=1) + 
  geom_line(aes(y=UCYNA_ASV_2, color="ASV2"), size=1, linetype=2) + 
  geom_line(aes(y=UCYNA_ASV_3, color="ASV3"), size=1, linetype=3) + 
  geom_line(aes(y=UCYNA_ASV_4, color="ASV4"), size=1, linetype=4) + 
  geom_line(aes(y=UCYNA_ASV_5, color="ASV5"), size=1, linetype=5) + 
  geom_line(aes(y=UCYNA_ASV_6, color="ASV6"), size=1, linetype=6) +

  geom_point(aes(y=UCYNA_ASV_1, color="ASV1"), size=3, shape=15) + 
  geom_point(aes(y=UCYNA_ASV_2, color="ASV2"), size=3, shape=0) + 
  geom_point(aes(y=UCYNA_ASV_3, color="ASV3"), size=3, shape=1) + 
  geom_point(aes(y=UCYNA_ASV_4, color="ASV4"), size=3, shape=2) +
  geom_point(aes(y=UCYNA_ASV_5, color="ASV5"), size=3, shape=4) + 
  geom_point(aes(y=UCYNA_ASV_6, color="ASV6"), size=3, shape=5) +
  
  scale_color_manual(name="ASV Number", values=c(ASV1="red", ASV2="orange", ASV3="yellow", ASV4="green", ASV5="blue", ASV6="purple")) + 
  
  scale_x_reverse("Depth", breaks=c(5, 22, 150, 500, 890)) + 
  scale_y_continuous("Relative abundance out of 16S sequences \n in 1-80um size fraction", position="right") + 
  annotate("label", y = 0.003, x = 890, label = "July 2009") +
  coord_flip() + 
  theme_bw()

#II. Plot July 2008, when UCYN-A ASV1 was at its maximum abundance-----
#Figure out what date UCYN-A ASV1 has maximum abundance and plot that
dim(data)
head(colnames(data))
summary(as.numeric(data[grep(UCYNA_ix[1], data$OTU_ID), 2:ncol(data)]))
max(as.numeric(data[grep(UCYNA_ix[1], data$OTU_ID), 2:ncol(data)]))
which(as.numeric(data[grep(UCYNA_ix[1], data$OTU_ID), 2:ncol(data)])==max(as.numeric(data[grep(UCYNA_ix[1], data$OTU_ID), 2:ncol(data)])))
colnames(data)[380+1] #2008.07.09.5m.AE 
data[grep(UCYNA_ix[1], data$OTU_ID), 381] #Yep, that's it! 

#Subset data, format it 
head(colnames(data))
colnames(data)[grep("2008.07", colnames(data))]

include_rows=c()
for(i in UCYNA_ix){
  include_rows=c(include_rows, grep(i, data$OTU_ID))
}
include_rows

data_subs=data[include_rows, grep("2008.07", colnames(data))]
dim(data_subs) #6 9, looks good
data_subs=as.data.frame(t(data_subs))
colnames(data_subs)=paste("UCYNA", "ASV", c(1:6), sep="_")

#Figure out how to get depth 
data_subs$ID=row.names(data_subs)
data_subs$depth_category="deep"

for(i in c(1:nrow(data_subs))){
  a=data_subs$ID[i]
  data_subs$depth_category[i]=unlist(strsplit(a, split="\\."))[5]
}
data_subs$depth_category #Ok so there is no 890m AE filter on this date 

#Wait, check DCM depth
dim(abio_data) #860 18
names(abio_data)
which(abio_data$Year==2008 & abio_data$Month==07)
abio_data[which(abio_data$Year==2008 & abio_data$Month==07), c(1:5)]
#DCM was 23.5m

#Write in depth 
data_subs$depth=5
data_subs$depth[which(data_subs$depth_category=="5m")]
data_subs$depth[which(data_subs$depth_category=="DCM")]=23.5
data_subs$depth[which(data_subs$depth_category=="150m")]=150
data_subs$depth[which(data_subs$depth_category=="500m")]=500
data_subs$depth[which(data_subs$depth_category=="890m")]=890

#Write in filter type
data_subs$Filter="AE"
for(i in c(1:nrow(data_subs))){
  a=data_subs$ID[i]
  data_subs$Filter[i]=unlist(strsplit(a, split="\\."))[6]
  #print(b)
}

#Actually, just remove all Durapore samples bc UCYN-A abundances are zero anyway
data_AE=data_subs[which(data_subs$Filter=="AE"),]
dim(data_AE) #4 10
dim(data_subs) #9 10 

#Now PLOT! 
ggplot(data_AE, aes(x=depth)) +
  geom_line(aes(y=UCYNA_ASV_1, color="ASV1"), size=1, linetype=1) + 
  geom_line(aes(y=UCYNA_ASV_2, color="ASV2"), size=1, linetype=2) + 
  geom_line(aes(y=UCYNA_ASV_3, color="ASV3"), size=1, linetype=3) + 
  geom_line(aes(y=UCYNA_ASV_4, color="ASV4"), size=1, linetype=4) + 
  geom_line(aes(y=UCYNA_ASV_5, color="ASV5"), size=1, linetype=5) + 
  geom_line(aes(y=UCYNA_ASV_6, color="ASV6"), size=1, linetype=6) +
  
  geom_point(aes(y=UCYNA_ASV_1, color="ASV1"), size=3, shape=15) + 
  geom_point(aes(y=UCYNA_ASV_2, color="ASV2"), size=3, shape=0) + 
  geom_point(aes(y=UCYNA_ASV_3, color="ASV3"), size=3, shape=1) + 
  geom_point(aes(y=UCYNA_ASV_4, color="ASV4"), size=3, shape=2) +
  geom_point(aes(y=UCYNA_ASV_5, color="ASV5"), size=3, shape=4) + 
  geom_point(aes(y=UCYNA_ASV_6, color="ASV6"), size=3, shape=5) +
  
  scale_color_manual(name="ASV Number", values=c(ASV1="red", ASV2="orange", ASV3="yellow", ASV4="green", ASV5="blue", ASV6="purple")) + 
  
  scale_x_reverse("Depth", breaks=c(5, 23.5, 150, 500)) + 
  scale_y_continuous("Relative abundance out of 16S sequences \n in 1-80um size fraction", position="right") + 
  annotate("label", y = 0.03, x = 500, label = "July 2008") +
  coord_flip() + 
  theme_bw()

