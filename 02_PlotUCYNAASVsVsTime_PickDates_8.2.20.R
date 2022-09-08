#Plotting up relative abundances of UCYNA ASVs!!!
#8.2.2020

rm(list=ls())
getwd()
setwd("/Users/cfletcherh1/Desktop/DataAnalysis/SPOT_16S_18S_eLSA/")

#Read in relative abundance of Proks and euks + taxonomy files-----
#Want to graph relative abundance of 16S with no chloro and 18S with no metazoans
proks_data <- read.table("ModifiedFiles/4.SPOT_16S_no_chloro_proportions.tsv", header=T, stringsAsFactors = F, sep="\t")
dim(proks_data) #72284 1206
head(names(proks_data))

proks_tax <- read.table("ModifiedFiles/Proks_tax_classified_23072020_SILVA132.tsv", header=T, stringsAsFactors = F, sep="\t")
dim(proks_tax) #78953 3
head(proks_tax)

length(grep("Chloroplast", proks_tax$Taxon))+nrow(proks_data) ==nrow(proks_tax) #TRUE

euks_data <- read.table("ModifiedFiles/8.SPOT_18S_no_metaz_proportions.tsv", header=T, stringsAsFactors = F, sep="\t")
dim(euks_data) #28402 1098
head(names(euks_data))

euks_tax <- read.table("ModifiedFiles/Euks_tax_classified_19052020_SILVA132.tsv", header=T, stringsAsFactors = F, sep="\t")
dim(euks_tax) #30537 3

length(grep("Metaz", euks_tax$Taxon))
nrow(euks_data)+length(grep("Metaz", euks_tax$Taxon)) == nrow(euks_tax) #TRUE

#Figure out which ASVs correspond to UCYNA, Brad, and Prymnesiophytes-----
#UCYN-A: 
length(grep("UCYN-A", proks_tax$Taxon))

UCYNA_ASVs=proks_tax$Feature.ID[grep("UCYN-A", proks_tax$Taxon)]
length(UCYNA_ASVs)
length(unique(UCYNA_ASVs))

UCYNA_ix <- c()
for(i in UCYNA_ASVs){
  UCYNA_ix <- c(UCYNA_ix, grep(i, proks_data$OTU_ID))
}
UCYNA_ix

#Write this out into a table
UCYNA_table=proks_tax[grep("UCYN-A", proks_tax$Taxon),]
dim(UCYNA_table)
UCYNA_table$ASV_num=c(1:6)
UCYNA_table$RowNumProks_data=UCYNA_ix
View(UCYNA_table)
write.table(x=UCYNA_table, file="UCYNA_ASVs_Numb_conf_hashes_8.2.2020.csv", quote=F, row.names=F, sep=",")

#Figure out which ASVs correspond to Brad
length(grep("Braarudosphaera", euks_tax$Taxon)) #7
length(grep("species_Braarudosphaeraceae_X_sp", euks_tax$Taxon)) #3
euks_tax$Taxon[grep("Braarudosphaera", euks_tax$Taxon)]

Brad_ASVs=euks_tax$Feature.ID[grep("Braarudosphaera_bigelowii", euks_tax$Taxon)]
Brad_ASVs #Only corresponds to B Bigelowii ASVs

Brad_ix <- c()
for(i in Brad_ASVs){
  Brad_ix <- c(Brad_ix, grep(i, euks_data$OTU_ID))
}
Brad_ix

genus_brad_ASVs=euks_tax$Feature.ID[grep("Braarudosphaera", euks_tax$Taxon)]
genus_brad_ASVs #Corresponds to all 7 ASVs from Braarudospharea

genus_brad_ix <- c()
for(i in genus_brad_ASVs){
  genus_brad_ix <- c(genus_brad_ix, grep(i, euks_data$OTU_ID))
}
genus_brad_ix

#Write this into a table
Brad_table=euks_tax[grep("Braarudosphaera", euks_tax$Taxon),]
dim(Brad_table)
Brad_table$ASV_num=c(1, "5x", 2, 3, 4, "6x", "7x")
Brad_table$RowNumEuks_data=genus_brad_ix

write.table(x=Brad_table, file="Bradarodusphaera_genus_ASVs_Numb_conf_hashes_8.2.2020.csv", quote=F, row.names=F, sep=",")

#Write out prymnesiophyte ASVs into a table
#Grep for "class_Prymn..."
length(grep("class_Prymnesiophyceae", euks_tax$Taxon)) #643 prymnesiophytes
prym_table=euks_tax[grep("class_Prymnesiophyceae", euks_tax$Taxon),]
dim(prym_table)

#Add in row numbers in euks_data (but not going to number the ASVs)
prym_ix <- c()
for(i in euks_tax$Feature.ID[grep("class_Prymnesiophyceae", euks_tax$Taxon)]){
  prym_ix <- c(prym_ix, grep(i, euks_data$OTU_ID))
}
length(prym_ix)

prym_table$RowNumEuks_data=prym_ix
dim(prym_table)

write.table(x=prym_table, file="Prymnesiophyte_ASVs_Numb_conf_hashes_8.2.2020.csv", quote=F, row.names=F, sep=",")

#Set up several vectors of SPOT sampling dates -----
#First all dates, so that it corresponds exactly to indices of SPOT data
all_dates <- c()
for(i in c(2:ncol(proks_data))){
  a <- strsplit(colnames(proks_data)[i], split=".", fixed=T)
  b <- paste(a[[1]][2],a[[1]][3],a[[1]][4], sep="-")
  all_dates <- c(all_dates, b)
}
head(all_dates)
length(all_dates) 
length(unique(all_dates)) #We only have samples from 213 dates?? NOOOOOOOOØ
#But I mean, ideally, we would have 240 samples
all_dates <- as.Date(all_dates, format="%Y-%m-%d") #Set to be dates
#all_dates <- c("OTU_ID", all_dates) #Add in "OTU ID" so indices match those of proks data perfectly
#crap that's not gonna work, just converts it to a seemingly random string of numbers
#Gonna have to remove 1 from indexing below
length(all_dates) #1205
head(all_dates)

#Now set up all dates but for euks data
euks_dates <- c()
for(i in c(2:ncol(euks_data))){
  a <- (strsplit(colnames(euks_data)[i], split=".", fixed=T))
  b <- paste(a[[1]][2],a[[1]][3],a[[1]][4], sep="-")
  euks_dates <- c(euks_dates, b)
}
length(colnames(euks_data)) #1098
length(euks_dates) #1097
length(unique(euks_dates)) #213 #So presumably all the sampling dates have proks and euks data
euks_dates <- as.Date(euks_dates)

#Set up dates where we have AE samples
AE_dates <- c()
for(i in (colnames(proks_data)[grep("5m.AE", colnames(proks_data))])){
  a <- strsplit(i, split=".", fixed=T)
  b <- paste(a[[1]][2],a[[1]][3],a[[1]][4], sep="-")
  AE_dates <- c(AE_dates, b)
}
length(AE_dates) #125
length(unique(AE_dates)) #125, good
AE_dates=as.Date(AE_dates, format="%Y-%m-%d")
AE_dates

#Durapore dates
D_dates <- c()
for(i in colnames(proks_data)[grep("5m.D", colnames(proks_data))]){
  a <- strsplit(i, split=".", fixed=T)
  b <- paste(a[[1]][2],a[[1]][3],a[[1]][4], sep="-")
  D_dates <- c(D_dates, b)
}
length(unique(D_dates)) #188
length(D_dates) #188
D_dates <- as.Date(D_dates, format="%Y-%m-%d")
D_dates

#Subset to take dates we have AE AND Durapore samples present
include_dates <- c()
for(i in c(1:length(AE_dates))){
  AE <- AE_dates[i]
  include_dates <- c(include_dates, grep(AE, D_dates))
}
length(include_dates) #107 #We only have 107 dates with 5m AE AND Dura samples?
AE_D_dates <- D_dates[include_dates]
length(AE_D_dates) #107

#now plot!
graphics.off()
par("mar")
par(mar=c(1,1,1,1))

UCYNA_ix

#Plot up UCYNA ASVs/ time, all depths, save it in a .PDF------
#Save it in a .PDF
pdf("Six_UCYNA_ASVs_over_time_all_depths_08042020.PDF")

#UCYNA ASV1: Row 223 in Proks data, 3d852410f44d21c92c9c55fbbb25187e--------
#AE from 5m depth first, color= red
op <- par(mar=c(5, 6, 4, 2) + 0.1)
plot(as.Date(all_dates[(grep("5m.AE", colnames(proks_data))+1)]),as.numeric(proks_data[223, grep("5m.AE", colnames(proks_data))]), col="red", pch=16, type="b", xlab="SPOT Sampling Date", ylab="Proportion of 16S sequences \n in 1-80um size class", main= "UCYN-A ASV1 AE 5m (3d852410f44d21c92c9c55fbbb25187e)")
par(op)

#Durapore from 5m depth, color=red, pch= something else
op <- par(mar=c(5, 6, 4, 2) + 0.1)
plot(as.Date(all_dates[(grep("5m.D", colnames(proks_data))+1)]),as.numeric(proks_data[223, grep("5m.D", colnames(proks_data))]), col="red", pch=4, type="b", lwd=2, xlab="SPOT Sampling Date", ylab="Proportion of 16S sequences \n in 0.22-1um size class", main= "UCYN-A ASV1 Dura. 5m (3d852410f44d21c92c9c55fbbb25187e)")
par(op)

#AE from DCM depth, color= green
op <- par(mar=c(5, 6, 4, 2) + 0.1)
plot(as.Date(all_dates[(grep("DCM.AE", colnames(proks_data))+1)]),as.numeric(proks_data[223, grep("DCM.AE", colnames(proks_data))]), col="forestgreen", pch=16, type="b", xlab="SPOT Sampling Date", ylab="Proportion of 16S sequences \n in 1-80um size class", main= "UCYN-A ASV1 AE DCM (3d852410f44d21c92c9c55fbbb25187e)")
par(op)

#Durapore from DCM depth, color= green
op <- par(mar=c(5, 6, 4, 2) + 0.1)
plot(as.Date(all_dates[(grep("DCM.D", colnames(proks_data))+1)]),as.numeric(proks_data[223, grep("DCM.D", colnames(proks_data))]), lwd=2, col="forestgreen", pch=4, type="b", xlab="SPOT Sampling Date", ylab="Proportion of 16S sequences \n in 0.22-1um size class", main= "UCYN-A ASV1 Dura. DCM (3d852410f44d21c92c9c55fbbb25187e)")
par(op)

#AE from 150m depth, color= blue
op <- par(mar=c(5, 6, 4, 2) + 0.1)
plot(as.Date(all_dates[(grep("150m.AE", colnames(proks_data))+1)]),as.numeric(proks_data[223, grep("150m.AE", colnames(proks_data))]), col="blue", pch=16, type="b", xlab="SPOT Sampling Date", ylab="Proportion of 16S sequences \n in 1-80um size class", main= "UCYN-A ASV1 AE 150m (3d852410f44d21c92c9c55fbbb25187e)")
par(op)

#Durapore from 150m depth, color= blue
op <- par(mar=c(5, 6, 4, 2) + 0.1)
plot(as.Date(all_dates[(grep("150m.D", colnames(proks_data))+1)]),as.numeric(proks_data[223, grep("150m.D", colnames(proks_data))]), lwd=2, col="blue", pch=4, type="b", xlab="SPOT Sampling Date", ylab="Proportion of 16S sequences \n in 0.22-1um size class", main= "UCYN-A ASV1 Dura. 150m (3d852410f44d21c92c9c55fbbb25187e)")
par(op)

#AE from 500m depth, color= purple
op <- par(mar=c(5, 6, 4, 2) + 0.1)
plot(as.Date(all_dates[(grep("500m.AE", colnames(proks_data))+1)]),as.numeric(proks_data[223, grep("500m.AE", colnames(proks_data))]), col="purple", pch=16, type="b", xlab="SPOT Sampling Date", ylab="Proportion of 16S sequences \n in 1-80um size class", main= "UCYN-A ASV1 AE 500m (3d852410f44d21c92c9c55fbbb25187e)")
par(op)

#Durapore from 500m depth, color= blue
op <- par(mar=c(5, 6, 4, 2) + 0.1)
plot(as.Date(all_dates[(grep("500m.D", colnames(proks_data))+1)]),as.numeric(proks_data[223, grep("500m.D", colnames(proks_data))]), lwd=2, col="purple", pch=4, type="b", xlab="SPOT Sampling Date", ylab="Proportion of 16S sequences \n in 0.22-1um size class", main= "UCYN-A ASV1 Dura. 500m (3d852410f44d21c92c9c55fbbb25187e)")
par(op)

#AE from 890m depth, color= purple
op <- par(mar=c(5, 6, 4, 2) + 0.1)
plot(as.Date(all_dates[(grep("890m.AE", colnames(proks_data))+1)]),as.numeric(proks_data[223, grep("890m.AE", colnames(proks_data))]), pch=16, type="b", xlab="SPOT Sampling Date", ylab="Proportion of 16S sequences \n in 1-80um size class", main= "UCYN-A ASV1 AE 890m (3d852410f44d21c92c9c55fbbb25187e)")
par(op)

#Durapore from 890m depth, color= black
op <- par(mar=c(5, 6, 4, 2) + 0.1)
plot(as.Date(all_dates[(grep("890m.D", colnames(proks_data))+1)]),as.numeric(proks_data[223, grep("890m.D", colnames(proks_data))]), lwd=2, pch=4, type="b", xlab="SPOT Sampling Date", ylab="Proportion of 16S sequences \n in 0.22-1um size class", main= "UCYN-A ASV1 Dura. 890m (3d852410f44d21c92c9c55fbbb25187e)")
par(op)


#UCYNA ASV2: Row 13635 in Proks data, 42c9bf8576acc275f3c9281e6b24f5a3--------
#AE from 5m depth first, color= red
op <- par(mar=c(5, 6, 4, 2) + 0.1)
plot(as.Date(all_dates[(grep("5m.AE", colnames(proks_data))+1)]),as.numeric(proks_data[13635, grep("5m.AE", colnames(proks_data))]), col="red", pch=16, type="b", xlab="SPOT Sampling Date", ylab="Proportion of 16S sequences \n in 1-80um size class", main= "UCYN-A ASV2 AE 5m (42c9bf8576acc275f3c9281e6b24f5a3)")
par(op)

#Durapore from 5m depth, color=red, pch= something else
op <- par(mar=c(5, 6, 4, 2) + 0.1)
plot(as.Date(all_dates[(grep("5m.D", colnames(proks_data))+1)]),as.numeric(proks_data[13635, grep("5m.D", colnames(proks_data))]), col="red", pch=4, type="b", lwd=2, xlab="SPOT Sampling Date", ylab="Proportion of 16S sequences \n in 0.22-1um size class", main= "UCYN-A ASV2 Dura. 5m (42c9bf8576acc275f3c9281e6b24f5a3)")
par(op)

#AE from DCM depth, color= green
op <- par(mar=c(5, 6, 4, 2) + 0.1)
plot(as.Date(all_dates[(grep("DCM.AE", colnames(proks_data))+1)]),as.numeric(proks_data[13635, grep("DCM.AE", colnames(proks_data))]), col="forestgreen", pch=16, type="b", xlab="SPOT Sampling Date", ylab="Proportion of 16S sequences \n in 1-80um size class", main= "UCYN-A ASV2 AE DCM (42c9bf8576acc275f3c9281e6b24f5a3)")
par(op)

#Durapore from DCM depth, color= green
op <- par(mar=c(5, 6, 4, 2) + 0.1)
plot(as.Date(all_dates[(grep("DCM.D", colnames(proks_data))+1)]),as.numeric(proks_data[13635, grep("DCM.D", colnames(proks_data))]), lwd=2, col="forestgreen", pch=4, type="b", xlab="SPOT Sampling Date", ylab="Proportion of 16S sequences \n in 0.22-1um size class", main= "UCYN-A ASV2 Dura. DCM (42c9bf8576acc275f3c9281e6b24f5a3)")
par(op)

#AE from 150m depth, color= blue
op <- par(mar=c(5, 6, 4, 2) + 0.1)
plot(as.Date(all_dates[(grep("150m.AE", colnames(proks_data))+1)]),as.numeric(proks_data[13635, grep("150m.AE", colnames(proks_data))]), col="blue", pch=16, type="b", xlab="SPOT Sampling Date", ylab="Proportion of 16S sequences \n in 1-80um size class", main= "UCYN-A ASV2 AE 150m (42c9bf8576acc275f3c9281e6b24f5a3)")
par(op)

#Durapore from 150m depth, color= blue
op <- par(mar=c(5, 6, 4, 2) + 0.1)
plot(as.Date(all_dates[(grep("150m.D", colnames(proks_data))+1)]),as.numeric(proks_data[13635, grep("150m.D", colnames(proks_data))]), lwd=2, col="blue", pch=4, type="b", xlab="SPOT Sampling Date", ylab="Proportion of 16S sequences \n in 0.22-1um size class", main= "UCYN-A ASV2 Dura. 150m (42c9bf8576acc275f3c9281e6b24f5a3)")
par(op)

#AE from 500m depth, color= purple
op <- par(mar=c(5, 6, 4, 2) + 0.1)
plot(as.Date(all_dates[(grep("500m.AE", colnames(proks_data))+1)]),as.numeric(proks_data[13635, grep("500m.AE", colnames(proks_data))]), col="purple", pch=16, type="b", xlab="SPOT Sampling Date", ylab="Proportion of 16S sequences \n in 1-80um size class", main= "UCYN-A ASV2 AE 500m (42c9bf8576acc275f3c9281e6b24f5a3)")
par(op)

#Durapore from 500m depth, color= blue
op <- par(mar=c(5, 6, 4, 2) + 0.1)
plot(as.Date(all_dates[(grep("500m.D", colnames(proks_data))+1)]),as.numeric(proks_data[13635, grep("500m.D", colnames(proks_data))]), lwd=2, col="purple", pch=4, type="b", xlab="SPOT Sampling Date", ylab="Proportion of 16S sequences \n in 0.22-1um size class", main= "UCYN-A ASV2 Dura. 500m (42c9bf8576acc275f3c9281e6b24f5a3)")
par(op)

#AE from 890m depth, color= purple
op <- par(mar=c(5, 6, 4, 2) + 0.1)
plot(as.Date(all_dates[(grep("890m.AE", colnames(proks_data))+1)]),as.numeric(proks_data[13635, grep("890m.AE", colnames(proks_data))]), pch=16, type="b", xlab="SPOT Sampling Date", ylab="Proportion of 16S sequences \n in 1-80um size class", main= "UCYN-A ASV2 AE 890m (42c9bf8576acc275f3c9281e6b24f5a3)")
par(op)

#Durapore from 890m depth, color= black
op <- par(mar=c(5, 6, 4, 2) + 0.1)
plot(as.Date(all_dates[(grep("890m.D", colnames(proks_data))+1)]),as.numeric(proks_data[13635, grep("890m.D", colnames(proks_data))]), lwd=2, pch=4, type="b", xlab="SPOT Sampling Date", ylab="Proportion of 16S sequences \n in 0.22-1um size class", main= "UCYN-A ASV2 Dura. 890m (42c9bf8576acc275f3c9281e6b24f5a3)")
par(op)

#UCYNA ASV3: Row 54043 in Proks data, 6115eab19c52bc45c6ba11d72ec88031--------
#AE from 5m depth first, color= red
op <- par(mar=c(5, 6, 4, 2) + 0.1)
plot(as.Date(all_dates[(grep("5m.AE", colnames(proks_data))+1)]),as.numeric(proks_data[54043, grep("5m.AE", colnames(proks_data))]), col="red", pch=16, type="b", xlab="SPOT Sampling Date", ylab="Proportion of 16S sequences \n in 1-80um size class", main= "UCYN-A ASV3 AE 5m (6115eab19c52bc45c6ba11d72ec88031)")
par(op)

#Durapore from 5m depth, color=red, pch= something else
op <- par(mar=c(5, 6, 4, 2) + 0.1)
plot(as.Date(all_dates[(grep("5m.D", colnames(proks_data))+1)]),as.numeric(proks_data[54043, grep("5m.D", colnames(proks_data))]), col="red", pch=4, type="b", lwd=2, xlab="SPOT Sampling Date", ylab="Proportion of 16S sequences \n in 0.22-1um size class", main= "UCYN-A ASV3 Dura. 5m (6115eab19c52bc45c6ba11d72ec88031)")
par(op)

#AE from DCM depth, color= green
op <- par(mar=c(5, 6, 4, 2) + 0.1)
plot(as.Date(all_dates[(grep("DCM.AE", colnames(proks_data))+1)]),as.numeric(proks_data[54043, grep("DCM.AE", colnames(proks_data))]), col="forestgreen", pch=16, type="b", xlab="SPOT Sampling Date", ylab="Proportion of 16S sequences \n in 1-80um size class", main= "UCYN-A ASV3 AE DCM (6115eab19c52bc45c6ba11d72ec88031)")
par(op)

#Durapore from DCM depth, color= green
op <- par(mar=c(5, 6, 4, 2) + 0.1)
plot(as.Date(all_dates[(grep("DCM.D", colnames(proks_data))+1)]),as.numeric(proks_data[54043, grep("DCM.D", colnames(proks_data))]), lwd=2, col="forestgreen", pch=4, type="b", xlab="SPOT Sampling Date", ylab="Proportion of 16S sequences \n in 0.22-1um size class", main= "UCYN-A ASV3 Dura. DCM (6115eab19c52bc45c6ba11d72ec88031)")
par(op)

#AE from 150m depth, color= blue
op <- par(mar=c(5, 6, 4, 2) + 0.1)
plot(as.Date(all_dates[(grep("150m.AE", colnames(proks_data))+1)]),as.numeric(proks_data[54043, grep("150m.AE", colnames(proks_data))]), col="blue", pch=16, type="b", xlab="SPOT Sampling Date", ylab="Proportion of 16S sequences \n in 1-80um size class", main= "UCYN-A ASV3 AE 150m (6115eab19c52bc45c6ba11d72ec88031)")
par(op)

#Durapore from 150m depth, color= blue
op <- par(mar=c(5, 6, 4, 2) + 0.1)
plot(as.Date(all_dates[(grep("150m.D", colnames(proks_data))+1)]),as.numeric(proks_data[54043, grep("150m.D", colnames(proks_data))]), lwd=2, col="blue", pch=4, type="b", xlab="SPOT Sampling Date", ylab="Proportion of 16S sequences \n in 0.22-1um size class", main= "UCYN-A ASV3 Dura. 150m (6115eab19c52bc45c6ba11d72ec88031)")
par(op)

#AE from 500m depth, color= purple
op <- par(mar=c(5, 6, 4, 2) + 0.1)
plot(as.Date(all_dates[(grep("500m.AE", colnames(proks_data))+1)]),as.numeric(proks_data[54043, grep("500m.AE", colnames(proks_data))]), col="purple", pch=16, type="b", xlab="SPOT Sampling Date", ylab="Proportion of 16S sequences \n in 1-80um size class", main= "UCYN-A ASV3 AE 500m (6115eab19c52bc45c6ba11d72ec88031)")
par(op)

#Durapore from 500m depth, color= blue
op <- par(mar=c(5, 6, 4, 2) + 0.1)
plot(as.Date(all_dates[(grep("500m.D", colnames(proks_data))+1)]),as.numeric(proks_data[54043, grep("500m.D", colnames(proks_data))]), lwd=2, col="purple", pch=4, type="b", xlab="SPOT Sampling Date", ylab="Proportion of 16S sequences \n in 0.22-1um size class", main= "UCYN-A ASV3 Dura. 500m (6115eab19c52bc45c6ba11d72ec88031)")
par(op)

#AE from 890m depth, color= purple
op <- par(mar=c(5, 6, 4, 2) + 0.1)
plot(as.Date(all_dates[(grep("890m.AE", colnames(proks_data))+1)]),as.numeric(proks_data[54043, grep("890m.AE", colnames(proks_data))]), pch=16, type="b", xlab="SPOT Sampling Date", ylab="Proportion of 16S sequences \n in 1-80um size class", main= "UCYN-A ASV3 AE 890m (6115eab19c52bc45c6ba11d72ec88031)")
par(op)

#Durapore from 890m depth, color= black
op <- par(mar=c(5, 6, 4, 2) + 0.1)
plot(as.Date(all_dates[(grep("890m.D", colnames(proks_data))+1)]),as.numeric(proks_data[54043, grep("890m.D", colnames(proks_data))]), lwd=2, pch=4, type="b", xlab="SPOT Sampling Date", ylab="Proportion of 16S sequences \n in 0.22-1um size class", main= "UCYN-A ASV3 Dura. 890m (6115eab19c52bc45c6ba11d72ec88031)")
par(op)

#UCYNA ASV4: Row 3880 in Proks data, a641110da9fb0da8f68143b5a79ba5d1--------
#AE from 5m depth first, color= red
op <- par(mar=c(5, 6, 4, 2) + 0.1)
plot(as.Date(all_dates[(grep("5m.AE", colnames(proks_data))+1)]),as.numeric(proks_data[3880, grep("5m.AE", colnames(proks_data))]), col="red", pch=16, type="b", xlab="SPOT Sampling Date", ylab="Proportion of 16S sequences \n in 1-80um size class", main= "UCYN-A ASV4 AE 5m (a641110da9fb0da8f68143b5a79ba5d1)")
par(op)

#Durapore from 5m depth, color=red, pch= something else
op <- par(mar=c(5, 6, 4, 2) + 0.1)
plot(as.Date(all_dates[(grep("5m.D", colnames(proks_data))+1)]),as.numeric(proks_data[3880, grep("5m.D", colnames(proks_data))]), col="red", pch=4, type="b", lwd=2, xlab="SPOT Sampling Date", ylab="Proportion of 16S sequences \n in 0.22-1um size class", main= "UCYN-A ASV4 Dura. 5m (a641110da9fb0da8f68143b5a79ba5d1)")
par(op)

#AE from DCM depth, color= green
op <- par(mar=c(5, 6, 4, 2) + 0.1)
plot(as.Date(all_dates[(grep("DCM.AE", colnames(proks_data))+1)]),as.numeric(proks_data[3880, grep("DCM.AE", colnames(proks_data))]), col="forestgreen", pch=16, type="b", xlab="SPOT Sampling Date", ylab="Proportion of 16S sequences \n in 1-80um size class", main= "UCYN-A ASV4 AE DCM (a641110da9fb0da8f68143b5a79ba5d1)")
par(op)

#Durapore from DCM depth, color= green
op <- par(mar=c(5, 6, 4, 2) + 0.1)
plot(as.Date(all_dates[(grep("DCM.D", colnames(proks_data))+1)]),as.numeric(proks_data[3880, grep("DCM.D", colnames(proks_data))]), lwd=2, col="forestgreen", pch=4, type="b", xlab="SPOT Sampling Date", ylab="Proportion of 16S sequences \n in 0.22-1um size class", main= "UCYN-A ASV4 Dura. DCM (a641110da9fb0da8f68143b5a79ba5d1)")
par(op)

#AE from 150m depth, color= blue
op <- par(mar=c(5, 6, 4, 2) + 0.1)
plot(as.Date(all_dates[(grep("150m.AE", colnames(proks_data))+1)]),as.numeric(proks_data[3880, grep("150m.AE", colnames(proks_data))]), col="blue", pch=16, type="b", xlab="SPOT Sampling Date", ylab="Proportion of 16S sequences \n in 1-80um size class", main= "UCYN-A ASV4 AE 150m (a641110da9fb0da8f68143b5a79ba5d1)")
par(op)

#Durapore from 150m depth, color= blue
op <- par(mar=c(5, 6, 4, 2) + 0.1)
plot(as.Date(all_dates[(grep("150m.D", colnames(proks_data))+1)]),as.numeric(proks_data[3880, grep("150m.D", colnames(proks_data))]), lwd=2, col="blue", pch=4, type="b", xlab="SPOT Sampling Date", ylab="Proportion of 16S sequences \n in 0.22-1um size class", main= "UCYN-A ASV4 Dura. 150m (a641110da9fb0da8f68143b5a79ba5d1)")
par(op)

#AE from 500m depth, color= purple
op <- par(mar=c(5, 6, 4, 2) + 0.1)
plot(as.Date(all_dates[(grep("500m.AE", colnames(proks_data))+1)]),as.numeric(proks_data[3880, grep("500m.AE", colnames(proks_data))]), col="purple", pch=16, type="b", xlab="SPOT Sampling Date", ylab="Proportion of 16S sequences \n in 1-80um size class", main= "UCYN-A ASV4 AE 500m (a641110da9fb0da8f68143b5a79ba5d1)")
par(op)

#Durapore from 500m depth, color= blue
op <- par(mar=c(5, 6, 4, 2) + 0.1)
plot(as.Date(all_dates[(grep("500m.D", colnames(proks_data))+1)]),as.numeric(proks_data[3880, grep("500m.D", colnames(proks_data))]), lwd=2, col="purple", pch=4, type="b", xlab="SPOT Sampling Date", ylab="Proportion of 16S sequences \n in 0.22-1um size class", main= "UCYN-A ASV4 Dura. 500m (a641110da9fb0da8f68143b5a79ba5d1)")
par(op)

#AE from 890m depth, color= purple
op <- par(mar=c(5, 6, 4, 2) + 0.1)
plot(as.Date(all_dates[(grep("890m.AE", colnames(proks_data))+1)]),as.numeric(proks_data[3880, grep("890m.AE", colnames(proks_data))]), pch=16, type="b", xlab="SPOT Sampling Date", ylab="Proportion of 16S sequences \n in 1-80um size class", main= "UCYN-A ASV4 AE 890m (a641110da9fb0da8f68143b5a79ba5d1)")
par(op)

#Durapore from 890m depth, color= black
op <- par(mar=c(5, 6, 4, 2) + 0.1)
plot(as.Date(all_dates[(grep("890m.D", colnames(proks_data))+1)]),as.numeric(proks_data[3880, grep("890m.D", colnames(proks_data))]), lwd=2, pch=4, type="b", xlab="SPOT Sampling Date", ylab="Proportion of 16S sequences \n in 0.22-1um size class", main= "UCYN-A ASV4 Dura. 890m (a641110da9fb0da8f68143b5a79ba5d1)")
par(op)

#UCYNA ASV5: Row 1820 in Proks data, af1bb1f9fb1c3f3d18571e711df407bb--------
#AE from 5m depth first, color= red
op <- par(mar=c(5, 6, 4, 2) + 0.1)
plot(as.Date(all_dates[(grep("5m.AE", colnames(proks_data))+1)]),as.numeric(proks_data[1820, grep("5m.AE", colnames(proks_data))]), col="red", pch=16, type="b", xlab="SPOT Sampling Date", ylab="Proportion of 16S sequences \n in 1-80um size class", main= "UCYN-A ASV5 AE 5m (af1bb1f9fb1c3f3d18571e711df407bb)")
par(op)

#Durapore from 5m depth, color=red, pch= something else
op <- par(mar=c(5, 6, 4, 2) + 0.1)
plot(as.Date(all_dates[(grep("5m.D", colnames(proks_data))+1)]),as.numeric(proks_data[1820, grep("5m.D", colnames(proks_data))]), col="red", pch=4, type="b", lwd=2, xlab="SPOT Sampling Date", ylab="Proportion of 16S sequences \n in 0.22-1um size class", main= "UCYN-A ASV5 Dura. 5m (af1bb1f9fb1c3f3d18571e711df407bb)")
par(op)

#AE from DCM depth, color= green
op <- par(mar=c(5, 6, 4, 2) + 0.1)
plot(as.Date(all_dates[(grep("DCM.AE", colnames(proks_data))+1)]),as.numeric(proks_data[1820, grep("DCM.AE", colnames(proks_data))]), col="forestgreen", pch=16, type="b", xlab="SPOT Sampling Date", ylab="Proportion of 16S sequences \n in 1-80um size class", main= "UCYN-A ASV5 AE DCM (af1bb1f9fb1c3f3d18571e711df407bb)")
par(op)

#Durapore from DCM depth, color= green
op <- par(mar=c(5, 6, 4, 2) + 0.1)
plot(as.Date(all_dates[(grep("DCM.D", colnames(proks_data))+1)]),as.numeric(proks_data[1820, grep("DCM.D", colnames(proks_data))]), lwd=2, col="forestgreen", pch=4, type="b", xlab="SPOT Sampling Date", ylab="Proportion of 16S sequences \n in 0.22-1um size class", main= "UCYN-A ASV5 Dura. DCM (af1bb1f9fb1c3f3d18571e711df407bb)")
par(op)

#AE from 150m depth, color= blue
op <- par(mar=c(5, 6, 4, 2) + 0.1)
plot(as.Date(all_dates[(grep("150m.AE", colnames(proks_data))+1)]),as.numeric(proks_data[1820, grep("150m.AE", colnames(proks_data))]), col="blue", pch=16, type="b", xlab="SPOT Sampling Date", ylab="Proportion of 16S sequences \n in 1-80um size class", main= "UCYN-A ASV5 AE 150m (af1bb1f9fb1c3f3d18571e711df407bb)")
par(op)

#Durapore from 150m depth, color= blue
op <- par(mar=c(5, 6, 4, 2) + 0.1)
plot(as.Date(all_dates[(grep("150m.D", colnames(proks_data))+1)]),as.numeric(proks_data[1820, grep("150m.D", colnames(proks_data))]), lwd=2, col="blue", pch=4, type="b", xlab="SPOT Sampling Date", ylab="Proportion of 16S sequences \n in 0.22-1um size class", main= "UCYN-A ASV5 Dura. 150m (af1bb1f9fb1c3f3d18571e711df407bb)")
par(op)

#AE from 500m depth, color= purple
op <- par(mar=c(5, 6, 4, 2) + 0.1)
plot(as.Date(all_dates[(grep("500m.AE", colnames(proks_data))+1)]),as.numeric(proks_data[1820, grep("500m.AE", colnames(proks_data))]), col="purple", pch=16, type="b", xlab="SPOT Sampling Date", ylab="Proportion of 16S sequences \n in 1-80um size class", main= "UCYN-A ASV5 AE 500m (af1bb1f9fb1c3f3d18571e711df407bb)")
par(op)

#Durapore from 500m depth, color= blue
op <- par(mar=c(5, 6, 4, 2) + 0.1)
plot(as.Date(all_dates[(grep("500m.D", colnames(proks_data))+1)]),as.numeric(proks_data[1820, grep("500m.D", colnames(proks_data))]), lwd=2, col="purple", pch=4, type="b", xlab="SPOT Sampling Date", ylab="Proportion of 16S sequences \n in 0.22-1um size class", main= "UCYN-A ASV5 Dura. 500m (af1bb1f9fb1c3f3d18571e711df407bb)")
par(op)

#AE from 890m depth, color= purple
op <- par(mar=c(5, 6, 4, 2) + 0.1)
plot(as.Date(all_dates[(grep("890m.AE", colnames(proks_data))+1)]),as.numeric(proks_data[1820, grep("890m.AE", colnames(proks_data))]), pch=16, type="b", xlab="SPOT Sampling Date", ylab="Proportion of 16S sequences \n in 1-80um size class", main= "UCYN-A ASV5 AE 890m (af1bb1f9fb1c3f3d18571e711df407bb)")
par(op)

#Durapore from 890m depth, color= black
op <- par(mar=c(5, 6, 4, 2) + 0.1)
plot(as.Date(all_dates[(grep("890m.D", colnames(proks_data))+1)]),as.numeric(proks_data[1820, grep("890m.D", colnames(proks_data))]), lwd=2, pch=4, type="b", xlab="SPOT Sampling Date", ylab="Proportion of 16S sequences \n in 0.22-1um size class", main= "UCYN-A ASV5 Dura. 890m (af1bb1f9fb1c3f3d18571e711df407bb)")
par(op)


#UCYNA ASV6: Row 5286 in Proks data, e6f42c535cf3849e1f1e12e7575561b7--------
#AE from 5m depth first, color= red
op <- par(mar=c(5, 6, 4, 2) + 0.1)
plot(as.Date(all_dates[(grep("5m.AE", colnames(proks_data))+1)]),as.numeric(proks_data[5286, grep("5m.AE", colnames(proks_data))]), col="red", pch=16, type="b", xlab="SPOT Sampling Date", ylab="Proportion of 16S sequences \n in 1-80um size class", main= "UCYN-A ASV6 AE 5m (e6f42c535cf3849e1f1e12e7575561b7)")
par(op)

#Durapore from 5m depth, color=red, pch= something else
op <- par(mar=c(5, 6, 4, 2) + 0.1)
plot(as.Date(all_dates[(grep("5m.D", colnames(proks_data))+1)]),as.numeric(proks_data[5286, grep("5m.D", colnames(proks_data))]), col="red", pch=4, type="b", lwd=2, xlab="SPOT Sampling Date", ylab="Proportion of 16S sequences \n in 0.22-1um size class", main= "UCYN-A ASV6 Dura. 5m (e6f42c535cf3849e1f1e12e7575561b7)")
par(op)

#AE from DCM depth, color= green
op <- par(mar=c(5, 6, 4, 2) + 0.1)
plot(as.Date(all_dates[(grep("DCM.AE", colnames(proks_data))+1)]),as.numeric(proks_data[5286, grep("DCM.AE", colnames(proks_data))]), col="forestgreen", pch=16, type="b", xlab="SPOT Sampling Date", ylab="Proportion of 16S sequences \n in 1-80um size class", main= "UCYN-A ASV6 AE DCM (e6f42c535cf3849e1f1e12e7575561b7)")
par(op)

#Durapore from DCM depth, color= green
op <- par(mar=c(5, 6, 4, 2) + 0.1)
plot(as.Date(all_dates[(grep("DCM.D", colnames(proks_data))+1)]),as.numeric(proks_data[5286, grep("DCM.D", colnames(proks_data))]), lwd=2, col="forestgreen", pch=4, type="b", xlab="SPOT Sampling Date", ylab="Proportion of 16S sequences \n in 0.22-1um size class", main= "UCYN-A ASV6 Dura. DCM (e6f42c535cf3849e1f1e12e7575561b7)")
par(op)

#AE from 150m depth, color= blue
op <- par(mar=c(5, 6, 4, 2) + 0.1)
plot(as.Date(all_dates[(grep("150m.AE", colnames(proks_data))+1)]),as.numeric(proks_data[5286, grep("150m.AE", colnames(proks_data))]), col="blue", pch=16, type="b", xlab="SPOT Sampling Date", ylab="Proportion of 16S sequences \n in 1-80um size class", main= "UCYN-A ASV6 AE 150m (e6f42c535cf3849e1f1e12e7575561b7)")
par(op)

#Durapore from 150m depth, color= blue
op <- par(mar=c(5, 6, 4, 2) + 0.1)
plot(as.Date(all_dates[(grep("150m.D", colnames(proks_data))+1)]),as.numeric(proks_data[5286, grep("150m.D", colnames(proks_data))]), lwd=2, col="blue", pch=4, type="b", xlab="SPOT Sampling Date", ylab="Proportion of 16S sequences \n in 0.22-1um size class", main= "UCYN-A ASV6 Dura. 150m (e6f42c535cf3849e1f1e12e7575561b7)")
par(op)

#AE from 500m depth, color= purple
op <- par(mar=c(5, 6, 4, 2) + 0.1)
plot(as.Date(all_dates[(grep("500m.AE", colnames(proks_data))+1)]),as.numeric(proks_data[5286, grep("500m.AE", colnames(proks_data))]), col="purple", pch=16, type="b", xlab="SPOT Sampling Date", ylab="Proportion of 16S sequences \n in 1-80um size class", main= "UCYN-A ASV6 AE 500m (e6f42c535cf3849e1f1e12e7575561b7)")
par(op)

#Durapore from 500m depth, color= blue
op <- par(mar=c(5, 6, 4, 2) + 0.1)
plot(as.Date(all_dates[(grep("500m.D", colnames(proks_data))+1)]),as.numeric(proks_data[5286, grep("500m.D", colnames(proks_data))]), lwd=2, col="purple", pch=4, type="b", xlab="SPOT Sampling Date", ylab="Proportion of 16S sequences \n in 0.22-1um size class", main= "UCYN-A ASV6 Dura. 500m (e6f42c535cf3849e1f1e12e7575561b7)")
par(op)

#AE from 890m depth, color= purple
op <- par(mar=c(5, 6, 4, 2) + 0.1)
plot(as.Date(all_dates[(grep("890m.AE", colnames(proks_data))+1)]),as.numeric(proks_data[5286, grep("890m.AE", colnames(proks_data))]), pch=16, type="b", xlab="SPOT Sampling Date", ylab="Proportion of 16S sequences \n in 1-80um size class", main= "UCYN-A ASV6 AE 890m (e6f42c535cf3849e1f1e12e7575561b7)")
par(op)

#Durapore from 890m depth, color= black
op <- par(mar=c(5, 6, 4, 2) + 0.1)
plot(as.Date(all_dates[(grep("890m.D", colnames(proks_data))+1)]),as.numeric(proks_data[5286, grep("890m.D", colnames(proks_data))]), lwd=2, pch=4, type="b", xlab="SPOT Sampling Date", ylab="Proportion of 16S sequences \n in 0.22-1um size class", main= "UCYN-A ASV6 Dura. 890m (e6f42c535cf3849e1f1e12e7575561b7)")
par(op)

#Stop sending graphs to .PDF
dev.off()

#At which dates does UCYNA ASV1 appear at 890m depth??------ 
#I must know! Is it the same as Hammerlsey's dates?
which(proks_data[223, grep("890m.AE", colnames(proks_data))]>0)
colnames(proks_data)[grep("890m.AE", colnames(proks_data))][which(proks_data[223, grep("890m.AE", colnames(proks_data))]>0)]
#"SPOT.2008.11.20.890m.AE" "SPOT.2009.07.09.890m.AE"
proks_data[223, grep("SPOT.2008.11.20.890m.AE|SPOT.2009.07.09.890m.AE", colnames(proks_data))] #Okay the 2008 date doesn't really count because abundance was 10^-5
#So it was 7.9.2009
#Hamersley didn't sample in 2009, darn it :( 

#Add relative abundances and dates of UCYN-A ASVs 1, 5, 6 to a table------

#Figure out number rows
length(which(as.numeric(proks_data[223,grep("5m.AE", colnames(proks_data))])>0.01))  #4
length(which(as.numeric(proks_data[223,grep("DCM.AE", colnames(proks_data))])>0.0005)) #4
length(which(as.numeric(proks_data[1820,grep("5m.AE", colnames(proks_data))])>0.002)) #6 
length(which(as.numeric(proks_data[1820,grep("DCM.AE", colnames(proks_data))])>0)) #5
length(which(as.numeric(proks_data[5286,grep("5m.AE", colnames(proks_data))])>0.0005)) #5
length(which(as.numeric(proks_data[5286,grep("DCM.AE", colnames(proks_data))])>0)) #5
4+4+6+5+5+5

#Set up table
dates_for_eLSA <- as.data.frame(matrix(nrow=29, ncol=6))
names(dates_for_eLSA)=c("SampleID", "Relabun_16S", "ASV_No", "ASV_Hash", "RowNum.", "Depth")
 
#1. Add in UCYNA ASV1, 5m
dates_for_eLSA$SampleID[1:4]=colnames(proks_data[223,grep("5m.AE", colnames(proks_data))])[which(as.numeric(proks_data[223,grep("5m.AE", colnames(proks_data))])>0.01)]

dates_for_eLSA$Relabun_16S[1:4]=as.numeric(proks_data[223,grep("5m.AE", colnames(proks_data))])[which(as.numeric(proks_data[223,grep("5m.AE", colnames(proks_data))])>0.01)]

dates_for_eLSA$ASV_No[1:4]=1
dates_for_eLSA$ASV_Hash[1:4]=UCYNA_ASVs[1]
dates_for_eLSA$RowNum.[1:4]=UCYNA_ix[1]
dates_for_eLSA$Depth[1:4]="5m"
View(dates_for_eLSA)

#2. UCYNA ASV1, DCM
dates_for_eLSA$SampleID[5:8]=colnames(proks_data[223,grep("DCM.AE", colnames(proks_data))])[which(as.numeric(proks_data[223,grep("DCM.AE", colnames(proks_data))])>0.0005)]

dates_for_eLSA$Relabun_16S[5:8]=as.numeric(proks_data[223,grep("DCM.AE", colnames(proks_data))])[which(as.numeric(proks_data[223,grep("DCM.AE", colnames(proks_data))])>0.0005)]

dates_for_eLSA$ASV_No[5:8]=1
dates_for_eLSA$ASV_Hash[5:8]=UCYNA_ASVs[1]
dates_for_eLSA$RowNum.[5:8]=UCYNA_ix[1]
dates_for_eLSA$Depth[5:8]="DCM"

#3. Add in UCYNA ASV5, 5m
dates_for_eLSA$SampleID[9:14]=colnames(proks_data[1820,grep("5m.AE", colnames(proks_data))])[which(as.numeric(proks_data[1820,grep("5m.AE", colnames(proks_data))])>0.002)]

dates_for_eLSA$Relabun_16S[9:14]=as.numeric(proks_data[1820,grep("5m.AE", colnames(proks_data))])[which(as.numeric(proks_data[1820,grep("5m.AE", colnames(proks_data))])>0.002)]

dates_for_eLSA$ASV_No[9:14]=5
dates_for_eLSA$ASV_Hash[9:14]=UCYNA_ASVs[5]
dates_for_eLSA$RowNum.[9:14]=UCYNA_ix[5]
dates_for_eLSA$Depth[9:14]="5m"
View(dates_for_eLSA)

#4. Add in UCYNA ASV5, DCM
dates_for_eLSA$SampleID[15:19]=colnames(proks_data[1820,grep("DCM.AE", colnames(proks_data))])[which(as.numeric(proks_data[1820,grep("DCM.AE", colnames(proks_data))])>0)]

dates_for_eLSA$Relabun_16S[15:19]=as.numeric(proks_data[1820,grep("DCM.AE", colnames(proks_data))])[which(as.numeric(proks_data[1820,grep("DCM.AE", colnames(proks_data))])>0)]

dates_for_eLSA$ASV_No[15:19]=5
dates_for_eLSA$ASV_Hash[15:19]=UCYNA_ASVs[5]
dates_for_eLSA$RowNum.[15:19]=UCYNA_ix[5]
dates_for_eLSA$Depth[15:19]="DCM"
View(dates_for_eLSA)

#5. Add in UCYNA ASV6, 5m
dates_for_eLSA$SampleID[20:23]=colnames(proks_data[5286,grep("5m.AE", colnames(proks_data))])[which(as.numeric(proks_data[5286,grep("5m.AE", colnames(proks_data))])>0.0005)]

dates_for_eLSA$Relabun_16S[20:23]=as.numeric(proks_data[5286,grep("5m.AE", colnames(proks_data))])[which(as.numeric(proks_data[5286,grep("5m.AE", colnames(proks_data))])>0.0005)]

dates_for_eLSA$ASV_No[20:23]=6
dates_for_eLSA$ASV_Hash[20:23]=UCYNA_ASVs[6]
dates_for_eLSA$RowNum.[20:23]=UCYNA_ix[6]
dates_for_eLSA$Depth[20:23]="5m"
View(dates_for_eLSA)

#6. Add in UCYNA ASV6, DCM
dates_for_eLSA$SampleID[24:28]=colnames(proks_data[5286,grep("DCM.AE", colnames(proks_data))])[which(as.numeric(proks_data[5286,grep("DCM.AE", colnames(proks_data))])>0)]

dates_for_eLSA$Relabun_16S[24:28]=as.numeric(proks_data[5286,grep("DCM.AE", colnames(proks_data))])[which(as.numeric(proks_data[5286,grep("DCM.AE", colnames(proks_data))])>0)]

dates_for_eLSA$ASV_No[24:28]=6
dates_for_eLSA$ASV_Hash[24:28]=UCYNA_ASVs[6]
dates_for_eLSA$RowNum.[24:28]=UCYNA_ix[6]
dates_for_eLSA$Depth[24:28]="DCM"
View(dates_for_eLSA)

#Write out the file! :) 
getwd()
setwd("/Users/cfletcherh1/Desktop/DataAnalysis/SPOT_16S_18S_eLSA/")
write.csv(x=dates_for_eLSA, file="SampleIDs_w.high_UCYN-A_relabun_for_eLSA_08182020.csv", row.names=F, quote=F)

length(dates_for_eLSA$SampleID) #29
length(unique(dates_for_eLSA$SampleID)) #27

#Plot up ASVs of genus Brad over time, save it in a .PDF------
#All 7 ASVs that correspond to the genus Braarudosphaera
#ASV numbers correspond to table written out on 8.3.20
#IS BRAD PRESENT IN THE DURAPORE SAMPLES?? 
Brad_ix
for(i in genus_brad_ix){
  print(which(euks_data[i, grep("5m.D", colnames(euks_data))]>0))
}
#Shit, 6x and 7 show up twice in durapore samples
for(i in genus_brad_ix){
  print(length(which(euks_data[i, grep("D.DCM", colnames(euks_data))]>0)))
}
#Also all zero, phew! 

#Set up .PDF
getwd()
setwd("/Users/cfletcherh1/Desktop/DataAnalysis/SPOT_16S_18S_eLSA/")
pdf("Seven_Braarudosphaera_ASVs_over_time_all_depths_08192020.pdf")

#Brad bigelowii ASV1: Row 8179 in Euks data, 04926e2fd1b8706b4866c02650f702dd-----
#AE from 5m depth first, color= red
op <- par(mar=c(5, 6, 4, 2) + 0.1)
plot(as.Date(euks_dates[grep("5m.AE", colnames(euks_data))+1]), as.numeric(euks_data[genus_brad_ix[1], grep("5m.AE", colnames(euks_data))]), type="b", col="red", pch=16, xlab="SPOT Sampling Date", ylab="Proportion of 18S sequences \n in 1-80um size class", main="Braarudosphaera bigelowii ASV1 \n AE 5m (04926e2fd1b8706b4866c02650f702dd)")
par(op)

#AE from DCM next, color=green
op <- par(mar=c(5, 6, 4, 2) + 0.1)
plot(as.Date(euks_dates[grep("DCM.AE", colnames(euks_data))+1]), as.numeric(euks_data[genus_brad_ix[1], grep("DCM.AE", colnames(euks_data))]), type="b", col="forestgreen", pch=16, xlab="SPOT Sampling Date", ylab="Proportion of 18S sequences \n in 1-80um size class", main="Braarudosphaera bigelowii ASV1 \n AE DCM (04926e2fd1b8706b4866c02650f702dd)")
par(op)

#AE from 150m, color=blue
op <- par(mar=c(5, 6, 4, 2) + 0.1)
plot(as.Date(euks_dates[grep("150m.AE", colnames(euks_data))+1]), as.numeric(euks_data[genus_brad_ix[1], grep("150m.AE", colnames(euks_data))]), type="b", col="blue", pch=16, xlab="SPOT Sampling Date", ylab="Proportion of 18S sequences \n in 1-80um size class", main="Braarudosphaera bigelowii ASV1 \n AE 150m (04926e2fd1b8706b4866c02650f702dd)")
par(op)

#AE from 500m, col=purple
op <- par(mar=c(5, 6, 4, 2) + 0.1)
plot(as.Date(euks_dates[grep("500m.AE", colnames(euks_data))+1]), as.numeric(euks_data[genus_brad_ix[1], grep("500m.AE", colnames(euks_data))]), type="b", col="purple", pch=16, xlab="SPOT Sampling Date", ylab="Proportion of 18S sequences \n in 1-80um size class", main="Braarudosphaera bigelowii ASV1 \n AE 500m (04926e2fd1b8706b4866c02650f702dd)")
par(op)

#AE from 890m, no color 
op <- par(mar=c(5, 6, 4, 2) + 0.1)
plot(as.Date(euks_dates[grep("890m.AE", colnames(euks_data))+1]), as.numeric(euks_data[genus_brad_ix[1], grep("890m.AE", colnames(euks_data))]), type="b", pch=16, xlab="SPOT Sampling Date", ylab="Proportion of 18S sequences \n in 1-80um size class", main="Braarudosphaera bigelowii ASV1 \n AE 890m (04926e2fd1b8706b4866c02650f702dd)")
par(op)

#Brad genus ASV5x: Row 1228 in Euks data, 324627f7f367298bbb5692fc5038e680-----
#AE from 5m depth first, color= red
op <- par(mar=c(5, 6, 4, 2) + 0.1)
plot(as.Date(euks_dates[grep("5m.AE", colnames(euks_data))+1]), as.numeric(euks_data[genus_brad_ix[2], grep("5m.AE", colnames(euks_data))]), type="b", col="red", pch=16, xlab="SPOT Sampling Date", ylab="Proportion of 18S sequences \n in 1-80um size class", main="Braarudosphaera sp. ASV5x \n AE 5m (324627f7f367298bbb5692fc5038e680)")
par(op)

#AE from DCM next, color=green
op <- par(mar=c(5, 6, 4, 2) + 0.1)
plot(as.Date(euks_dates[grep("DCM.AE", colnames(euks_data))+1]), as.numeric(euks_data[genus_brad_ix[2], grep("DCM.AE", colnames(euks_data))]), type="b", col="forestgreen", pch=16, xlab="SPOT Sampling Date", ylab="Proportion of 18S sequences \n in 1-80um size class", main="Braarudosphaera sp. ASV5x \n AE DCM (324627f7f367298bbb5692fc5038e680)")
par(op)

#AE from 150m, color=blue
op <- par(mar=c(5, 6, 4, 2) + 0.1)
plot(as.Date(euks_dates[grep("150m.AE", colnames(euks_data))+1]), as.numeric(euks_data[genus_brad_ix[2], grep("150m.AE", colnames(euks_data))]), type="b", col="blue", pch=16, xlab="SPOT Sampling Date", ylab="Proportion of 18S sequences \n in 1-80um size class", main="Braarudosphaera sp. ASV5x \n AE 150m (324627f7f367298bbb5692fc5038e680)")
par(op)

#AE from 500m, col=purple
op <- par(mar=c(5, 6, 4, 2) + 0.1)
plot(as.Date(euks_dates[grep("500m.AE", colnames(euks_data))+1]), as.numeric(euks_data[genus_brad_ix[2], grep("500m.AE", colnames(euks_data))]), type="b", col="purple", pch=16, xlab="SPOT Sampling Date", ylab="Proportion of 18S sequences \n in 1-80um size class", main="Braarudosphaera sp. ASV5x \n AE 500m (324627f7f367298bbb5692fc5038e680)")
par(op)

#AE from 890m, no color 
op <- par(mar=c(5, 6, 4, 2) + 0.1)
plot(as.Date(euks_dates[grep("890m.AE", colnames(euks_data))+1]), as.numeric(euks_data[genus_brad_ix[2], grep("890m.AE", colnames(euks_data))]), type="b", pch=16, xlab="SPOT Sampling Date", ylab="Proportion of 18S sequences \n in 1-80um size class", main="Braarudosphaera sp. ASV5x \n AE 890m (324627f7f367298bbb5692fc5038e680)")
par(op)

#Brad bigelowii ASV2: Row 21645 in Euks data, 529269deeb5fb7fbf0d0ebda989d9d82-----
#AE from 5m depth first, color= red
op <- par(mar=c(5, 6, 4, 2) + 0.1)
plot(as.Date(euks_dates[grep("5m.AE", colnames(euks_data))+1]), as.numeric(euks_data[genus_brad_ix[3], grep("5m.AE", colnames(euks_data))]), type="b", col="red", pch=16, xlab="SPOT Sampling Date", ylab="Proportion of 18S sequences \n in 1-80um size class", main="Braarudosphaera bigelowii ASV2 \n AE 5m (529269deeb5fb7fbf0d0ebda989d9d82)")
par(op)

#AE from DCM next, color=green
op <- par(mar=c(5, 6, 4, 2) + 0.1)
plot(as.Date(euks_dates[grep("DCM.AE", colnames(euks_data))+1]), as.numeric(euks_data[genus_brad_ix[3], grep("DCM.AE", colnames(euks_data))]), type="b", col="forestgreen", pch=16, xlab="SPOT Sampling Date", ylab="Proportion of 18S sequences \n in 1-80um size class", main="Braarudosphaera bigelowii ASV2 \n AE DCM (529269deeb5fb7fbf0d0ebda989d9d82)")
par(op)

#AE from 150m, color=blue
op <- par(mar=c(5, 6, 4, 2) + 0.1)
plot(as.Date(euks_dates[grep("150m.AE", colnames(euks_data))+1]), as.numeric(euks_data[genus_brad_ix[3], grep("150m.AE", colnames(euks_data))]), type="b", col="blue", pch=16, xlab="SPOT Sampling Date", ylab="Proportion of 18S sequences \n in 1-80um size class", main="Braarudosphaera bigelowii ASV2 \n AE 150m (529269deeb5fb7fbf0d0ebda989d9d82)")
par(op)

#AE from 500m, col=purple
op <- par(mar=c(5, 6, 4, 2) + 0.1)
plot(as.Date(euks_dates[grep("500m.AE", colnames(euks_data))+1]), as.numeric(euks_data[genus_brad_ix[3], grep("500m.AE", colnames(euks_data))]), type="b", col="purple", pch=16, xlab="SPOT Sampling Date", ylab="Proportion of 18S sequences \n in 1-80um size class", main="Braarudosphaera bigelowii ASV2 \n AE 500m (529269deeb5fb7fbf0d0ebda989d9d82)")
par(op)

#AE from 890m, no color 
op <- par(mar=c(5, 6, 4, 2) + 0.1)
plot(as.Date(euks_dates[grep("890m.AE", colnames(euks_data))+1]), as.numeric(euks_data[genus_brad_ix[3], grep("890m.AE", colnames(euks_data))]), type="b", pch=16, xlab="SPOT Sampling Date", ylab="Proportion of 18S sequences \n in 1-80um size class", main="Braarudosphaera bigelowii ASV2 \n AE 890m (529269deeb5fb7fbf0d0ebda989d9d82)")
par(op)

#Brad bigelowii ASV3: Row 991 in Euks data, 70a5283da28db501a349c5beb22881e7-----
#AE from 5m depth first, color= red
op <- par(mar=c(5, 6, 4, 2) + 0.1)
plot(as.Date(euks_dates[grep("5m.AE", colnames(euks_data))+1]), as.numeric(euks_data[genus_brad_ix[4], grep("5m.AE", colnames(euks_data))]), type="b", col="red", pch=16, xlab="SPOT Sampling Date", ylab="Proportion of 18S sequences \n in 1-80um size class", main="Braarudosphaera bigelowii ASV3 \n AE 5m (70a5283da28db501a349c5beb22881e7)")
par(op)

#AE from DCM next, color=green
op <- par(mar=c(5, 6, 4, 2) + 0.1)
plot(as.Date(euks_dates[grep("DCM.AE", colnames(euks_data))+1]), as.numeric(euks_data[genus_brad_ix[4], grep("DCM.AE", colnames(euks_data))]), type="b", col="forestgreen", pch=16, xlab="SPOT Sampling Date", ylab="Proportion of 18S sequences \n in 1-80um size class", main="Braarudosphaera bigelowii ASV3 \n AE DCM (70a5283da28db501a349c5beb22881e7)")
par(op)

#AE from 150m, color=blue
op <- par(mar=c(5, 6, 4, 2) + 0.1)
plot(as.Date(euks_dates[grep("150m.AE", colnames(euks_data))+1]), as.numeric(euks_data[genus_brad_ix[4], grep("150m.AE", colnames(euks_data))]), type="b", col="blue", pch=16, xlab="SPOT Sampling Date", ylab="Proportion of 18S sequences \n in 1-80um size class", main="Braarudosphaera bigelowii ASV3 \n AE 150m (70a5283da28db501a349c5beb22881e7)")
par(op)

#AE from 500m, col=purple
op <- par(mar=c(5, 6, 4, 2) + 0.1)
plot(as.Date(euks_dates[grep("500m.AE", colnames(euks_data))+1]), as.numeric(euks_data[genus_brad_ix[4], grep("500m.AE", colnames(euks_data))]), type="b", col="purple", pch=16, xlab="SPOT Sampling Date", ylab="Proportion of 18S sequences \n in 1-80um size class", main="Braarudosphaera bigelowii ASV3 \n AE 500m (70a5283da28db501a349c5beb22881e7)")
par(op)

#AE from 890m, no color 
op <- par(mar=c(5, 6, 4, 2) + 0.1)
plot(as.Date(euks_dates[grep("890m.AE", colnames(euks_data))+1]), as.numeric(euks_data[genus_brad_ix[4], grep("890m.AE", colnames(euks_data))]), type="b", pch=16, xlab="SPOT Sampling Date", ylab="Proportion of 18S sequences \n in 1-80um size class", main="Braarudosphaera bigelowii ASV3 \n AE 890m (70a5283da28db501a349c5beb22881e7)")
par(op)

#Brad bigelowii ASV4: Row 984 in Euks data, 8c144683114fbb1ad2e9425f7dcd1b02-----
#AE from 5m depth first, color= red
op <- par(mar=c(5, 6, 4, 2) + 0.1)
plot(as.Date(euks_dates[grep("5m.AE", colnames(euks_data))+1]), as.numeric(euks_data[genus_brad_ix[5], grep("5m.AE", colnames(euks_data))]), type="b", col="red", pch=16, xlab="SPOT Sampling Date", ylab="Proportion of 18S sequences \n in 1-80um size class", main="Braarudosphaera bigelowii ASV4 \n AE 5m (8c144683114fbb1ad2e9425f7dcd1b02)")
par(op)

#AE from DCM next, color=green
op <- par(mar=c(5, 6, 4, 2) + 0.1)
plot(as.Date(euks_dates[grep("DCM.AE", colnames(euks_data))+1]), as.numeric(euks_data[genus_brad_ix[5], grep("DCM.AE", colnames(euks_data))]), type="b", col="forestgreen", pch=16, xlab="SPOT Sampling Date", ylab="Proportion of 18S sequences \n in 1-80um size class", main="Braarudosphaera bigelowii ASV4 \n AE DCM (8c144683114fbb1ad2e9425f7dcd1b02)")
par(op)

#AE from 150m, color=blue
op <- par(mar=c(5, 6, 4, 2) + 0.1)
plot(as.Date(euks_dates[grep("150m.AE", colnames(euks_data))+1]), as.numeric(euks_data[genus_brad_ix[5], grep("150m.AE", colnames(euks_data))]), type="b", col="blue", pch=16, xlab="SPOT Sampling Date", ylab="Proportion of 18S sequences \n in 1-80um size class", main="Braarudosphaera bigelowii ASV4 \n AE 150m (8c144683114fbb1ad2e9425f7dcd1b02)")
par(op)

#AE from 500m, col=purple
op <- par(mar=c(5, 6, 4, 2) + 0.1)
plot(as.Date(euks_dates[grep("500m.AE", colnames(euks_data))+1]), as.numeric(euks_data[genus_brad_ix[5], grep("500m.AE", colnames(euks_data))]), type="b", col="purple", pch=16, xlab="SPOT Sampling Date", ylab="Proportion of 18S sequences \n in 1-80um size class", main="Braarudosphaera bigelowii ASV4 \n AE 500m (8c144683114fbb1ad2e9425f7dcd1b02)")
par(op)

#AE from 890m, no color 
op <- par(mar=c(5, 6, 4, 2) + 0.1)
plot(as.Date(euks_dates[grep("890m.AE", colnames(euks_data))+1]), as.numeric(euks_data[genus_brad_ix[5], grep("890m.AE", colnames(euks_data))]), type="b", pch=16, xlab="SPOT Sampling Date", ylab="Proportion of 18S sequences \n in 1-80um size class", main="Braarudosphaera bigelowii ASV4 \n AE 890m (8c144683114fbb1ad2e9425f7dcd1b02)")
par(op)

#Brad genus ASV6x: Row 1493 in Euks data, ab5338a49f7e9307027c50b3256a7f59-----
#AE from 5m depth first, color= red
op <- par(mar=c(5, 6, 4, 2) + 0.1)
plot(as.Date(euks_dates[grep("5m.AE", colnames(euks_data))+1]), as.numeric(euks_data[genus_brad_ix[6], grep("5m.AE", colnames(euks_data))]), type="b", col="red", pch=16, xlab="SPOT Sampling Date", ylab="Proportion of 18S sequences \n in 1-80um size class", main="Braarudosphaera sp. ASV6x \n AE 5m (ab5338a49f7e9307027c50b3256a7f59)")
par(op)

#AE from DCM next, color=green
op <- par(mar=c(5, 6, 4, 2) + 0.1)
plot(as.Date(euks_dates[grep("DCM.AE", colnames(euks_data))+1]), as.numeric(euks_data[genus_brad_ix[6], grep("DCM.AE", colnames(euks_data))]), type="b", col="forestgreen", pch=16, xlab="SPOT Sampling Date", ylab="Proportion of 18S sequences \n in 1-80um size class", main="Braarudosphaera sp. ASV6x \n AE DCM (ab5338a49f7e9307027c50b3256a7f59)")
par(op)

#AE from 150m, color=blue
op <- par(mar=c(5, 6, 4, 2) + 0.1)
plot(as.Date(euks_dates[grep("150m.AE", colnames(euks_data))+1]), as.numeric(euks_data[genus_brad_ix[6], grep("150m.AE", colnames(euks_data))]), type="b", col="blue", pch=16, xlab="SPOT Sampling Date", ylab="Proportion of 18S sequences \n in 1-80um size class", main="Braarudosphaera sp. ASV6x \n AE 150m (ab5338a49f7e9307027c50b3256a7f59)")
par(op)

#AE from 500m, col=purple
op <- par(mar=c(5, 6, 4, 2) + 0.1)
plot(as.Date(euks_dates[grep("500m.AE", colnames(euks_data))+1]), as.numeric(euks_data[genus_brad_ix[6], grep("500m.AE", colnames(euks_data))]), type="b", col="purple", pch=16, xlab="SPOT Sampling Date", ylab="Proportion of 18S sequences \n in 1-80um size class", main="Braarudosphaera sp. ASV6x \n AE 500m (ab5338a49f7e9307027c50b3256a7f59)")
par(op)

#AE from 890m, no color 
op <- par(mar=c(5, 6, 4, 2) + 0.1)
plot(as.Date(euks_dates[grep("890m.AE", colnames(euks_data))+1]), as.numeric(euks_data[genus_brad_ix[6], grep("890m.AE", colnames(euks_data))]), type="b", pch=16, xlab="SPOT Sampling Date", ylab="Proportion of 18S sequences \n in 1-80um size class", main="Braarudosphaera sp. ASV6x \n AE 890m (ab5338a49f7e9307027c50b3256a7f59)")
par(op)

which(euks_data[1493,2:ncol(euks_data)]>0)
names(euks_data)[970] #Add one to indexing to account for "OTU_ID"
euks_data[1493,970]
#So this ASV was once 2% of 18S sequences in Durapore size fraction... once.

#Brad genus ASV7x: Row 730 in Euks data, be3cdecefbceb0d8b25a2e42ed058b50-----
#AE from 5m depth first, color= red
op <- par(mar=c(5, 6, 4, 2) + 0.1)
plot(as.Date(euks_dates[grep("5m.AE", colnames(euks_data))+1]), as.numeric(euks_data[genus_brad_ix[7], grep("5m.AE", colnames(euks_data))]), type="b", col="red", pch=16, xlab="SPOT Sampling Date", ylab="Proportion of 18S sequences \n in 1-80um size class", main="Braarudosphaera sp. ASV7x \n AE 5m (be3cdecefbceb0d8b25a2e42ed058b50)")
par(op)

#AE from DCM next, color=green
op <- par(mar=c(5, 6, 4, 2) + 0.1)
plot(as.Date(euks_dates[grep("DCM.AE", colnames(euks_data))+1]), as.numeric(euks_data[genus_brad_ix[7], grep("DCM.AE", colnames(euks_data))]), type="b", col="forestgreen", pch=16, xlab="SPOT Sampling Date", ylab="Proportion of 18S sequences \n in 1-80um size class", main="Braarudosphaera sp. ASV7x \n AE DCM (be3cdecefbceb0d8b25a2e42ed058b50)")
par(op)

#AE from 150m, color=blue
op <- par(mar=c(5, 6, 4, 2) + 0.1)
plot(as.Date(euks_dates[grep("150m.AE", colnames(euks_data))+1]), as.numeric(euks_data[genus_brad_ix[7], grep("150m.AE", colnames(euks_data))]), type="b", col="blue", pch=16, xlab="SPOT Sampling Date", ylab="Proportion of 18S sequences \n in 1-80um size class", main="Braarudosphaera sp. ASV7x \n AE 150m (be3cdecefbceb0d8b25a2e42ed058b50)")
par(op)

#AE from 500m, col=purple
op <- par(mar=c(5, 6, 4, 2) + 0.1)
plot(as.Date(euks_dates[grep("500m.AE", colnames(euks_data))+1]), as.numeric(euks_data[genus_brad_ix[7], grep("500m.AE", colnames(euks_data))]), type="b", col="purple", pch=16, xlab="SPOT Sampling Date", ylab="Proportion of 18S sequences \n in 1-80um size class", main="Braarudosphaera sp. ASV7x \n AE 500m (be3cdecefbceb0d8b25a2e42ed058b50)")
par(op)

#AE from 890m, no color 
op <- par(mar=c(5, 6, 4, 2) + 0.1)
plot(as.Date(euks_dates[grep("890m.AE", colnames(euks_data))+1]), as.numeric(euks_data[genus_brad_ix[7], grep("890m.AE", colnames(euks_data))]), type="b", pch=16, xlab="SPOT Sampling Date", ylab="Proportion of 18S sequences \n in 1-80um size class", main="Braarudosphaera sp. ASV7x \n AE 890m (be3cdecefbceb0d8b25a2e42ed058b50)")
par(op)

dev.off()

#Plot Brad genus ASV6 and 7 in 5m Durapore samples---
op <- par(mar=c(5, 6, 4, 2) + 0.1)
plot(as.Date(euks_dates[grep("5m.D", colnames(euks_data))+1]), as.numeric(euks_data[genus_brad_ix[6], grep("5m.D", colnames(euks_data))]), type="b", pch=16, xlab="SPOT Sampling Date", ylab="Proportion of 18S sequences \n in 0.22-1um size class", main="Braarudosphaera sp. ASV6x \n Dura. 5m (ab5338a49f7e9307027c50b3256a7f59)")
par(op)

op <- par(mar=c(5, 6, 4, 2) + 0.1)
plot(as.Date(euks_dates[grep("5m.D", colnames(euks_data))+1]), as.numeric(euks_data[genus_brad_ix[7], grep("5m.D", colnames(euks_data))]), type="b", pch=16, xlab="SPOT Sampling Date", ylab="Proportion of 18S sequences \n in 0.22-1um size class", main="Braarudosphaera sp. ASV7x \n Dura. 5m (be3cdecefbceb0d8b25a2e42ed058b50)")
par(op)
#This one is up to 6% of the 18S sequences in Durapore 5m 
