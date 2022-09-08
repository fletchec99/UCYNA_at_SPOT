#Format SPOT tag data for ELSA!!
#Goalz: 1) Subset to just get dates of interest (2008-2018) and AE filters, 2) Linearly interpolate, 3) CLR transform
#Ver. 10.07.20 <- updated to write out interpolated files and to include mclr() transform via SPRING package

rm(list=ls())
getwd()
setwd("/Users/cfletcherh1/Desktop/DataAnalysis/SPOT_16S_18S_eLSA/")
#1. First, read in files-----

#Proks with chloroplasts - proportions
#Set rownames equal to OTU IDs
proks_data <- read.table("ModifiedFiles/3.SPOT_16S_w.chloro_proportions.tsv", header=T, sep="\t", stringsAsFactors = F, row.names=1) 
#Apparently the first file (#1, raw counts) is not tab delimited anymore... will that be a problem?
#sep=" "
dim(proks_data) #78953 1207
head(names(proks_data))

#Proks taxonomy 
proks_tax <- read.table("ModifiedFiles/Proks_tax_classified_23072020_SILVA132.tsv", sep="\t", header=T, stringsAsFactors = F)
dim(proks_tax)
names(proks_tax)
head(proks_tax)

#Euks without metazoans - proportions
euks_data <- read.table("ModifiedFiles/8.SPOT_18S_no_metaz_proportions.tsv", header=T, stringsAsFactors = F, sep="\t", row.names=1) #Okay this one is separated by a tab, WTF? 
dim(euks_data) #24802 1097
head(names(euks_data))
head(row.names(euks_data))

#Euks taxonomy
euks_tax <- read.table("ModifiedFiles/Euks_tax_classified_19052020_SILVA132.tsv", header=T, stringsAsFactors = F, sep="\t")
dim(euks_tax)
head(euks_tax)

#I. PROKS DATA------

#2. Subset to only take a) 5m AE, b) DCM AE, c) 5m D, d) DCM D from  proks samples after 2008-----

#a) Subset to take only 5m AE filters-----
proks_subs_AE_5m <- proks_data[,grep("5m.AE", colnames(proks_data))]
dim(proks_subs_AE_5m) #78953 125
head(row.names(proks_subs_AE_5m))
head(names(proks_subs_AE_5m))

#Subset to take only filters after 2008
include_ix <- c()
for(i in colnames(proks_subs_AE_5m)){
  a <- strsplit(i, split=".", fixed=T)
  #print(a)
  b <- paste(a[[1]][2], a[[1]][3], a[[1]][4], sep=".")
  #b is is just the date. If the year is later than 2007, include this column
  if(as.numeric(a[[1]][2])>2007){
    include_ix <- c(include_ix, grep(b, colnames(proks_subs_AE_5m)))
  } 
}
length(include_ix) #110

proks_subs_AE_5m_2008 <- proks_subs_AE_5m[,include_ix]
dim(proks_subs_AE_5m_2008) #78953 110

#b) Subset to only take AE samples from the DCM-----
length(grep("DCM.AE", colnames(proks_data))) #118
proks_subs_AE_DCM <- proks_data[,grep("DCM.AE", colnames(proks_data))]
dim(proks_subs_AE_DCM) #78953 118 
head(row.names(proks_subs_AE_DCM))
head(colnames(proks_subs_AE_DCM))

#Now take only filters after 2008
include_ix <- c()
for(i in colnames(proks_subs_AE_DCM)){
  a <- strsplit(i, split=".", fixed=T)
  #print(a)
  b <- paste(a[[1]][2], a[[1]][3], a[[1]][4], sep=".")
  #b is is just the date. If the year is later than 2007, include this column
  if(as.numeric(a[[1]][2])>2007){
    include_ix <- c(include_ix, grep(b, colnames(proks_subs_AE_DCM)))
  } 
}
length(include_ix) #104

proks_subs_AE_DCM_2008 <- proks_subs_AE_DCM[,include_ix]
dim(proks_subs_AE_DCM_2008) #78953 104
head(colnames(proks_subs_AE_DCM_2008))

#c) Subset to take Durapore samples from 5m-----
length(grep("5m.D", colnames(proks_data))) #188 
proks_subs_D_5m <- proks_data[,grep("5m.D", colnames(proks_data))]
dim(proks_subs_D_5m) #78953 188

#Now let's find out how many of those 188 samples are from 2008 or later
include_ix <- c()
for(i in colnames(proks_subs_D_5m)){
  a <- strsplit(i, split=".", fixed=T)
  #print(a)
  b <- paste(a[[1]][2], a[[1]][3], a[[1]][4], sep=".")
  #b is is just the date. If the year is later than 2007, include this column
  if(as.numeric(a[[1]][2])>2007){
    include_ix <- c(include_ix, grep(b, colnames(proks_subs_D_5m)))
  } 
}
length(include_ix) #118, the most of any matrix so far #Highest possible is 120

proks_subs_D_5m_2008 <- proks_subs_D_5m[,include_ix]
dim(proks_subs_D_5m_2008) #78953 118 
head(colnames(proks_subs_D_5m_2008))

#d) Lastly, subset to take Durapore samples from DCM-----
length(grep("DCM.D", colnames(proks_data))) #174
proks_subs_D_DCM <- proks_data[,grep("DCM.D", colnames(proks_data))]
dim(proks_subs_D_DCM) #78953 174

#Now let's find out how many of those samples are from 2008 or later
include_ix <- c()
for(i in colnames(proks_subs_D_DCM)){
  a <- strsplit(i, split=".", fixed=T)
  #print(a)
  b <- paste(a[[1]][2], a[[1]][3], a[[1]][4], sep=".")
  #b is is just the date. If the year is later than 2007, include this column
  if(as.numeric(a[[1]][2])>2007){
    include_ix <- c(include_ix, grep(b, colnames(proks_subs_D_DCM)))
  } 
}
length(include_ix) #112

proks_subs_D_DCM_2008 <- proks_subs_D_DCM[,include_ix]
dim(proks_subs_D_DCM_2008) #78953 112

#3. Subset to exclude samples with few 18S ASVs (Not necessary for proks data, include all samples)

#4. Linearly interpolate missing dates and write out interpolated data-----
library(tidyverse)
library(tidyr)
library(zoo)

#a) 5m AE samples----- 

#Set up a new data frame
length(colnames(proks_subs_AE_5m_2008))
length(rep(2008:2018, each=12))
length(rep(c("01", "02", "03", "04", "05", "06", "07", "08", "09", "10", "11", "12"), times=11))

proks_int_AE_5m <- as.data.frame(matrix(nrow=nrow(proks_subs_AE_5m_2008), ncol=132))
dim(proks_int_AE_5m)
row.names(proks_int_AE_5m)=row.names(proks_subs_AE_5m_2008)
colnames(proks_int_AE_5m)=paste(rep("SPOT", times=132), rep(2008:2018, each=12), rep(c("01", "02", "03", "04", "05", "06", "07", "08", "09", "10", "11", "12"), times=11), "01", rep("5m.AE", times=132), sep=".")
#Data frame filled with NA's automatically 

#Change colnames on proks_subs_AE_5m_2008 to match those of new matrix
##Ie replace the actual day of sampling with "01" 
length(colnames(proks_subs_AE_5m_2008))
newcolnames <- c()
for(i in colnames(proks_subs_AE_5m_2008)){
  a <- (strsplit(i, split=".", fixed=T))
  b <- paste(a[[1]][1], a[[1]][2], a[[1]][3], "01", a[[1]][5], a[[1]][6], sep=".")
  newcolnames <- c(newcolnames, b)
}
head(newcolnames)
head(colnames(proks_subs_AE_5m_2008))
tail(newcolnames)
length(newcolnames)
length(colnames(proks_subs_AE_5m_2008))

colnames(proks_subs_AE_5m_2008)=newcolnames

#Fill in spreadsheet with data available in proks data
for(i in c(1:110)){
  a <- newcolnames[i]
  b <- grep(a, colnames(proks_int_AE_5m))
  proks_int_AE_5m[,b]=proks_subs_AE_5m_2008[,i]
}
View(proks_int_AE_5m)
dim(proks_int_AE_5m) #78953 132

#Transpose data frame and linearly interpolate
t_proks_int_AE_5m <- na.approx(t(proks_int_AE_5m))
dim(t_proks_int_AE_5m) #125 78953 -- why did we lose 8 sampling dates? And which 8?  
head(row.names(t_proks_int_AE_5m)) #Rownames (sampling dates) have disappeared. Hm. 
head(colnames(t_proks_int_AE_5m))

#Double check which samples we lost
for(i in c((ncol(proks_int_AE_5m)-5):ncol(proks_int_AE_5m))){
  print(sum(proks_int_AE_5m[,i]))
}
#Lost the first two and the last 5

#Transpose back and write it out!
writeout_proks_int_AE_5m=t(t_proks_int_AE_5m)
dim(writeout_proks_int_AE_5m) #78953 125 
head(rownames(writeout_proks_int_AE_5m))
head(colnames(writeout_proks_int_AE_5m)) #NULL 

#Add back in colnames and double check first and last several 
colnames(writeout_proks_int_AE_5m)=colnames(proks_int_AE_5m)[-c(1:2, 128:132)]
head(colnames(writeout_proks_int_AE_5m))
head(colnames(proks_int_AE_5m))
tail(colnames(writeout_proks_int_AE_5m))
tail(colnames(proks_int_AE_5m))

#Convert to DF and add "OTU ID" as a column
writeout_proks_int_AE_5m=as.data.frame(as.matrix(writeout_proks_int_AE_5m))
typeof(writeout_proks_int_AE_5m)
dim(writeout_proks_int_AE_5m) #78953 125

writeout_proks_int_AE_5m$OTU_ID=row.names(writeout_proks_int_AE_5m)
head(writeout_proks_int_AE_5m$OTU_ID)
tail(colnames(writeout_proks_int_AE_5m))

#Make OTU_ID come first
writeout_proks_int_AE_5m=writeout_proks_int_AE_5m[,order(colnames(writeout_proks_int_AE_5m))]
head(colnames(writeout_proks_int_AE_5m))
head(writeout_proks_int_AE_5m[,1])

#Write it into int_CLR_data/ folder
write.table(x=writeout_proks_int_AE_5m, file="int_CLR_data/proks_AE_5m_int_10.08.2020.tsv", quote=F, sep="\t", row.names=F)

#b) AE DCM samples-----
#Set up a new dataframe
proks_int_AE_DCM <- as.data.frame(matrix(nrow=nrow(proks_subs_AE_DCM_2008), ncol=132))
dim(proks_int_AE_DCM) #78953 132 
row.names(proks_int_AE_DCM)=row.names(proks_subs_AE_DCM_2008)
colnames(proks_int_AE_DCM)=paste(rep("SPOT", times=132), rep(2008:2018, each=12), rep(c("01", "02", "03", "04", "05", "06", "07", "08", "09", "10", "11", "12"), times=11), "01", rep("DCM.AE", times=132), sep=".")
head(colnames(proks_int_AE_DCM))

#Change the colnames on subs_2008 to match the colnames above
##Ie replace the actual day of sampling with "01" 
length(colnames(proks_subs_AE_DCM_2008)) #104 
newcolnames <- c()
for(i in colnames(proks_subs_AE_DCM_2008)){
  a <- (strsplit(i, split=".", fixed=T))
  b <- paste(a[[1]][1], a[[1]][2], a[[1]][3], "01", a[[1]][5], a[[1]][6], sep=".")
  newcolnames <- c(newcolnames, b)
}
length(newcolnames) #104
head(newcolnames)
head(colnames(proks_subs_AE_DCM_2008))

colnames(proks_subs_AE_DCM_2008)=newcolnames

#Fill in the spreadsheet with the data available
for(i in c(1:ncol(proks_subs_AE_DCM_2008))){
  a <- newcolnames[i]
  b <- grep(a, colnames(proks_int_AE_DCM))
  #print(b)
  proks_int_AE_DCM[,b]=proks_subs_AE_DCM[,i]
}
View(proks_int_AE_DCM)

#Transpose and linearly interpolate
t_proks_int_AE_DCM <- na.approx(t(proks_int_AE_DCM))
dim(t_proks_int_AE_DCM) #125 78953 #interpolated 21 columns 
dim(proks_int_AE_DCM) #78953 132
head(colnames(t_proks_int_AE_DCM))
head(row.names(t_proks_int_AE_DCM))

#Double check which columns we lost
for(i in c(127:ncol(proks_int_AE_DCM))){
  print(sum(proks_int_AE_DCM[,i]))
}
#Lost the first two and the last 5: include 3-127

#Transpose it back and add the colnames back in 
writeout_proks_int_AE_DCM=t(t_proks_int_AE_DCM)
dim(writeout_proks_int_AE_DCM) #78953 125
head(rownames(writeout_proks_int_AE_DCM))
head(colnames(writeout_proks_int_AE_DCM))

colnames(writeout_proks_int_AE_DCM)=colnames(proks_int_AE_DCM)[3:127]
head(colnames(writeout_proks_int_AE_DCM))
head(colnames(proks_int_AE_DCM))
tail(colnames(writeout_proks_int_AE_DCM))
tail(colnames(proks_int_AE_DCM))

#Convert to dataframe and add in OTU ID
writeout_proks_int_AE_DCM=as.data.frame(as.matrix(writeout_proks_int_AE_DCM))
dim(writeout_proks_int_AE_DCM) #78953 125
head(rownames(writeout_proks_int_AE_DCM))
head(colnames(writeout_proks_int_AE_DCM))

writeout_proks_int_AE_DCM$OTU_ID=row.names(writeout_proks_int_AE_DCM)
dim(writeout_proks_int_AE_DCM) #78953 126
tail(colnames(writeout_proks_int_AE_DCM))

writeout_proks_int_AE_DCM=writeout_proks_int_AE_DCM[,order(colnames(writeout_proks_int_AE_DCM))]
head(colnames(writeout_proks_int_AE_DCM))

#Write it into int_CLR_data/ folder
write.table(x=writeout_proks_int_AE_DCM, file="int_CLR_data/proks_AE_DCM_int_10.08.2020.tsv", quote=F, sep="\t", row.names=F)

#c) 5m Durapore samples----
dim(proks_subs_D_5m_2008) #78953 118
head(colnames(proks_subs_D_5m_2008)) #Missing a lot from late 2008

#Set up a new dataframe 
proks_int_D_5m <- as.data.frame(matrix(nrow=nrow(proks_subs_D_5m_2008), ncol=132))
dim(proks_int_D_5m) #78953 132
row.names(proks_int_D_5m)=row.names(proks_subs_D_5m_2008)
colnames(proks_int_D_5m)=paste(rep("SPOT", times=132), rep(2008:2018, each=12), rep(c("01", "02", "03", "04", "05", "06", "07", "08", "09", "10", "11", "12"), times=11), "01", rep("5m.D", times=132), sep=".")
head(row.names(proks_int_D_5m))
head(colnames(proks_int_D_5m))

#Change colnames on proks_subs_D_2008 to match 
length(colnames(proks_subs_D_5m_2008)) #118
newcolnames <- c()
for(i in colnames(proks_subs_D_5m_2008)){
  a <- strsplit(i, split=".", fixed=T)
  b <- paste(a[[1]][1], a[[1]][2], a[[1]][3], "01", "5m.D", sep=".")
  newcolnames <- c(newcolnames, b)
}
length(newcolnames)
head(newcolnames)
head(colnames(proks_subs_D_5m_2008))
tail(newcolnames)
tail(colnames(proks_subs_D_DCM_2008))

colnames(proks_subs_D_5m_2008)=newcolnames

#Fill in data frame with what samples we have
for(i in c(1:ncol(proks_subs_D_5m_2008))){
  a <- newcolnames[i]
  b <- grep(a, colnames(proks_int_D_5m))
  proks_int_D_5m[,b]=proks_subs_D_5m_2008[,i]
}
View(proks_int_D_5m)

#Linearly interpolate and transpose all in one!
t_proks_int_D_5m <- na.approx(t(proks_int_D_5m))
dim(t_proks_int_D_5m) #124 78953

#I'm pretty sure we lost the first three dates (not two, i.e. no data till April 2008), and the last five #Which would be 8 dates lost 
for(i in c((ncol(proks_int_D_5m)-5):ncol(proks_int_D_5m))){
  print(colnames(proks_int_D_5m[i]))
  print(sum(proks_int_D_5m[,i]))
}
#Yeah we are missing until April 2008 
#And after July 2018

#Write it out! 
#Transpose it back
writeout_proks_int_D_5m <- t(t_proks_int_D_5m)
dim(writeout_proks_int_D_5m) #78953 124
head(rownames(writeout_proks_int_D_5m))
head(colnames(writeout_proks_int_D_5m))

#Add colnames back in 
colnames(writeout_proks_int_D_5m)=colnames(proks_int_D_5m)[4:127]
head(writeout_proks_int_D_5m[,1:5])
head(proks_int_D_5m[,1:5])
head(writeout_proks_int_D_5m[,120:124])
head(proks_int_D_5m[,127:132])

#Convert to DF and add "OTU_ID"
writeout_proks_int_D_5m=as.data.frame(as.matrix(writeout_proks_int_D_5m))
dim(writeout_proks_int_D_5m) #78953 124
writeout_proks_int_D_5m$OTU_ID=row.names(writeout_proks_int_D_5m)
dim(writeout_proks_int_D_5m) #78953 125
tail(colnames(writeout_proks_int_D_5m))

#Re-order
writeout_proks_int_D_5m=writeout_proks_int_D_5m[,order(colnames(writeout_proks_int_D_5m))]

#Write it out!
write.table(x=writeout_proks_int_D_5m, file="int_CLR_data/proks_D_5m_int_10.08.2020.tsv", sep="\t", quote=F, row.names=F)

#d) DCM Durapore samples-----
dim(proks_subs_D_DCM_2008) #78953 112 

#Set up a new dataframe
proks_int_D_DCM <- as.data.frame(matrix(ncol=132, nrow=nrow(proks_subs_D_DCM_2008)))
dim(proks_int_D_DCM) #78953 112
row.names(proks_int_D_DCM)=row.names(proks_subs_D_DCM_2008)
colnames(proks_int_D_DCM)=paste(rep("SPOT", times=132), rep(2008:2018, each=12), rep(c("01", "02", "03", "04", "05", "06", "07", "08", "09", "10", "11", "12"), times=11), "01", rep("DCM.D", times=132), sep=".")
head(row.names(proks_int_D_DCM))
head(colnames(proks_int_D_DCM))

#Change the colnames of proks_subs_D_DCM_2008
newcolnames <- c()
for(i in colnames(proks_subs_D_DCM_2008)){
  a <- strsplit(i, split=".", fixed=T)
  b <- paste(a[[1]][1], a[[1]][2], a[[1]][3], "01", "DCM.D", sep=".")
  newcolnames <- c(newcolnames, b)
}
length(newcolnames)==length(colnames(proks_subs_D_DCM_2008))
head(newcolnames)
head(colnames(proks_subs_D_DCM_2008))
tail(newcolnames)
tail(colnames(proks_subs_D_DCM_2008))

colnames(proks_subs_D_DCM_2008)=newcolnames

#Fill in data for the dates on which we have data
for(i in c(1:ncol(proks_subs_D_DCM_2008))){
  a <- newcolnames[i]
  b <- grep(a, colnames(proks_int_D_DCM))
  proks_int_D_DCM[,b]=proks_subs_D_DCM_2008[,i]
}
View(proks_int_D_DCM)
dim(proks_int_D_DCM) #78953 132

#Transpose and linearly interpolate
t_proks_int_D_DCM <- na.approx(t(proks_int_D_DCM))
dim(t_proks_int_D_DCM) #125 78953
head(colnames(t_proks_int_D_DCM))
head(row.names(t_proks_int_D_DCM))

#Which samples are missing?? 
for(i in c((ncol(proks_int_D_DCM)-5):ncol(proks_int_D_DCM))){
  print(colnames(proks_int_D_DCM[i]))
  print(sum(proks_int_D_DCM[,i]))
}
#Yeah lose the first two 
#And the last 5 columns

#Write it out!
#Transpose it back
writeout_proks_int_D_DCM <- t(t_proks_int_D_DCM)
dim(writeout_proks_int_D_DCM) #78953 125
head(rownames(writeout_proks_int_D_DCM))
head(colnames(writeout_proks_int_D_DCM))

#Add colnames back in
length(colnames(proks_int_D_DCM)[3:127])
ncol(writeout_proks_int_D_DCM)
colnames(writeout_proks_int_D_DCM)=colnames(proks_int_D_DCM)[3:127]

#Convert to DF and add "OTU_ID"
writeout_proks_int_D_DCM=as.data.frame(as.matrix(writeout_proks_int_D_DCM))
typeof(writeout_proks_int_D_DCM)
dim(writeout_proks_int_D_DCM) #78953 125

writeout_proks_int_D_DCM$OTU_ID=rownames(writeout_proks_int_D_DCM)
dim(writeout_proks_int_D_DCM) #78953 126

#Re-order
tail(colnames(writeout_proks_int_D_DCM))
writeout_proks_int_D_DCM=writeout_proks_int_D_DCM[,order(colnames(writeout_proks_int_D_DCM))]
head(colnames(writeout_proks_int_D_DCM))
head(writeout_proks_int_D_DCM$OTU_ID)
head(rownames(writeout_proks_int_D_DCM))

#Write it out!!
write.table(x=writeout_proks_int_D_DCM, file="int_CLR_data/proks_D_DCM_int_10.08.2020.tsv", sep="\t", quote=F, row.names=F)

#5. CLR transform data-----
install.packages("devtools")
library(devtools)
install_github("GraceYoon/SPRING")
library(SPRING)

#SPRING mclr() requires data (count or relabun) in n by p format: # samples(n) as columns, ASVs (p) as rows
#Transpose data back or use writeout data

#a) Proks 5m AE
#CLR it!! And remove OTU_ID col
proks_CLR_AE_5m <- mclr(dat=writeout_proks_int_AE_5m[,-1])
dim(proks_CLR_AE_5m) #78953 125
head(row.names(proks_CLR_AE_5m))
head(colnames(proks_CLR_AE_5m))
#No need to transpose it back now :)

#b) Proks DCM AE samples
proks_CLR_AE_DCM <- mclr(dat=writeout_proks_int_AE_DCM[,-1])
dim(proks_CLR_AE_DCM) #78953 125
head(rownames(proks_CLR_AE_DCM))
head(colnames(proks_CLR_AE_DCM))

#c) Proks 5m Durapore samples
proks_CLR_D_5m <- mclr(dat=writeout_proks_int_D_5m[,-1])
dim(proks_CLR_D_5m) #78953 124
head(rownames(proks_CLR_D_5m))
head(colnames(proks_CLR_D_5m))

#d) Proks DCM Durapore samples
proks_CLR_D_DCM <- mclr(dat=writeout_proks_int_D_DCM[,-1])
dim(proks_CLR_D_DCM) #78953 124
head(rownames(proks_CLR_D_DCM))
head(colnames(proks_CLR_D_DCM))

#6. Write out data-----

#a) Proks AE 5m data
#Convert to data frame
proks_CLR_AE_5m=as.data.frame(as.matrix(proks_CLR_AE_5m))
dim(proks_CLR_AE_5m)
head(proks_CLR_AE_5m$SPOT.2008.03.01.5m.AE)

#Write a column, OTU_ID, and set it equal to row names
head(row.names(proks_CLR_AE_5m))
proks_CLR_AE_5m$OTU_ID=row.names(proks_CLR_AE_5m)
head(proks_CLR_AE_5m$OTU_ID)
tail(colnames(proks_CLR_AE_5m))

#Then sort, so that OTU ID comes first (all other colnames start with "SPOT")
writeout_proks_CLR_AE_5m <- proks_CLR_AE_5m[,order(colnames(proks_CLR_AE_5m))]
dim(writeout_proks_CLR_AE_5m) #78953 126
head(colnames(writeout_proks_CLR_AE_5m)) 

#Then write it out!! 
write.table(x=writeout_proks_CLR_AE_5m, file="int_CLR_data/proks_AE_5m_int_CLR_10.08.2020.tsv", sep="\t", quote=F, row.names=F)

#b) Proks AE DCM data
#Convert to data frame
proks_CLR_AE_DCM=as.data.frame(as.matrix(proks_CLR_AE_DCM))
dim(proks_CLR_AE_DCM) #78953 125 
head(proks_CLR_AE_DCM$SPOT.2008.03.01.DCM.AE)

#Write in a column for OTU_ID and set it equal to rownames
head(row.names(proks_CLR_AE_DCM))
proks_CLR_AE_DCM$OTU_ID=row.names(proks_CLR_AE_DCM)
head(proks_CLR_AE_DCM$OTU_ID)

#Sort so that OTU ID comes first
writeout_proks_CLR_AE_DCM=proks_CLR_AE_DCM[, order(colnames(proks_CLR_AE_DCM))]
dim(writeout_proks_AE_DCM) #78953 126
head(colnames(writeout_proks_CLR_AE_DCM))

#Write it out!
write.table(x=writeout_proks_CLR_AE_DCM, file="int_CLR_data/proks_AE_DCM_int_CLR_10.08.20.tsv", sep="\t", quote=F, row.names=F)

#c) Proks Durapore 5m data
#Convert to data frame
proks_CLR_D_5m=as.data.frame(as.matrix(proks_CLR_D_5m))
dim(proks_CLR_D_5m) #78953 124
head(proks_CLR_D_5m$SPOT.2008.04.01.5m.D)

#Write in OTU IDs as a column 
head(row.names(proks_CLR_D_5m))
proks_CLR_D_5m$OTU_ID=row.names(proks_CLR_D_5m)
dim(proks_CLR_D_5m) #78953 125
head(proks_CLR_D_5m$OTU_ID)

#Re order so that OTU ID comes first
order(colnames(proks_CLR_D_5m))
writeout_proks_CLR_D_5m=proks_CLR_D_5m[, order(colnames(proks_CLR_D_5m))]
head(colnames(writeout_proks_CLR_D_5m))
head(writeout_proks_CLR_D_5m[,1])

#Write it out!
write.table(x=writeout_proks_CLR_D_5m, file="int_CLR_data/proks_D_5m_int_CLR_10.08.2020.tsv", sep="\t", quote=F, row.names=F)

#d) Proks Durapore DCM 
#Convert to data frame 
proks_CLR_D_DCM=as.data.frame(as.matrix(proks_CLR_D_DCM))
dim(proks_CLR_D_DCM) #78953 125
head(proks_CLR_D_DCM$SPOT.2008.03.01.DCM.D)

#write in OTU ID as a column 
head(row.names(proks_CLR_D_DCM))
proks_CLR_D_DCM$OTU_ID=row.names(proks_CLR_D_DCM)
head(proks_CLR_D_DCM$OTU_ID)

#Re order so that OTU ID comes first
order(colnames(proks_CLR_D_DCM))
writeout_proks_CLR_D_DCM=proks_CLR_D_DCM[, order(colnames(proks_CLR_D_DCM))]
dim(writeout_proks_CLR_D_DCM)
head(colnames(writeout_proks_CLR_D_DCM))

#Write it out! 
write.table(x=writeout_proks_CLR_D_DCM, file="int_CLR_data/proks_D_DCM_int_CLR_10.08.2020.tsv", sep="\t", quote=F, row.names=F)

#II. EUKS DATA-----

#2. Subset to take only a) AE samples from 5m, b) AE samples from DCM from euks samples after 2008----
#AE only; we don't care about the Micromonas in Durapore samples :) 

#a) Subset to take AE samples from 5m-----
length(grep("5m.AE", colnames(euks_data))) #124
euks_subs_AE_5m <- euks_data[,grep("5m.AE", colnames(euks_data))]
dim(euks_subs_AE_5m) #24802 124

#How many of those samples are from 2008 or later? 
include_ix <- c()
for(i in colnames(euks_subs_AE_5m)){
  a <- strsplit(i, split=".", fixed=T)
  #print(a)
  b <- paste(a[[1]][2], a[[1]][3], a[[1]][4], sep=".")
  if(as.numeric(a[[1]][2]>2007)){
    include_ix <- c(include_ix, grep(b, colnames(euks_subs_AE_5m)))
  }
}
length(include_ix) #109, not bad

#Subset
euks_subs_AE_5m_2008 <- euks_subs_AE_5m[,include_ix]
dim(euks_subs_AE_5m_2008) #28402 109

#b) Subset to take AE samples from DCM----
length(grep("DCM.AE", colnames(euks_data))) #116
euks_subs_AE_DCM <- euks_data[,grep("DCM.AE", colnames(euks_data))]
dim(euks_subs_AE_DCM) #28402 116

#Take samples after 2008
include_ix <- c()
for(i in colnames(euks_subs_AE_DCM)){
  a <- strsplit(i, split=".", fixed=T)
  #print(a)
  b <- paste(a[[1]][2], a[[1]][3], a[[1]][4], sep=".")
  if(as.numeric(a[[1]][2]>2007)){
    include_ix <- c(include_ix, grep(b, colnames(euks_subs_AE_DCM)))
  }
}
length(include_ix) #102 #oh dear

euks_subs_AE_DCM_2008 <- euks_subs_AE_DCM[,include_ix]
dim(euks_subs_AE_DCM_2008) #28402 102

#3. Exclude samples with too few sequences to capture all 18S diversity-----
#5m: 2015.04.22, 2013.11.13, 2013.08.14, 2010.08.11, 2010.04.27
#DCM: 2010.4.27 and 2010.8.11 

#a) 5m AE samples
length(grep("2015.04.22|2013.11.13|2013.08.14|2010.08.11|2010.04.27", colnames(euks_subs_AE_5m_2008)))
dim(euks_subs_AE_5m_2008) #28402 109 #Need to remove 5 cols
exclude_cols <- grep("2015.04.22|2013.11.13|2013.08.14|2010.08.11|2010.04.27", colnames(euks_subs_AE_5m_2008))
exclude_cols
euks_subs_AE_5m_2008 <- euks_subs_AE_5m_2008[,-c(exclude_cols)]
dim(euks_subs_AE_5m_2008) #24802 104

#b) DCM AE samples
length(grep("2010.04.27|2010.08.11", colnames(euks_subs_AE_DCM_2008)))
dim(euks_subs_AE_DCM_2008) #28402 102
euks_subs_AE_DCM_2008 <- euks_subs_AE_DCM_2008[,-grep("2010.04.27|2010.08.11", colnames(euks_subs_AE_DCM_2008))]
dim(euks_subs_AE_DCM_2008) #28402 100

#4. Linearly interpolate missing samples and write out interpolated data-------
library(tidyverse)
library(tidyr)
library(zoo)

#a) 5m AE samples-----
ncol(euks_subs_AE_5m_2008) #104
length(rep(2008:2018, each=12))
length(rep(c("01", "02", "03", "04", "05", "06", "07", "08", "09", "10", "11", "12"), times=11))

#Set up a new dataframe with dates for 2008-2018 evenly spaced (pretend we sampled on the first day of the month every month)
euks_int_AE_5m <- as.data.frame(matrix(nrow=nrow(euks_subs_AE_5m_2008), ncol=132))
dim(euks_int_AE_5m) #24802 132
row.names(euks_int_AE_5m)=row.names(euks_subs_AE_5m_2008)
colnames(euks_int_AE_5m)=paste(rep("SPOT", times=132), rep(2008:2018, each=12), rep(c("01", "02", "03", "04", "05", "06", "07", "08", "09", "10", "11", "12"), times=11), "01", rep("5m.AE", times=132), sep=".")

#Change colnames of euks_subs_AE_5m_2008 to be the same as colnames of euks_int_AE_5m
#i.e. sampling date is "01"
newcolnames <- c()
for(i in colnames(euks_subs_AE_5m_2008)){
  a <- strsplit(i, split=".", fixed=T)
  b <- paste(a[[1]][1],a[[1]][2],a[[1]][3],"01", "5m.AE", sep=".")
  newcolnames <- c(newcolnames, b)
}
length(newcolnames)==length(colnames(euks_subs_AE_5m_2008)) #TRUE
head(newcolnames)
head(colnames(euks_subs_AE_5m_2008))
tail(newcolnames)
tail(colnames(euks_subs_AE_5m_2008))
#Everything looks good!

#Set em up!
colnames(euks_subs_AE_5m_2008)=newcolnames

#Fill in data where we have it
for(i in c(1:ncol(euks_subs_AE_5m_2008))){
  a <- newcolnames[i]
  b <- grep(a, colnames(euks_int_AE_5m))
  euks_int_AE_5m[,b]=euks_subs_AE_5m_2008[,i]
}
dim(euks_int_AE_5m) #28402 132
View(euks_int_AE_5m)

#Check the samples that should be missing are actually missing
#Dates have changed, so just check month
for(i in c("2015.04", "2013.11", "2013.08", "2010.08", "2010.04")){
  print(sum(euks_int_AE_5m[,grep(i, colnames(euks_int_AE_5m))]))
} #all NAs means it worked! 

#Transpose and linearly interpolate! 
t_euks_int_AE_5m <- na.approx(t(euks_int_AE_5m))
dim(t_euks_int_AE_5m) #125 28402
head(row.names(t_euks_int_AE_5m)) #rownames have vanished
head(colnames(t_euks_int_AE_5m)) #colnames are ASV names

#Figure out which dates are missing data
for(i in c(122:132)){
  print(colnames(euks_int_AE_5m)[i])
  print(sum(euks_int_AE_5m[,i]))
}
#Lose the first two columns - NAs
grep("SPOT.2018.08.01.5m.AE", colnames(euks_int_AE_5m)) #128
length(c(128:132))
#And from column 128:132 - lose the last 5 columns 
length(c(3:127)) #125

#Transpose it back
writeout_euks_int_AE_5m=t(t_euks_int_AE_5m)
dim(writeout_euks_int_AE_5m) #28402 125
head(rownames(writeout_euks_int_AE_5m))
head(colnames(writeout_euks_int_AE_5m))

#Add colnames back in 
colnames(writeout_euks_int_AE_5m)=colnames(euks_int_AE_5m)[3:127]
head(colnames(writeout_euks_int_AE_5m))
head(colnames(euks_int_AE_5m))
tail(colnames(writeout_euks_int_AE_5m))
tail(colnames(euks_int_AE_5m))

#Convert to DF
typeof(writeout_euks_int_AE_5m)
writeout_euks_int_AE_5m=as.data.frame(as.matrix(writeout_euks_int_AE_5m))
typeof(writeout_euks_int_AE_5m)
dim(writeout_euks_int_AE_5m)

#Add OTU ID in as a column
head(rownames(writeout_euks_int_AE_5m))
writeout_euks_int_AE_5m$OTU_ID=rownames(writeout_euks_int_AE_5m)
head(writeout_euks_int_AE_5m$OTU_ID)
dim(writeout_euks_int_AE_5m) #28402 126

#Sort it
order(colnames(writeout_euks_int_AE_5m))
tail(colnames(writeout_euks_int_AE_5m))
writeout_euks_int_AE_5m=writeout_euks_int_AE_5m[,order(colnames(writeout_euks_int_AE_5m))]
tail(colnames(writeout_euks_int_AE_5m))
head(writeout_euks_int_AE_5m[,1])

#Write it out!
write.table(x=writeout_euks_int_AE_5m, file="int_CLR_data/euks_AE_5m_int_10.09.2020.tsv", sep="\t", quote=F, row.names=F)

#b) AE DCM samples-----
length(colnames(euks_subs_AE_DCM_2008)) #100
#We are missing 25 dates #whoooo

#Set up a new dataframe
euks_int_AE_DCM <- as.data.frame(matrix(nrow=nrow(euks_subs_AE_DCM_2008), ncol=132))
dim(euks_int_AE_DCM) #28402 132
row.names(euks_int_AE_DCM)=row.names(euks_subs_AE_DCM_2008)
colnames(euks_int_AE_DCM)=paste(rep("SPOT", times=132), rep(2008:2018, each=12), rep(c("01", "02", "03", "04", "05", "06", "07", "08", "09", "10", "11", "12"), times=11), "01", rep("DCM.AE", times=132), sep=".")
tail(colnames(euks_int_AE_DCM))

#Change colnames on euks_subs_AE_DCM_2008 to match interpolated data frame 
newcolnames <- c()
for(i in colnames(euks_subs_AE_DCM_2008)){
  a <- strsplit(i, split=".", fixed=T)
  b <- paste(a[[1]][1],a[[1]][2],a[[1]][3],"01", "DCM.AE", sep=".")
  newcolnames <- c(newcolnames, b)
}
length(unique(newcolnames)) #98 #Why are there duplicates? 
which(duplicated(newcolnames))
colnames(euks_subs_AE_DCM_2008)[which(duplicated(newcolnames))] #Because we collected samples in the same month, same year, two different days 2x #But not from 5m? 
length(newcolnames)==length(euks_subs_AE_DCM_2008)
head(newcolnames)
head(colnames(euks_subs_AE_DCM_2008))
tail(newcolnames)
tail(colnames(euks_subs_AE_DCM_2008))

#Looks good, set 'em!
colnames(euks_subs_AE_DCM_2008)=newcolnames

#Fill in data on the days that we do have samples
for(i in c(1:ncol(euks_subs_AE_DCM_2008))){
  a <- newcolnames[i]
  b <- grep(a, colnames(euks_int_AE_DCM))
  euks_int_AE_DCM[,b]=euks_subs_AE_DCM_2008[,i]
}
dim(euks_int_AE_DCM) #28402 132
View(euks_int_AE_DCM)

#Transpose and linearly interpolate
t_euks_int_AE_DCM <- na.approx(t(euks_int_AE_DCM))
View(t_euks_int_AE_DCM)
dim(t_euks_int_AE_DCM) #125 28402
head(row.names(t_euks_int_AE_DCM))
head(colnames(t_euks_int_AE_DCM))

#Figure out which samples are missing 
for(i in c((ncol(euks_int_AE_DCM)-5):ncol(euks_int_AE_DCM))){
  print(sum(euks_int_AE_DCM[,i]))
}
#Lose the first two
#and the last five
sum(euks_int_AE_DCM[,128])

#Transpose it back
writeout_euks_int_AE_DCM <- t(t_euks_int_AE_DCM)
dim(writeout_euks_int_AE_DCM) #28402 125
head(rownames(writeout_euks_int_AE_DCM))
head(colnames(writeout_euks_int_AE_DCM))

#Add in colnames
colnames(writeout_euks_int_AE_DCM)=colnames(euks_int_AE_DCM)[3:127]
head(colnames(writeout_euks_int_AE_DCM))
head(colnames(euks_int_AE_DCM))
tail(colnames(writeout_euks_int_AE_DCM))
tail(colnames(euks_int_AE_DCM))

#Convert to DF
writeout_euks_int_AE_DCM <- as.data.frame(as.matrix(writeout_euks_int_AE_DCM))
dim(writeout_euks_int_AE_DCM) #28402 125
head(rownames(writeout_euks_int_AE_DCM))
head(colnames(writeout_euks_int_AE_DCM))

#Add in OTU ID
writeout_euks_int_AE_DCM$OTU_ID=rownames(writeout_euks_int_AE_DCM)
head(writeout_euks_int_AE_DCM$OTU_ID)
dim(writeout_euks_int_AE_DCM) #28402 126

#Sort
order(colnames(writeout_euks_int_AE_DCM))
writeout_euks_int_AE_DCM=writeout_euks_int_AE_DCM[,order(colnames(writeout_euks_int_AE_DCM))]
head(colnames(writeout_euks_int_AE_DCM))

#Write it out!
write.table(x=writeout_euks_int_AE_DCM, file="int_CLR_data/euks_AE_DCM_int_10.09.20.tsv", sep="\t", quote=F, row.names=F)

#5. CLR transform samples and write out CLR data!-----
install.packages("devtools")
library(devtools)
install_github("GraceYoon/SPRING")
library(SPRING)

#a) AE 5m samples
#CLR transform
euks_CLR_AE_5m <- mclr(dat=writeout_euks_int_AE_5m[,-1])
dim(euks_CLR_AE_5m) #28402 125
head(rownames(euks_CLR_AE_5m))
head(colnames(euks_CLR_AE_5m))
grep("OTU", colnames(euks_CLR_AE_5m)) #Not found

#Convert to DF
typeof(euks_CLR_AE_5m) #"double"
euks_CLR_AE_5m=as.data.frame(as.matrix(euks_CLR_AE_5m))
typeof(euks_CLR_AE_5m) #"list"
dim(euks_CLR_AE_5m) #28502 125
head(rownames(euks_CLR_AE_5m))
head(colnames(euks_CLR_AE_5m))

#Add in "OTU_ID"
euks_CLR_AE_5m$OTU_ID=row.names(euks_CLR_AE_5m)
dim(euks_CLR_AE_5m) #28402 126

#Order
order(colnames(euks_CLR_AE_5m))
writeout_euks_CLR_AE_5m=euks_CLR_AE_5m[,order(colnames(euks_CLR_AE_5m))]
head(writeout_euks_CLR_AE_5m$OTU_ID)

#Write it out!
write.table(x=writeout_euks_CLR_AE_5m, file="int_CLR_data/euks_AE_5m_int_CLR_10.09.2020.tsv", quote=F, row.names=F, sep="\t")

#b) AE DCM samples
#First, CLR transform
euks_CLR_AE_DCM <- mclr(dat=writeout_euks_int_AE_DCM[,-1])
dim(euks_CLR_AE_DCM) #28402 125
head(rownames(euks_CLR_AE_DCM))
head(colnames(euks_CLR_AE_DCM))

#Convert to DF
euks_CLR_AE_DCM <- as.data.frame(as.matrix(euks_CLR_AE_DCM))
typeof(euks_CLR_AE_DCM) #"list"
head(rownames(euks_CLR_AE_DCM))
head(colnames(euks_CLR_AE_DCM))
tail(colnames(euks_CLR_AE_DCM))

#Add in OTU ID as a column
euks_CLR_AE_DCM$OTU_ID=row.names(euks_CLR_AE_DCM)
dim(euks_CLR_AE_DCM) #28402 126
tail(colnames(euks_CLR_AE_DCM))

#Sort
order(colnames(euks_CLR_AE_DCM))
writeout_euks_CLR_AE_DCM=euks_CLR_AE_DCM[,order(colnames(euks_CLR_AE_DCM))]
head(colnames(writeout_euks_CLR_AE_DCM))
dim(writeout_euks_CLR_AE_DCM)

#Write it out!
write.table(x=writeout_euks_CLR_AE_DCM, file="int_CLR_data/euks_AE_DCM_int_CLR_10.09.2020.tsv", sep="\t", quote=F, row.names=F)
