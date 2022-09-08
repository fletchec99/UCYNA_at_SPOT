#Setting up eLSA (again): take only 200 taxa at a time, and run eLSA on ada piecewise
#This time exluding UCYN-A
#Ver. 11.2.2020

getwd()
setwd("Desktop/DataAnalysis/SPOT_16S_18S_eLSA/")

#1. CLR transformed from 5m----- 

#Read in the data using my trick for the colnames
temp <- "eLSA_input/eLSA_input_5m_CLR_abund_and_lit_10.12.2020.tsv"
colnames <- scan(text=readLines(temp, 2), what="", quiet=T)
length(colnames) #253, yikes
head(colnames)
which(grepl("SPOT", colnames)==FALSE)
colnames[125:130] #127 is the last good one
colnames <- c("OTU_ID", colnames[3:127])

#Now read in data
CLR_5m <- read.table(temp, col.names=colnames, stringsAsFactors = F, sep="\t")
dim(CLR_5m) #580 126

#Okay, now divide this into three chunks 
#DON"T include UCYN-A in each one

#A: First 189 rows
subs_A <- as.data.frame(matrix(nrow=(189), ncol=ncol(CLR_5m)))
dim(subs_A) #189 126
colnames(subs_A)=colnames(CLR_5m)
#Add in data
subs_A=CLR_5m[1:189,]

#Re-do the colname "OTU ID" to put a hash first, then write it out!
colnames(subs_A)[1]="#OTU ID"
write.table(x=subs_A, file="eLSA_input/subs_A_eLSA_input_5m_CLR_11.08.2020.tsv", sep="\t", quote=F, row.names=F)

#B: Next third (rows 190-378)
subs_B <- as.data.frame(matrix(nrow=(189), ncol=ncol(CLR_5m)))
dim(subs_B) #189 126
colnames(subs_B)=colnames(CLR_5m)
#Add in data
subs_B[1:189,]=CLR_5m[190:378,]

#Format "OTU ID" and write out!
colnames(subs_B)[1]="#OTU ID"
write.table(x=subs_A, file="eLSA_input/subs_B_eLSA_input_5m_CLR_11.08.2020.tsv", sep="\t", quote=F, row.names=F)

#C: Aaaaand final third! 
subs_C <-as.data.frame(matrix(nrow=nrow(CLR_5m)-378, ncol=ncol(CLR_5m))) 
dim(subs_C) #202 126
(nrow(subs_C))+(nrow(subs_B))+nrow(subs_A)==nrow(CLR_5m) #TRUE
colnames(subs_C)=colnames(CLR_5m)
#Add in data
subs_C=CLR_5m[(nrow(subs_C))+(nrow(subs_B)):nrow(CLR_5m),]
#No need to add in UCYN-A data! :) 

#Fix "OTU_ID" and then write out!
colnames(subs_C)[1]="#OTU ID"
write.table(x=subs_C, file="eLSA_input/subs_C_eLSA_input_5m_CLR_11.08.2020.tsv", sep="\t", quote=F, row.names=F)

#2. CLR transformed from DCM----- 
#Read in data, using my trick for colnames
temp <- "eLSA_input/eLSA_input_DCM_CLR_abund-lit_10.13.2020.tsv"
colnames <- scan(text=readLines(temp, 1), what="", quiet=T)
colnames
length(colnames) #Hint, use temp, 1 to only scan the first line!!
colnames[1]="OTU_ID" #Get rid of hash (#) 
CLR_DCM <- read.table(temp, col.names=colnames, stringsAsFactors = F, sep="\t")
dim(CLR_DCM) #733 126
tail(CLR_DCM$OTU_ID, n=20)
733/3 #245 2x and 244 1x

#Divide into three chunks, set colnames, write out with a "#OTU ID"
#First chunk
subs_A=CLR_DCM[1:245,]
dim(subs_A) #245 126
colnames(subs_A) #wow, kept its colnames!
head(subs_A$SPOT.2008.03.01.DCM.AE) #Looks good!
colnames(subs_A)[1]="#OTU ID"
#Write it out!
write.table(x=subs_A, file="eLSA_input/subs_A_eLSA_input_DCM_CLR_11.16.2020.tsv", sep="\t", quote=F, row.names=F)

#Second chunk 
subs_B=CLR_DCM[(nrow(subs_A)+1):(nrow(subs_A)+245),]
dim(subs_B) #245 126 
colnames(subs_B)
head(subs_B$SPOT.2008.03.01.DCM.AE)
colnames(subs_B)[1]="#OTU ID"
#Write it out!
write.table(x=subs_B, file="eLSA_input/subs_B_eLSA_input_DCM_CLR_11.16.2020.tsv", sep="\t", quote=F, row.names=F)

#Third chunk
subs_C=CLR_DCM[(nrow(subs_A)+nrow(subs_B)+1):nrow(CLR_DCM),]
dim(subs_C) #243 126
nrow(subs_A)+nrow(subs_B)+nrow(subs_C)==nrow(CLR_DCM) #TRUE 
colnames(subs_C)
head(subs_C$SPOT.2008.03.01.DCM.AE)
tail(subs_C$OTU_ID, n=20)
colnames(subs_C)[1]="#OTU ID"
#Write out
write.table(x=subs_C, file="eLSA_input/subs_C_eLSA_input_DCM_CLR_11.16.2020.tsv", sep="\t", quote=F, row.names=F)

#3. Interpolated, non CLR transformed data from 5m------
#Read in the file
temp <- "eLSA_input/eLSA_input_5m_int_abund_and_lit_10.12.2020.tsv"
colnames <- scan(text=readLines(temp, 1), what="", quiet=T)
colnames #For some reason, it is splitting "#OTU" and "ID"
colnames=c("OTU_ID", colnames[3:length(colnames)])
colnames
length(colnames) #126 #That's better
int_5m=read.table(temp, col.names=colnames, stringsAsFactors = F, sep="\t")
dim(int_5m) #579 126 #Why is this not 580? 

#Divide into thirds
#First third
subs_A=int_5m[1:(nrow(int_5m)/3),]
dim(subs_A) #193 126
colnames(subs_A)
head(subs_A$SPOT.2008.03.01.5m.AE)
colnames(subs_A)[1]="#OTU ID"
write.table(x=subs_A, file="eLSA_input/subs_A_input_5m_int_11.16.2020.tsv", sep="\t", quote=F, row.names=F) #Names should include "eLSA_input"

#Second third
subs_B=int_5m[(nrow(subs_A)+1):(2*nrow(subs_A)),]
dim(subs_B) #193 126
colnames(subs_B)
head(subs_B$SPOT.2008.03.01.5m.AE, n=20)
colnames(subs_B)[1]="#OTU ID"
write.table(x=subs_B, file="eLSA_input/subs_B_input_5m_int_11.16.2020.tsv", sep="\t", quote=F, row.names=F) #Names should include "eLSA_input"

#Third third
subs_C=int_5m[(nrow(subs_A)+nrow(subs_B)+1):nrow(int_5m),]
dim(subs_C) #193 126
nrow(subs_A)+nrow(subs_B)+nrow(subs_C)==nrow(int_5m) 
head(subs_C$SPOT.2008.03.01.5m.AE)
colnames(subs_C)
colnames(subs_C)[1]="#OTU ID"
write.table(x=subs_C, file="eLSA_input/subs_C_input_5m_int_11.16.2020.tsv", sep="\t", quote=F, row.names=F) #Names should include "eLSA_input"

#4. Interpolated, non CLR-transformed data from DCM-----
#Read in data
temp <- "eLSA_input/eLSA_input_DCM_int_abund-lit_10.13.2020.tsv"
colnames=scan(text=readLines(temp, 1), what="", quiet=T)
length(colnames) #126 
head(colnames)
tail(colnames)
colnames[1]="OTU_ID"
int_DCM=read.table(temp, col.names=colnames, stringsAsFactors = F, sep="\t")
dim(int_DCM) #733 126 
head(int_DCM$SPOT.2008.03.01.DCM.AE)
round(nrow(int_DCM)/3)

#First third
subs_A=int_DCM[1:round(nrow(int_DCM)/3),]
dim(subs_A) #244 126
colnames(subs_A)
colnames(subs_A)[1]="#OTU ID"
#write it out!
write.table(subs_A, file="eLSA_input/subs_A_eLSA_input_DCM_int_11.24.2020.tsv", sep="\t", quote=F, row.names=F)

#Second third
subs_B=int_DCM[(1+round(nrow(int_DCM)/3)):(2*round(nrow(int_DCM)/3)),]
dim(subs_B) #244 126
colnames(subs_B)
colnames(subs_B)[1]="#OTU ID"
#Write out
write.table(subs_B, file="eLSA_input/subs_B_eLSA_input_DCM_int_11.24.2020.tsv", sep="\t", quote=F, row.names=F)

#And third third
#length(c(2*round(nrow(int_DCM)/3):nrow(int_DCM)))
nrow(int_DCM)-(nrow(subs_A)+nrow(subs_B)) #245
length(c((nrow(subs_A)+nrow(subs_B)+1):nrow(int_DCM))) #245
subs_C=int_DCM[(nrow(subs_A)+nrow(subs_B)+1):nrow(int_DCM),]
dim(subs_C) #245 126
#Check: 
nrow(subs_A)+nrow(subs_B)+nrow(subs_C)==nrow(int_DCM) #YEAH!

colnames(subs_C)
colnames(subs_C)[1]="#OTU ID"
#Write out!
write.table(subs_C, file="eLSA_input/subs_C_eLSA_input_DCM_int_11.24.2020.tsv", sep="\t", quote=F, row.names=F)
