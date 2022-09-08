#Setting up eLSA!! (ASV data only)
#NOT CLR-transformed data, just interpolated data
#ver. 10.13.2020

rm(list=ls())
getwd()

#I. Read in files-------
#Read in interpolated, non-CLR transformed data from a) euks in the AE filter, b) proks in the AE filter, c) proks in the durapore filters for both 1. 5m (surface), 2. DCM

#0. Files relevant to both 5m and DCM: hashes and taxonomy----
#Read in hashes of abundant 18S data to include
#ASVs abundant at 5m first
abund_hashes_5m <- read.table("eLSA_input/abundant_18S_hashes_5m_all_UCYNA_ASVs_08202020.txt", header=T,  stringsAsFactors = F)
dim(abund_hashes_5m) #543 rows
names(abund_hashes_5m) #x
head(abund_hashes_5m) 

#Then ASVs abundant at DCM
abund_hashes_DCM <- read.table("eLSA_input/abundant_18S_hashes_DCM_all_UCYNA_ASVs_08202020.txt", header=T, stringsAsFactors = F)
dim(abund_hashes_DCM) #697 1 #Is 18S more diverse at the DCM? 
names(abund_hashes_DCM)

#Hashes associated in other publications 
lit_hashes <- read.table("eLSA_input/18S_assoc_ASVs.txt", header=F, stringsAsFactors = F)
dim(lit_hashes) #79 rows #Does not have a header
names(lit_hashes) #V1
head(lit_hashes)

#Can I combine them?
length(abund_hashes_5m$x) #543
length(lit_hashes$V1) #79
length(unique(c(abund_hashes_5m$x, lit_hashes$V1))) #568 #54 ASVs overlaop

#Yes I totally can 
hashes_5m <- unique(c(abund_hashes_5m$x, lit_hashes$V1))
dim(hashes_5m)
length(hashes_5m) #568
#Remove that troublesome bivalve (see "TESTING" below)
grep("ff51151c2a1e7919ecb5b8f1f9e4262a", hashes_5m)
hashes_5m=hashes_5m[-grep("ff51151c2a1e7919ecb5b8f1f9e4262a", hashes_5m)]
length(hashes_5m) #567

#Repeat for DCM
hashes_DCM <- unique(c(abund_hashes_DCM$x, lit_hashes$V1))
length(hashes_DCM) #722 
length(abund_hashes_DCM$x)+length(lit_hashes$V1) #Again, overlap of 54 ASVs
#Remove that troublesome bivalve (see "TESTING" below)
grep("ff51151c2a1e7919ecb5b8f1f9e4262a", hashes_DCM)
hashes_DCM <- hashes_DCM[-grep("ff51151c2a1e7919ecb5b8f1f9e4262a", hashes_DCM)]
length(hashes_DCM) #721

#Taxonomy
proks_tax <- read.table("ModifiedFiles/Proks_tax_classified_23072020_SILVA132.tsv", sep="\t", header=T, stringsAsFactors = F)
dim(proks_tax) #78953 3
head(proks_tax)

#Euks tazonomy
euks_tax <- read.table("ModifiedFiles/Euks_tax_classified_19052020_SILVA132.tsv", sep="\t", header=T, stringsAsFactors = F)

#Set up UCYNA variables
proks_tax$Taxon[grep("UCYN-A", proks_tax$Taxon)]
UCYNA_ASVs <- proks_tax$Feature.ID[grep("UCYN-A", proks_tax$Taxon)]
UCYNA_ASVs

#1. 5m (surface) samples ----

#Euks first
euks_AE_5m <- read.table("int_CLR_data/euks_AE_5m_int_10.09.2020.tsv", sep="\t", header=T, stringsAsFactors = F, row.names=1) #Make ASV hashes equal to rownames
dim(euks_AE_5m) #28402 125
head(colnames(euks_AE_5m))

#Then proks AE
proks_AE_5m <- read.table("int_CLR_data/proks_AE_5m_int_10.08.2020.tsv", header=T, stringsAsFactors = F, row.names=1)
dim(proks_AE_5m) #78953 125
head(colnames(proks_AE_5m))

#Last, proks durapore
proks_D_5m <- read.table("int_CLR_data/proks_D_5m_int_10.08.2020.tsv", header=T, stringsAsFactors = F, row.names=1)
dim(proks_D_5m) #78953 124 <- one column less, missing the first date (March 2008, I think)
head(colnames(proks_D_5m)) #Yes

#Are the colnames all the same? (Except for one missing durapore date)
grep("FALSE", colnames(euks_AE_5m)==colnames(proks_AE_5m)) #All the same
colnames(euks_AE_5m)[-1]==colnames(proks_D_5m) #All different, but only because one has "D" the other "AE"


#2. DCM samples-----
#Euks first
euks_AE_DCM <- read.table("int_CLR_data/euks_AE_DCM_int_10.09.20.tsv", sep="\t", stringsAsFactors = F, header=T, row.names=1)
dim(euks_AE_DCM) #28402 125
head(colnames(euks_AE_DCM))

#Proks AE
proks_AE_DCM <- read.table("int_CLR_data/proks_AE_DCM_int_10.08.2020.tsv", sep="\t", stringsAsFactors = F, header=T, row.names=1)
dim(proks_AE_DCM) #78953 125
head(colnames(proks_AE_DCM))

#Proks durapore
proks_D_DCM <- read.table("int_CLR_data/proks_D_DCM_int_10.08.2020.tsv", sep="\t", stringsAsFactors = F, header=T, row.names=1)
dim(proks_D_DCM) #78953 125
head(colnames(proks_D_DCM))

#II. Setting up matrix for 5m (surface)-----

#1. Pull out relevant 18S hashes and put them into a new matrix----
#Create a dataframe for new data
eLSA_in_5m=as.data.frame(matrix(ncol=ncol(euks_AE_5m), nrow=(length(hashes_5m)+(length(UCYNA_ASVs)*2)))) #Add 12 to reflect we'll be adding in UCYN-A ASVs from both AE and durapore
dim(eLSA_in_5m) #579 125

#Set colnames
colnames(eLSA_in_5m)=colnames(euks_AE_5m)

#Set rownames - remember to include each UCYNA ASV + "AE" or "D"
paste(UCYNA_ASVs, "AE", sep="_") #Probably better to leave these IDs as hashes rather than numbers
length(c(hashes_5m, (paste(UCYNA_ASVs, "AE", sep="_")), (paste(UCYNA_ASVs, "D", sep="_"))))==nrow(eLSA_in_5m) #TRUE!
rownames(eLSA_in_5m)=c(hashes_5m, (paste(UCYNA_ASVs, "AE", sep="_")), (paste(UCYNA_ASVs, "D", sep="_")))
#Double check that
tail(rownames(eLSA_in_5m), n=20)
length(grep("D", rownames(eLSA_in_5m)))
length(grep("AE", rownames(eLSA_in_5m)))
#All looks good! 

#Fill in the matrix with the 18S data we have
for(i in c(1:length(hashes_5m))){
  eLSA_in_5m[grep(hashes_5m[i], rownames(eLSA_in_5m)),]=euks_AE_5m[grep(hashes_5m[i], rownames(euks_AE_5m)),]
}

#Double check that
for(i in sample(x=c(1:length(hashes_5m)), size=20)){
  print(sum(eLSA_in_5m[grep(hashes_5m[i], rownames(eLSA_in_5m)),]))
  print(sum(euks_AE_5m[grep(hashes_5m[i], rownames(euks_AE_5m)),]))
  print(sum(eLSA_in_5m[grep(hashes_5m[i], rownames(eLSA_in_5m)),])==sum(euks_AE_5m[grep(hashes_5m[i], rownames(euks_AE_5m)),]))
} #All true, all the same!

#2. Write in hashes corresponding to UCYN-A-----
#AE first
for(i in UCYNA_ASVs){
  eLSA_in_5m[grep(paste(i, "AE", sep="_"), row.names(eLSA_in_5m)),]=proks_AE_5m[grep(i, row.names(proks_AE_5m)),]
  #Double check it worked: 
  print(sum(eLSA_in_5m[grep(paste(i, "AE", sep="_"), row.names(eLSA_in_5m)),]))
  print(sum(proks_AE_5m[grep(i, row.names(proks_AE_5m)),]))
} #Great, all the sums match! It worked! 

#Now add in relative abundances from durapore
for(i in UCYNA_ASVs){
  eLSA_in_5m[grep(paste(i, "D", sep="_"), row.names(eLSA_in_5m)),]=proks_D_5m[grep(i, row.names(proks_D_5m)),]
  #Double check it worked: 
  print(sum(eLSA_in_5m[grep(paste(i, "D", sep="_"), row.names(eLSA_in_5m)),]))
  print(sum(proks_D_5m[grep(i, row.names(proks_D_5m)),]))
} #They all match, but a lot of the CLR transformed abundances are 0... should I also do an eLSA without the CLR transform, i.e. just using interpolated data?

#What happened to the missing data on 2008.03?
tail(eLSA_in_5m$SPOT.2008.03.01.5m.AE, n=12)
#It just got set to zero
#Should actually be "na"'s 
for(i in UCYNA_ASVs){
  eLSA_in_5m$SPOT.2008.03.01.5m.AE[grep(paste(i, "D", sep="_"), row.names(eLSA_in_5m))]="na"
}
tail(eLSA_in_5m$SPOT.2008.03.01.5m.AE) #Okay great, now just have to remove the quotes when I write out

#3. Format matrix a little and write it out! -----
#Add in OTU ID as a column 
tail(rownames(eLSA_in_5m), n=20)
eLSA_in_5m$OTU_ID=rownames(eLSA_in_5m)
dim(eLSA_in_5m)

#Re-order it
order(colnames(eLSA_in_5m))
writeout_eLSA_in_5m=eLSA_in_5m[,order(colnames(eLSA_in_5m))]
head(colnames(writeout_eLSA_in_5m))
head(writeout_eLSA_in_5m$OTU_ID)

#Change OTU_ID to be "#OTU ID"
colnames(writeout_eLSA_in_5m)[1]="#OTU ID"

#Write it out!
write.table(x=writeout_eLSA_in_5m, file="eLSA_input/eLSA_input_5m_int_abund_and_lit_10.13.2020.tsv", sep="\t", quote=F, row.names=F)

#III. Setting up matrix for DCM samples-----

#1. Pull out 18S hashes and put them into a new matrix-----
#Create new dataframe
eLSA_in_DCM=as.data.frame(matrix(nrow=(length(hashes_DCM)+(length(UCYNA_ASVs)*2)), ncol=ncol(euks_AE_DCM)))
dim(eLSA_in_DCM)
length(hashes_DCM)

#Set colnames
colnames(eLSA_in_DCM)=colnames(euks_AE_DCM)

#Set rownames -- for UCYNA, AE first, then durapore
paste(UCYNA_ASVs, "AE", sep="_")
row.names(eLSA_in_DCM)=c(hashes_DCM, paste(UCYNA_ASVs, "AE", sep="_"), paste(UCYNA_ASVs, "D", sep="_"))

#Check rownames
tail(rownames(eLSA_in_DCM), n=20)
length(grep("AE", rownames(eLSA_in_DCM)))
length(grep("D", rownames(eLSA_in_DCM)))

#Fill in spreadsheet with 18S data
for(i in c(1:length(hashes_DCM))){
  eLSA_in_DCM[grep(hashes_DCM[i], row.names(eLSA_in_DCM)),]=euks_AE_DCM[grep(hashes_DCM[i], row.names(euks_AE_DCM)),]
}
View(eLSA_in_DCM)

#Double check this
for(i in sample(1:length(hashes_DCM), size=20)){
  print(sum(eLSA_in_DCM[grep(hashes_DCM[i], rownames(eLSA_in_DCM)),]))
  print(sum(euks_AE_DCM[grep(hashes_DCM[i],rownames(euks_AE_DCM)),]))
  print(sum(eLSA_in_DCM[grep(hashes_DCM[i], rownames(eLSA_in_DCM)),])==sum(euks_AE_DCM[grep(hashes_DCM[i],rownames(euks_AE_DCM)),]))
} #All true, all the same!


#2. Add in UCYN-A data from AE and Durapore----

#AE first
for(i in UCYNA_ASVs){
  eLSA_in_DCM[grep(paste(i, "AE", sep="_"), row.names(eLSA_in_DCM)),]=proks_AE_DCM[grep(i, row.names(proks_AE_DCM)),]
  print(sum(eLSA_in_DCM[grep(paste(i, "AE", sep="_"), row.names(eLSA_in_DCM)),]))
  print(sum(proks_AE_DCM[grep(i, row.names(proks_AE_DCM)),]))
} #YUP all the sums match! 

#Durapore second
for(i in UCYNA_ASVs){
  eLSA_in_DCM[grep(paste(i, "D", sep="_"), row.names(eLSA_in_DCM)),]=proks_D_DCM[grep(i, row.names(proks_D_DCM)),]
  print(sum(eLSA_in_DCM[grep(paste(i, "D", sep="_"), row.names(eLSA_in_DCM)),]))
  print(sum(proks_D_DCM[grep(i, row.names(proks_D_DCM)),]))
} #OK much better! #A lot of the abundances are zero...

#3. Format a little and write it out!----
#Add in "OTU_ID" as a column
dim(eLSA_in_DCM) #733 125
eLSA_in_DCM$OTU_ID=row.names(eLSA_in_DCM)
dim(eLSA_in_DCM) #733 126
tail(colnames(eLSA_in_DCM))
order(colnames(eLSA_in_DCM))

#Re-order matrix
writeout_eLSA_in_DCM=eLSA_in_DCM[,order(colnames(eLSA_in_DCM))]
head(colnames(writeout_eLSA_in_DCM))

#Add in hash before "OTU_ID"
colnames(writeout_eLSA_in_DCM)[1]="#OTU_ID"
head(colnames(writeout_eLSA_in_DCM))

#Write it out! 
write.table(x=writeout_eLSA_in_DCM, file="eLSA_input/eLSA_input_DCM_int_abund-lit_10.13.2020.tsv", sep="\t", quote=F, row.names=F)

