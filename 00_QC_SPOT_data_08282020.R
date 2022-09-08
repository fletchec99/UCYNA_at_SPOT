#Trying to figure out what is going on in samples with fewer than 100 18S ASVs
#8.24.20

#Read in data----
setwd("/Users/cfletcherh1/Desktop/DataAnalysis/SPOT_16S_18S_eLSA/")

#Read in proks data: 
proks_data <- read.table("ModifiedFiles/3.SPOT_16S_w.chloro_proportions.tsv", header=T, stringsAsFactors = F, sep="\t")
dim(proks_data) #78953 1208
head(names(proks_data))
length(grep("5m.AE", colnames(proks_data))) #125 
length(grep("DCM.AE", colnames(proks_data))) #118

#Proks data with raw counts + w chloroplasts
proks_data <- read.table("ModifiedFiles/1.SPOT_16S_w.chloro_counts.tsv", header=T, stringsAsFactors = F, sep=" ")
dim(proks_data)
head(colnames(proks_data))
length(grep("5m.AE", colnames(proks_data))) #125
length(grep("DCM.AE", colnames(proks_data))) #118

proks_tax <- read.table("ModifiedFiles/Proks_tax_classified_23072020_SILVA132.tsv", header=T, stringsAsFactors = F, sep="\t")

#Read in euks data: 
euks_data <- read.table("ModifiedFiles/8.SPOT_18S_no_metaz_proportions.tsv", header=T, stringsAsFactors = F, sep="\t")
dim(euks_data) #28402 1098
length(grep("5m.AE", colnames(euks_data))) #124
length(grep("DCM.AE", colnames(euks_data))) #116

#Read in data with raw counts, not proportions (without metazoans)
euks_data <- read.table("ModifiedFiles/6.SPOT_18S_no_metaz_counts.tsv", header=T, stringsAsFactors = F, sep="\t")
dim(euks_data) #28402 1098

#Read in euks data with metazoans, proportions: 
euks_data_w_metaz <- read.table("ModifiedFiles/7.SPOT_18S_w.metaz_proportions.tsv", header=T, stringsAsFactors = F, sep = "\t")
dim(euks_data_w_metaz)

euks_tax <- read.table("ModifiedFiles/Euks_tax_classified_19052020_SILVA132.tsv", header=T, stringsAsFactors = F, sep="\t")

#Figure out which AE samples have fewer than 100 18S taxa----
names <- c()
no_taxa <- c()
for(i in grep("AE", colnames(euks_data))){
  if(length(which(euks_data[,i]>0))<100){
    names <- c(names, colnames(euks_data)[i])
    no_taxa <- c(no_taxa, length(which(euks_data[,i]>0)))
  }
}
length(names) #195
length(no_taxa) #195
length(grep("AE", colnames(euks_data))) #473

#Write this out in a data frame
troublesome_samples <- as.data.frame(matrix(ncol=2, nrow=195))
names(troublesome_samples) <- c("SampleID", "No_18S_taxa")
troublesome_samples$SampleID=names
troublesome_samples$No_18S_taxa=no_taxa

grep("5m", names)
troublesome_samples[grep("5m", names),]
troublesome_samples[grep("DCM", names),]

#Add in depth
for(i in c(1:195)){
  a <- strsplit(troublesome_samples$SampleID[i], split=".", fixed=T)
  troublesome_samples$Depth[i]=(a[[1]][5])
}
View(troublesome_samples)

#which run did these samples come from? 
list.files("OriginalFiles/")

euks_test <- read.table("OriginalFiles/18S_mock_SPOT_Run40_Run58Lane2_dada2q2_feature-table.biom.tsv", header=T, stringsAsFactors = F, sep="\t")
dim(euks_test) #30537 1251
head(colnames(euks_test))
colnames(euks_test)[10:20]

#Convert back to old cols format
oldcolnames <- c()
for(i in troublesome_samples$SampleID){
  a <- strsplit(i, split=".", fixed=T)
  #print(a)
  if(as.numeric(a[[1]][3])<10){
    b <- (strsplit(a[[1]][3], split="0"))
    oldcolnames <- c(oldcolnames, paste("AE", b[[1]][2], a[[1]][4], a[[1]][2], "SPOT", a[[1]][5], "Monthlies.515Y.926R", sep="."))
  } else{
    oldcolnames <- c(oldcolnames, paste("AE", b[[1]][2], a[[1]][4], a[[1]][2], "SPOT", a[[1]][5], "Monthlies.515Y.926R", sep="."))
  }
}
length(oldcolnames)
length(unique(oldcolnames))
length(troublesome_samples$SampleID)
head(oldcolnames)

#Figure out run number and old sample ID name
old_col_names <- c()
run_no <- c()
for(i in oldcolnames){
  if(length(grep(i, colnames(euks_test)))>0){
    a <- colnames(euks_test)[(grep(i, colnames(euks_test)))]
    old_col_names <- c(old_col_names, a)
    b <- strsplit(a, split=".", fixed=T)
    run_no <- c(run_no, b[[1]][1])
  } else {
    d <- "NoMatch"
    old_col_names <- c(old_col_names, d)
    run_no <- c(run_no, "NA")
  }
}
old_col_names
run_no
length(old_col_names)
length(run_no)

#Add into table
troublesome_samples$Old_col_name=old_col_names[1:195]
troublesome_samples$RunNo=run_no

#Write out file
write.csv(x=troublesome_samples, file="Samples_w_fewer_100_18S_taxa_08212020.csv", row.names=F, quote=F)

#How many 18S and 16S ASVs does each AE sample actually have?----- 
colnames <- c()
numb_ASVs <- c()
for(i in grep("AE", colnames(euks_data))){
  colnames <- c(colnames, colnames(euks_data)[i])
  numb_ASVs <- c(numb_ASVs, length(which(euks_data[,i]>0)))
}
length(grep("AE", colnames(euks_data)))
length(colnames)
length(numb_ASVs)

#Figure out the number of 16S ASVs these AEs have
numb_16S_ASVs <- c()
for(i in colnames){
  a <- grep(i, colnames(proks_data))
  numb_16S_ASVs <- c(numb_16S_ASVs, length(which(proks_data[,a]>0)))
}

#Write this into a table
ASVs_per_sample <- as.data.frame(matrix(ncol=3, nrow=length(numb_ASVs)))
names(ASVs_per_sample)=c("SampleID", "Number_18S_ASVs", "Number_16S_ASVs")
ASVs_per_sample$SampleID=colnames
ASVs_per_sample$Number_18S_ASVs=numb_ASVs
ASVs_per_sample$Number_16S_ASVs=numb_16S_ASVs

#Wait, hold on, add in 18S sequences WITH metazoans
euks_data_w_metaz <- read.table("ModifiedFiles/7.SPOT_18S_w.metaz_proportions.tsv", header=T, stringsAsFactors = F, sep = "\t")
dim(euks_data_w_metaz)

numb_ASVs_w_metaz <- c()
for(i in ASVs_per_sample$SampleID){
  a <- grep(i, colnames(euks_data_w_metaz))
  numb_ASVs_w_metaz <- c(numb_ASVs_w_metaz, length(which(euks_data_w_metaz[,a]>0)))
}

ASVs_per_sample$Number_18S_ASVS_w_metaz=numb_ASVs_w_metaz

write.table(x=ASVs_per_sample, file="Data_QC/Number_18S_16S_ASVs_present_per_AE_sample_08242020.csv", row.names=F, quote=F)

#Plot it up-------
#Histogram of number of 18S ASVs, all depths
hist(ASVs_per_sample$Number_18S_ASVs, breaks=22, xlab="Number of 18S ASVs", main="Number of 18S ASVs per sample \n excluding metazoans, all depths")
summary(ASVs_per_sample$Number_18S_ASVs)

hist(ASVs_per_sample$Number_18S_ASVs, breaks=220, xlab="Number of 18S ASVs", main="Number of 18S ASVs per sample \n excluding metazoans, all depths")
summary(ASVs_per_sample$Number_18S_ASVs)

#Histogram of number of 18S ASVs, 5m only
length(grep("5m", ASVs_per_sample$SampleID))
hist(ASVs_per_sample$Number_18S_ASVs[grep("5m", ASVs_per_sample$SampleID)], breaks=11, xlab="Number of 18S ASVs", main="Number of 18S ASVs per sample \n excluding metazoans, 5m depth")

#Histogram of number of 18S ASVs, DCM only 
hist(ASVs_per_sample$Number_18S_ASVs[grep("DCM", ASVs_per_sample$SampleID)], breaks=11, xlab="Number of 18S ASVs", main="Number of 18S ASVs per sample \n excluding metazoans, DCM depth")

#Plot number of 18S ASVs vs. 16S ASVs
plot(ASVs_per_sample$Number_18S_ASVs, ASVs_per_sample$Number_16S_ASVs, pch=16, xlab="Number of 18S ASVs present (rel. abundance > 0) in \n AE samples without metazoans", ylab="Number of 16S ASVs present (rel. abundance > 0) in \n AE samples with chloroplasts")

#Histogram of number of 16S ASVs all depths
hist(ASVs_per_sample$Number_16S_ASVs, breaks=44, xlab="Number of 16S ASVs", main="Number of 16S ASVs per sample \n with chloroplasts, all depths")

#Histogram of number of 16S ASVs 5m depth only
hist(ASVs_per_sample$Number_16S_ASVs[grep("5m", ASVs_per_sample$SampleID)], breaks=22, xlab="Number of 16S ASVs", main="Number of 16S ASVs per sample \n with chloroplasts, 5m")

#Same but DCM depth only 
hist(ASVs_per_sample$Number_16S_ASVs[grep("DCM", ASVs_per_sample$SampleID)], breaks=22, xlab="Number of 16S ASVs", main="Number of 16S ASVs per sample \n with chloroplasts, DCM")

#Histogram of number of 18S ASVs per sample with metazoans included
hist(ASVs_per_sample$Number_18S_ASVS_w_metaz, breaks=22, xlab="Number of 18S ASVs", main="Number of 18S ASVs per sample \n including metazoans, all depths")

#Plot up one vs other
plot(ASVs_per_sample$Number_18S_ASVS, ASVs_per_sample$Number_18S_ASVS_w_metaz, pch=16, xlab="Number of 18S ASVs \n excluding metazoans", ylab="Number 18S ASVs \n including metazoans")
#Okay that is a STRAIGHT line, so CLEARLY it doesn't matter if metazoans are included or not! 

#Which samples are in the euks data with metazoans but not in the euks data without metazoans? -----
#(Which samples have only metazoans)
problems <- c()
for(i in colnames(euks_data_w_metaz)){
  a <- grep(i, colnames(euks_data))
  if(length(a)<1){
    problems <- c(problems, i)
  }
}
length(problems)

#Which run are they from? 
oldcolnames <- c()
for(i in problems){
  a <- strsplit(i, split=".", fixed=T)
  if(as.numeric(a[[1]][3])<10){
    b <- strsplit(a[[1]][3], split="0")
    oldcolnames <- c(oldcolnames, paste("D", b[[1]][2], a[[1]][4], a[[1]][2], "SPOT", a[[1]][5], "Monthlies.515Y.926R", sep="."))
  } else{
    oldcolnames <- c(oldcolnames, paste("D", b[[1]][2], a[[1]][4], a[[1]][2], "SPOT", a[[1]][5], "Monthlies.515Y.926R", sep="."))
  }
}
length(oldcolnames)

for(i in oldcolnames){
  #print(grep(i, colnames(euks_test)))
  print(colnames(euks_test)[grep(i, colnames(euks_test))])
}
#Manually go in and fix colnames with the wrong date
oldcolnames[4]="D.9.22.2004.SPOT.500m.Monthlies.515Y.926R.a"
oldcolnames[6]="D.6.14.2005.SPOT.890m.Monthlies.515Y.926R"
oldcolnames[7]="D.6.14.2005.SPOT.500m.Monthlies.515Y.926R"
oldcolnames[8]="D.6.14.2005.SPOT.890m.Monthlies.515Y.926R"
oldcolnames[11]="D.4.29.2011.SPOT.500m.Monthlies.515Y.926R"

#How many samples have <50 ASVs and how many of those are from 5m/ DCM? -----
length(which(troublesome_samples$No_18S_taxa<50)) #101
length(grep("5m", troublesome_samples$SampleID[which(troublesome_samples$No_18S_taxa<50)])) #19
length(grep("DCM", troublesome_samples$SampleID[which(troublesome_samples$No_18S_taxa<50)])) #15

#How does number of 18S ASVs vary with sequencing depth? ----
#Read in the spreadsheet that Liv sent me with the number of sequences
num_sequences <- read.table("Data_QC/Samples_w_fewer_100_18S_taxa_08212020_Liv.csv", header=T, stringsAsFactors = F, sep =",")
dim(num_sequences)

length(grep("NA", num_sequences$Liv_NumAllReads))
length(which(num_sequences$Liv_NumAllReads>0)) #179

plot(num_sequences$Liv_NumAllReads[which(num_sequences$Liv_NumAllReads>0)], num_sequences$No_18S_ASVs[which(num_sequences$Liv_NumAllReads>0)], xlab="Number of reads (16S & 18S)", ylab="Number 18S ASVs", main ="Number 18S ASVs vs Total Number Reads", pch=16)
abline(v=51770, col="red", lwd=2)
text(4E+05, 10, " <- Median number reads= 51,770", col="red")

summary(num_sequences$Liv_NumAllReads[which(num_sequences$Liv_NumAllReads>0)]) #Median number of reads=51,770

#Histogram of number of reads
hist(num_sequences$Liv_NumAllReads[which(num_sequences$Liv_NumAllReads>0)], breaks=20, xlab="Number of reads per sample", main="Histogram of number of reads (16S & 18S) of \n samples with <100 18S ASVs")

#Plot number 18S ASVs by number 18S reads
plot(num_sequences$Liv_Num18SReads[which(num_sequences$Liv_Num18SReads>0)], num_sequences$No_18S_ASVs[which(num_sequences$Liv_Num18SReads>0)], pch=16, xlab="Number of reads (18S only)", ylab="Number of 18S ASVs", main="Number of reads (18S) vs. number of 18S ASVs\n in weird AE samples")

#Histogram of number of 18S reads
hist(num_sequences$Liv_Num18SReads[which(num_sequences$Liv_Num18SReads>0)], breaks=20, xlab=" Number 18S Reads/ Sample", main="Histogram of number of reads (18S only) \n of samples with <100 18S ASVs", pch=16 )

#Number 18S reads vs number all reads
plot(num_sequences$Liv_NumAllReads[which(num_sequences$Liv_NumAllReads>0)], num_sequences$Liv_Num18SReads[which(num_sequences$Liv_Num18SReads>0)], pch=16, xlab="Number All Reads", ylab="Number 18S Reads", main="Number 18S reads vs. Number all reads")

#Rarefaction curves with vegan-----
install.packages("vegan")
library(vegan)

#First, an example
data(BCI)
S <- specnumber(BCI) # observed number of species
(raremax <- min(rowSums(BCI))) #raremax=340
Srare <- rarefy(BCI, raremax)
plot(S, Srare, xlab = "Observed No. of Species", ylab = "Rarefied No. of Species")
abline(0, 1)
rarecurve(BCI, step = 20, sample = raremax, col = "blue", cex = 0.6)

#Then need to create a matrix of SPOT data (without metazoans)
dim(euks_data)
#Let's do one for all AE euks samples
AE_all <- as.matrix(euks_data[,grep("AE", colnames(euks_data))])
dim(AE_all)
#Remove rows that sum to zero
length(which(rowSums(AE_all)==0)) #840 rows sum to zero
AE_all <- AE_all[-which(rowSums(AE_all)==0),]
dim(AE_all)
#Transform the matrix so that species (ASVs) are columns and AE samples are rows
AE_all <- t(AE_all)
dim(AE_all)
head(AE_all)

S <- specnumber(AE_all)
raremin <- min(rowSums(AE_all)) #still 1? 
Srare <- rarefy(AE_all, raremin) 
plot(S, Srare, xlab = "Observed No. of Species", ylab = "Rarefied No. of Species")
rarecurve(AE_all, col = "blue", step = 20, sample = raremax, cex = 0.6, main="Rarefaction curves of 18S ASVs \n from AE samples from all depths")
#This is taking way too long

#Rarefy AE samples from 5m only-Euks---- 
AE_5m <- as.matrix(euks_data[,grep("5m.AE", colnames(euks_data))])
dim(AE_5m)
#Again, remove ASVs with 0 abundance across all samples
length(which(rowSums(AE_5m)==0))
length(which(rowSums(AE_5m)>0))
length(which(rowSums(AE_5m)==0))+length(which(rowSums(AE_5m)>0)) == nrow(AE_5m)

AE_5m <- AE_5m[which(rowSums(AE_5m)>0),]
dim(AE_5m)
#Transpose
AE_5m <- t(AE_5m)
dim(AE_5m) #124 8657

#Now rarefy
S <- specnumber(AE_5m)
raremin <- min(rowSums(AE_5m)) #Raremin=4 
Srare <- rarefy(AE_5m, raremin)
plot(S, Srare, xlab = "Observed No. of Species", ylab = "Rarefied No. of Species")
rarecurve(AE_5m, col = "blue", step = 20, sample = raremin, cex = 0.6, main="Rarefaction curves of 18S ASVs \n from AE samples from 5m", label=F)

#Subset to just take samples with <200 ASVs? 
length(which(rowSums(AE_5m)<100))
length(which(rowSums(AE_5m)<200))
length(which(AE_5m))
sort(rowSums(AE_5m)) 

#Okay, let's do a for-loop
no_ASVs_5m <- c()
for(i in c(1:nrow(AE_5m))){
  no_ASVs_5m <- c(no_ASVs_5m, length(which(AE_5m[i,]>0)))
}
length(no_ASVs_5m) #124
length(which(no_ASVs_5m<100)) #43, perfect! 
#Indexing should be the same between the vector and the matrix

#Now subset to only include AE samples from 5m with <100 18S ASVs
AE_5m_subs <- AE_5m[which(no_ASVs_5m<100),]
dim(AE_5m_subs)

S <- specnumber(AE_5m_subs)
raremin <- min(rowSums(AE_5m_subs)) #Raremin=4 
Srare <- rarefy(AE_5m_subs, raremin)
plot(S, Srare, xlab = "Observed No. of Species", ylab = "Rarefied No. of Species", pch=16)
rarecurve(AE_5m_subs, col = "blue", step = 20, sample = raremax, cex = 0.6, main="Rarefaction curves of 18S ASVs \n AE samples with <100 18S ASVs from 5m", label=F)

#So it looks like 20 18S ASVs is a reasonable cut-off
#How many 5m AE samples have <20 ASVs? 11
#Okay, let's repeat with samples with <10 18S ASVs
temp=troublesome_samples[which(troublesome_samples$No_18S_taxa<10),]
dim(temp)
length(grep("5m", temp$SampleID))
temp$SampleID[grep("5m", temp$SampleID)]

#Let's try to rarefy just samples with fewer than 20 ASVs
AE_5m_subs <- AE_5m[which(no_ASVs_5m<20),]
dim(AE_5m_subs)

S <- specnumber(AE_5m_subs)
raremin <- min(rowSums(AE_5m_subs)) #Raremin=4 
Srare <- rarefy(AE_5m_subs, raremin)
plot(S, Srare, xlab = "Observed No. of Species", ylab = "Rarefied No. of Species", pch=16)
rarecurve(AE_5m_subs, col = "blue", step = 20, sample = raremax, cex = 0.6, main="Rarefaction curves of 18S ASVs \n AE samples with <20 18S ASVs from 5m", label=F)

#Wow, ok. Actually 4 ASVs looks like a more reasonable cut-off.
temp=troublesome_samples[which(troublesome_samples$No_18S_taxa<4),]
dim(temp)
length(grep("5m", temp$SampleID)) #So that's only three 5m AE samples! 
temp$SampleID[grep("5m", temp$SampleID)]

#9.4.20: Now subset to only include AE samples from 5m with <15 18S ASVs
AE_5m_subs <- AE_5m[which(no_ASVs_5m<15),]
dim(AE_5m_subs)  #11 8657

S <- specnumber(AE_5m_subs)
raremin <- min(rowSums(AE_5m_subs)) #Raremin=4 
Srare <- rarefy(AE_5m_subs, raremin)
raremin
rarecurve(AE_5m_subs, col = "blue", step = 20, sample = raremin, cex = 0.6, main="Rarefaction curves of 18S ASVs \n AE samples with <15 18S ASVs from 5m", label=T)

#Rarefy AE samples from DCM only-Euks----
AE_DCM <- as.matrix(euks_data[,grep("DCM.AE", colnames(euks_data))])
dim(AE_DCM)
#Remove ASVs that are not present
which(rowSums(AE_DCM)==0)
AE_DCM <- AE_DCM[-which(rowSums(AE_DCM)==0),]
dim(AE_DCM) #10191 116
#Transpose #I don't think this makes a difference? 
AE_DCM <- t(AE_DCM)
dim(AE_DCM) #116 10191

#Now rarefy
S <- specnumber(AE_DCM)
raremin <- min(rowSums(AE_DCM)) #Raremin=2? 
Srare <- rarefy(AE_DCM, raremin)
plot(S, Srare, xlab = "Observed No. of Species", ylab = "Rarefied No. of Species")
rarecurve(AE_DCM, col = "blue", step = 20, sample = raremin, cex = 0.6, main="Rarefaction curves of 18S ASVs \n from AE samples from DCM", label=F)

#Subset to take fewer than 100 18S ASVs 
no_ASVs_DCM <- c()
for(i in c(1:nrow(AE_DCM))){
  no_ASVs_DCM <- c(no_ASVs_DCM, length(which(AE_DCM[i,]>0)))
}
length(no_ASVs_DCM) #116
nrow(AE_DCM) #116

#Subset
length(which(no_ASVs_DCM<100)) #33
AE_DCM_subs <- AE_DCM[which(no_ASVs_DCM<100),]
dim(AE_DCM_subs)

#Now rarefy
S <- specnumber(AE_DCM_subs)
raremin <- min(rowSums(AE_DCM_subs)) #Raremin=8 
Srare <- rarefy(AE_DCM_subs, raremin)
plot(S, Srare, xlab = "Observed No. of Species", ylab = "Rarefied No. of Species")
rarecurve(AE_DCM_subs, col = "blue", step = 20, sample = raremax, cex = 0.6, main="Rarefaction curves of 18S ASVs \n from AE samples with <100 18S taxa from DCM", label=F)

#Subset to samples with fewer than 20 ASVs
length(which(no_ASVs_DCM<20)) #33
AE_DCM_subs <- AE_DCM[which(no_ASVs_DCM<20),]
dim(AE_DCM_subs)

#Now rarefy
S <- specnumber(AE_DCM_subs)
raremin <- min(rowSums(AE_DCM_subs)) #Raremin=8 
Srare <- rarefy(AE_DCM_subs, raremin)
plot(S, Srare, xlab = "Observed No. of Species", ylab = "Rarefied No. of Species")
rarecurve(AE_DCM_subs, col = "blue", step = 20, sample = raremax, cex = 0.6, main="Rarefaction curves of 18S ASVs \n from AE samples with <20 18S taxa from DCM", label=F)

#Rarefy to take samples with <15 taxa
#Subset
length(which(no_ASVs_DCM<15)) #7
AE_DCM_subs <- AE_DCM[which(no_ASVs_DCM<15),]
dim(AE_DCM_subs) #4 10191

#Now rarefy
S <- specnumber(AE_DCM_subs)
raremin <- min(rowSums(AE_DCM_subs)) #Raremin=8
raremin
Srare <- rarefy(AE_DCM_subs, raremin)
rarecurve(AE_DCM_subs, col = "blue", step = 20, sample = raremin, cex = 0.6, main="Rarefaction curves of 18S ASVs \n from AE samples with <15 18S taxa from DCM", label=T)

#Rarefy AE samples from 5m only -Proks-----
AE_5m <- as.matrix(proks_data[,grep("5m.AE", colnames(proks_data))])
dim(AE_5m) #78953 125
#Remove ASVs that are not present
length(which(rowSums(AE_5m)==0))
AE_5m <- AE_5m[-which(rowSums(AE_5m)==0),]
dim(AE_5m) #17104 125
#Transpose 
AE_5m <- t(AE_5m)
dim(AE_5m) #125 17104

#Now rarefy
S <- specnumber(AE_5m)
raremin <- min(rowSums(AE_5m)) #Raremin=8
raremin #8029
Srare <- rarefy(AE_5m, raremin)
rarecurve(AE_5m, col = "blue", step = 20, sample = raremin, cex = 0.6, main="Rarefaction curves of 16S ASVs \n from AE samples from 5m", label=F)

#Subset to take samples with <1000 ASVs
no_ASVs_proks_5m <- c()
for(i in c(1:nrow(AE_5m))){
  no_ASVs_proks_5m <- c(no_ASVs_proks_5m, length(which(AE_5m[i,]>0)))
}
length(no_ASVs_proks_5m) #125 

length(which(no_ASVs_proks_5m>1000)) #33

AE_5m_subs <- AE_5m[which(no_ASVs_proks_5m>1000),]
dim(AE_5m_subs) #33 17104

#Now rarefy
S <- specnumber(AE_5m_subs)
raremin <- min(rowSums(AE_5m_subs))
raremin #82304
Srare <- rarefy(AE_5m_subs, raremin)
rarecurve(AE_5m_subs, col = "blue", step = 20, sample = raremin, cex = 0.6, main="Rarefaction curves of 16S ASVs \n from AE samples from 5m with <1000 16S ASVs", label=F)

#Which euks taxa are present in the troublesome samples (with >100 18S ASVs)?-----
taxonomy_data <- c()
for(i in troublesome_samples$SampleID[which(troublesome_samples$No_18S_taxa<20)]){
  taxonomy_data <- c(taxonomy_data, i)
  a <- grep(i, colnames(euks_data))
  for(i in euks_data$OTU_ID[which(euks_data[,a]>0)]){
    taxonomy_data <- c(taxonomy_data, euks_tax$Taxon[grep(i, euks_tax$Feature.ID)])
  }
}

#Write this into a CSV
write.table(x=taxonomy_data, file="Data_QC/18S_tax_of_AE_samples<20_18S_ASVs_08282020.csv", sep=",", quote=F, row.names=F)

#Plot up % Metazoans vs. number of 18S ASVs-------

#Figure out which rows correspond to metazoans
length(euks_tax$Taxon[grep("division_Metazoa;", euks_tax$Taxon)]) #2052

metaz_rows <- c()
for(i in euks_tax$Feature.ID[grep("division_Metazoa;", euks_tax$Taxon)]){
  metaz_rows <- c(metaz_rows, grep(i, euks_data_w_metaz$OTU_ID))
}
length(metaz_rows)

#Now sum across those rows for each sample
rel_metaz <- c()
for(i in ASVs_per_sample$SampleID){
  a <- grep(i, colnames(euks_data_w_metaz))
  rel_metaz <- c(rel_metaz, sum(euks_data_w_metaz[metaz_rows, a]))
}
length(rel_metaz)
length(ASVs_per_sample$SampleID) #473

#Write it into the data frame
ASVs_per_sample$Relabun_all_metaz=rel_metaz

#Write this file out 
write.table(x=ASVs_per_sample, file="Data_QC/Number_18S_16S_ASVs_relabun_metaz_present_per_AE_sample_08312020.csv", row.names=F, quote=F)

#Now plot number of 18S ASVs without metazoans vs. % metazoans
plot(ASVs_per_sample$Number_18S_ASVs, ASVs_per_sample$Relabun_all_metaz, pch=16, ylab="Relative abundance all metazoan ASVs \n (out of all 18S sequences)", xlab="Number of 18S ASVs \n excluding metazoans", main="Relative abundance metazoans vs. Number 18S ASVs \n SPOT AE Samples, all depths")

#Take the log of 18S taxa
ASVs_per_sample$Log_number_18S=log(x=ASVs_per_sample$Number_18S_ASVs, base=10)

#Plot that
plot(ASVs_per_sample$Log_number_18S, ASVs_per_sample$Relabun_all_metaz, pch=16, ylab="Relative abundance all metazoan ASVs \n (out of all 18S sequences)", xlab="Log10 number of 18S ASVs \n excluding metazoans", main="Relative abundance metazoans vs. Log # 18S ASVs \n SPOT AE Samples, all depths")
