##Formatting SPOT Data for eLSA - remove duplicates, etc
#Ver 7.31.20

#Program versions: 
sessionInfo()
#R version 3.5.1 (2018-07-02)
#Platform: x86_64-apple-darwin15.6.0 (64-bit)
#Running under: macOS  10.15.4

rm(list=ls())
getwd()
setwd("/Users/cfletcherh1/Desktop/DataAnalysis/SPOT_16S_18S_eLSA/")
list.files()
list.dirs()

#1. Raw counts of 16S data with chloroplasts included: 1.SPOT_w.chloro_counts.tsv --------
#1.1) Read in files
proks_data <- read.table("OriginalFiles/16S_mock_SPOT_Run03_Run58Lane2_dada2q2_feature-table.biom.tsv", sep="\t", header=T, stringsAsFactors = F)
dim(proks_data) #78953 1899
head(names(proks_data))

proks_tax <- read.table("ModifiedFiles/Proks_tax_classified_23072020_SILVA132.tsv", sep="\t", header=T, stringsAsFactors = F)
#Use taxonomy that I classified on 7.23.20 against SILVA132
dim(proks_tax) #78953 3 

#Read in taxonomy of chloroplast sequences
proks_tax_PR <- read.table("ModifiedFiles/Chloro_tax_classified_31072020_PhytoRef.tsv", sep="\t", header=T, stringsAsFactors = F)
#Use taxonomy that I classified on 7.31.20 against PR2
dim(proks_tax_PR) #6669 3 
#Exactly the same number as those classified by SILVA, yay :) 
names(proks_tax_PR)
head(proks_tax_PR)

#1.2) Remove blanks, mocks, old samples, some duplicates, etc. #Keep replicates that say "mb"
length(grep("Blank|Lab|515C|DOBD|mock|Mock|926R.b|926R.c|926R.d|926R.e|926R.f|ma|Ma|Rep2|Original|unknown|temp.notmerged",colnames(proks_data))) #675 #So some will still be duplicates after removing these
#Subset to exclude these
#Keep replicates that say "mb", but will have to rename them to just say 890m or 500m, etc.
subsetted_proks_data <- proks_data[-grep("Blank|Lab|515C|DOBD|mock|926R.b|926R.c|926R.d|926R.e|926R.f|ma|Ma|Rep2|Original|unknown|temp.notmerged",colnames(proks_data))]
dim(subsetted_proks_data) #78953 1224
colnames(subsetted_proks_data)

#1.3) Re-format colnames and figure out which samples are duplicates
#Re-format colnames thusly: word "SPOT" comes first, then year, month (but reformat month), date(but reformat date), depth, filter type
newcolnames <- c()
for(i in c(2:ncol(subsetted_proks_data))){
  a <- strsplit(colnames(subsetted_proks_data)[i], split=".", fixed=T)
  b <- as.numeric(a[[1]][3])
  d <- as.numeric(a[[1]][4])
  if(b<10){
    tempmonth <- paste("0", a[[1]][3], sep="")
  } else{
    tempmonth <- a[[1]][3]
  }
  if(d<10){
    tempdate <- paste("0", a[[1]][4], sep="")
  } else{
    tempdate <- a[[1]][4]
  }
  e <- paste(a[[1]][6],a[[1]][5],tempmonth,tempdate,a[[1]][7],a[[1]][2],sep=".")
  newcolnames <- c(newcolnames, e)
}
length(newcolnames) #1223
length(unique(newcolnames)) #1207 #so we gotta get rid of 16 columns #Or is it 1208? 
head(newcolnames)
tail(newcolnames)

#Add back in "OTU.ID" so that newcolnames exactly match original colnames of subsetted data
newcolnames <- c("OTU_ID",newcolnames)
length(newcolnames)
length(unique(newcolnames))

#Rename samples with "ma"
grep("[a-l,n-z]", newcolnames)
replace <- grep("[a-l,n-z]", newcolnames)
newcolnames[grep("[a-l,n-z]", newcolnames)]

newcolnames <- gsub(pattern="mb", replacement="m", newcolnames, fixed=T)
newcolnames <- gsub(pattern="Mb", replacement="M", newcolnames, fixed=T)
newcolnames[replace] #perfect! 

#Figure out which ones are duplicates
which(duplicated(newcolnames)==TRUE)
newcolnames[which(duplicated(newcolnames)==TRUE)]
du=colnames(subsetted_proks_data)[which(duplicated(newcolnames)==TRUE)] #This is only giving me the second duplicate out of two, but I want the first
#Create a list of column numbers to exclude 
badcols <- c()
for(i in du){
  #print(i)
  a <- strsplit(i, split=".", fixed=T)
  b <- paste(a[[1]][2],a[[1]][3],a[[1]][4],a[[1]][5],a[[1]][6],a[[1]][7],a[[1]][8],a[[1]][9],a[[1]][10],a[[1]][11],sep=".")
  d <- grep(b, colnames(subsetted_proks_data))
  #print(d)
  badcols <- c(badcols, d[1])
}
colnames(subsetted_proks_data)[badcols]
length(badcols) #16
#For some reason, the sample from 2009.6.18.5m is duplicated, but only reported once in the loop above #Because the primer combination was a typo - #729 lists 926R twice #Just on this one sample
#Want to get rid of column 622, not 729 is the most recent version of this sample
colnames(subsetted_proks_data)[728]
colnames(subsetted_proks_data)[621]

badcols[grep("728", badcols)]=621
badcols

#Make sure to keep "Modified" sample dates, which are: 1.14.2013, 10.17.2012, 5.24.2011, 6.22.2011, 7.21.2012 ##Not duplicated above, CHA-ching! 

#1.4a) Remove duplicated samples and keep old colnames - save this variable for now
subs_proks_data_oldcols=subsetted_proks_data[,-badcols]
dim(subs_proks_data_oldcols) #78953 1208
dim(subsetted_proks_data) #78953 1224
length(badcols) #16
ncol(subsetted_proks_data)-length(badcols)==ncol(subs_proks_data_oldcols) #TRUE!
head(names(subs_proks_data_oldcols))

#1.4b) Re-name columns and remove duplicated samples
subs_proks_data_newcols=subsetted_proks_data
dim(subs_proks_data_newcols) #78953 1224
dim(subsetted_proks_data) #78953 1224 
length(newcolnames)
colnames(subs_proks_data_newcols)=newcolnames
head(names(subs_proks_data_newcols))
head(subs_proks_data_newcols$OTU_ID)

#Testing whether old column names and new column names match
test=as.data.frame(matrix(ncol=2, nrow=length(newcolnames)))
names(test)=c("oldcolnames", "newcolnames")
test$oldcolnames=colnames(subsetted_proks_data)
test$newcolnames=newcolnames
View(test)

for(i in (sample(x=c(1:length(newcolnames)), size=10, replace=F))){
  print(colnames(subsetted_proks_data)[i])
  print(colnames(subs_proks_data_newcols)[i])
}

#Remove badcols (colnames won't match any more)
subs_proks_data_newcols=subs_proks_data_newcols[,-badcols]
dim(subs_proks_data_newcols) #78953 1208
which(duplicated(colnames(subs_proks_data_newcols))==TRUE) #integer(0), !! No duplicates
grep("[a-l,n-z]", newcolnames) #ALSO integer(0) 
dim(subsetted_proks_data) #78953 1224
dim(subs_proks_data_newcols) #78953 1208

#1.5) Re-order matrix to sort by date
order(colnames(subs_proks_data_newcols)) #OTU_ID will actually be first because all others start with "SPOT" :) 
ordered_subs_proks_data=subs_proks_data_newcols[,order(colnames(subs_proks_data_newcols))]
dim(ordered_subs_proks_data) #78953 1208

#1.6) Double check a few dates and taxa
UCYNA_ASVs=proks_tax$Feature.ID[grep("UCYN-A", proks_tax$Taxon)]
for(i in UCYNA_ASVs){
  print(grep(i, proks_data$OTU.ID))
}

UCYNA_ix <- c()
for(i in UCYNA_ASVs){
  UCYNA_ix <- c(UCYNA_ix, grep(i, ordered_subs_proks_data$OTU_ID))
}
UCYNA_ix
#UCYNA indices match between proks data and ordered_subs_proks_data

#Double check that the number of times UCYNA appears is consistent in proks data and subs ordered proks data
for(i in UCYNA_ix){
  print(rowSums(proks_data[i,] >0))
  print(rowSums(ordered_subs_proks_data[i,] >0))
}
#Lost ASV1 (3d85...) on 7 dates and ASV4 (af1bb...) on 8 dates

#Which samples? Any of them duplicates? 
colnames(proks_data)[which(proks_data[2131,]>0)] #I think a lot of these are duplicates of each other, i.e. samples sequenced 2x
colnames(ordered_subs_proks_data)[which(ordered_subs_proks_data[2131,]>0)]
grep("d", colnames(proks_data)[which(proks_data[2131,]>0)])
colnames(proks_data)[which(proks_data[2131,]>0)][12] #Okay, did lose one "modified" sample

for(i in UCYNA_ix){
  print(sum(proks_data[i,2:ncol(proks_data)]))
  print(sum(ordered_subs_proks_data[i,2:ncol(ordered_subs_proks_data)]))
}
#Lost ASV1 a total of 850 times and A4 100 times - why the big discrepancy?? 
#Probably by removing replicates, not blanks and mocks, etc.

for(i in (grep("2008.10.20", colnames(ordered_subs_proks_data)))){
  print(colnames(ordered_subs_proks_data)[i])
  print(sum(ordered_subs_proks_data[,i]))
}

for(i in (grep("10.20.2008", colnames(proks_data)))){
  print(colnames(proks_data)[i])
  print(sum(proks_data[,i]))
}
#Match, but they are out of order 

#1.7) Write out the file (counts of 16S data with chloroplasts)
write.table(x=ordered_subs_proks_data, file="ModifiedFiles/1.SPOT_16S_w.chloro_counts.tsv", row.names=F, quote=F)

#2. Raw counts of 16S data with no chloroplasts: 2.SPOT_16S_no_chloro_counts.tsv-----
#2.1) Figure out which rows correspond to chloroplasts
proks_tax$Taxon[grep("D_3__Chloroplast", proks_tax$Taxon)]
chloro_ASVs=proks_tax$Feature.ID[grep("D_0__Bacteria;D_1__Cyanobacteria;D_2__Oxyphotobacteria;D_3__Chloroplast", proks_tax$Taxon)]
length(chloro_ASVs)
length(unique(chloro_ASVs)) #6669

#Which of these correspond to PR chloroplast ASVs? 
grep("kingdom", proks_tax$Taxon) #integer(0)
#So chloroplast taxonomy does NOT start with the word "kingdom" 
dim(proks_tax_PR) #6669 3
head(proks_tax_PR$Taxon)
test_ASVs=proks_tax_PR$Feature.ID

head(sort(test_ASVs))
head(sort(chloro_ASVs))
#They look exactly the same...
sort(test_ASVs)==sort(chloro_ASVs)
grep("FALSE", sort(test_ASVs)==sort(chloro_ASVs)) #None! They are exactly the same!
#Can use the SILVA classification of taxonomy in place of the PR2 classification 

#Figure out which rows to exclude
exclude_rows <- c()
for(i in chloro_ASVs){
  exclude_rows <- c(exclude_rows, grep(i, ordered_subs_proks_data$OTU_ID))
}
length(exclude_rows) #6669
length(unique(exclude_rows))
head(exclude_rows, n=20)
head(sort(exclude_rows), n=20) #"sort" is different than "order"? 

#2.2) Exclude those rows
ordered_subs_proks_data_nochloro <- ordered_subs_proks_data[-exclude_rows,]
dim(ordered_subs_proks_data_nochloro) #72284 1208
nrow(ordered_subs_proks_data)-nrow(ordered_subs_proks_data_nochloro)==length(exclude_rows) #TRUE
head(names(ordered_subs_proks_data_nochloro))
head(names(ordered_subs_proks_data))

#Double check that no samples are missing data
check_sums <- c()
for(i in c(2:ncol(ordered_subs_proks_data_nochloro))){
  check_sums <- c(check_sums, sum(ordered_subs_proks_data_nochloro[,i]))
}
which(check_sums<1) #638 and 944
colnames(ordered_subs_proks_data_nochloro)[which(check_sums<1)+1]
sum(ordered_subs_proks_data_nochloro[,639])
sum(ordered_subs_proks_data_nochloro[,945])
#These are indeed both zero

#So remove columns 639 and 945
dim(ordered_subs_proks_data_nochloro) #72284 1208
ordered_subs_proks_data_nochloro <- ordered_subs_proks_data_nochloro[,-c(639, 945)]
dim(ordered_subs_proks_data_nochloro) #72284 1206

#2.3) Write out new dataset
write.table(x=ordered_subs_proks_data_nochloro, file="ModifiedFiles/2.SPOT_16S_no_chloro_counts.tsv", row.names=F, quote=F)

#2.4) Write out 16S data for normalization (i.e. subs_proks_data_oldcols) but without chloro's
#Exclude chloroplasts
length(exclude_rows) #6669
subs_proks_data_oldcols_nochloro <- subs_proks_data_oldcols[-exclude_rows,-c(639, 945)]
dim(subs_proks_data_oldcols_nochloro) #72284 1206
dim(ordered_subs_proks_data_nochloro) 
dim(subsetted_proks_data) #78953 1224 
head(names(subs_proks_data_oldcols_nochloro))

#Write out
write.table(x=subs_proks_data_oldcols_nochloro, file="ModifiedFiles/temp_normalize_16S_nochloro_counts.tsv", sep="\t", row.names=F, quote=F)

#3. Relative abundances of 16S data with chloroplasts: 3.SPOT_16S_w.chloro_proportions.tsv-----
#3.1) Write out a temporary file in the correct format for Jesse's Python script 

#Add in fake taxonomy 
tail(colnames(ordered_subs_proks_data))
dim(ordered_subs_proks_data) #Still 78953 1208
ordered_subs_proks_data$taxonomy="SAR11" #Must be lowercase "taxonomy"
#colnames(ordered_subs_proks_data)[grep("Taxonomy", colnames(ordered_subs_proks_data))]="taxonomy"
dim(ordered_subs_proks_data) #Now 78953 1209
head(ordered_subs_proks_data$taxonomy)

#Write this file out
write.table(x=ordered_subs_proks_data, file="ModifiedFiles/temp_SPOT_16S_w.chloro_counts.tsv", row.names=F, sep="\t", quote=F) #Apparently you really do need all these flags! 

#Read it back in without colnames
temp_proks <- read.table("ModifiedFiles/temp_SPOT_16S_w.chloro_counts.tsv", header=F, stringsAsFactors = F, sep="\t", row.names=NULL)
dim(temp_proks) #78954 1209

#Insert "#Constructed from biom file" and "#OTU ID"
colnames(temp_proks)[1]="#Constructed from biom file"
colnames(temp_proks)[2:ncol(temp_proks)]=""
temp_proks[1,1]="#OTU ID"

#Write this file out
write.table(x=temp_proks, file="ModifiedFiles/temp_SPOT_16S_w.chloro_counts_final.tsv", row.names=F, sep="\t", quote=F)

#Remove old file
file.remove("ModifiedFiles/temp_SPOT_16S_w.chloro_counts.tsv")

#3.2) SCP the temp file to Kraken and transform to proportions
#Don't forget to source activate qiime! 
#Download the output (with proportions) to a folder called ProportionsFiles/

#3.3) Read back in and re-format (delete "taxonomy" column)
#Make sure "#OTU ID" does not totally wipe out colnames
temp <- "ProportionsFiles/output_temp_SPOT_16S_w.chloro_proportions.tsv"
colnames <- scan(text=readLines(temp, 2), what="", quiet=T)
colnames <- colnames[-c(grep("OTU", colnames))]
colnames[1]="OTU_ID"
head(colnames) #That looks better! 
length(colnames) #why is this 2416?
#First column to exclude will be the third colname that does NOT have the word "SPOT" in it (the first two are "OTU_ID" and "taxonomy")
which(grepl("SPOT", colnames)==FALSE)[3] #1210
colnames=colnames[1:(which(grepl("SPOT", colnames)==FALSE)[3]-1)]
tail(colnames)
length(colnames) #THERE we go! 

#NOW read in file
temp_proks <- read.table(temp, col.names=colnames, stringsAsFactors = F, sep="\t")
dim(temp_proks)
head(temp_proks$OTU_ID)
head(temp_proks$taxonomy)
head(temp_proks[,1:5])

#Remove "taxonomy" column
grep("taxonomy", colnames(temp_proks))
temp_proks <- temp_proks[,-grep("taxonomy", colnames(temp_proks))]
dim(temp_proks) #78953 1208
head(colnames(temp_proks))

#Now check sums to make sure all the columns sum to 1
check_sums <- c()
for(i in c(2:ncol(temp_proks))){
  check_sums <- c(check_sums, sum(temp_proks[,i]))
}
which(check_sums!=1) #NONE! :) 

#Write out new file and delete temporary files
list.files("ProportionsFiles/")
file.remove("ProportionsFiles/output_temp_SPOT_16S_w.chloro_proportions.tsv")
list.files("ModifiedFiles/")
file.remove("ModifiedFiles/temp_SPOT_16S_w.chloro_counts_final.tsv")

write.table(x=temp_proks, file="ModifiedFiles/3.SPOT_16S_w.chloro_proportions.tsv", sep="\t", row.names=F, quote=F) 

#4. Relative abundance of 16S data without chloroplasts: 4.SPOT_16S_no_chloro_proportions.tsv-----

#4.1) write out a temporary file in the correct format for Jesse's script
#Add in "taxonomy"
dim(ordered_subs_proks_data_nochloro) #72284 1206
ordered_subs_proks_data_nochloro$taxonomy="SAR11"
dim(ordered_subs_proks_data_nochloro) #72284 1207
head(ordered_subs_proks_data_nochloro$taxonomy)

#Write this file out 
write.table(x=ordered_subs_proks_data_nochloro, file="ModifiedFiles/temp_SPOT_16S_no_chloro_counts.tsv", row.names=F, sep="\t", quote=F)

#Read it back in with no colnames
temp_proks <- read.table("ModifiedFiles/temp_SPOT_16S_no_chloro_counts.tsv", header=F, stringsAsFactors = F, row.names=NULL, sep="\t")
dim(temp_proks)

#Insert "#Constructed from biom file" and "#OTU ID"
colnames(temp_proks)[1]="#Constructed from biom file"
colnames(temp_proks)[2:ncol(temp_proks)]=""
temp_proks[1,1]="#OTU ID"

#Write this file out
write.table(x=temp_proks, file="ModifiedFiles/temp_SPOT_16S_no_chloro_counts_final.tsv", row.names=F, sep="\t", quote=F)

#Remove old file
file.remove("ModifiedFiles/temp_SPOT_16S_no_chloro_counts.tsv")

#4.2) SCP the temp file to Kraken and transform to proportions
#Don't forget to source activate qiime! 
#Download the output (with proportions) to a folder called ProportionsFiles/

#4.3) Read back in and remove the taxonomy column, then write out formatted file
#Make sure colnames are included
temp <- "ProportionsFiles/output_temp_SPOT_16S_no_chloro_proportions.tsv"
colnames <- scan(text=readLines(temp, 2), what="", quiet=T)
head(colnames)
length(colnames)
colnames <- colnames[3:(which(grepl("SPOT", colnames)==FALSE)[4]-1)]
length(colnames) #1206
tail(colnames)
colnames <- c("OTU_ID", colnames)
length(colnames) #1207 
length(unique(colnames)) #1207
head(colnames)

#NOW read in file
temp_proks <- read.table(temp, col.names=colnames, stringsAsFactors = F, sep="\t")
dim(temp_proks) #72284 1207 
head(temp_proks$OTU_ID)
head(temp_proks$taxonomy)
head(temp_proks[,1:5])

#Exclude taxonomy
temp_proks <- temp_proks[,-grep("taxonomy", colnames(temp_proks))]
dim(temp_proks) #72284 1206
grep("taxonomy", colnames(temp_proks)) #none, good! 

#Check sums
check_sums <- c()
for(i in c(2:ncol(temp_proks))){
  check_sums <- c(check_sums, sum(temp_proks[,i]))
}
length(check_sums)
which(check_sums!=1) #integer(0)

#Write out new file and delete temporary file
list.files("ProportionsFiles/")
file.remove("ProportionsFiles/output_temp_SPOT_16S_no_chloro_proportions.tsv")
list.files("ModifiedFiles/")
file.remove("ModifiedFiles/temp_SPOT_16S_no_chloro_counts_final.tsv")

write.table(x=temp_proks, file="ModifiedFiles/4.SPOT_16S_no_chloro_proportions.tsv", sep="\t", row.names=F, quote=F) 

#5. Raw counts of 18S data with metazoans: 5.SPOT_18S_w.metaz_counts.tsv-------

#5.1) Read in files
euks_data <- read.table("OriginalFiles/18S_mock_SPOT_Run40_Run58Lane2_dada2q2_feature-table.biom.tsv", header=T, stringsAsFactors = F, sep="\t")
dim(euks_data) #30537 1251
names(euks_data)

euks_tax <- read.table("OriginalFiles/18S_mock_SPOT_Run40_Run58Lane2_dada2q2_taxonomy.tsv", header=T, stringsAsFactors = F, sep="\t")
dim(euks_tax)
names(euks_tax)
head(euks_tax)

#5.2) Remove blanks, mocks, duplicate samples, etc. 
length(grep("Blank|Lab|515C|DOBD|mock|Mock|926R.b|926R.c|926R.d|926R.e|926R.f|ma|Ma|Rep2|Original|unknown|temp.notmerged",colnames(euks_data))) #120 
subsetted_euks_data <- euks_data[,-grep("Blank|Lab|515C|DOBD|mock|Mock|926R.b|926R.c|926R.d|926R.e|926R.f|ma|Ma|Rep2|Original|unknown|temp.notmerged",colnames(euks_data))]
dim(subsetted_euks_data) #30537 1131
ncol(euks_data)-ncol(subsetted_euks_data) #120

#5.3) Re-format colnames and figure out which samples are duplicates
#Re-format colnames thusly: word "SPOT" comes first, then year, month (but reformat month), date(but reformat date), depth, filter type
newcolnames <- c()
for(i in c(2:ncol(subsetted_euks_data))){
  a <- strsplit(colnames(subsetted_euks_data)[i], split=".", fixed=T)
  b <- as.numeric(a[[1]][3])
  d <- as.numeric(a[[1]][4])
  if(b<10){
    tempmonth <- paste("0", a[[1]][3], sep="")
  } else{
    tempmonth <- a[[1]][3]
  }
  if(d<10){
    tempdate <- paste("0", a[[1]][4], sep="")
  } else{
    tempdate <- a[[1]][4]
  }
  e <- paste(a[[1]][6],a[[1]][5],tempmonth,tempdate,a[[1]][7],a[[1]][2],sep=".")
  newcolnames <- c(newcolnames, e)
}
length(newcolnames) #1130
length(unique(newcolnames)) #1119 - so there are 11 duplicates
head(newcolnames)
tail(newcolnames) #Wow the year is really out of order here! 2008 samples come last, 2013 first

#Add back in "OTU.ID" so that newcolnames exactly match original colnames of subsetted data
newcolnames <- c("OTU_ID",newcolnames)
length(newcolnames) #1131
length(unique(newcolnames)) #1119 #Or 1120?
head(newcolnames)
head(colnames(subsetted_euks_data))

#Rename samples with "ma"
grep("[a-l,n-z]", newcolnames)
replace <- grep("[a-l,n-z]", newcolnames)
newcolnames[grep("[a-l,n-z]", newcolnames)]

newcolnames <- gsub(pattern="mb", replacement="m", newcolnames, fixed=T)
newcolnames <- gsub(pattern="Mb", replacement="M", newcolnames, fixed=T)
newcolnames[replace] #perfect!

#Find duplicates
which(duplicated(newcolnames)==TRUE)
length(which(duplicated(newcolnames)==TRUE))
newcolnames[which(duplicated(newcolnames)==TRUE)]
du=colnames(subsetted_euks_data)[which(duplicated(newcolnames)==TRUE)]
#Again, this is giving me the second duplicate, but I want the first
badcols <- c()
for(i in du){
  a <- strsplit(i, split=".", fixed=T)
  b <- paste(a[[1]][2],a[[1]][3],a[[1]][4],a[[1]][5],a[[1]][6],a[[1]][7],a[[1]][8],a[[1]][9],a[[1]][10],a[[1]][11], sep=".") #Paste together everything except run name
  #grep to find first duplicate
  d <- grep(b, colnames(subsetted_euks_data))
  badcols <- c(badcols, d)
}
#6.18.2009 5m AE sample still has primers 926R.926R, but is not duplicated in euks, so not to worry
#Most of the duplicates are samples from June, July, September of 2009 run on Run 53 and Run 55 - keep run 55
#The first three duplicates don't actually seem to be duplicates - double check this
colnames(subsetted_euks_data)[282]
grep("AE.5.24.2011.SPOT.890m", colnames(subsetted_euks_data))
colnames(subsetted_euks_data)[grep("AE.5.24.2011.SPOT.890m", colnames(subsetted_euks_data))]
##This is because they are both "Modified" samples
#282 is a duplicate of col #4 ("Run40.AE.5.24.2011.SPOT.890m.Monthlies.515Y.926R.Modifieda" vs. "Run47.AE.5.24.2011.SPOT.890m.Monthlies.515Y.926R.a")- get rid of col 4

colnames(subsetted_euks_data)[439]
grep("AE.10.17.2012.SPOT.890m", colnames(subsetted_euks_data))
#439 is a duplicate of col #3 ("Run40.AE.10.17.2012.SPOT.890m.Monthlies.515Y.926R.Modifieda" vs. "Run53.AE.10.17.2012.SPOT.890m.Monthlies.515Y.926R.a")- get rid of col 3
#Later run number trumps "modified" - probably the modified samples were run 2x, keep the later run

#Also running into trouble with "DCMa" vs "DCMb" which I gsubbed above
grep("AE.5.24.2011.SPOT.DCM", colnames(subsetted_euks_data)) #283 284
colnames(subsetted_euks_data)[grep("AE.5.24.2011.SPOT.DCM", colnames(subsetted_euks_data))]
#Get rid of 283 - that's DCMa

badcols[1]=3
badcols[2]=283
badcols[3]=4
badcols

#5.4a) Remove duplicated samples and keep old colnames 
subs_euks_data_oldcols <- subsetted_euks_data[,-badcols]
dim(subs_euks_data_oldcols) #30537, 1110
head(colnames(subs_euks_data_oldcols))
dim(subsetted_euks_data) #30537 1131
length(badcols) #21

#5.4b) Rename columns with newcolnames and then remove duplicated samples
subs_euks_data_newcols <- subsetted_euks_data
dim(subs_euks_data_newcols) #30537 1131
length(newcolnames)
colnames(subs_euks_data_newcols)=newcolnames
head(colnames(subs_euks_data_newcols))

#Test whether colnames actually match
test=as.data.frame(matrix(ncol=2, nrow=length(newcolnames)))
names(test)=c("oldcolnames", "newcolnames")
test$oldcolnames=colnames(subsetted_euks_data)
test$newcolnames=newcolnames
View(test)

#Randomly sample 10  colnames
for(i in (sample(2:ncol(subsetted_euks_data), size=10, replace=F))){
  print(colnames(subsetted_euks_data)[i])
  print(colnames(subs_euks_data_newcols)[i])
}

#Remove duplicates (now colnames won't match)
subs_euks_data_newcols=subs_euks_data_newcols[,-badcols]
dim(subs_euks_data_newcols) #30537 1110
ncol(subsetted_euks_data)-ncol(subs_euks_data_newcols)==length(badcols) #TRUE

#Double check for duplicates and lowercase letters
which(duplicated(colnames(subs_euks_data_newcols))==TRUE) #integer(0) :) 
length(which(duplicated(colnames(subs_euks_data_newcols))==TRUE))
length(grep("[a-l,n-z]", colnames(subs_euks_data_newcols))) #0

#5.5) Order by sampling date
order(colnames(subs_euks_data_newcols))==sort(colnames(subs_euks_data_newcols)) #FALSE?
#One spits out the index, the other spits out the names 
head(colnames(subs_euks_data_newcols)[order(colnames(subs_euks_data_newcols))])
head(colnames(subs_euks_data_newcols)[sort(colnames(subs_euks_data_newcols))])
colnames(subs_euks_data_newcols)[order(colnames(subs_euks_data_newcols))]==sort(colnames(subs_euks_data_newcols)) #THAT is true 

ordered_subs_euks_data <- subs_euks_data_newcols[,order(colnames(subs_euks_data_newcols))]
dim(ordered_subs_euks_data) #30537 1110
head(colnames(ordered_subs_euks_data))
tail(colnames(ordered_subs_euks_data))

#5.6) Double check Brad ASVs made it in, and in the correct order too
euks_tax$Taxon[grep("Braarudosphaera_bigelowii", euks_tax$Taxon)] #Why are there only 4? 
euks_tax$Taxon[grep("Braarudosphaera", euks_tax$Taxon)] #7 #Okay so there are 3 ASVs in Braarodsphaeara that are not B. bigelowii (Braarudosphaeraceae_X_sp.)

#Check that row numbers are the same
for(i in Brad_ASVs){
  print(grep(i, euks_data$OTU.ID))
}

Brad_ix <- c()
for(i in Brad_ASVs){
  Brad_ix <- c(Brad_ix, (grep(i, ordered_subs_euks_data_newcols$OTU_ID)))
}
Brad_ix
#They match!

#Check number of times they appear
for(i in Brad_ix){
  print(rowSums(euks_data[i,]>0))
}

for(i in Brad_ix){
  print(rowSums(ordered_subs_euks_data[i,]>0))
}
#Lost #1 on one date, #4 on two dates, #7 only once

#Check total sum
for(i in Brad_ix){
  print(rowSums(euks_data[i,2:ncol(euks_data)]))
}

for(i in Brad_ix){
  print(rowSums(ordered_subs_euks_data[i,2:ncol(ordered_subs_euks_data)]))
}
#Lost #1 17 times, #4 86 times, #7 106 times  #Wow, super abundant on the dates we lost!

#5.7) write out file!
write.table(x=ordered_subs_euks_data, file="ModifiedFiles/5.SPOT_18S_w.metaz_counts.tsv", sep="\t", row.names=F, quote=F)

#6. Raw counts of 18S data without metazoans: 6.SPOT_18S_no_metaz_counts.tsv------

#6.1) Figure out which rows correspond to metazoans
euks_tax$Taxon[grep("kingdom_Eukaryota;supergroup_Opisthokonta;division_Metazoa", euks_tax$Taxon)]
length(euks_tax$Taxon[grep("kingdom_Eukaryota;supergroup_Opisthokonta;division_Metazoa", euks_tax$Taxon)]) #2135
metaz_ASVs <- euks_tax$Feature.ID[grep("kingdom_Eukaryota;supergroup_Opisthokonta;division_Metazoa", euks_tax$Taxon)]
length(metaz_ASVs)
length(unique(metaz_ASVs))

#6.2) Remove 'em!
exclude_rows <- c()
for(i in metaz_ASVs){
  exclude_rows <- c(exclude_rows, grep(i, ordered_subs_euks_data$OTU_ID))
}
length(exclude_rows)
length(unique(exclude_rows)) #2135 

ordered_subs_euks_data_nometaz <- ordered_subs_euks_data[-exclude_rows,]
dim(ordered_subs_euks_data_nometaz) #24802 1110
dim(ordered_subs_euks_data) #30537 1110

nrow(ordered_subs_euks_data)-nrow(ordered_subs_euks_data_nometaz)==length(exclude_rows)

#6.3) double check that none of the columns sum to zero before writing out
check_sums <- c()
for(i in c(2:ncol(ordered_subs_euks_data_nometaz))){
  check_sums <- c(check_sums, sum(ordered_subs_euks_data_nometaz[,i]))
}
length(check_sums) #1109
which(check_sums<1)
check_sums[which(check_sums<1)]
length(which(check_sums<1)) #12 #uh oh... 12 columns have no data without metazoans :(

#Add 1 to each of these numbers to get the colnumber (have to take "OTU_ID" into account)
colnames(ordered_subs_euks_data_nometaz)[which(check_sums<1)+1] #Yep, these are all Durapore samples, and all are from the depths (500m and 890m)

#Exclude these columns
filt_ord_subs_euks_data_nometaz <- ordered_subs_euks_data_nometaz[,-(which(check_sums<1)+1)]
dim(filt_ord_subs_euks_data_nometaz)
dim(ordered_subs_euks_data_nometaz)

#Then write out the file
write.table(x=filt_ord_subs_euks_data_nometaz, file="ModifiedFiles/6.SPOT_18S_no_metaz_counts.tsv", sep="\t", row.names=F, quote=F)

#6.4) Write out 18S data for normalization (i.e. subs_euks_data_oldcols) but with no metazoans
dim(subs_euks_data_oldcols) #30537  1110
#Remove metazoans
subs_euks_data_oldcols_nometaz <- subs_euks_data_oldcols[-exclude_rows,]
dim(subs_euks_data_oldcols_nometaz) #28402 1110 
dim(subs_euks_data_oldcols) #30537 1110
nrow(subs_euks_data_oldcols)-nrow(subs_euks_data_oldcols_nometaz)==length(exclude_rows) #TRUE

write.table(x=subs_euks_data_oldcols_nometaz, file="ModifiedFiles/temp_normalize_18S_nometaz_counts.tsv", sep="\t", row.names=F)

#7. Relative abundance of 18S data with metazoans: 7.SPOT_18S_w.metaz_proportions.tsv-----
#7.1) Write out a temporary file in the correct format for Jesse's Python script 

#Figure out if any dates have zero abundances and remove them
dim(ordered_subs_euks_data)
head(colnames(ordered_subs_euks_data))
tail(colnames(ordered_subs_euks_data))

check_sums <- c()
for(i in c(2:ncol(ordered_subs_euks_data))){
  check_sums <- c(check_sums, sum(ordered_subs_euks_data[,i]))
}
which(check_sums<1) #integer(0) - ok to proceed

#Add in fake taxonomy 
ordered_subs_euks_data$taxonomy="supergroup_SAR"
dim(ordered_subs_euks_data) #30537 1111

#Write this temporary file out to "ModifiedFiles"
write.table(x=ordered_subs_euks_data, file="ModifiedFiles/temp_SPOT_18S_w.metaz_counts.tsv", row.names=F, sep="\t", quote=F)

#Read back in without colnames
temp_euks <- read.table("ModifiedFiles/temp_SPOT_18S_w.metaz_counts.tsv", header=F, stringsAsFactors = F, row.names=NULL, sep="\t")
dim(temp_euks)

#Insert the phrases "#Constructed from biom file" and "#OTU ID" 
colnames(temp_euks)[1]="#Constructed from biom file"
colnames(temp_euks)[2:ncol(temp_euks)]=""
temp_euks[1,1]="#OTU ID"

#Write this temporary file out to "ModifiedFiles"
write.table(x=temp_euks, file="ModifiedFiles/temp_SPOT_18S_w.metaz_counts_final.tsv", row.names=F, sep="\t", quote=F) 

#Remove old temporary file
file.remove("ModifiedFiles/temp_SPOT_18S_w.metaz_counts.tsv")

#7.2) SCP the temp file to Kraken and transform to proportions
#Don't forget to source activate qiime! 
#Download the output (with proportions) to a folder called ProportionsFiles/

#7.3) Read back in and remove the taxonomy column, then write out formatted file
#Make sure colnames are included
temp <- "ProportionsFiles/output_temp_SPOT_18S_w.metaz_proportions.tsv"
colnames <- scan(text=readLines(temp, 2), what="", quiet=T)
head(colnames)
length(colnames)

which(grepl("SPOT", colnames)==FALSE)
colnames <- c("OTU_ID", colnames[grep("taxonomy", colnames):(which(grepl("SPOT", colnames)==FALSE)[4]-1)])
head(colnames)
tail(colnames)
length(colnames)

#NOW read in data
temp_euks <- read.table(temp, col.names=colnames, stringsAsFactors = F, sep="\t")
dim(temp_euks) #30537 1111
head(temp_euks$taxonomy)
head(temp_euks$OTU_ID)

#Exclude taxonomy
temp_euks <- temp_euks[,-grep("taxonomy", colnames(temp_euks))]
dim(temp_euks) #30537 1110

#Check sums
check_sums <- c()
for(i in 2:ncol(temp_euks)){
  check_sums <- c(check_sums, sum(temp_euks[,i]))
}
which(check_sums!=1)
check_sums[which(check_sums!=1)] #nvm, these columns do sum to 1 

#Write out new file, delete temporary files
write.table(x=temp_euks, file="ModifiedFiles/7.SPOT_18S_w.metaz_proportions.tsv", sep="\t", row.names=F, quote=F)

list.files("ProportionsFiles/")
file.remove("ProportionsFiles/output_temp_SPOT_18S_w.metaz_proportions.tsv")
list.files("ModifiedFiles/")
file.remove("ModifiedFiles/temp_SPOT_18S_w.metaz_counts_final.tsv")

#8. Relative abundance of 18S data with no metazoans: 8.SPOT_18S_no_metaz_proportions.tsv------
#8.1) Write out a temporary file formatted so that Jesse's script can recognize it 
dim(filt_ord_subs_euks_data_nometaz) #28402 1098

#Add in fake taxonomy
head(colnames(filt_ord_subs_euks_data_nometaz))
tail(colnames(filt_ord_subs_euks_data_nometaz))
filt_ord_subs_euks_data_nometaz$taxonomy="supergroup_SAR"
dim(filt_ord_subs_euks_data_nometaz) #28402 1099
head(filt_ord_subs_euks_data_nometaz$taxonomy)

#Write this temporary file out to "ModifiedFiles"
write.table(x=filt_ord_subs_euks_data_nometaz, file="ModifiedFiles/temp_SPOT_18S_no_metaz_counts.tsv", row.names=F, sep="\t", quote=F)

#Read back in without colnames
temp_euks <- read.table("ModifiedFiles/temp_SPOT_18S_no_metaz_counts.tsv", header=F, stringsAsFactors = F, row.names=NULL, sep="\t")
dim(temp_euks) #28403 1099

#Insert the phrases "#Constructed from biom file" and "#OTU ID" 
colnames(temp_euks)[1]="#Constructed from biom file"
colnames(temp_euks)[2:ncol(temp_euks)]=""
temp_euks[1,1]="#OTU ID"

#Write this temporary file out to "ModifiedFiles" 
write.table(x=temp_euks, file="ModifiedFiles/temp_SPOT_18S_no_metaz_counts_final.tsv", row.names=F, sep="\t", quote=F)

#Remove old temp file
file.remove("ModifiedFiles/temp_SPOT_18S_no_metaz_counts.tsv")

#8.2) SCP the temp file to Kraken and transform to proportions
#Don't forget to source activate qiime! 
#Download the output (with proportions) to a folder called ProportionsFiles/

#8.3) Read back in and remove the taxonomy column, then write out formatted file
#Make sure colnames are included
temp <- "ProportionsFiles/output_temp_SPOT_18S_no_metaz_proportions.tsv"
colnames <- scan(text=readLines(temp, 2), what="", quiet=T)
head(colnames)
length(colnames)

which(grepl("SPOT", colnames)==FALSE)[4]
colnames <- colnames[1:(which(grepl("SPOT", colnames)==FALSE)[4]-1)]
tail(colnames)
head(colnames)
colnames <- c("OTU_ID", colnames[grep("taxonomy", colnames):length(colnames)])
head(colnames)
length(colnames) #1099
length(unique(colnames)) #1099

#NOW read in file
temp_euks <- read.table(temp, col.names=colnames, stringsAsFactors = F, sep="\t")
dim(temp_euks) #20842 1099
head(temp_euks$OTU_ID)
head(temp_euks$taxonomy)
head(colnames(temp_euks))
tail(colnames(temp_euks))

#Exclude taxonomy
temp_euks <- temp_euks[,-grep("taxonomy", colnames(temp_euks))]
dim(temp_euks) #28402 1098
grep("taxonomy", colnames(temp_euks)) #integer(0)

#Check all the cols sum to 1
check_sums <- c()
for(i in c(2:ncol(temp_euks))){
  check_sums <- c(check_sums, sum(temp_euks[,i]))
}
which(check_sums!=1)
check_sums[which(check_sums!=1)] #LIES! all the sums are 1! :) 

#Write out new file, delete temporary files
write.table(x=temp_euks, file="ModifiedFiles/8.SPOT_18S_no_metaz_proportions.tsv", sep="\t", row.names=F, quote=F)

list.files("ProportionsFiles/")
file.remove("ProportionsFiles/output_temp_SPOT_18S_no_metaz_proportions.tsv")
list.files("ModifiedFiles/")
file.remove("ModifiedFiles//temp_SPOT_18S_no_metaz_counts_final.tsv")
