#Purpose: normalize raw counts of 16S and 18S ASVs in SPOT Data by doubling counts of euks and dividing counts of proks and euks by the % passing DADA2
#Then need to upload to Kraken and convert to proportions
#Output: This returns a file with ASV relative abundances out of (16S + 18S).
#Note: we are assuming a 2-fold bias against 18S sequences, which has been found with Illumina HiSeq or MiSeq data (Yeh et al. 2018)
#Ver 12.1.20

rm(list=ls())
getwd()
setwd("normalize_16S_18S/") #Set wd to normalizing folder 

#1. Calculate percent of reads that passed DADA2 denoising for both proks and euks -----
#Read in the data in a for-loop
directory <- getwd()
filenames<-system(paste('ls ', directory, '/DADA2_stats', sep=''), intern=TRUE)
filenames
#We want only files 1:19, not 20, which is the copy/ paste file

#Read in files - PROKS and EUKS
for(i in c(1:19)){
  a <- strsplit(filenames[i], split="_", fixed=T) 
  print(paste(a[[1]][1], a[[1]][2], sep="_"))
  assign(paste(a[[1]][1], a[[1]][2], sep="_"), value=read.table(file=paste('DADA2_stats/', filenames[i], sep=""), sep="\t", stringsAsFactors = F, header=T))
}

#Rbind the proks and euks, then calculate the percent passing DADA2 for each 
proks_stats=rbind(PROKS_run031, PROKS_run040, PROKS_run046, PROKS_run047, PROKS_run053, PROKS_run055, PROKS_run056Lane1, PROKS_run056Lane2, PROKS_run058Lane1, PROKS_run058Lane2)
dim(proks_stats) #1451 6
colnames(proks_stats)
proks_stats$perc.passed=proks_stats$non.chimeric/proks_stats$input
dim(proks_stats) #1451 7 

euks_stats <- rbind(EUKS_run040, EUKS_run046, EUKS_run047, EUKS_run053, EUKS_run055, EUKS_run056Lane1, EUKS_run056Lane2, EUKS_run058Lane1, EUKS_run058Lane2)
dim(euks_stats) #1283 5 
colnames(euks_stats) 
euks_stats$perc.passed=euks_stats$non.chimeric/euks_stats$input
colnames(euks_stats)
dim(euks_stats) #1285 6


#2. Normalize ASV counts (divide counts of ASVs/ percent passed for each sample, multiply euks ASV counts by the bias you specified)------

#a. Proks--------
#Read in ASV counts data
proks_data=read.table(file="temp_normalize_16S_nochloro_absolute.tsv", sep="\t", stringsAsFactors = F, header=T, row.names=1) #Make the OTU IDs be the row names
dim(proks_data) #722284 1205
head(colnames(proks_data))

#Set up a new matrix for normalized data
proks_norm <- as.data.frame(matrix(nrow=nrow(proks_data), ncol=ncol(proks_data)))
row.names(proks_norm)=row.names(proks_data)
colnames(proks_norm)=colnames(proks_data)

ncol(proks_data) #1205
nrow(proks_stats) #1451 
#There are >200 samples in proks_stats that don't match rows in proks_data

#Divide ASV count for each sample by percent passing, write into the new matrix
missing_samples=c()
for(i in proks_stats$sample.id){
  #print(i)
  if(length(grep(i, colnames(proks_data)))<1){
    missing_samples=c(missing_samples, i)
  } else{
  #print(grep(i, colnames(proks_data)))
  proks_norm[,grep(i, colnames(proks_norm))]=proks_data[,grep(i, colnames(proks_data))]/proks_stats$perc.passed[grep(i, proks_stats$sample.id)]
  }
}
#Expecting to have 246 missing samples, mostly from blanks and replicates (b-f, etc.)
nrow(proks_stats)-ncol(proks_data) #246
length(missing_samples)
length(missing_samples)==nrow(proks_stats)-ncol(proks_data)
missing_samples
length(grep("Blank", missing_samples)) #80
length(grep("Mock", missing_samples, ignore.case = T)) #115 - did any of the mocks get included? 
grep("Mock", colnames(proks_data)) #Nope!! 
missing_samples[-grep("Mock|Blank", missing_samples, ignore.case = T)] #51 "samples"
#Here's what should be missing: "Blank|Lab|515C|DOBD|mock|Mock|926R.b|926R.c|926R.d|926R.e|926R.f|ma|Ma|Rep2|Original|unknown|temp.notmerged"
length(grep("Blank|Lab|515C|DOBD|mock|Mock|926R.b|926R.c|926R.d|926R.e|926R.f|ma|Ma|Rep2|Original|unknown|temp.notmerged", missing_samples)) #228
#And that's what's in "missing samples" - looks good! 
missing_samples[-grep("Blank|Lab|515C|DOBD|mock|Mock|926R.b|926R.c|926R.d|926R.e|926R.f|ma|Ma|Rep2|Original|unknown|temp.notmerged", missing_samples)]
#These are probably the duplicates 

#b. Now repeat normalization for the euks, and also multiply by 2x the bias-----
#Read in files (use file with PR2 taxonomy)
euks_data=read.table("temp_normalize_18S_nometaz_absolute.tsv", sep="\t", stringsAsFactors = F, header=T, row.names=1) #OTU ID is the row names
dim(euks_data) #28402 1109 
head(colnames(euks_data))
grep("taxonomy", colnames(euks_data))

#Set up a new matrix for normalized data
euks_norm <- as.data.frame(matrix(nrow=nrow(euks_data), ncol=ncol(euks_data)))
row.names(euks_norm)=row.names(euks_data)
colnames(euks_norm)=colnames(euks_data)

#Divide ASV count by percent passing DADA2 for each sample and multiply by the given bias to normalize ASV counts
#Figure out what samples are missing
missing_samples_euks=c()
for(i in euks_stats$sample.id){
  if(length(grep(i, colnames(euks_data)))<1){
    missing_samples_euks=c(missing_samples_euks, i)
  }else{
  euks_norm[,grep(i, colnames(euks_norm))]=2*euks_data[,grep(i, colnames(euks_data))]/euks_stats$perc.passed[grep(i, euks_stats$sample.id)]    
  }
}
#Expect 174 missing samples 
length(missing_samples_euks)==nrow(euks_stats)-ncol(euks_data) #TRUE
length(grep("Blank|Mock", missing_samples_euks, ignore.case = T)) #131
missing_samples_euks[-grep("Blank|Mock", missing_samples_euks, ignore.case = T)]
missing_samples_euks[-grep("Blank|Lab|515C|DOBD|mock|Mock|926R.b|926R.c|926R.d|926R.e|926R.f|ma|Ma|Rep2|Original|unknown|temp.notmerged", missing_samples_euks)] #25 real samples left out

missing_samples[-grep("Blank|Lab|515C|DOBD|mock|Mock|926R.b|926R.c|926R.d|926R.e|926R.f|ma|Ma|Rep2|Original|unknown|temp.notmerged", missing_samples)]==missing_samples_euks[-grep("Blank|Lab|515C|DOBD|mock|Mock|926R.b|926R.c|926R.d|926R.e|926R.f|ma|Ma|Rep2|Original|unknown|temp.notmerged", missing_samples_euks)] #Some of which are also left out of proks data

#3. Combine proks and euks tables of normalized sequencing counts-------
#First, we need to format a little, because this next part only works if colnames for proks and euks data are the same.

#a. Add in dummy columns for samples missing from one matrix or the other (i.e. the mocks) ------

ncol(proks_norm)-ncol(euks_norm) #One has 96 more columns more than the other - this is probably 

#1. First, check if there are any samples missing from proks spreadsheet, and if so, write them in 
missing_from_proks <- c()
for(i in colnames(euks_norm)){
  if(i %in% colnames(proks_norm)){
    #print(i)
  } else{
    missing_from_proks <- c(missing_from_proks, i)
  }
}
missing_from_proks #5 samples -- really?? One of them is a durapore sample from DCM, which has euks in it but no proks 

#Add in dummy columns for these samples--Must create a new matrix to do so
norm_proks <- as.data.frame(matrix(nrow=nrow(proks_norm), ncol=(ncol(proks_norm)+length(missing_from_proks))))
rownames(norm_proks)=rownames(proks_norm)
colnames(norm_proks)=c(colnames(proks_norm), missing_from_proks)
norm_proks[,1:ncol(proks_norm)]=proks_norm
#Now add in dummy columns
for(i in c(1:length(missing_from_proks))){
  norm_proks[,i+ncol(proks_norm)]=0
}

#2. Then, check if there are columns missing in euks spreadsheet and if so, write them in
missing_from_euks <- c()
for(i in colnames(proks_norm)){
  if(i %in% colnames(euks_norm)){
  } else {
    missing_from_euks <- c(missing_from_euks, i)
  }
}
length(missing_from_euks) #101 columns actually 

#Add dummy columns for these samples into normalized euks data
norm_euks <- as.data.frame(matrix(nrow=nrow(euks_norm), ncol=(ncol(euks_norm)+length(missing_from_euks))))
rownames(norm_euks)=rownames(euks_norm)
colnames(norm_euks)=c(colnames(euks_norm), missing_from_euks)
norm_euks[,1:ncol(euks_norm)]=euks_norm
#Now add in dummy columns 
for(i in c(1:length(missing_from_euks))){
  norm_euks[,i+ncol(euks_norm)]=0
}

#b. Make sure normalized proks and euks data are in the same order, but OTU ID comes first------
print("3b. Re-ordering samples")
#Proks first
proks_ordered <- as.data.frame(matrix(nrow=nrow(norm_proks), ncol=ncol(norm_proks)+1))
colnames(proks_ordered)[1]="OTU_ID"
proks_ordered$OTU_ID=row.names(norm_proks)
proks_ordered[,2:ncol(proks_ordered)]=norm_proks[,order(colnames(norm_proks))]
colnames(proks_ordered)[2:ncol(proks_ordered)]=colnames(norm_proks)[order(colnames(norm_proks))]

#Now euks second
euks_ordered <- as.data.frame(matrix(nrow=nrow(norm_euks), ncol=ncol(norm_euks)+1))
colnames(euks_ordered)[1]="OTU_ID"
euks_ordered$OTU_ID=rownames(norm_euks)
euks_ordered[,2:ncol(euks_ordered)]=norm_euks[,order(colnames(norm_euks))]
colnames(euks_ordered)[2:ncol(euks_ordered)]=colnames(norm_euks)[order(colnames(norm_euks))]

#c. Merge normalized proks and euks data-----
print("3c. Merging proks and euks data")
norm_ordered_proks_euks <- rbind(proks_ordered, euks_ordered, stringsAsFactors=F)
dim(norm_ordered_proks_euks) #100686 1211

#4. Format so that Jesse's script can convert this monster file!! -----

#First write in a dummy taxonomy #Can't write in real taxonomy because we don't have that in the original files 
norm_ordered_proks_euks$taxonomy="SAR11"

#Add in "#OTU ID" #Not yet
colnames(norm_ordered_proks_euks)[1]="OTU_ID"

#Write out this temp file 
write.table(x=norm_ordered_proks_euks, file="temp_normalized_counts_proks_euks_12.7.2020.tsv", row.names=F, sep="\t", quote=F)

#Read this file back in without colnames -- will it work with the hash? #NO!
temp_data=read.table("temp_normalized_counts_proks_euks_12.7.2020.tsv", header=F, stringsAsFactors = F, sep="\t", row.names=NULL)
dim(temp_data) #100687 1212
head(colnames(temp_data))
temp_data[1,1:6] #Looking good! 

#Insert "#Constructed from biom file" and "#OTU ID"
colnames(temp_data)[1]="#Constructed from biom file"
colnames(temp_data)[2:ncol(temp_data)]=""
temp_data[1,1]="#OTU ID"

#Write out file
write.table("temp_formatted_normalized_counts_proks_euks_12.7.2020.tsv", x=temp_data, sep="\t", row.names=F, quote=F)

#5. SCP to Kraken and use Jesse's script to convert to proportions---

#6. Re-download file and format colnames properly-----
#Read in the data, using that trick to ignore #OTU ID, which is apparently the first colname
temp="normalize_16S_18S/normalized_proportions_proks_euks_12.7.2020.tsv" 
colnames <- scan(text=readLines(temp, 1), what="", quiet=T)
head(colnames, n=10) #OTU ID are the first 2, then "taxonomy"
tail(colnames)
length(colnames) #1213

colnames=c("OTU_ID", colnames[3:1213])
length(colnames) #1212 - more like it 

data=read.table(temp, col.names=colnames, stringsAsFactors = F, sep="\t")
head(colnames(data))
dim(data) #100686 1212 

#Get rid of fake taxonomy
data=data[,-grep("tax", colnames(data))]
dim(data) #100686 1211 

#Double check sums are all 1 (or close to it)
check_sums=c()
for(i in c(2:ncol(data))){
  check_sums=c(check_sums, sum(data[,i]))
}
length(check_sums) #1210
length(which(check_sums!=1)) #1125
length(which(round(check_sums)!=1)) #0
check_sums[which(check_sums!=1)] #It appears all of them are not quite 1, but close to 1
check_sums[which(check_sums>1)]

#2. Format data: shorten colnames, and remove duplicates (if necessary)
head(colnames(data))

#Re-format colnames thusly: word "SPOT" comes first, then year, month (but reformat month), date(but reformat date), depth, filter type
###Stolen straight outta the 01... script 
newcolnames <- c()
for(i in c(2:ncol(data))){
  a <- strsplit(colnames(data)[i], split=".", fixed=T)
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
length(newcolnames) #1210
head(newcolnames)
tail(newcolnames)

#Add in "OTU_ID" so that newcolnames and colnames have same indices 
newcolnames=c("OTU_ID", newcolnames)
length(newcolnames)

#Test whether oldcolnames match newcolnames 
for(i in (sample(x=c(1:length(newcolnames)), size=10, replace=F))){
  print(colnames(data)[i])
  print(newcolnames[i])
}

#Figure out which ones are duplicates
length(unique(newcolnames)) #1209 - only one duplicate?? Two? 
which(duplicated(newcolnames))
newcolnames[which(duplicated(newcolnames))] #"SPOT.2011.05.24.890m.AE" "SPOT.2012.10.17.890m.AE"
newcolnames[grep("2011.05.24.890m.AE", newcolnames)]
colnames[grep("AE.5.24.2011.SPOT.890m", colnames(data))] 
#Ok, so one of these samples is "Modified", from Run40, and the other is not, from Run47 
#Doesn't REALLY matter bc both sets of duplicates are from 890m
#Remove the one from run40, keep the one from Run47 #Remove 89 from newcolnames/ 90 from colnames
newcolnames[grep("SPOT.2012.10.17.890m.AE", newcolnames)]
colnames[grep("AE.10.17.2012.SPOT.890m", colnames(data))] 
#Run 40, modified, and Run53 <- keep Run53 #Remove 88 from newcolnames

#Set newcolnames
length(colnames(data))==length(newcolnames)
head(colnames(data))
head(newcolnames)
colnames(data)=newcolnames

#Remove duplicated columns
data=data[,-c(88,89)]

#Now order by date, not by Run#
colnames(data)[order(colnames(data))]
ordered_data=data[,order(colnames(data))]
head(colnames(ordered_data))
tail(colnames(ordered_data))
dim(ordered_data) #100686 1209

#And write it out into "ModifiedFiles!" 
write.table(x=ordered_data, file="ModifiedFiles/9.SPOT_16S_no_chloro_18S_no_metaz_normalized_proportions.tsv", sep="\t", row.names=F, quote=F)

