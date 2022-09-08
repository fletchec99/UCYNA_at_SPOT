#Identify hashes of 18S sequences with highest relative abundance on date that UCYNA ASVs also peaked in abundance
#Date: 08.20.2020

rm(list=ls())
getwd()
list.dirs()

#Read in data----
#Read in proks data: 
proks_data <- read.table("ModifiedFiles/3.SPOT_16S_w.chloro_proportions.tsv", header=T, stringsAsFactors = F, sep="\t")
dim(proks_data) #78953 1208
head(names(proks_data))
length(grep("5m.AE", colnames(proks_data))) #125 
length(grep("DCM.AE", colnames(proks_data))) #118

#Read in euks data: 
euks_data <- read.table("ModifiedFiles/8.SPOT_18S_no_metaz_proportions.tsv", header=T, stringsAsFactors = F, sep="\t")
dim(euks_data) #28402 1098
length(grep("5m.AE", colnames(euks_data))) #124
length(grep("DCM.AE", colnames(euks_data))) #116

#Get hashes for ASV1 at 5m and DCM----
#5m
include_hashes_ASV1_5m <- c()
for(i in c("SPOT.2008.07.09.5m.AE", "SPOT.2015.08.05.5m.AE", "SPOT.2016.09.14.5m.AE", "SPOT.2016.10.05.5m.AE")){
  print(i)
  for(j in c(1:100)){
    include_hashes_ASV1_5m <- c(include_hashes_ASV1_5m, euks_data$OTU_ID[grep(sort(euks_data[,grep(i, colnames(euks_data))], decreasing=T)[j], euks_data[,grep(i, colnames(euks_data))])])
  }
}
length(include_hashes_ASV1_5m) #919
length(unique(include_hashes_ASV1_5m)) #283
include_hashes_ASV1_5m <- unique(include_hashes_ASV1_5m)
length(include_hashes_ASV1_5m)

#Write it out: 
write.table(x=include_hashes_ASV1_5m, file="abundant_18S_hashes_ASV1_5m_08.20.2020.txt", row.names=F, quote=F)

#DCM
include_hashes_ASV1_DCM <- c()
for(i in c("SPOT.2012.02.16.DCM.AE", "SPOT.2013.01.14.DCM.AE", "SPOT.2016.06.15.DCM.AE", "SPOT.2016.11.30.DCM.AE")){
  print(i)
  for(j in c(1:100)){
    include_hashes_ASV1_DCM <- c(include_hashes_ASV1_5m, euks_data$OTU_ID[grep(sort(euks_data[,grep(i, colnames(euks_data))], decreasing=T)[j], euks_data[,grep(i, colnames(euks_data))])])
  }
}
length(include_hashes_ASV1_DCM) #284
length(unique(include_hashes_ASV1_DCM)) #283
include_hashes_ASV1_DCM <- unique(include_hashes_ASV1_DCM)

#Write it out
write.table(x=include_hashes_ASV1_DCM, file="abundant_18S_hashes_ASV1_DCM_08.20.2020.txt", row.names=F, quote=F)

#Get hashes for ASV5 at 5m and DCM------
#5m
include_hashes_ASV5_5m <- c()
for(i in c("SPOT.2005.09.15.5m.AE", "SPOT.2008.06.18.5m.AE", "SPOT.2014.05.21.5m.AE", "SPOT.2014.08.13.5m.AE", "SPOT.2015.07.14.5m.AE", "SPOT.2016.02.17.5m.AE")){
  print(i)
  L <- length(which(as.numeric(euks_data[,grep(i, colnames(euks_data))])>0))
  print(L)
  if(L<100){
    for(j in c(1:L)){
      include_hashes_ASV5_5m <- c(include_hashes_ASV5_5m, euks_data$OTU_ID[grep(sort(euks_data[,grep(i, colnames(euks_data))], decreasing=T)[j], euks_data[,grep(i, colnames(euks_data))])])
    }
  } else{
    for(k in c(1:100)){
      include_hashes_ASV5_5m <- c(include_hashes_ASV5_5m, euks_data$OTU_ID[grep(sort(euks_data[,grep(i, colnames(euks_data))], decreasing=T)[k], euks_data[,grep(i, colnames(euks_data))])]) 
    }
  }
}
length(include_hashes_ASV5_5m) #1091
length(unique(include_hashes_ASV5_5m)) #312

include_hashes_ASV5_5m <- unique(include_hashes_ASV5_5m)

#Write out the list, this time to the folder "eLSA input"
write.table(x=include_hashes_ASV5_5m, file="eLSA_input/abundant_18S_hashes_ASV5_5m_08202020.txt", quote=F, row.names=F)

#DCM
include_hashes_ASV5_DCM <- c()
for(i in c("SPOT.2005.09.15.DCM.AE", "SPOT.2009.02.18.DCM.AE", "SPOT.2009.07.09.DCM.AE", "SPOT.2013.01.14.DCM.AE", "SPOT.2013.01.16.DCM.AE")){
  print(i)
  L <- length(which(as.numeric(euks_data[,grep(i, colnames(euks_data))])>0))
  print(L)
  if(L<100){
    for(j in c(1:L)){
      include_hashes_ASV5_DCM <- c(include_hashes_ASV5_DCM, euks_data$OTU_ID[grep(sort(euks_data[,grep(i, colnames(euks_data))], decreasing=T)[j], euks_data[,grep(i, colnames(euks_data))])])
    }
  } else{
    for(k in c(1:100)){
      include_hashes_ASV5_DCM <- c(include_hashes_ASV5_DCM, euks_data$OTU_ID[grep(sort(euks_data[,grep(i, colnames(euks_data))], decreasing=T)[k], euks_data[,grep(i, colnames(euks_data))])]) 
    }
  }
}
length(include_hashes_ASV5_DCM) #934
length(unique(include_hashes_ASV5_DCM)) #351 #Much better!

include_hashes_ASV5_DCM <- unique(include_hashes_ASV5_DCM)

#Write out to eLSA input
write.table(x=include_hashes_ASV5_DCM, file="eLSA_input/abundant_18S_hashes_ASV5_DCM.txt", quote=F, row.names=F)

#Get hashes for ASV6 at 5m and DCM----- 
#5m 
include_hashes_ASV6_5m <- c()
for(i in c("SPOT.2009.03.11.5m.AE", "SPOT.2013.09.18.5m.AE", "SPOT.2013.11.13.5m.AE", "SPOT.2017.11.15.5m.AE")){
  print(i)
  L <- length(which(as.numeric(euks_data[,grep(i, colnames(euks_data))])>0))
  print(L)
  if(L<100){
    for(j in c(1:L)){
      include_hashes_ASV6_5m <- c(include_hashes_ASV6_5m, euks_data$OTU_ID[grep(sort(euks_data[,grep(i, colnames(euks_data))], decreasing=T)[j], euks_data[,grep(i, colnames(euks_data))])])
    }
  } else{
    for(k in c(1:100)){
      include_hashes_ASV6_5m <- c(include_hashes_ASV6_5m, euks_data$OTU_ID[grep(sort(euks_data[,grep(i, colnames(euks_data))], decreasing=T)[k], euks_data[,grep(i, colnames(euks_data))])]) 
    }
  }
}
length(include_hashes_ASV6_5m) #453
length(unique(include_hashes_ASV6_5m)) #153

include_hashes_ASV6_5m <- unique(include_hashes_ASV6_5m)
#Write out file
write.table(x=include_hashes_ASV6_5m, file="eLSA_input/abundant_18S_hashes_ASV6_5m_08.20.2020.txt", quote=F, row.names = F)

#DCM
include_hashes_ASV6_DCM <- c()
for(i in c("SPOT.2005.03.16.DCM.AE", "SPOT.2009.08.19.DCM.AE", "SPOT.2012.06.13.DCM.AE", "SPOT.2013.12.23.DCM.AE", "SPOT.2017.10.11.DCM.AE")){
  print(i)
  L <- length(which(as.numeric(euks_data[,grep(i, colnames(euks_data))])>0))
  print(L)
  if(L<100){
    for(j in c(1:L)){
      include_hashes_ASV6_DCM <- c(include_hashes_ASV6_DCM, euks_data$OTU_ID[grep(sort(euks_data[,grep(i, colnames(euks_data))], decreasing=T)[j], euks_data[,grep(i, colnames(euks_data))])])
    }
  } else{
    for(k in c(1:100)){
      include_hashes_ASV6_DCM <- c(include_hashes_ASV6_DCM, euks_data$OTU_ID[grep(sort(euks_data[,grep(i, colnames(euks_data))], decreasing=T)[k], euks_data[,grep(i, colnames(euks_data))])]) 
    }
  }
}
length(include_hashes_ASV6_DCM) #799
length(unique(include_hashes_ASV6_DCM)) #294

include_hashes_ASV6_DCM <- unique(include_hashes_ASV6_DCM)

write.table(x=include_hashes_ASV6_5m, file="eLSA_input/abundant_18S_hashes_ASV6_DCM_08.20.2020.txt", quote=F, row.names = F)

#Combine across ASVs and write out new list of hashes----
length(c(include_hashes_ASV1_5m, include_hashes_ASV5_5m, include_hashes_ASV6_5m)) #748
length(include_hashes_ASV1_5m) #283
length(include_hashes_ASV5_5m) #312
length(include_hashes_ASV6_5m) #153
length(include_hashes_ASV1_5m) + length(include_hashes_ASV5_5m) + length(include_hashes_ASV6_5m) #748
length(unique(c(include_hashes_ASV1_5m, include_hashes_ASV5_5m, include_hashes_ASV6_5m))) #543

include_hashes_5m_all <- unique(c(include_hashes_ASV1_5m, include_hashes_ASV5_5m, include_hashes_ASV6_5m))
length(include_hashes_5m_all) #543

write.table(x=include_hashes_5m_all, file="eLSA_input/abundant_18S_hashes_5m_all_UCYNA_ASVs_08202020.txt", quote=F, row.names = F)

length(c(include_hashes_ASV1_DCM, include_hashes_ASV5_DCM, include_hashes_ASV6_DCM)) #928
length(unique(c(include_hashes_ASV1_DCM, include_hashes_ASV5_DCM, include_hashes_ASV6_DCM))) #697

include_hashes_DCM_all <- unique(c(include_hashes_ASV1_DCM, include_hashes_ASV5_DCM, include_hashes_ASV6_DCM))
length(include_hashes_DCM_all) #697

write.table(x=include_hashes_DCM_all, file="eLSA_input/abundant_18S_hashes_DCM_all_UCYNA_ASVs_08202020.txt", quote=F, row.names = F)

