#Parse eLSA output
#ver. 12.18.2020

#Read in relevant files-----
setwd("/Users/cfletcherh1/Desktop/DataAnalysis/SPOT_16S_18S_eLSA/")

#Proks taxa
proks_tax <- read.table("ModifiedFiles/Proks_tax_classified_23072020_SILVA132.tsv", header=T, stringsAsFactors = F, sep="\t")
dim(proks_tax) #78953 3
head(proks_tax) 
names(proks_tax)

UCYNA_ix=proks_tax$Feature.ID[grep("UCYN-A", proks_tax$Taxon)]
UCYNA_ix

#Euks taxa
euks_tax <- read.table("ModifiedFiles/Euks_tax_classified_19052020_SILVA132.tsv", header=T, stringsAsFactors = F, sep="\t")
dim(euks_tax) #30537 3 
head(euks_tax)
names(euks_tax)

Brad_ix=euks_tax$Feature.ID[grep("Braarudosphaera", euks_tax$Taxon)]
Brad_ix

#1. CLR transformed data from 5m-----

#Read in the data in a for-loop
directory <- getwd()
filenames<-system(paste('ls ', directory, '/eLSA_output/piecewise_eLSA/CLR_5m', sep=''), intern=TRUE)

#Change i to be "5:10"
for(i in c(5:10)){
  a <- strsplit(filenames[i], split="_", fixed=T) 
  print(paste("eLSA_data", a[[1]][3], a[[1]][4], a[[1]][5], sep="_"))
  assign(x=paste("eLSA_data", a[[1]][3], a[[1]][4], a[[1]][5], sep="_"), value=read.table(file=paste("eLSA_output/piecewise_eLSA/CLR_5m/", filenames[i], sep=""), sep="\t", stringsAsFactors = F, header=T))
}

#How to remove data objects en masse: rm(list=ls()[grep("eLSA_data", ls())])

#Use rbind to combine them
eLSA_data_concat=rbind(eLSA_data_subs_A_B, eLSA_data_subs_A_C, eLSA_data_subs_A_CLR, eLSA_data_subs_B_C, eLSA_data_subs_B_CLR, eLSA_data_subs_C_CLR)
dim(eLSA_data_concat) #296065     27 
#WHOA

#Calculate q values from p values in the table-----
BiocManager::install("qvalue")
library(qvalue)

#Normal Q values: 
Q=qvalue(eLSA_data_concat$P)
length(Q$qvalues)==length(eLSA_data_concat$P) #TRUE
eLSA_data_concat$Q=Q$qvalues

#Pearson's: Qpcc and Qspcc 
pearsons_P=qvalue(eLSA_data_concat$Ppcc)
pearsons_Q=pearsons_P$qvalues
length(eLSA_data_concat$Ppcc)==length(pearsons_P$qvalues)
eLSA_data_concat$Qpcc=pearsons_Q

pearsons_Q_delay=qvalue(eLSA_data_concat$Pspcc)
length(pearsons_Q_delay$qvalues)==length(eLSA_data_concat$Qspcc)
eLSA_data_concat$Qspcc=pearsons_Q_delay$qvalues

#Then Spearman's: Qscc and Qsscc
Spearman_Q=qvalue(eLSA_data_concat$Pscc)
length(Spearman_Q$qvalues)==length(eLSA_data_concat$Pscc)
eLSA_data_concat$Qscc=Spearman_Q$qvalues

Spearman_Q_delay=qvalue(eLSA_data_concat$Psscc)
length(Spearman_Q_delay$qvalues)==length(eLSA_data_concat$Psscc)
eLSA_data_concat$Qsscc=Spearman_Q_delay$qvalues

dim(eLSA_data_concat)

#Shorten names of X and Y taxa 
eLSA_data_concat$X_abr=abbreviate(eLSA_data_concat$X, minlength=4)
head(eLSA_data_concat$X_abr)
length(unique(eLSA_data_concat$X)) #380
length(unique(eLSA_data_concat$X_abr)) #380

eLSA_data_concat$Y_abr=abbreviate(eLSA_data_concat$Y, minlength=4)
head(eLSA_data_concat$Y_abr)
length(unique(eLSA_data_concat$Y)) #380
length(unique(eLSA_data_concat$Y_abr)) #380
length(grep(eLSA_data_concat$Y[1], eLSA_data_concat$Y)) #189
length(grep(eLSA_data_concat$Y_abr[1], eLSA_data_concat$Y_abr)) #189

#Write out modified file!! 
#Not on 11.22.20, but: write.table(x=eLSA_data_concat, file="eLSA_output/piecewise_eLSA/CLR_5m/concat_eLSA_output_CLR_5m_11.20.2020.tsv", sep="\t", quote=F, row.names=F)

#Parse out only rows that are significant by Pearson's and Spearmans
length(which(eLSA_data_concat$Ppcc < 0.05 & eLSA_data_concat$Qpcc < 0.05 & eLSA_data_concat$Pscc < 0.05 & eLSA_data_concat$Qscc < 0.05)) #34108, or 11%
#Try more stringent criteria
length(which(eLSA_data_concat$P< 0.05 & eLSA_data_concat$Q<0.05 & eLSA_data_concat$Ppcc < 0.05 & eLSA_data_concat$Qpcc < 0.05 & eLSA_data_concat$Pscc < 0.05 & eLSA_data_concat$Qscc < 0.05)) #21811
length(which(eLSA_data_concat$P< 0.005 & eLSA_data_concat$Q<0.01 & eLSA_data_concat$Ppcc < 0.005 & eLSA_data_concat$Qpcc < 0.01 & eLSA_data_concat$Pscc < 0.005 & eLSA_data_concat$Qscc < 0.01)) #12313, or 4% of original rows

parsed_eLSA_CLR_5m=eLSA_data_concat[which(eLSA_data_concat$P< 0.005 & eLSA_data_concat$Q<0.01 & eLSA_data_concat$Ppcc < 0.005 & eLSA_data_concat$Qpcc < 0.01 & eLSA_data_concat$Pscc < 0.005 & eLSA_data_concat$Qscc < 0.01),]
dim(parsed_eLSA_CLR_5m) #12313 29

#Add in taxonomy: X
parsed_eLSA_CLR_5m$X_tax="SAR"
length(unique(parsed_eLSA_CLR_5m$X)) #304 unique taxa 

#Add in euks taxonomy, get indices for UCYNA
ix <- c()
for(i in c(1:nrow(parsed_eLSA_CLR_5m))){
  a <- parsed_eLSA_CLR_5m$X[i]
  b <- grep(a, euks_tax$Feature.ID)
  if(as.numeric(length(b)<1)){
    print(b)
    ix <- c(ix, i)
  }else{
    parsed_eLSA_CLR_5m$X_tax[i]=euks_tax$Taxon[b]
  }
}

#Add in UCYN-A
parsed_eLSA_CLR_5m$X_tax[grep("_AE", parsed_eLSA_CLR_5m$X)]="UCYN-A_1_AE"

#add in taxonomy for Y
parsed_eLSA_CLR_5m$Y_tax="SAR"
length(unique(parsed_eLSA_CLR_5m$Y)) #323

ix <- c()
for(i in c(1:nrow(parsed_eLSA_CLR_5m))){
  a <- parsed_eLSA_CLR_5m$Y[i]
  b <- grep(a, euks_tax$Feature.ID)
  if(as.numeric(length(b)<1)){
    print(i)
    ix <- c(ix, i)
  }else{
    parsed_eLSA_CLR_5m$Y_tax[i]=euks_tax$Taxon[b]
  }
}

#Add in UCYN-A
length(ix)==length(grep("_AE|_D", parsed_eLSA_CLR_5m$Y)) #TRUE #Both 179
head(parsed_eLSA_CLR_5m$Y[ix])

for(i in ix){
  a <- parsed_eLSA_CLR_5m$Y[i]
  b <- strsplit(a, split="_", fixed=T)
  d <- grep(b[[1]][1], UCYNA_ix)
  parsed_eLSA_CLR_5m$Y_tax[i]=(paste("UCYN-A", d, b[[1]][2], sep="_"))
}

grep("SAR", parsed_eLSA_CLR_5m$X_tax)
grep("SAR", parsed_eLSA_CLR_5m$Y_tax)
#Both integer(0) #YAY! 


#Abbreviated X taxonomy----- 
#Start by labelling everything "Eukaryote"
parsed_eLSA_CLR_5m$X_tax_abr="Eukaryote"

#Manually change taxa, oh god :( #At the "supergroup" level 

##Start with the Stramenopiles: All stramenopiles except Pseudo-nitzschia and Chaetoceros set to "Stramenopiles"

length(grep("kingdom_Eukaryota;supergroup_Stramenopiles", parsed_eLSA_CLR_5m$X_tax)) #That's like 10% of the X taxa
parsed_eLSA_CLR_5m$X_tax_abr[grep("kingdom_Eukaryota;supergroup_Stramenopiles;", parsed_eLSA_CLR_5m$X_tax)]="Stramenopile"
grep("Stramenopile", temp)
parsed_eLSA_CLR_5m$X_tax_abr[grep(temp[grep("Stramenopile", temp)], parsed_eLSA_CLR_5m$X_tax)]="Stramenopile"


length(grep("kingdom_Eukaryota;supergroup_Stramenopiles;division_Ochrophyta;class_Bacillariophyta;order_Bacillariophyta_X;family_Raphid-pennate;genus_Pseudo-nitzschia", parsed_eLSA_CLR_5m$X_tax))
parsed_eLSA_CLR_5m$X_tax_abr[grep("kingdom_Eukaryota;supergroup_Stramenopiles;division_Ochrophyta;class_Bacillariophyta;order_Bacillariophyta_X;family_Raphid-pennate;genus_Pseudo-nitzschia", parsed_eLSA_CLR_5m$X_tax)]="Pseudo-nitzchia"

length(grep("kingdom_Eukaryota;supergroup_Stramenopiles;division_Ochrophyta;class_Bacillariophyta;order_Bacillariophyta_X;family_Polar-centric-Mediophyceae;genus_Chaetoceros", parsed_eLSA_CLR_5m$X_tax))
parsed_eLSA_CLR_5m$X_tax_abr[grep("kingdom_Eukaryota;supergroup_Stramenopiles;division_Ochrophyta;class_Bacillariophyta;order_Bacillariophyta_X;family_Polar-centric-Mediophyceae;genus_Chaetoceros", parsed_eLSA_CLR_5m$X_tax)]="Chaetoceros"

length(grep("kingdom_Eukaryota;supergroup_Rhizaria", parsed_eLSA_CLR_5m$X_tax))
parsed_eLSA_CLR_5m$X_tax_abr[grep("Rhizaria", parsed_eLSA_CLR_5m$X_tax)]="Rhizaria"

length(grep("kingdom_Eukaryota;supergroup_Opisthokonta", parsed_eLSA_CLR_5m$X_tax))
parsed_eLSA_CLR_5m$X_tax_abr[grep("kingdom_Eukaryota;supergroup_Opisthokonta", parsed_eLSA_CLR_5m$X_tax)]="Opisthokonta"

##All hacrobia except prymnesiophytes, Haptophytes, and Phaeocystis set to "Hacrobia" 

length(grep("kingdom_Eukaryota;supergroup_Hacrobia", parsed_eLSA_CLR_5m$X_tax))
parsed_eLSA_CLR_5m$X_tax_abr[grep("kingdom_Eukaryota;supergroup_Hacrobia", parsed_eLSA_CLR_5m$X_tax)]="Hacrobia"

length(grep("kingdom_Eukaryota;supergroup_Hacrobia;division_Haptophyta;class_Prymnesiophyceae;order_Prymnesiales;family_Prymnesiaceae;genus_Prymnesium;species_Prymnesium_sp.", parsed_eLSA_CLR_5m$X_tax))
parsed_eLSA_CLR_5m$X_tax_abr[grep("kingdom_Eukaryota;supergroup_Hacrobia;division_Haptophyta;class_Prymnesiophyceae;order_Prymnesiales;family_Prymnesiaceae;genus_Prymnesium;species_Prymnesium_sp.", parsed_eLSA_CLR_5m$X_tax)]="Prymnesiophyte"

length(grep("kingdom_Eukaryota;supergroup_Hacrobia;division_Haptophyta;class_Prymnesiophyceae;order_Prymnesiales;family_Chrysochromulinaceae;genus_Chrysochromulina", parsed_eLSA_CLR_5m$X_tax)) #287
parsed_eLSA_CLR_5m$X_tax_abr[grep("kingdom_Eukaryota;supergroup_Hacrobia;division_Haptophyta;class_Prymnesiophyceae;order_Prymnesiales;family_Chrysochromulinaceae;genus_Chrysochromulina", parsed_eLSA_CLR_5m$X_tax)]="Chrysochromulina"

#Phaeocystis
length(grep("kingdom_Eukaryota;supergroup_Hacrobia;division_Haptophyta;class_Prymnesiophyceae;order_Phaeocystales;family_Phaeocystaceae;genus_Phaeocystis", parsed_eLSA_CLR_5m$X_tax))
parsed_eLSA_CLR_5m$X_tax_abr[grep("kingdom_Eukaryota;supergroup_Hacrobia;division_Haptophyta;class_Prymnesiophyceae;order_Phaeocystales;family_Phaeocystaceae;genus_Phaeocystis", parsed_eLSA_CLR_5m$X_tax)]="Phaeocystis"

#Haptophyte
length(grep("kingdom_Eukaryota;supergroup_Hacrobia;division_Haptophyta;class_Haptophyta", parsed_eLSA_CLR_5m$X_tax))
parsed_eLSA_CLR_5m$X_tax_abr[grep("kingdom_Eukaryota;supergroup_Hacrobia;division_Haptophyta;class_Haptophyta", parsed_eLSA_CLR_5m$X_tax)]="Haptophyte"

#Insert Brad ix's
length(grep("kingdom_Eukaryota;supergroup_Hacrobia;division_Haptophyta;class_Prymnesiophyceae;order_Prymnesiophyceae_X;family_Braarudosphaeraceae;genus_Braarudosphaeraceae_X", parsed_eLSA_CLR_5m$X_tax)) #169
parsed_eLSA_CLR_5m$X[grep("kingdom_Eukaryota;supergroup_Hacrobia;division_Haptophyta;class_Prymnesiophyceae;order_Prymnesiophyceae_X;family_Braarudosphaeraceae;genus_Braarudosphaeraceae_X", parsed_eLSA_CLR_5m$X_tax)]

for(i in Brad_ix){
  if(length(grep(i, parsed_eLSA_CLR_5m$X))>1){
    print(i)
    print(grep(i, parsed_eLSA_CLR_5m$X)) 
  }
} #Hmm ok, only 3 Brad ASVs showed up: #3 (formerly 4), 4 (formerly 5), 7x (formerly 7) 

parsed_eLSA_CLR_5m$X_tax_abr[grep(Brad_ix[4], parsed_eLSA_CLR_5m$X)]="Braarudospharea_bigelowii_3"
parsed_eLSA_CLR_5m$X_tax_abr[grep(Brad_ix[5], parsed_eLSA_CLR_5m$X)]="Braarudospharea_bigelowii_4"
parsed_eLSA_CLR_5m$X_tax_abr[grep(Brad_ix[7], parsed_eLSA_CLR_5m$X)]="Braarudosphaera_sp_7x"

#Moving on!
length(grep("kingdom_Eukaryota;supergroup_Archaeplastida;division_Rhodophyta", parsed_eLSA_CLR_5m$X_tax))
parsed_eLSA_CLR_5m$X_tax_abr[grep("kingdom_Eukaryota;supergroup_Archaeplastida;division_Rhodophyta", parsed_eLSA_CLR_5m$X_tax)]="Archeplastida"

length(grep("kingdom_Eukaryota;supergroup_Archaeplastida;division_Chlorophyta", parsed_eLSA_CLR_5m$X_tax))
parsed_eLSA_CLR_5m$X_tax_abr[grep("kingdom_Eukaryota;supergroup_Archaeplastida;division_Chlorophyta", parsed_eLSA_CLR_5m$X_tax)]="Chlorophyta"

length(grep("kingdom_Eukaryota;supergroup_Archaeplastida;division_Chlorophyta;class_Mamiellophyceae;order_Mamiellales;family_Bathycoccaceae", parsed_eLSA_CLR_5m$X_tax))
parsed_eLSA_CLR_5m$X_tax_abr[grep("kingdom_Eukaryota;supergroup_Archaeplastida;division_Chlorophyta;class_Mamiellophyceae;order_Mamiellales;family_Bathycoccaceae", parsed_eLSA_CLR_5m$X_tax)]="Bathycoccus"

#Dinoflagellates
length(grep("kingdom_Eukaryota;supergroup_Alveolata;division_Dinoflagellata", parsed_eLSA_CLR_5m$X_tax))
parsed_eLSA_CLR_5m$X_tax_abr[grep("kingdom_Eukaryota;supergroup_Alveolata;division_Dinoflagellata", parsed_eLSA_CLR_5m$X_tax)]="Dinoflagellate"

length(grep("kingdom_Eukaryota;supergroup_Alveolata;division_Dinoflagellata;class_Syndiniales", parsed_eLSA_CLR_5m$X_tax))
parsed_eLSA_CLR_5m$X_tax_abr[grep("kingdom_Eukaryota;supergroup_Alveolata;division_Dinoflagellata;class_Syndiniales", parsed_eLSA_CLR_5m$X_tax)]="Syndiniales"

length(grep("kingdom_Eukaryota;supergroup_Alveolata;division_Dinoflagellata;class_Dinophyta", parsed_eLSA_CLR_5m$X_tax))
parsed_eLSA_CLR_5m$X_tax_abr[grep("kingdom_Eukaryota;supergroup_Alveolata;division_Dinoflagellata;class_Dinophyta", parsed_eLSA_CLR_5m$X_tax)]="Dinophyta"

length(grep("kingdom_Eukaryota;supergroup_Alveolata;division_Dinoflagellata;class_Dinophyceae", parsed_eLSA_CLR_5m$X_tax))
parsed_eLSA_CLR_5m$X_tax_abr[grep("kingdom_Eukaryota;supergroup_Alveolata;division_Dinoflagellata;class_Dinophyceae", parsed_eLSA_CLR_5m$X_tax)]="Dinophyceae"

length(grep("kingdom_Eukaryota;supergroup_Alveolata;division_Dinoflagellata;class_Dinophyceae;order_Gymnodiniales;family_Gymnodiniaceae;genus_Lepidodinium", parsed_eLSA_CLR_5m$X_tax)) #174
parsed_eLSA_CLR_5m$X_tax_abr[grep("kingdom_Eukaryota;supergroup_Alveolata;division_Dinoflagellata;class_Dinophyceae;order_Gymnodiniales;family_Gymnodiniaceae;genus_Lepidodinium", parsed_eLSA_CLR_5m$X_tax)]="Lepidodinium"

length(grep("kingdom_Eukaryota;supergroup_Alveolata;division_Dinoflagellata;class_Dinophyceae;order_Gymnodiniales", parsed_eLSA_CLR_5m$X_tax)) #1062
parsed_eLSA_CLR_5m$X_tax_abr[grep("kingdom_Eukaryota;supergroup_Alveolata;division_Dinoflagellata;class_Dinophyceae;order_Gymnodiniales", parsed_eLSA_CLR_5m$X_tax)]="Gymnodinium"

length(grep("kingdom_Eukaryota;supergroup_Alveolata;division_Ciliophora;class_Spirotrichea;order_Tintinnida", parsed_eLSA_CLR_5m$X_tax))
parsed_eLSA_CLR_5m$X_tax_abr[grep("kingdom_Eukaryota;supergroup_Alveolata;division_Ciliophora;class_Spirotrichea;order_Tintinnida", parsed_eLSA_CLR_5m$X_tax)]="Tintinnida"
  
length(grep(" kingdom_Eukaryota;supergroup_Alveolata;division_Ciliophora", parsed_eLSA_CLR_5m$X_tax)) #754
parsed_eLSA_CLR_5m$X_tax_abr[grep(" kingdom_Eukaryota;supergroup_Alveolata;division_Ciliophora;class_Spirotrichea;order_Strombidiida", parsed_eLSA_CLR_5m$X_tax)]="Ciliophora"

for(i in grep("Ciliophora", temp)){
  parsed_eLSA_CLR_5m$X_tax_abr[grep(temp[i], parsed_eLSA_CLR_5m$X_tax)]="Ciliophora"
}

length(grep("kingdom_Eukaryota;supergroup_Alveolata;division_Apicomplexa", parsed_eLSA_CLR_5m$X_tax))
parsed_eLSA_CLR_5m$X_tax_abr[grep("kingdom_Eukaryota;supergroup_Alveolata;division_Apicomplexa", parsed_eLSA_CLR_5m$X_tax)]="Ampicomplexa"

#OKAY DOUBLE CHECK IT 
length(grep("Eukaryote", parsed_eLSA_CLR_5m$X_tax_abr)) #270 only!! 
unique(parsed_eLSA_CLR_5m$X_tax[grep("Eukaryote", parsed_eLSA_CLR_5m$X_tax_abr)]) #FUCK YEAH!!!! 
parsed_eLSA_CLR_5m$X_tax_abr[grep("UCYN-A", parsed_eLSA_CLR_5m$X_tax)]="UCYN-A_1_AE"

#Abbreviated Y taxonomy-----
#Start by labelling everything "SAR"
parsed_eLSA_CLR_5m$Y_tax_abr="SAR"

#Manually change taxa, oh god :( #At the "supergroup" level 

##Start with the Stramenopiles: All stramenopiles except Pseudo-nitzschia and Chaetoceros set to "Stramenopiles"

length(grep("kingdom_Eukaryota;supergroup_Stramenopiles", parsed_eLSA_CLR_5m$Y_tax))
parsed_eLSA_CLR_5m$Y_tax_abr[grep("kingdom_Eukaryota;supergroup_Stramenopiles;", parsed_eLSA_CLR_5m$Y_tax)]="Stramenopile"
parsed_eLSA_CLR_5m$Y_tax_abr[grep(temp[17], parsed_eLSA_CLR_5m$Y_tax)]="Stramenopile"

length(grep("kingdom_Eukaryota;supergroup_Stramenopiles;division_Ochrophyta;class_Bacillariophyta;order_Bacillariophyta_X;family_Raphid-pennate;genus_Pseudo-nitzschia", parsed_eLSA_CLR_5m$Y_tax))
parsed_eLSA_CLR_5m$Y_tax_abr[grep("kingdom_Eukaryota;supergroup_Stramenopiles;division_Ochrophyta;class_Bacillariophyta;order_Bacillariophyta_X;family_Raphid-pennate;genus_Pseudo-nitzschia", parsed_eLSA_CLR_5m$Y_tax)]="Pseudo-nitzchia"

length(grep("kingdom_Eukaryota;supergroup_Stramenopiles;division_Ochrophyta;class_Bacillariophyta;order_Bacillariophyta_X;family_Polar-centric-Mediophyceae;genus_Chaetoceros", parsed_eLSA_CLR_5m$Y_tax))
parsed_eLSA_CLR_5m$Y_tax_abr[grep("kingdom_Eukaryota;supergroup_Stramenopiles;division_Ochrophyta;class_Bacillariophyta;order_Bacillariophyta_X;family_Polar-centric-Mediophyceae;genus_Chaetoceros", parsed_eLSA_CLR_5m$Y_tax)]="Chaetoceros"

length(grep("kingdom_Eukaryota;supergroup_Rhizaria", parsed_eLSA_CLR_5m$Y_tax))
parsed_eLSA_CLR_5m$Y_tax_abr[grep("Rhizaria", parsed_eLSA_CLR_5m$Y_tax)]="Rhizaria"

length(grep("kingdom_Eukaryota;supergroup_Opisthokonta", parsed_eLSA_CLR_5m$Y_tax))
parsed_eLSA_CLR_5m$Y_tax_abr[grep("kingdom_Eukaryota;supergroup_Opisthokonta", parsed_eLSA_CLR_5m$Y_tax)]="Opisthokonta"

##All hacrobia except prymnesiophytes, Haptophytes, and Phaeocystis set to "Hacrobia" 

length(grep("kingdom_Eukaryota;supergroup_Hacrobia", parsed_eLSA_CLR_5m$Y_tax))
parsed_eLSA_CLR_5m$Y_tax_abr[grep("kingdom_Eukaryota;supergroup_Hacrobia", parsed_eLSA_CLR_5m$Y_tax)]="Hacrobia"

length(grep("kingdom_Eukaryota;supergroup_Hacrobia;division_Haptophyta;class_Prymnesiophyceae;order_Prymnesiales;family_Prymnesiaceae;genus_Prymnesium;species_Prymnesium_sp.", parsed_eLSA_CLR_5m$Y_tax))
parsed_eLSA_CLR_5m$Y_tax_abr[grep("kingdom_Eukaryota;supergroup_Hacrobia;division_Haptophyta;class_Prymnesiophyceae;order_Prymnesiales;family_Prymnesiaceae;genus_Prymnesium;species_Prymnesium_sp.", parsed_eLSA_CLR_5m$Y_tax)]="Prymnesiophyte"

length(grep("kingdom_Eukaryota;supergroup_Hacrobia;division_Haptophyta;class_Prymnesiophyceae;order_Prymnesiales;family_Chrysochromulinaceae;genus_Chrysochromulina", parsed_eLSA_CLR_5m$Y_tax)) #281
parsed_eLSA_CLR_5m$Y_tax_abr[grep("kingdom_Eukaryota;supergroup_Hacrobia;division_Haptophyta;class_Prymnesiophyceae;order_Prymnesiales;family_Chrysochromulinaceae;genus_Chrysochromulina", parsed_eLSA_CLR_5m$Y_tax)]="Chrysochromulina"

#Phaeocystis
length(grep("kingdom_Eukaryota;supergroup_Hacrobia;division_Haptophyta;class_Prymnesiophyceae;order_Phaeocystales;family_Phaeocystaceae;genus_Phaeocystis", parsed_eLSA_CLR_5m$Y_tax))
parsed_eLSA_CLR_5m$Y_tax_abr[grep("kingdom_Eukaryota;supergroup_Hacrobia;division_Haptophyta;class_Prymnesiophyceae;order_Phaeocystales;family_Phaeocystaceae;genus_Phaeocystis", parsed_eLSA_CLR_5m$Y_tax)]="Phaeocystis"

#Haptophyte
length(grep("kingdom_Eukaryota;supergroup_Hacrobia;division_Haptophyta;class_Haptophyta", parsed_eLSA_CLR_5m$Y_tax))
parsed_eLSA_CLR_5m$Y_tax_abr[grep("kingdom_Eukaryota;supergroup_Hacrobia;division_Haptophyta;class_Haptophyta", parsed_eLSA_CLR_5m$Y_tax)]="Haptophyte"

#Insert Brad ix's
length(grep("kingdom_Eukaryota;supergroup_Hacrobia;division_Haptophyta;class_Prymnesiophyceae;order_Prymnesiophyceae_X;family_Braarudosphaeraceae;genus_Braarudosphaeraceae", parsed_eLSA_CLR_5m$Y_tax)) #47
parsed_eLSA_CLR_5m$Y[grep("kingdom_Eukaryota;supergroup_Hacrobia;division_Haptophyta;class_Prymnesiophyceae;order_Prymnesiophyceae_X;family_Braarudosphaeraceae;genus_Braarudosphaeraceae_X", parsed_eLSA_CLR_5m$Y_tax)]

for(i in Brad_ix){
  if(length(grep(i, parsed_eLSA_CLR_5m$Y))>1){
    print(i)
    print(grep(i, parsed_eLSA_CLR_5m$Y)) 
  }
} #Hmm ok, only 3 Brad ASVs showed up: #3 (formerly 4), 4 (formerly 5), 7x (formerly 7) 

parsed_eLSA_CLR_5m$Y_tax_abr[grep(Brad_ix[4], parsed_eLSA_CLR_5m$Y)]="Braarudospharea_bigelowii_3"
parsed_eLSA_CLR_5m$Y_tax_abr[grep(Brad_ix[5], parsed_eLSA_CLR_5m$Y)]="Braarudospharea_bigelowii_4"
parsed_eLSA_CLR_5m$Y_tax_abr[grep(Brad_ix[7], parsed_eLSA_CLR_5m$Y)]="Braarudosphaera_sp_7x"

#Moving on!
length(grep("kingdom_Eukaryota;supergroup_Archaeplastida;division_Rhodophyta", parsed_eLSA_CLR_5m$Y_tax))
parsed_eLSA_CLR_5m$Y_tax_abr[grep("kingdom_Eukaryota;supergroup_Archaeplastida;division_Rhodophyta", parsed_eLSA_CLR_5m$Y_tax)]="Archeplastida"

length(grep("kingdom_Eukaryota;supergroup_Archaeplastida;division_Chlorophyta", parsed_eLSA_CLR_5m$Y_tax))
parsed_eLSA_CLR_5m$Y_tax_abr[grep("kingdom_Eukaryota;supergroup_Archaeplastida;division_Chlorophyta", parsed_eLSA_CLR_5m$Y_tax)]="Chlorophyta"

length(grep("kingdom_Eukaryota;supergroup_Archaeplastida;division_Chlorophyta;class_Mamiellophyceae;order_Mamiellales;family_Bathycoccaceae", parsed_eLSA_CLR_5m$Y_tax))
parsed_eLSA_CLR_5m$Y_tax_abr[grep("kingdom_Eukaryota;supergroup_Archaeplastida;division_Chlorophyta;class_Mamiellophyceae;order_Mamiellales;family_Bathycoccaceae", parsed_eLSA_CLR_5m$Y_tax)]="Bathycoccus"

#Dinoflagellates
length(grep("kingdom_Eukaryota;supergroup_Alveolata;division_Dinoflagellata", parsed_eLSA_CLR_5m$Y_tax))
parsed_eLSA_CLR_5m$Y_tax_abr[grep("kingdom_Eukaryota;supergroup_Alveolata;division_Dinoflagellata", parsed_eLSA_CLR_5m$Y_tax)]="Dinoflagellate"

length(grep("kingdom_Eukaryota;supergroup_Alveolata;division_Dinoflagellata;class_Syndiniales", parsed_eLSA_CLR_5m$Y_tax))
parsed_eLSA_CLR_5m$Y_tax_abr[grep("kingdom_Eukaryota;supergroup_Alveolata;division_Dinoflagellata;class_Syndiniales", parsed_eLSA_CLR_5m$Y_tax)]="Syndiniales"

length(grep("kingdom_Eukaryota;supergroup_Alveolata;division_Dinoflagellata;class_Dinophyta", parsed_eLSA_CLR_5m$Y_tax))
parsed_eLSA_CLR_5m$Y_tax_abr[grep("kingdom_Eukaryota;supergroup_Alveolata;division_Dinoflagellata;class_Dinophyta", parsed_eLSA_CLR_5m$Y_tax)]="Dinophyta"

length(grep("kingdom_Eukaryota;supergroup_Alveolata;division_Dinoflagellata;class_Dinophyceae", parsed_eLSA_CLR_5m$Y_tax))
parsed_eLSA_CLR_5m$Y_tax_abr[grep("kingdom_Eukaryota;supergroup_Alveolata;division_Dinoflagellata;class_Dinophyceae", parsed_eLSA_CLR_5m$Y_tax)]="Dinophyceae"

length(grep("kingdom_Eukaryota;supergroup_Alveolata;division_Dinoflagellata;class_Dinophyceae;order_Gymnodiniales;family_Gymnodiniaceae;genus_Lepidodinium", parsed_eLSA_CLR_5m$Y_tax)) #105
parsed_eLSA_CLR_5m$Y_tax_abr[grep("kingdom_Eukaryota;supergroup_Alveolata;division_Dinoflagellata;class_Dinophyceae;order_Gymnodiniales;family_Gymnodiniaceae;genus_Lepidodinium", parsed_eLSA_CLR_5m$Y_tax)]="Lepidodinium"

length(grep("kingdom_Eukaryota;supergroup_Alveolata;division_Dinoflagellata;class_Dinophyceae;order_Gymnodiniales", parsed_eLSA_CLR_5m$Y_tax)) #657
parsed_eLSA_CLR_5m$Y_tax_abr[grep("kingdom_Eukaryota;supergroup_Alveolata;division_Dinoflagellata;class_Dinophyceae;order_Gymnodiniales", parsed_eLSA_CLR_5m$Y_tax)]="Gymnodinium"

length(grep("kingdom_Eukaryota;supergroup_Alveolata;division_Ciliophora;class_Spirotrichea;order_Tintinnida", parsed_eLSA_CLR_5m$Y_tax))
parsed_eLSA_CLR_5m$Y_tax_abr[grep("kingdom_Eukaryota;supergroup_Alveolata;division_Ciliophora;class_Spirotrichea;order_Tintinnida", parsed_eLSA_CLR_5m$Y_tax)]="Tintinnida"

length(grep(" kingdom_Eukaryota;supergroup_Alveolata;division_Ciliophora", parsed_eLSA_CLR_5m$Y_tax)) #1504
parsed_eLSA_CLR_5m$Y_tax_abr[grep(" kingdom_Eukaryota;supergroup_Alveolata;division_Ciliophora;class_Spirotrichea;order_Strombidiida", parsed_eLSA_CLR_5m$Y_tax)]="Ciliophora"

for(i in grep("Ciliophora", temp)){
  parsed_eLSA_CLR_5m$Y_tax_abr[grep(temp[i], parsed_eLSA_CLR_5m$Y_tax)]="Ciliophora"
}

length(grep("kingdom_Eukaryota;supergroup_Alveolata;division_Apicomplexa", parsed_eLSA_CLR_5m$Y_tax))
parsed_eLSA_CLR_5m$Y_tax_abr[grep("kingdom_Eukaryota;supergroup_Alveolata;division_Apicomplexa", parsed_eLSA_CLR_5m$Y_tax)]="Ampicomplexa"

#Okay, now the last of them
length(grep("SAR", parsed_eLSA_CLR_5m$Y_tax_abr)) #1121 
sort(unique(parsed_eLSA_CLR_5m$Y_tax[grep("SAR", parsed_eLSA_CLR_5m$Y_tax_abr)]))
#EFF
temp=sort(unique(parsed_eLSA_CLR_5m$Y_tax[grep("SAR", parsed_eLSA_CLR_5m$Y_tax_abr)]))
length(grep("Ciliophora", temp)) #14 + 6 UCYN-A ASVs
#Go back and assign Ciliophora in a for-loop

#Assign taxonomy for UCYN-A 
grep("UCYN-A", parsed_eLSA_CLR_5m$Y_tax)
parsed_eLSA_CLR_5m$Y_tax_abr[grep("UCYN-A", parsed_eLSA_CLR_5m$Y_tax)]=parsed_eLSA_CLR_5m$Y_tax[grep("UCYN-A", parsed_eLSA_CLR_5m$Y_tax)]

#Now all the rest should just be "Eukaryote"
length(grep("SAR", parsed_eLSA_CLR_5m$Y_tax_abr)) #236
parsed_eLSA_CLR_5m$Y_tax[grep("SAR", parsed_eLSA_CLR_5m$Y_tax_abr)]
parsed_eLSA_CLR_5m$Y_tax_abr[grep("SAR", parsed_eLSA_CLR_5m$Y_tax_abr)]="Eukaryote"

length(grep("SAR", parsed_eLSA_CLR_5m$Y_tax_abr)) #0 #YAY! 

#Sort by either SSCC or LS + length #Try to do this in Cytoskape first 

#Write it out!
write.table(x=parsed_eLSA_CLR_5m, file="eLSA_output/piecewise_eLSA/CLR_5m/parsed_eLSA_output_CLR_5m_11.24.2020.tsv", sep="\t", quote=F, row.names=F)

#Parse 5m CLR transformed data for supplementary table----
setwd("/Users/colettef/Desktop/DataAnalysis/SPOT_16S_18S_eLSA/")

eLSA_data=read.table("eLSA_output/piecewise_eLSA/CLR_5m/parsed_eLSA_output_CLR_5m_11.24.2020.tsv", header=T, stringsAsFactors = F, sep="\t")
dim(eLSA_data) #12313 33
names(eLSA_data)

#Export interaction data from Cytoskape session
inter_data=read.csv("eLSA_output/piecewise_eLSA/CLR_5m/UCYNA-all3ASVs_default_edge_datafrom_11.24.2020.csv", header=T, stringsAsFactors = F)
dim(inter_data)
names(inter_data)

eLSA_data$X_abr[grep("UCYN", eLSA_data$X_tax_abr)] #once
length(grep("UCYN", eLSA_data$Y_tax_abr)) #179 times
unique(eLSA_data$Y_abr[grep("UCYN", eLSA_data$Y_tax_abr)])

grep("e6425353849111275755617_A", inter_data$name) #That'x ix 6
unique(eLSA_data$Y_abr[grep("UCYN", eLSA_data$Y_tax_abr)])

#For each of the unique UCYN-A abbreviated names in eLSA data (X_abr or Y_abr)
#Pull those rows
euks=c(inter_data$name[grep(unique(eLSA_data$Y_abr[grep("UCYN", eLSA_data$Y_tax_abr)])[1], inter_data$name)],
      inter_data$name[grep(unique(eLSA_data$Y_abr[grep("UCYN", eLSA_data$Y_tax_abr)])[2], inter_data$name)], 
      inter_data$name[grep(unique(eLSA_data$Y_abr[grep("UCYN", eLSA_data$Y_tax_abr)])[3], inter_data$name)])
euks=as.data.frame(euks)
head(euks)

#Split column in two
library(tidyr)
library(dplyr)
DF_interact=separate(euks, col="euks", into=c("abr_hash", "UCYNA_ASV"), sep=" \\(interacts with\\) ")
head(DF_interact)
dim(DF_interact) #163 2
head(DF_interact$abr_hash) #YES!

#Figure out the taxa in them, comparing "names" column to abbreviated X or Y 
DF_interact$taxon_X="SAR11"
DF_interact$taxon_Y="SAR11"
DF_interact$Number=0
DF_interact$hash_X="hash"
DF_interact$hash_Y="hash"

for(i in c(1:nrow(DF_interact))){
  if(length(unique(eLSA_data$X_tax_abr[grep(DF_interact$abr_hash[i], eLSA_data$X_abr)]))>1){
    print(i)
    print(unique(eLSA_data$X_tax_abr[grep(DF_interact$abr_hash[i], eLSA_data$X_abr)]))
    print(unique(eLSA_data$X[grep(DF_interact$abr_hash[i], eLSA_data$X_abr)]))
  } else{
    DF_interact$taxon_X[i]=unique(eLSA_data$X_tax_abr[grep(DF_interact$abr_hash[i], eLSA_data$X_abr)])
    DF_interact$hash_X[i]=unique(eLSA_data$X[grep(DF_interact$abr_hash[i], eLSA_data$X_abr)])
  }
}
DF_interact$taxon_X[144]="Hacrobia"
DF_interact$hash_X[144]="55bc25d102de6079cc48ba48515e72e2"
grep("SAR11", DF_interact$taxon_X) #YUSS
sort(unique(DF_interact$taxon_X))

#Compare X and Y taxonomy columns 
for(i in c(1:nrow(DF_interact))){
  if(length(unique(eLSA_data$Y_tax_abr[grep(DF_interact$abr_hash[i], eLSA_data$Y_abr)]))>1){
    print(i)
    print(unique(eLSA_data$Y_tax_abr[grep(DF_interact$abr_hash[i], eLSA_data$Y_abr)]))
    print(unique(eLSA_data$Y[grep(DF_interact$abr_hash[i], eLSA_data$Y_abr)]))
  } else{
    DF_interact$taxon_Y[i]=unique(eLSA_data$Y_tax_abr[grep(DF_interact$abr_hash[i], eLSA_data$Y_abr)])
    DF_interact$hash_Y[i]=unique(eLSA_data$Y[grep(DF_interact$abr_hash[i], eLSA_data$Y_abr)])
  }
}
DF_interact$taxon_Y[144]="Hacrobia"
DF_interact$hash_Y[144]="55bc25d102de6079cc48ba48515e72e2"
grep("SAR11", DF_interact$taxon_Y) #YUSS
sort(unique(DF_interact$taxon_Y))

sort(unique(c(DF_interact$taxon_X, DF_interact$taxon_Y)))
DF_interact$taxon_X==DF_interact$taxon_Y #Not all
DF_interact$hash_X==DF_interact$hash_Y #all true?
#Where do discrepancies between X tax and Y tax come from?
View(DF_interact[which(DF_interact$taxon_X!=DF_interact$taxon_Y),])
#3 Lepidodinium, one dinophyceae and a pseudonitzchia 
which(DF_interact$abr_hash=="c114")
DF_interact$taxon_X[which(DF_interact$taxon_X!=DF_interact$taxon_Y)]=DF_interact$taxon_Y[which(DF_interact$taxon_X!=DF_interact$taxon_Y)]

#From here on out, use X taxonomy 

#Collapse taxonomy levels to just be what they are on figure
unique(DF_interact$taxon_X)
#Get rid of: Ciliophora (-> Alveolata), Dinophyceae (-> Dinoflagellate?), Hacrobia (-> Prymnesiophyte), 
#Missing: Alveolata, Dinoflagellate

#Okay, make these changes 
DF_interact$taxon_X[which(DF_interact$taxon_X=="Ciliophora")]="Alveolata"
DF_interact$taxon_X[which(DF_interact$taxon_X=="Dinophyceae")]="Dinoflagellate"
DF_interact$taxon_X[which(DF_interact$taxon_X=="Hacrobia")]="Prymnesiophyte"
sort(unique(DF_interact$taxon_X))

#Then put a column if it interacts with UCYNA1 in AE/ Durapore or ASV5
DF_interact$UCYNA1_larger=""
DF_interact$UCYNA1_larger[which(DF_interact$UCYNA_ASV=="385241044219295525187_A")]="x"
DF_interact$UCYNA1_smaller=""
DF_interact$UCYNA1_smaller[which(DF_interact$UCYNA_ASV=="385241044219295525187_D")]="y"
DF_interact$UCYNA2=""
DF_interact$UCYNA2[which(DF_interact$UCYNA_ASV=="a11913318571711407_A")]="z"

#Paste them
DF_interact$Interaction=paste(DF_interact$UCYNA1_larger, DF_interact$UCYNA1_smaller, DF_interact$UCYNA2, sep="")

#Subset to just take the columns I want
subs_interact=DF_interact[-grep("UCYN", DF_interact$taxon_X), grep("_X|Interact", names(DF_interact))]
dim(subs_interact) #161 3
names(subs_interact)[c(1,2)]=c("Taxonomy", "ASV_hash")
head(subs_interact)

subs_interact=group_by(subs_interact, ASV_hash) %>% summarize(Interaction=toString(Interaction), Taxonomy=unique(Taxonomy))
sort(unique(subs_interact$Interaction))
head(subs_interact)
dim(subs_interact) #82 2
length(unique(subs_interact$ASV_hash)) #82

# #subs_interact=group_by(subs_interact, ASV_hash) %>% summarize(Interaction=toString(Interaction))
# subs_interact=group_by(DF_interact, hash_X) %>% select(grep("_X|Interact", names(DF_interact))) %>% mutate(Interaction=toString(Interaction))
# names(subs_interact)[c(1,2)]=c("Taxonomy", "ASV_hash")

#Change X's and Y's to be "UCYN-A1_AE" etc.
subs_interact$Interacts_with=paste("UCYN-A1_AE", "UCYN-A1_Durapore", "UCYN-A2", sep=", ")
head(subs_interact$Interacts_with)
sort(unique(subs_interact$Interacts_with))
sort(unique(subs_interact$Interaction))

subs_interact$Interacts_with[which(subs_interact$Interaction %in% c("x", "x, x"))]="UCYN-A1_1-80um"
subs_interact$Interacts_with[which(subs_interact$Interaction %in% c("x, x, y", "x, x, y, y", "x, y", "x, y, y", "y, y, x"))]=paste("UCYN-A1_1-80um", "UCYN-A1_0.22-1um", sep=", ")
subs_interact$Interacts_with[which(subs_interact$Interaction %in% c("x, x, z", "x, x, z, z", "x, z", "x, z, z"))]=paste("UCYN-A1_1-80um", "UCYN-A2", sep=", ")
subs_interact$Interacts_with[which(subs_interact$Interaction %in% c("y", "y, y"))]="UCYN-A1_0.22-1um"
subs_interact$Interacts_with[which(subs_interact$Interaction %in% c("y, z"))] #None #No yz
subs_interact$Interacts_with[which(subs_interact$Interaction %in% c("z", "z, z"))]="UCYN-A2"
subs_interact$Interacts_with[which(subs_interact$Interaction %in% c("x, y, z, z", "x, x, y, y, z", "x, x, y, y, z", "x, y, y, z, z"))]=paste("UCYN-A1_1-80um", "UCYN-A1_0.22-1um", "UCYN-A2", sep=", ")

subs_interact$Interaction[grep("AE|Durapore", subs_interact$Interacts_with)]
sort(unique(subs_interact$Interacts_with))
dim(subs_interact) #161 4

names(subs_interact)
subs_interact=subs_interact[,c(3, 1, 4)]
head(subs_interact)

#WRITE IT OUT!
write.table(subs_interact, file="eLSA_output/piecewise_eLSA/CLR_5m/supp_table_18S_hashes_tax_assoc_09.13.2022_B.tsv", row.names=F, quote=F, sep="\t")

#2. CLR transformed data from DCM----

#Read in the data
directory <- getwd()
filenames<-system(paste('ls ', directory, '/eLSA_output/piecewise_eLSA/CLR_DCM', sep=''), intern=TRUE)
filenames

#How to remove data objects en masse: 
rm(list=ls()[grep("eLSA_data", ls())])

for(i in c(1:6)){
  a <- strsplit(filenames[i], split="_", fixed=T) 
  print(paste("eLSA_data", a[[1]][3], a[[1]][4], a[[1]][5], sep="_"))
  assign(x=paste("eLSA_data", a[[1]][3], a[[1]][4], a[[1]][5], sep="_"), value=read.table(file=paste("eLSA_output/piecewise_eLSA/CLR_DCM/", filenames[i], sep=""), sep="\t", stringsAsFactors = F, header=T))
}

#Rbind em together!
eLSA_data_concat=rbind(eLSA_data_subs_A_B, eLSA_data_subs_A_C, eLSA_data_subs_A_CLR, eLSA_data_subs_B_C, eLSA_data_subs_B_CLR, eLSA_data_subs_C_CLR)
dim(eLSA_data_concat) #268278 27 

#Assign q values----
library(qvalue)

#normal Q values
Q=qvalue(eLSA_data_concat$P)
length(eLSA_data_concat$Q)==length(Q$qvalues) #TRUE
eLSA_data_concat$Q=Q$qvalues

#Pearson's: Qpcc and Qspcc
pearsons_Q=qvalue(eLSA_data_concat$Ppcc)
length(pearsons_Q$qvalues)==nrow(eLSA_data_concat) #TRUE
eLSA_data_concat$Qpcc=pearsons_Q$qvalues

delay_pearsons_Q=qvalue(eLSA_data_concat$Pspcc)
length(delay_pearsons_Q$qvalues)==length(eLSA_data_concat$Pspcc) #TRUE
eLSA_data_concat$Qspcc=delay_pearsons_Q$qvalues

#Spearman's: Qscc and Qsscc
spearmans_Q=qvalue(eLSA_data_concat$Pscc)
length(spearmans_Q$qvalues)==length(eLSA_data_concat$Pscc)
eLSA_data_concat$Qscc=spearmans_Q$qvalues

delay_spearmans_Q=qvalue(eLSA_data_concat$Psscc)
length(delay_spearmans_Q$qvalues)==length(eLSA_data_concat$Psscc)
eLSA_data_concat$Qsscc=delay_spearmans_Q$qvalues

#Write out the file
write.table(eLSA_data_concat, file="eLSA_output/piecewise_eLSA/CLR_DCM/concat_eLSA_output_CLR_DCM_11.24.2020.tsv", sep="\t", quote=F, row.names=F)

#Parse down this massive df!
nrow(eLSA_data_concat) #268278
length(which(eLSA_data_concat$P<0.005 & eLSA_data_concat$Q<0.01)) #14692
length(which(eLSA_data_concat$P< 0.005 & eLSA_data_concat$Q<0.01 & eLSA_data_concat$Ppcc < 0.005 & eLSA_data_concat$Qpcc < 0.01 & eLSA_data_concat$Pscc < 0.005 & eLSA_data_concat$Qscc < 0.01)) 
#12397, or ~4% of the original data frame 
parsed_eLSA_CLR_DCM=eLSA_data_concat[which(eLSA_data_concat$P< 0.005 & eLSA_data_concat$Q<0.01 & eLSA_data_concat$Ppcc < 0.005 & eLSA_data_concat$Qpcc < 0.01 & eLSA_data_concat$Pscc < 0.005 & eLSA_data_concat$Qscc < 0.01),]
dim(parsed_eLSA_CLR_DCM)

#Abbreviated X and Y hashes
length(unique(parsed_eLSA_CLR_DCM$X)) #599
length(unique(abbreviate(parsed_eLSA_CLR_DCM$X, minlength=4))) #599
parsed_eLSA_CLR_DCM$X_abr=abbreviate(parsed_eLSA_CLR_DCM$X, minlength=4)
head(unique(parsed_eLSA_CLR_DCM$X_abr), n=15)

length(unique(parsed_eLSA_CLR_DCM$Y)) #571
length(unique(abbreviate(parsed_eLSA_CLR_DCM$Y))) #571
parsed_eLSA_CLR_DCM$Y_abr=abbreviate(parsed_eLSA_CLR_DCM$Y, minlength=4)
head(unique(parsed_eLSA_CLR_DCM$Y_abr), n=15)

#Add in taxonomy
#X taxa: 
parsed_eLSA_CLR_DCM$X_tax="SAR"

#Add in euks taxonomy, get indices for UCYNA
ix <- c()
for(i in c(1:nrow(parsed_eLSA_CLR_DCM))){
  a <- parsed_eLSA_CLR_DCM$X[i]
  b <- grep(a, euks_tax$Feature.ID)
  if(as.numeric(length(b)<1)){
    print(b)
    ix <- c(ix, i)
  }else{
    parsed_eLSA_CLR_DCM$X_tax[i]=euks_tax$Taxon[b]
  }
}
#Why is ix empty?? 

grep("_AE|_D", parsed_eLSA_CLR_DCM$X) #UCYN-A is not in the X column 
length(grep("_AE|_D", parsed_eLSA_CLR_DCM$Y)) #But it appears in the Y column 12 times #Is any of that statistically significant?? 
parsed_eLSA_CLR_DCM$LS[grep("_AE|_D", parsed_eLSA_CLR_DCM$Y)]
parsed_eLSA_CLR_DCM$Len[grep("_AE|_D", parsed_eLSA_CLR_DCM$Y)] #Maybe!! 

length(grep("SAR", parsed_eLSA_CLR_DCM$X_tax)) #0

#Add in Y taxonomy, get ix's for UCYN-A
parsed_eLSA_CLR_DCM$Y_tax="SAR"

ix <- c()
for(i in c(1:nrow(parsed_eLSA_CLR_DCM))){
  a <- parsed_eLSA_CLR_DCM$Y[i]
  b <- grep(a, euks_tax$Feature.ID)
  if(as.numeric(length(b)<1)){
    ix <- c(ix, i)
  }else{
    parsed_eLSA_CLR_DCM$Y_tax[i]=euks_tax$Taxon[b]
  }
}  

length(ix)
parsed_eLSA_CLR_DCM$Y[ix] #Only in the AE size fraction #And it's mostly 6115, the new one #HMM

#Now add in UCYN-A
for(i in ix){
  a <- parsed_eLSA_CLR_DCM$Y[i]
  b <- strsplit(a, split="_", fixed=T)
  parsed_eLSA_CLR_DCM$Y_tax[i]=paste("UCYN-A", grep(b[[1]][1], UCYNA_ix), b[[1]][2], sep="_")
}

#Abbreviate X taxonomy-----
#Start by labelling everything "SAR"
parsed_eLSA_CLR_DCM$X_tax_abr="SAR"

#Manually change taxa, oh god :( #At the "supergroup" level 

##Start with the Stramenopiles: All stramenopiles except Pseudo-nitzschia and Chaetoceros set to "Stramenopiles"

length(grep("kingdom_Eukaryota;supergroup_Stramenopiles", parsed_eLSA_CLR_DCM$X_tax)) #That's like 10% of the X taxa
parsed_eLSA_CLR_DCM$X_tax_abr[grep("kingdom_Eukaryota;supergroup_Stramenopiles;", parsed_eLSA_CLR_DCM$X_tax)]="Stramenopile"
parsed_eLSA_CLR_DCM$X_tax_abr[grep(temp[grep("Stramenopile", temp)], parsed_eLSA_CLR_DCM$X_tax)]="Stramenopile"

length(grep("kingdom_Eukaryota;supergroup_Stramenopiles;division_Ochrophyta;class_Bacillariophyta;order_Bacillariophyta_X;family_Raphid-pennate;genus_Pseudo-nitzschia", parsed_eLSA_CLR_DCM$X_tax))
parsed_eLSA_CLR_DCM$X_tax_abr[grep("kingdom_Eukaryota;supergroup_Stramenopiles;division_Ochrophyta;class_Bacillariophyta;order_Bacillariophyta_X;family_Raphid-pennate;genus_Pseudo-nitzschia", parsed_eLSA_CLR_DCM$X_tax)]="Pseudo-nitzchia"

length(grep("kingdom_Eukaryota;supergroup_Stramenopiles;division_Ochrophyta;class_Bacillariophyta;order_Bacillariophyta_X;family_Polar-centric-Mediophyceae;genus_Chaetoceros", parsed_eLSA_CLR_DCM$X_tax))
parsed_eLSA_CLR_DCM$X_tax_abr[grep("kingdom_Eukaryota;supergroup_Stramenopiles;division_Ochrophyta;class_Bacillariophyta;order_Bacillariophyta_X;family_Polar-centric-Mediophyceae;genus_Chaetoceros", parsed_eLSA_CLR_DCM$X_tax)]="Chaetoceros"

length(grep("kingdom_Eukaryota;supergroup_Rhizaria", parsed_eLSA_CLR_DCM$X_tax))
parsed_eLSA_CLR_DCM$X_tax_abr[grep("Rhizaria", parsed_eLSA_CLR_DCM$X_tax)]="Rhizaria"

length(grep("kingdom_Eukaryota;supergroup_Opisthokonta", parsed_eLSA_CLR_DCM$X_tax))
parsed_eLSA_CLR_DCM$X_tax_abr[grep("kingdom_Eukaryota;supergroup_Opisthokonta", parsed_eLSA_CLR_DCM$X_tax)]="Opisthokonta"

##All hacrobia except prymnesiophytes, Haptophytes, and Phaeocystis set to "Hacrobia" 

length(grep("kingdom_Eukaryota;supergroup_Hacrobia", parsed_eLSA_CLR_DCM$X_tax))
parsed_eLSA_CLR_DCM$X_tax_abr[grep("kingdom_Eukaryota;supergroup_Hacrobia", parsed_eLSA_CLR_DCM$X_tax)]="Hacrobia"

length(grep("kingdom_Eukaryota;supergroup_Hacrobia;division_Haptophyta;class_Prymnesiophyceae;order_Prymnesiales;family_Prymnesiaceae;genus_Prymnesium;species_Prymnesium_sp.", parsed_eLSA_CLR_DCM$X_tax))
parsed_eLSA_CLR_DCM$X_tax_abr[grep("kingdom_Eukaryota;supergroup_Hacrobia;division_Haptophyta;class_Prymnesiophyceae;order_Prymnesiales;family_Prymnesiaceae;genus_Prymnesium;species_Prymnesium_sp.", parsed_eLSA_CLR_DCM$X_tax)]="Prymnesiophyte"

length(grep("kingdom_Eukaryota;supergroup_Hacrobia;division_Haptophyta;class_Prymnesiophyceae;order_Prymnesiales;family_Chrysochromulinaceae;genus_Chrysochromulina", parsed_eLSA_CLR_DCM$X_tax)) #388
parsed_eLSA_CLR_DCM$X_tax_abr[grep("kingdom_Eukaryota;supergroup_Hacrobia;division_Haptophyta;class_Prymnesiophyceae;order_Prymnesiales;family_Chrysochromulinaceae;genus_Chrysochromulina", parsed_eLSA_CLR_DCM$X_tax)]="Chrysochromulina"

#Phaeocystis
length(grep("kingdom_Eukaryota;supergroup_Hacrobia;division_Haptophyta;class_Prymnesiophyceae;order_Phaeocystales;family_Phaeocystaceae;genus_Phaeocystis", parsed_eLSA_CLR_DCM$X_tax))
parsed_eLSA_CLR_DCM$X_tax_abr[grep("kingdom_Eukaryota;supergroup_Hacrobia;division_Haptophyta;class_Prymnesiophyceae;order_Phaeocystales;family_Phaeocystaceae;genus_Phaeocystis", parsed_eLSA_CLR_DCM$X_tax)]="Phaeocystis"

#Haptophyte
length(grep("kingdom_Eukaryota;supergroup_Hacrobia;division_Haptophyta;class_Haptophyta", parsed_eLSA_CLR_DCM$X_tax))
parsed_eLSA_CLR_DCM$X_tax_abr[grep("kingdom_Eukaryota;supergroup_Hacrobia;division_Haptophyta;class_Haptophyta", parsed_eLSA_CLR_DCM$X_tax)]="Haptophyte"

#Insert Brad ix's
length(grep("kingdom_Eukaryota;supergroup_Hacrobia;division_Haptophyta;class_Prymnesiophyceae;order_Prymnesiophyceae_X;family_Braarudosphaeraceae;genus_Braarudosphaeraceae_X", parsed_eLSA_CLR_DCM$X_tax)) #106
parsed_eLSA_CLR_DCM$X[grep("kingdom_Eukaryota;supergroup_Hacrobia;division_Haptophyta;class_Prymnesiophyceae;order_Prymnesiophyceae_X;family_Braarudosphaeraceae;genus_Braarudosphaeraceae_X", parsed_eLSA_CLR_DCM$X_tax)] #They are all #7 and #4

for(i in Brad_ix){
  if(length(grep(i, parsed_eLSA_CLR_DCM$X))>1){
    print(i)
    print(grep(i, parsed_eLSA_CLR_DCM$X)) 
  }
} #Hmm ok, only 3 Brad ASVs showed up: #3 (formerly 4), 4 (formerly 5), 7x (formerly 7) 

parsed_eLSA_CLR_DCM$X_tax_abr[grep(Brad_ix[4], parsed_eLSA_CLR_DCM$X)]="Braarudospharea_bigelowii_3"
parsed_eLSA_CLR_DCM$X_tax_abr[grep(Brad_ix[7], parsed_eLSA_CLR_DCM$X)]="Braarudosphaera_sp_7x"

#Moving on!
length(grep("kingdom_Eukaryota;supergroup_Archaeplastida;division_Rhodophyta", parsed_eLSA_CLR_DCM$X_tax))
parsed_eLSA_CLR_DCM$X_tax_abr[grep("kingdom_Eukaryota;supergroup_Archaeplastida;division_Rhodophyta", parsed_eLSA_CLR_DCM$X_tax)]="Archeplastida"

length(grep("kingdom_Eukaryota;supergroup_Archaeplastida;division_Chlorophyta", parsed_eLSA_CLR_DCM$X_tax))
parsed_eLSA_CLR_DCM$X_tax_abr[grep("kingdom_Eukaryota;supergroup_Archaeplastida;division_Chlorophyta", parsed_eLSA_CLR_DCM$X_tax)]="Chlorophyta"

length(grep("kingdom_Eukaryota;supergroup_Archaeplastida;division_Chlorophyta;class_Mamiellophyceae;order_Mamiellales;family_Bathycoccaceae", parsed_eLSA_CLR_DCM$X_tax))
parsed_eLSA_CLR_DCM$X_tax_abr[grep("kingdom_Eukaryota;supergroup_Archaeplastida;division_Chlorophyta;class_Mamiellophyceae;order_Mamiellales;family_Bathycoccaceae", parsed_eLSA_CLR_DCM$X_tax)]="Bathycoccus"

#Dinoflagellates
length(grep("kingdom_Eukaryota;supergroup_Alveolata;division_Dinoflagellata", parsed_eLSA_CLR_DCM$X_tax))
parsed_eLSA_CLR_DCM$X_tax_abr[grep("kingdom_Eukaryota;supergroup_Alveolata;division_Dinoflagellata", parsed_eLSA_CLR_DCM$X_tax)]="Dinoflagellate"

length(grep("kingdom_Eukaryota;supergroup_Alveolata;division_Dinoflagellata;class_Syndiniales", parsed_eLSA_CLR_DCM$X_tax))
parsed_eLSA_CLR_DCM$X_tax_abr[grep("kingdom_Eukaryota;supergroup_Alveolata;division_Dinoflagellata;class_Syndiniales", parsed_eLSA_CLR_DCM$X_tax)]="Syndiniales"

length(grep("kingdom_Eukaryota;supergroup_Alveolata;division_Dinoflagellata;class_Dinophyta", parsed_eLSA_CLR_DCM$X_tax)) #0
#parsed_eLSA_CLR_DCM$X_tax_abr[grep("kingdom_Eukaryota;supergroup_Alveolata;division_Dinoflagellata;class_Dinophyta", parsed_eLSA_CLR_DCM$X_tax)]="Dinophyta"

length(grep("kingdom_Eukaryota;supergroup_Alveolata;division_Dinoflagellata;class_Dinophyceae", parsed_eLSA_CLR_DCM$X_tax))
parsed_eLSA_CLR_DCM$X_tax_abr[grep("kingdom_Eukaryota;supergroup_Alveolata;division_Dinoflagellata;class_Dinophyceae", parsed_eLSA_CLR_DCM$X_tax)]="Dinophyceae"

length(grep("kingdom_Eukaryota;supergroup_Alveolata;division_Dinoflagellata;class_Dinophyceae;order_Gymnodiniales;family_Gymnodiniaceae;genus_Lepidodinium", parsed_eLSA_CLR_DCM$X_tax)) #65
parsed_eLSA_CLR_DCM$X_tax_abr[grep("kingdom_Eukaryota;supergroup_Alveolata;division_Dinoflagellata;class_Dinophyceae;order_Gymnodiniales;family_Gymnodiniaceae;genus_Lepidodinium", parsed_eLSA_CLR_DCM$X_tax)]="Lepidodinium"

length(grep("kingdom_Eukaryota;supergroup_Alveolata;division_Dinoflagellata;class_Dinophyceae;order_Gymnodiniales", parsed_eLSA_CLR_DCM$X_tax)) #606
parsed_eLSA_CLR_DCM$X_tax_abr[grep("kingdom_Eukaryota;supergroup_Alveolata;division_Dinoflagellata;class_Dinophyceae;order_Gymnodiniales", parsed_eLSA_CLR_DCM$X_tax)]="Gymnodinium"

length(grep("kingdom_Eukaryota;supergroup_Alveolata;division_Ciliophora;class_Spirotrichea;order_Tintinnida", parsed_eLSA_CLR_DCM$X_tax))
parsed_eLSA_CLR_DCM$X_tax_abr[grep("kingdom_Eukaryota;supergroup_Alveolata;division_Ciliophora;class_Spirotrichea;order_Tintinnida", parsed_eLSA_CLR_DCM$X_tax)]="Tintinnida"

length(grep(" kingdom_Eukaryota;supergroup_Alveolata;division_Ciliophora", parsed_eLSA_CLR_DCM$X_tax)) #754
parsed_eLSA_CLR_DCM$X_tax_abr[grep(" kingdom_Eukaryota;supergroup_Alveolata;division_Ciliophora;class_Spirotrichea;order_Strombidiida", parsed_eLSA_CLR_DCM$X_tax)]="Ciliophora"

for(i in grep("Ciliophora", temp)){
  parsed_eLSA_CLR_DCM$X_tax_abr[grep(temp[i], parsed_eLSA_CLR_DCM$X_tax)]="Ciliophora"
}

length(grep("kingdom_Eukaryota;supergroup_Alveolata;division_Apicomplexa", parsed_eLSA_CLR_DCM$X_tax))
parsed_eLSA_CLR_DCM$X_tax_abr[grep("kingdom_Eukaryota;supergroup_Alveolata;division_Apicomplexa", parsed_eLSA_CLR_DCM$X_tax)]="Ampicomplexa"

#one more!
length(grep("kingdom_Eukaryota;supergroup_Apusozoa", parsed_eLSA_CLR_DCM$X_tax))
parsed_eLSA_CLR_DCM$X_tax_abr[grep("kingdom_Eukaryota;supergroup_Apusozoa", parsed_eLSA_CLR_DCM$X_tax)]="Apusozoa"

#OKAY DOUBLE CHECK IT 
length(grep("SAR", parsed_eLSA_CLR_DCM$X_tax_abr)) #851 
temp=sort(unique(parsed_eLSA_CLR_DCM$X_tax[grep("SAR", parsed_eLSA_CLR_DCM$X_tax_abr)]))
length(temp) #18
length(grep("Ciliophora", temp)) #16
temp[grep("Ciliophora", temp, invert=T)]

#FINALLY, all the remainder are just "Eukaryotes" 
parsed_eLSA_CLR_DCM$X_tax_abr[grep("SAR", parsed_eLSA_CLR_DCM$X_tax_abr)]="Eukaryote"
unique(parsed_eLSA_CLR_DCM$X_tax_abr)

#parsed_eLSA_CLR_DCM$X_tax_abr[grep("UCYN-A", parsed_eLSA_CLR_DCM$X_tax)]="UCYN-A_1_AE"

#Abbreviate Y taxonomy-----
#Start by labelling everything "Eukaryote"
parsed_eLSA_CLR_DCM$Y_tax_abr="SAR"

#Manually change taxa, oh god :( #At the "supergroup" level 

##Start with the Stramenopiles: All stramenopiles except Pseudo-nitzschia and Chaetoceros set to "Stramenopiles"

length(grep("kingdom_Eukaryota;supergroup_Stramenopiles", parsed_eLSA_CLR_DCM$Y_tax)) #That's like 10% of the X taxa
parsed_eLSA_CLR_DCM$Y_tax_abr[grep("kingdom_Eukaryota;supergroup_Stramenopiles;", parsed_eLSA_CLR_DCM$Y_tax)]="Stramenopile"
parsed_eLSA_CLR_DCM$Y_tax_abr[grep(temp[grep("Stramenopile", temp)], parsed_eLSA_CLR_DCM$Y_tax)]="Stramenopile"

length(grep("kingdom_Eukaryota;supergroup_Stramenopiles;division_Ochrophyta;class_Bacillariophyta;order_Bacillariophyta_X;family_Raphid-pennate;genus_Pseudo-nitzschia", parsed_eLSA_CLR_DCM$Y_tax))
parsed_eLSA_CLR_DCM$Y_tax_abr[grep("kingdom_Eukaryota;supergroup_Stramenopiles;division_Ochrophyta;class_Bacillariophyta;order_Bacillariophyta_X;family_Raphid-pennate;genus_Pseudo-nitzschia", parsed_eLSA_CLR_DCM$Y_tax)]="Pseudo-nitzchia"

length(grep("kingdom_Eukaryota;supergroup_Stramenopiles;division_Ochrophyta;class_Bacillariophyta;order_Bacillariophyta_X;family_Polar-centric-Mediophyceae;genus_Chaetoceros", parsed_eLSA_CLR_DCM$Y_tax))
parsed_eLSA_CLR_DCM$Y_tax_abr[grep("kingdom_Eukaryota;supergroup_Stramenopiles;division_Ochrophyta;class_Bacillariophyta;order_Bacillariophyta_X;family_Polar-centric-Mediophyceae;genus_Chaetoceros", parsed_eLSA_CLR_DCM$Y_tax)]="Chaetoceros"

length(grep("kingdom_Eukaryota;supergroup_Rhizaria", parsed_eLSA_CLR_DCM$Y_tax))
parsed_eLSA_CLR_DCM$Y_tax_abr[grep("Rhizaria", parsed_eLSA_CLR_DCM$Y_tax)]="Rhizaria"

length(grep("kingdom_Eukaryota;supergroup_Opisthokonta", parsed_eLSA_CLR_DCM$Y_tax))
parsed_eLSA_CLR_DCM$Y_tax_abr[grep("kingdom_Eukaryota;supergroup_Opisthokonta", parsed_eLSA_CLR_DCM$Y_tax)]="Opisthokonta"

##All hacrobia except prymnesiophytes, Haptophytes, and Phaeocystis set to "Hacrobia" 

length(grep("kingdom_Eukaryota;supergroup_Hacrobia", parsed_eLSA_CLR_DCM$Y_tax))
parsed_eLSA_CLR_DCM$Y_tax_abr[grep("kingdom_Eukaryota;supergroup_Hacrobia", parsed_eLSA_CLR_DCM$Y_tax)]="Hacrobia"

length(grep("kingdom_Eukaryota;supergroup_Hacrobia;division_Haptophyta;class_Prymnesiophyceae;order_Prymnesiales;family_Prymnesiaceae;genus_Prymnesium;species_Prymnesium_sp.", parsed_eLSA_CLR_DCM$Y_tax))
parsed_eLSA_CLR_DCM$Y_tax_abr[grep("kingdom_Eukaryota;supergroup_Hacrobia;division_Haptophyta;class_Prymnesiophyceae;order_Prymnesiales;family_Prymnesiaceae;genus_Prymnesium;species_Prymnesium_sp.", parsed_eLSA_CLR_DCM$Y_tax)]="Prymnesiophyte"

length(grep("kingdom_Eukaryota;supergroup_Hacrobia;division_Haptophyta;class_Prymnesiophyceae;order_Prymnesiales;family_Chrysochromulinaceae;genus_Chrysochromulina", parsed_eLSA_CLR_DCM$Y_tax)) #434
parsed_eLSA_CLR_DCM$Y_tax_abr[grep("kingdom_Eukaryota;supergroup_Hacrobia;division_Haptophyta;class_Prymnesiophyceae;order_Prymnesiales;family_Chrysochromulinaceae;genus_Chrysochromulina", parsed_eLSA_CLR_DCM$Y_tax)]="Chrysochromulina"

#Phaeocystis
length(grep("kingdom_Eukaryota;supergroup_Hacrobia;division_Haptophyta;class_Prymnesiophyceae;order_Phaeocystales;family_Phaeocystaceae;genus_Phaeocystis", parsed_eLSA_CLR_DCM$Y_tax))
parsed_eLSA_CLR_DCM$Y_tax_abr[grep("kingdom_Eukaryota;supergroup_Hacrobia;division_Haptophyta;class_Prymnesiophyceae;order_Phaeocystales;family_Phaeocystaceae;genus_Phaeocystis", parsed_eLSA_CLR_DCM$Y_tax)]="Phaeocystis"

#Haptophyte
length(grep("kingdom_Eukaryota;supergroup_Hacrobia;division_Haptophyta;class_Haptophyta", parsed_eLSA_CLR_DCM$Y_tax)) #only 15...
parsed_eLSA_CLR_DCM$Y_tax_abr[grep("kingdom_Eukaryota;supergroup_Hacrobia;division_Haptophyta;class_Haptophyta", parsed_eLSA_CLR_DCM$Y_tax)]="Haptophyte"

#Insert Brad ix's
length(grep("kingdom_Eukaryota;supergroup_Hacrobia;division_Haptophyta;class_Prymnesiophyceae;order_Prymnesiophyceae_X;family_Braarudosphaeraceae;genus_Braarudosphaeraceae", parsed_eLSA_CLR_DCM$Y_tax)) #1
parsed_eLSA_CLR_DCM$Y[grep("kingdom_Eukaryota;supergroup_Hacrobia;division_Haptophyta;class_Prymnesiophyceae;order_Prymnesiophyceae_X;family_Braarudosphaeraceae;genus_Braarudosphaeraceae_X", parsed_eLSA_CLR_DCM$Y_tax)] #There is only 1 #7 #0r 7? 

for(i in Brad_ix){
  if(length(grep(i, parsed_eLSA_CLR_DCM$Y))>1){
    print(i)
    print(grep(i, parsed_eLSA_CLR_DCM$Y)) 
  }
}  

#parsed_eLSA_CLR_DCM$Y_tax_abr[grep(Brad_ix[4], parsed_eLSA_CLR_DCM$X)]="Braarudospharea_bigelowii_3"
parsed_eLSA_CLR_DCM$Y_tax_abr[grep(Brad_ix[7], parsed_eLSA_CLR_DCM$Y)]="Braarudosphaera_sp_7x"

#Moving on!
length(grep("kingdom_Eukaryota;supergroup_Archaeplastida;division_Rhodophyta", parsed_eLSA_CLR_DCM$Y_tax))
parsed_eLSA_CLR_DCM$Y_tax_abr[grep("kingdom_Eukaryota;supergroup_Archaeplastida;division_Rhodophyta", parsed_eLSA_CLR_DCM$Y_tax)]="Archeplastida"

length(grep("kingdom_Eukaryota;supergroup_Archaeplastida;division_Chlorophyta", parsed_eLSA_CLR_DCM$Y_tax))
parsed_eLSA_CLR_DCM$Y_tax_abr[grep("kingdom_Eukaryota;supergroup_Archaeplastida;division_Chlorophyta", parsed_eLSA_CLR_DCM$Y_tax)]="Chlorophyta"

length(grep("kingdom_Eukaryota;supergroup_Archaeplastida;division_Chlorophyta;class_Mamiellophyceae;order_Mamiellales;family_Bathycoccaceae", parsed_eLSA_CLR_DCM$Y_tax))
parsed_eLSA_CLR_DCM$Y_tax_abr[grep("kingdom_Eukaryota;supergroup_Archaeplastida;division_Chlorophyta;class_Mamiellophyceae;order_Mamiellales;family_Bathycoccaceae", parsed_eLSA_CLR_DCM$Y_tax)]="Bathycoccus"

#Dinoflagellates
length(grep("kingdom_Eukaryota;supergroup_Alveolata;division_Dinoflagellata", parsed_eLSA_CLR_DCM$Y_tax))
parsed_eLSA_CLR_DCM$Y_tax_abr[grep("kingdom_Eukaryota;supergroup_Alveolata;division_Dinoflagellata", parsed_eLSA_CLR_DCM$Y_tax)]="Dinoflagellate"

length(grep("kingdom_Eukaryota;supergroup_Alveolata;division_Dinoflagellata;class_Syndiniales", parsed_eLSA_CLR_DCM$Y_tax))
parsed_eLSA_CLR_DCM$Y_tax_abr[grep("kingdom_Eukaryota;supergroup_Alveolata;division_Dinoflagellata;class_Syndiniales", parsed_eLSA_CLR_DCM$Y_tax)]="Syndiniales"

length(grep("kingdom_Eukaryota;supergroup_Alveolata;division_Dinoflagellata;class_Dinophyta", parsed_eLSA_CLR_DCM$Y_tax)) #0
#parsed_eLSA_CLR_DCM$Y_tax_abr[grep("kingdom_Eukaryota;supergroup_Alveolata;division_Dinoflagellata;class_Dinophyta", parsed_eLSA_CLR_DCM$Y_tax)]="Dinophyta"

length(grep("kingdom_Eukaryota;supergroup_Alveolata;division_Dinoflagellata;class_Dinophyceae", parsed_eLSA_CLR_DCM$Y_tax))
parsed_eLSA_CLR_DCM$Y_tax_abr[grep("kingdom_Eukaryota;supergroup_Alveolata;division_Dinoflagellata;class_Dinophyceae", parsed_eLSA_CLR_DCM$Y_tax)]="Dinophyceae"

length(grep("kingdom_Eukaryota;supergroup_Alveolata;division_Dinoflagellata;class_Dinophyceae;order_Gymnodiniales;family_Gymnodiniaceae;genus_Lepidodinium", parsed_eLSA_CLR_DCM$Y_tax)) #13
parsed_eLSA_CLR_DCM$Y_tax_abr[grep("kingdom_Eukaryota;supergroup_Alveolata;division_Dinoflagellata;class_Dinophyceae;order_Gymnodiniales;family_Gymnodiniaceae;genus_Lepidodinium", parsed_eLSA_CLR_DCM$Y_tax)]="Lepidodinium"

length(grep("kingdom_Eukaryota;supergroup_Alveolata;division_Dinoflagellata;class_Dinophyceae;order_Gymnodiniales", parsed_eLSA_CLR_DCM$Y_tax)) #384
parsed_eLSA_CLR_DCM$Y_tax_abr[grep("kingdom_Eukaryota;supergroup_Alveolata;division_Dinoflagellata;class_Dinophyceae;order_Gymnodiniales", parsed_eLSA_CLR_DCM$Y_tax)]="Gymnodinium"

length(grep("kingdom_Eukaryota;supergroup_Alveolata;division_Ciliophora;class_Spirotrichea;order_Tintinnida", parsed_eLSA_CLR_DCM$Y_tax))
parsed_eLSA_CLR_DCM$Y_tax_abr[grep("kingdom_Eukaryota;supergroup_Alveolata;division_Ciliophora;class_Spirotrichea;order_Tintinnida", parsed_eLSA_CLR_DCM$Y_tax)]="Tintinnida"

length(grep(" kingdom_Eukaryota;supergroup_Alveolata;division_Ciliophora", parsed_eLSA_CLR_DCM$Y_tax)) #1681
parsed_eLSA_CLR_DCM$Y_tax_abr[grep(" kingdom_Eukaryota;supergroup_Alveolata;division_Ciliophora;class_Spirotrichea;order_Strombidiida", parsed_eLSA_CLR_DCM$Y_tax)]="Ciliophora"

for(i in grep("Ciliophora", temp)){
  parsed_eLSA_CLR_DCM$Y_tax_abr[grep(temp[i], parsed_eLSA_CLR_DCM$Y_tax)]="Ciliophora"
}

length(grep("kingdom_Eukaryota;supergroup_Alveolata;division_Apicomplexa", parsed_eLSA_CLR_DCM$Y_tax))
parsed_eLSA_CLR_DCM$Y_tax_abr[grep("kingdom_Eukaryota;supergroup_Alveolata;division_Apicomplexa", parsed_eLSA_CLR_DCM$Y_tax)]="Ampicomplexa"

#one more!
length(grep("kingdom_Eukaryota;supergroup_Apusozoa", parsed_eLSA_CLR_DCM$Y_tax))
parsed_eLSA_CLR_DCM$Y_tax_abr[grep("kingdom_Eukaryota;supergroup_Apusozoa", parsed_eLSA_CLR_DCM$Y_tax)]="Apusozoa"

#Another one more!
length(grep("kingdom_Eukaryota;supergroup_Amoebozoa", parsed_eLSA_CLR_DCM$Y_tax)) #Just one...
parsed_eLSA_CLR_DCM$Y_tax_abr[grep("kingdom_Eukaryota;supergroup_Amoebozoa", parsed_eLSA_CLR_DCM$Y_tax)]="Amebozoa"

#Check if we're all done 
length(grep("SAR", parsed_eLSA_CLR_DCM$Y_tax_abr)) #790
unique(parsed_eLSA_CLR_DCM$Y_tax[grep("SAR", parsed_eLSA_CLR_DCM$Y_tax_abr)])
temp=sort(unique(parsed_eLSA_CLR_DCM$Y_tax[grep("SAR", parsed_eLSA_CLR_DCM$Y_tax_abr)]))
length(grep("Ciliophora", temp)) #13
temp[grep("Ciliophora", temp, invert=T)]

#Ok great, add in UCYN-A
parsed_eLSA_CLR_DCM$Y_tax[grep("UCYN-A", parsed_eLSA_CLR_DCM$Y_tax)]
parsed_eLSA_CLR_DCM$Y_tax_abr[grep("UCYN-A", parsed_eLSA_CLR_DCM$Y_tax)]=parsed_eLSA_CLR_DCM$Y_tax[grep("UCYN-A", parsed_eLSA_CLR_DCM$Y_tax)]

#And the rest are "Eukaryote"
length(grep("SAR", parsed_eLSA_CLR_DCM$Y_tax_abr))
unique(parsed_eLSA_CLR_DCM$Y_tax_abr[grep("SAR", parsed_eLSA_CLR_DCM$Y_tax_abr)])
parsed_eLSA_CLR_DCM$Y_tax_abr[grep("SAR", parsed_eLSA_CLR_DCM$Y_tax_abr)]="Eukaryote"

unique(parsed_eLSA_CLR_DCM$Y_tax_abr)
dim(parsed_eLSA_CLR_DCM)

#WRITE IT OUT!
write.table(x=parsed_eLSA_CLR_DCM, file="eLSA_output/piecewise_eLSA/CLR_DCM/parsed_eLSA_output_CLR_DCM_11.25.2020.tsv", sep="\t", quote=F, row.names = F)

#3. Non-CLR transformed data from 5m-----

#Read in the data in a for-loop
setwd("../") #SPOT_16S_18S_eLSA folder
directory <- getwd()
filenames<-system(paste('ls ', directory, '/eLSA_output/piecewise_eLSA/int_5m', sep=''), intern=TRUE)

for(i in c(1:6)){
  a <- strsplit(filenames[i], split="_", fixed=T) 
  print(paste("eLSA_data_int", a[[1]][3], a[[1]][4], a[[1]][5], sep="_"))
  assign(x=paste("eLSA_data_int", a[[1]][3], a[[1]][4], a[[1]][5], sep="_"), value=read.table(file=paste("eLSA_output/piecewise_eLSA/int_5m/", filenames[i], sep=""), sep="\t", stringsAsFactors = F, header=T))
}

#Rbind em together! 
eLSA_int_5m_concat=rbind(eLSA_data_int_subs_A_B, eLSA_data_int_subs_A_C, eLSA_data_int_subs_A_int, eLSA_data_int_subs_B_C, eLSA_data_int_subs_B_int, eLSA_data_int_subs_C_int)
dim(eLSA_int_5m_concat) #167331 27 #Smaller file than non-CLR transformed 

#Calculate Q values----
library(qvalue)

#First, normal Q:
Q=qvalue(eLSA_int_5m_concat$P)
length(Q$qvalues)==length(eLSA_int_5m_concat$P)
eLSA_int_5m_concat$Q=Q$qvalues

#Then Q values for Pearsons's: Qpcc and Qspcc
pearsons_Q=qvalue(eLSA_int_5m_concat$Ppcc)
length(pearsons_Q$qvalues)==length(eLSA_int_5m_concat$Ppcc)
eLSA_int_5m_concat$Qpcc=pearsons_Q$qvalues

pearsons_Q_delay=qvalue(eLSA_int_5m_concat$Pspcc)
length(pearsons_Q_delay$qvalues)==length(eLSA_int_5m_concat$Pspcc)
eLSA_int_5m_concat$Qspcc=pearsons_Q_delay$qvalues

#Next Q values for Spearman's: Qscc and Qsscc
spearmans_Q=qvalue(eLSA_int_5m_concat$Pscc)
length(spearmans_Q$qvalues)==length(eLSA_int_5m_concat$Pscc)
eLSA_int_5m_concat$Qscc=spearmans_Q$qvalues

spearmans_Q_delay=qvalue(eLSA_int_5m_concat$Psscc)
length(spearmans_Q_delay$qvalues)==length(eLSA_int_5m_concat$Psscc)
eLSA_int_5m_concat$Qsscc=spearmans_Q_delay$qvalues

#Write this concatenated file out
write.table(x=eLSA_int_5m_concat, file="eLSA_output/piecewise_eLSA/int_5m/concat_eLSA_output_int_5m_12.08.2020.tsv", sep="\t", quote=F, row.names=F)

#Parse the data: P<0.005 and Q<0.01
nrow(eLSA_int_5m_concat) #167331
length(which(eLSA_int_5m_concat$P<0.005 & eLSA_int_5m_concat$Q<0.01)) #7635
length(which(eLSA_int_5m_concat$P<0.005 & eLSA_int_5m_concat$Q<0.01 & eLSA_int_5m_concat$Ppcc<0.005 & eLSA_int_5m_concat$Qpcc< 0.01)) #only 2252
length(which(eLSA_int_5m_concat$P<0.005 & eLSA_int_5m_concat$Q<0.01 & eLSA_int_5m_concat$Ppcc<0.005 & eLSA_int_5m_concat$Qpcc< 0.01 & eLSA_int_5m_concat$Pscc<0.005 & eLSA_int_5m_concat$Qscc < 0.01)) #2227

dim(eLSA_int_5m_concat) #167331     27
parsed_eLSA_int_5m=eLSA_int_5m_concat[which(eLSA_int_5m_concat$P<0.005 & eLSA_int_5m_concat$Q<0.01 & eLSA_int_5m_concat$Ppcc<0.005 & eLSA_int_5m_concat$Qpcc< 0.01 & eLSA_int_5m_concat$Pscc<0.005 & eLSA_int_5m_concat$Qscc < 0.01),]
dim(parsed_eLSA_int_5m) #2227 27

#Abbreviated X and Y hashes
length(unique(parsed_eLSA_int_5m$X)) #397
length(unique(abbreviate(parsed_eLSA_int_5m$X, minlength=4))) #397
parsed_eLSA_int_5m$X_abr=abbreviate(parsed_eLSA_int_5m$X, minlength=4)
head(unique(parsed_eLSA_int_5m$X_abr), n=15)

length(unique(parsed_eLSA_int_5m$Y)) #425
length(unique(abbreviate(parsed_eLSA_int_5m$Y))) #425
parsed_eLSA_int_5m$Y_abr=abbreviate(parsed_eLSA_int_5m$Y, minlength=4)
head(unique(parsed_eLSA_int_5m$Y_abr), n=15)

#Add in taxonomy and get UCYN-A indices 
#X:
parsed_eLSA_int_5m$X_taxonomy="SAR"
ix=c()
for(i in c(1:nrow(parsed_eLSA_int_5m))){
  if(length(as.numeric(grep(parsed_eLSA_int_5m$X[i], euks_tax$Feature.ID)))<1){
    ix=c(ix, i)
  }else{
    parsed_eLSA_int_5m$X_taxonomy[i]=euks_tax$Taxon[grep(parsed_eLSA_int_5m$X[i], euks_tax$Feature.ID)]    
  }
}
length(ix) #0
length(grep("_AE|_D", parsed_eLSA_int_5m$X)) #0 #That's why 
length(grep("SAR", parsed_eLSA_int_5m$X_taxonomy))

#Y: 
parsed_eLSA_int_5m$Y_taxonomy="SAR"
ix=c()
for(i in c(1:nrow(parsed_eLSA_int_5m))){
  if(length(as.numeric(grep(parsed_eLSA_int_5m$Y[i], euks_tax$Feature.ID)))<1){
    ix=c(ix, i)
  }else{
  parsed_eLSA_int_5m$Y_taxonomy[i]=euks_tax$Taxon[grep(parsed_eLSA_int_5m$Y[i], euks_tax$Feature.ID)]    
  }
}
length(ix)
length(grep("_AE|_D", parsed_eLSA_int_5m$Y)) #42
parsed_eLSA_int_5m$Y[ix]

#Write in UCYN-A indices
for(i in ix){
  a=strsplit(parsed_eLSA_int_5m$Y[i], split="_", fixed=T)
  parsed_eLSA_int_5m$Y_taxonomy[i]=paste("UCYN-A", grep(a[[1]][1], UCYNA_ix), a[[1]][2], sep="_")
}
length(grep("UCYN", parsed_eLSA_int_5m$Y_taxonomy)) #42

#Abbreviate X taxonomy-----
#Start by labeling everything "SAR"
parsed_eLSA_int_5m$X_tax_abr="SAR"

#Manually change taxa, oh god :( #At the "supergroup" level 

##Start with the Stramenopiles: All stramenopiles except Pseudo-nitzschia and Chaetoceros set to "Stramenopiles"

length(grep("kingdom_Eukaryota;supergroup_Stramenopiles", parsed_eLSA_int_5m$X_taxonomy)) #287 #That's like 10% of the X taxa
parsed_eLSA_int_5m$X_tax_abr[grep("kingdom_Eukaryota;supergroup_Stramenopiles;", parsed_eLSA_int_5m$X_taxonomy)]="Stramenopile"
parsed_eLSA_int_5m$X_tax_abr[grep(temp[grep("Stramenopile", temp)], parsed_eLSA_int_5m$X_taxonomy)]="Stramenopile"

length(grep("kingdom_Eukaryota;supergroup_Stramenopiles;division_Ochrophyta;class_Bacillariophyta;order_Bacillariophyta_X;family_Raphid-pennate;genus_Pseudo-nitzschia", parsed_eLSA_int_5m$X_taxonomy))
parsed_eLSA_int_5m$X_tax_abr[grep("kingdom_Eukaryota;supergroup_Stramenopiles;division_Ochrophyta;class_Bacillariophyta;order_Bacillariophyta_X;family_Raphid-pennate;genus_Pseudo-nitzschia", parsed_eLSA_int_5m$X_taxonomy)]="Pseudo-nitzchia"

length(grep("kingdom_Eukaryota;supergroup_Stramenopiles;division_Ochrophyta;class_Bacillariophyta;order_Bacillariophyta_X;family_Polar-centric-Mediophyceae;genus_Chaetoceros", parsed_eLSA_int_5m$X_taxonomy))
parsed_eLSA_int_5m$X_tax_abr[grep("kingdom_Eukaryota;supergroup_Stramenopiles;division_Ochrophyta;class_Bacillariophyta;order_Bacillariophyta_X;family_Polar-centric-Mediophyceae;genus_Chaetoceros", parsed_eLSA_int_5m$X_taxonomy)]="Chaetoceros"

length(grep("kingdom_Eukaryota;supergroup_Rhizaria", parsed_eLSA_int_5m$X_taxonomy))
parsed_eLSA_int_5m$X_tax_abr[grep("Rhizaria", parsed_eLSA_int_5m$X_taxonomy)]="Rhizaria"

length(grep("kingdom_Eukaryota;supergroup_Opisthokonta", parsed_eLSA_int_5m$X_taxonomy))
parsed_eLSA_int_5m$X_tax_abr[grep("kingdom_Eukaryota;supergroup_Opisthokonta", parsed_eLSA_int_5m$X_taxonomy)]="Opisthokonta"

##All hacrobia except prymnesiophytes, Haptophytes, and Phaeocystis set to "Hacrobia" 

length(grep("kingdom_Eukaryota;supergroup_Hacrobia", parsed_eLSA_int_5m$X_taxonomy))
parsed_eLSA_int_5m$X_tax_abr[grep("kingdom_Eukaryota;supergroup_Hacrobia", parsed_eLSA_int_5m$X_taxonomy)]="Hacrobia"

length(grep("kingdom_Eukaryota;supergroup_Hacrobia;division_Haptophyta;class_Prymnesiophyceae;order_Prymnesiales;family_Prymnesiaceae;genus_Prymnesium;species_Prymnesium_sp.", parsed_eLSA_int_5m$X_taxonomy))
parsed_eLSA_int_5m$X_tax_abr[grep("kingdom_Eukaryota;supergroup_Hacrobia;division_Haptophyta;class_Prymnesiophyceae;order_Prymnesiales;family_Prymnesiaceae;genus_Prymnesium;species_Prymnesium_sp.", parsed_eLSA_int_5m$X_taxonomy)]="Prymnesiophyte"

length(grep("kingdom_Eukaryota;supergroup_Hacrobia;division_Haptophyta;class_Prymnesiophyceae;order_Prymnesiales;family_Chrysochromulinaceae;genus_Chrysochromulina", parsed_eLSA_int_5m$X_taxonomy)) #62
parsed_eLSA_int_5m$X_tax_abr[grep("kingdom_Eukaryota;supergroup_Hacrobia;division_Haptophyta;class_Prymnesiophyceae;order_Prymnesiales;family_Chrysochromulinaceae;genus_Chrysochromulina", parsed_eLSA_int_5m$X_taxonomy)]="Chrysochromulina"

#Phaeocystis
length(grep("kingdom_Eukaryota;supergroup_Hacrobia;division_Haptophyta;class_Prymnesiophyceae;order_Phaeocystales;family_Phaeocystaceae;genus_Phaeocystis", parsed_eLSA_int_5m$X_taxonomy))
parsed_eLSA_int_5m$X_tax_abr[grep("kingdom_Eukaryota;supergroup_Hacrobia;division_Haptophyta;class_Prymnesiophyceae;order_Phaeocystales;family_Phaeocystaceae;genus_Phaeocystis", parsed_eLSA_int_5m$X_taxonomy)]="Phaeocystis"

#Haptophyte
length(grep("kingdom_Eukaryota;supergroup_Hacrobia;division_Haptophyta;class_Haptophyta", parsed_eLSA_int_5m$X_taxonomy))
parsed_eLSA_int_5m$X_tax_abr[grep("kingdom_Eukaryota;supergroup_Hacrobia;division_Haptophyta;class_Haptophyta", parsed_eLSA_int_5m$X_taxonomy)]="Haptophyte"

#Insert Brad ix's
length(grep("kingdom_Eukaryota;supergroup_Hacrobia;division_Haptophyta;class_Prymnesiophyceae;order_Prymnesiophyceae_X;family_Braarudosphaeraceae;genus_Braarudosphaeraceae_X", parsed_eLSA_int_5m$X_taxonomy)) #24
parsed_eLSA_int_5m$X[grep("kingdom_Eukaryota;supergroup_Hacrobia;division_Haptophyta;class_Prymnesiophyceae;order_Prymnesiophyceae_X;family_Braarudosphaeraceae;genus_Braarudosphaeraceae_X", parsed_eLSA_int_5m$X_taxonomy)] #They are all #7

for(i in Brad_ix){
  if(length(grep(i, parsed_eLSA_int_5m$X))>1){
    print(i)
    print(grep(i, parsed_eLSA_int_5m$X)) 
  }
} #Hmm ok, only 3 Brad ASVs showed up: #3 (formerly 4), 4 (formerly 5), 7x (formerly 7) 

parsed_eLSA_int_5m$X_tax_abr[grep(Brad_ix[4], parsed_eLSA_int_5m$X)]="Braarudospharea_bigelowii_3"
parsed_eLSA_int_5m$X_tax_abr[grep(Brad_ix[7], parsed_eLSA_int_5m$X)]="Braarudosphaera_sp_7x"
parsed_eLSA_int_5m$X_tax_abr[grep(Brad_ix[5], parsed_eLSA_int_5m$X)]="Braarudospharea_bigelowii_4"

#Moving on!
length(grep("kingdom_Eukaryota;supergroup_Archaeplastida;division_Rhodophyta", parsed_eLSA_int_5m$X_taxonomy))
parsed_eLSA_int_5m$X_tax_abr[grep("kingdom_Eukaryota;supergroup_Archaeplastida;division_Rhodophyta", parsed_eLSA_int_5m$X_taxonomy)]="Archeplastida"

length(grep("kingdom_Eukaryota;supergroup_Archaeplastida;division_Chlorophyta", parsed_eLSA_int_5m$X_taxonomy))
parsed_eLSA_int_5m$X_tax_abr[grep("kingdom_Eukaryota;supergroup_Archaeplastida;division_Chlorophyta", parsed_eLSA_int_5m$X_taxonomy)]="Chlorophyta"

length(grep("kingdom_Eukaryota;supergroup_Archaeplastida;division_Chlorophyta;class_Mamiellophyceae;order_Mamiellales;family_Bathycoccaceae", parsed_eLSA_int_5m$X_taxonomy))
parsed_eLSA_int_5m$X_tax_abr[grep("kingdom_Eukaryota;supergroup_Archaeplastida;division_Chlorophyta;class_Mamiellophyceae;order_Mamiellales;family_Bathycoccaceae", parsed_eLSA_int_5m$X_taxonomy)]="Bathycoccus"

#Dinoflagellates
parsed_eLSA_int_5m$X_tax_abr[grep(temp[grep("Alveolata", temp)][1], parsed_eLSA_int_5m$X_taxonomy)]="Alveolata"

length(grep("kingdom_Eukaryota;supergroup_Alveolata;division_Dinoflagellata", parsed_eLSA_int_5m$X_taxonomy)) #1089 #WOWZA! 
parsed_eLSA_int_5m$X_tax_abr[grep("kingdom_Eukaryota;supergroup_Alveolata;division_Dinoflagellata", parsed_eLSA_int_5m$X_taxonomy)]="Dinoflagellate"

length(grep("kingdom_Eukaryota;supergroup_Alveolata;division_Dinoflagellata;class_Syndiniales", parsed_eLSA_int_5m$X_taxonomy)) #629
parsed_eLSA_int_5m$X_tax_abr[grep("kingdom_Eukaryota;supergroup_Alveolata;division_Dinoflagellata;class_Syndiniales", parsed_eLSA_int_5m$X_taxonomy)]="Syndiniales"

length(grep("kingdom_Eukaryota;supergroup_Alveolata;division_Dinoflagellata;class_Dinophyta", parsed_eLSA_int_5m$X_taxonomy)) #9
parsed_eLSA_int_5m$X_tax_abr[grep("kingdom_Eukaryota;supergroup_Alveolata;division_Dinoflagellata;class_Dinophyta", parsed_eLSA_int_5m$X_taxonomy)]="Dinophyta"

length(grep("kingdom_Eukaryota;supergroup_Alveolata;division_Dinoflagellata;class_Dinophyceae", parsed_eLSA_int_5m$X_taxonomy)) #438
parsed_eLSA_int_5m$X_tax_abr[grep("kingdom_Eukaryota;supergroup_Alveolata;division_Dinoflagellata;class_Dinophyceae", parsed_eLSA_int_5m$X_taxonomy)]="Dinophyceae"

length(grep("kingdom_Eukaryota;supergroup_Alveolata;division_Dinoflagellata;class_Dinophyceae;order_Gymnodiniales;family_Gymnodiniaceae;genus_Lepidodinium", parsed_eLSA_int_5m$X_taxonomy)) #14
parsed_eLSA_int_5m$X_tax_abr[grep("kingdom_Eukaryota;supergroup_Alveolata;division_Dinoflagellata;class_Dinophyceae;order_Gymnodiniales;family_Gymnodiniaceae;genus_Lepidodinium", parsed_eLSA_int_5m$X_taxonomy)]="Lepidodinium"

length(grep("kingdom_Eukaryota;supergroup_Alveolata;division_Dinoflagellata;class_Dinophyceae;order_Gymnodiniales", parsed_eLSA_int_5m$X_taxonomy)) #135
parsed_eLSA_int_5m$X_tax_abr[grep("kingdom_Eukaryota;supergroup_Alveolata;division_Dinoflagellata;class_Dinophyceae;order_Gymnodiniales", parsed_eLSA_int_5m$X_taxonomy)]="Gymnodinium"

length(grep("kingdom_Eukaryota;supergroup_Alveolata;division_Ciliophora;class_Spirotrichea;order_Tintinnida", parsed_eLSA_int_5m$X_taxonomy))
parsed_eLSA_int_5m$X_tax_abr[grep("kingdom_Eukaryota;supergroup_Alveolata;division_Ciliophora;class_Spirotrichea;order_Tintinnida", parsed_eLSA_int_5m$X_taxonomy)]="Tintinnida"

length(grep(" kingdom_Eukaryota;supergroup_Alveolata;division_Ciliophora", parsed_eLSA_int_5m$X_taxonomy)) #260
parsed_eLSA_int_5m$X_tax_abr[grep(" kingdom_Eukaryota;supergroup_Alveolata;division_Ciliophora;class_Spirotrichea;order_Strombidiida", parsed_eLSA_int_5m$X_taxonomy)]="Ciliophora"

for(i in grep("Ciliophora", temp)){
  parsed_eLSA_int_5m$X_tax_abr[grep(temp[i], parsed_eLSA_int_5m$X_taxonomy)]="Ciliophora"
}

length(grep("kingdom_Eukaryota;supergroup_Alveolata;division_Apicomplexa", parsed_eLSA_int_5m$X_taxonomy)) #0
#parsed_eLSA_int_5m$X_tax_abr[grep("kingdom_Eukaryota;supergroup_Alveolata;division_Apicomplexa", parsed_eLSA_int_5m$X_tax)]="Ampicomplexa"


#one more!
parsed_eLSA_int_5m$X_tax_abr[grep(temp[grep("Amoebozoa", temp)], parsed_eLSA_int_5m$X_taxonomy)]="Amoebozoa"
length(grep("kingdom_Eukaryota;supergroup_Amoebozoa", parsed_eLSA_int_5m$X_taxonomy))
parsed_eLSA_int_5m$X_tax_abr[grep("kingdom_Eukaryota;supergroup_Apusozoa", parsed_eLSA_int_5m$X_taxonomy)]="Amoebozoa"

#OKAY DOUBLE CHECK IT 
length(grep("SAR", parsed_eLSA_int_5m$X_tax_abr)) #181 #64  
temp=sort(unique(parsed_eLSA_int_5m$X_taxonomy[grep("SAR", parsed_eLSA_int_5m$X_tax_abr)]))
length(temp) #15
length(grep("Ciliophora", temp)) #11
temp[grep("Ciliophora", temp, invert=T)]

#FINALLY, all the remainder are just "Eukaryotes" 
parsed_eLSA_int_5m$X_tax_abr[grep("SAR", parsed_eLSA_int_5m$X_tax_abr)]="Eukaryote"
unique(parsed_eLSA_int_5m$X_tax_abr) #25 long 

#parsed_eLSA_int_5m$X_tax_abr[grep("UCYN-A", parsed_eLSA_int_5m$X_tax)]="UCYN-A_1_AE"

#Abbreviate Y taxonomy-----
#Start by labelling everything "SAR"
parsed_eLSA_int_5m$Y_tax_abr="SAR"

#Manually change taxa, oh god :( #At the "supergroup" level 

##Start with the Stramenopiles: All stramenopiles except Pseudo-nitzschia and Chaetoceros set to "Stramenopiles"

length(grep("kingdom_Eukaryota;supergroup_Stramenopiles", parsed_eLSA_int_5m$Y_taxonomy)) #464 #That's like 10% of the X taxa
parsed_eLSA_int_5m$Y_tax_abr[grep("kingdom_Eukaryota;supergroup_Stramenopiles;", parsed_eLSA_int_5m$Y_taxonomy)]="Stramenopile"
parsed_eLSA_int_5m$Y_tax_abr[grep(temp[grep("Stramenopile", temp)], parsed_eLSA_int_5m$Y_taxonomy)]="Stramenopile"

length(grep("kingdom_Eukaryota;supergroup_Stramenopiles;division_Ochrophyta;class_Bacillariophyta;order_Bacillariophyta_X;family_Raphid-pennate;genus_Pseudo-nitzschia", parsed_eLSA_int_5m$Y_taxonomy))
parsed_eLSA_int_5m$Y_tax_abr[grep("kingdom_Eukaryota;supergroup_Stramenopiles;division_Ochrophyta;class_Bacillariophyta;order_Bacillariophyta_X;family_Raphid-pennate;genus_Pseudo-nitzschia", parsed_eLSA_int_5m$Y_taxonomy)]="Pseudo-nitzchia"

length(grep("kingdom_Eukaryota;supergroup_Stramenopiles;division_Ochrophyta;class_Bacillariophyta;order_Bacillariophyta_X;family_Polar-centric-Mediophyceae;genus_Chaetoceros", parsed_eLSA_int_5m$Y_taxonomy))
parsed_eLSA_int_5m$Y_tax_abr[grep("kingdom_Eukaryota;supergroup_Stramenopiles;division_Ochrophyta;class_Bacillariophyta;order_Bacillariophyta_X;family_Polar-centric-Mediophyceae;genus_Chaetoceros", parsed_eLSA_int_5m$Y_taxonomy)]="Chaetoceros"

length(grep("kingdom_Eukaryota;supergroup_Rhizaria", parsed_eLSA_int_5m$Y_taxonomy))
parsed_eLSA_int_5m$Y_tax_abr[grep("Rhizaria", parsed_eLSA_int_5m$Y_taxonomy)]="Rhizaria"

length(grep("kingdom_Eukaryota;supergroup_Opisthokonta", parsed_eLSA_int_5m$Y_taxonomy))
parsed_eLSA_int_5m$Y_tax_abr[grep("kingdom_Eukaryota;supergroup_Opisthokonta", parsed_eLSA_int_5m$Y_taxonomy)]="Opisthokonta"

##All hacrobia except prymnesiophytes, Haptophytes, and Phaeocystis set to "Hacrobia" 

length(grep("kingdom_Eukaryota;supergroup_Hacrobia", parsed_eLSA_int_5m$Y_taxonomy))
parsed_eLSA_int_5m$Y_tax_abr[grep("kingdom_Eukaryota;supergroup_Hacrobia", parsed_eLSA_int_5m$Y_taxonomy)]="Hacrobia"

length(grep("kingdom_Eukaryota;supergroup_Hacrobia;division_Haptophyta;class_Prymnesiophyceae;order_Prymnesiales;family_Prymnesiaceae;genus_Prymnesium;species_Prymnesium_sp.", parsed_eLSA_int_5m$Y_taxonomy))
parsed_eLSA_int_5m$Y_tax_abr[grep("kingdom_Eukaryota;supergroup_Hacrobia;division_Haptophyta;class_Prymnesiophyceae;order_Prymnesiales;family_Prymnesiaceae;genus_Prymnesium;species_Prymnesium_sp.", parsed_eLSA_int_5m$Y_taxonomy)]="Prymnesiophyte"

length(grep("kingdom_Eukaryota;supergroup_Hacrobia;division_Haptophyta;class_Prymnesiophyceae;order_Prymnesiales;family_Chrysochromulinaceae;genus_Chrysochromulina", parsed_eLSA_int_5m$Y_taxonomy)) #62
parsed_eLSA_int_5m$Y_tax_abr[grep("kingdom_Eukaryota;supergroup_Hacrobia;division_Haptophyta;class_Prymnesiophyceae;order_Prymnesiales;family_Chrysochromulinaceae;genus_Chrysochromulina", parsed_eLSA_int_5m$Y_taxonomy)]="Chrysochromulina"

#Phaeocystis
length(grep("kingdom_Eukaryota;supergroup_Hacrobia;division_Haptophyta;class_Prymnesiophyceae;order_Phaeocystales;family_Phaeocystaceae;genus_Phaeocystis", parsed_eLSA_int_5m$Y_taxonomy))
parsed_eLSA_int_5m$Y_tax_abr[grep("kingdom_Eukaryota;supergroup_Hacrobia;division_Haptophyta;class_Prymnesiophyceae;order_Phaeocystales;family_Phaeocystaceae;genus_Phaeocystis", parsed_eLSA_int_5m$Y_taxonomy)]="Phaeocystis"

#Haptophyte
length(grep("kingdom_Eukaryota;supergroup_Hacrobia;division_Haptophyta;class_Haptophyta", parsed_eLSA_int_5m$Y_taxonomy))
parsed_eLSA_int_5m$Y_tax_abr[grep("kingdom_Eukaryota;supergroup_Hacrobia;division_Haptophyta;class_Haptophyta", parsed_eLSA_int_5m$Y_taxonomy)]="Haptophyte"

#Insert Brad ix's
length(grep("kingdom_Eukaryota;supergroup_Hacrobia;division_Haptophyta;class_Prymnesiophyceae;order_Prymnesiophyceae_X;family_Braarudosphaeraceae;genus_Braarudosphaeraceae_X", parsed_eLSA_int_5m$Y_taxonomy)) #0
#parsed_eLSA_int_5m$Y_taxonomy[grep("kingdom_Eukaryota;supergroup_Hacrobia;division_Haptophyta;class_Prymnesiophyceae;order_Prymnesiophyceae_X;family_Braarudosphaeraceae", parsed_eLSA_int_5m$Y_taxonomy)]
#Nada

#parsed_eLSA_int_5m$X[grep("kingdom_Eukaryota;supergroup_Hacrobia;division_Haptophyta;class_Prymnesiophyceae;order_Prymnesiophyceae_X;family_Braarudosphaeraceae;genus_Braarudosphaeraceae_X", parsed_eLSA_int_5m$Y_taxonomy)] #They are all #7

#for(i in Brad_ix){
  if(length(grep(i, parsed_eLSA_int_5m$X))>1){
    print(i)
    print(grep(i, parsed_eLSA_int_5m$X)) 
  }
} #Hmm ok, only 3 Brad ASVs showed up: #3 (formerly 4), 4 (formerly 5), 7x (formerly 7) 

#parsed_eLSA_int_5m$Y_tax_abr[grep(Brad_ix[4], parsed_eLSA_int_5m$Y)]="Braarudospharea_bigelowii_3"
#parsed_eLSA_int_5m$Y_tax_abr[grep(Brad_ix[7], parsed_eLSA_int_5m$Y)]="Braarudosphaera_sp_7x"
#parsed_eLSA_int_5m$Y_tax_abr[grep(Brad_ix[5], parsed_eLSA_int_5m$Y)]="Braarudospharea_bigelowii_4"

#Moving on!
length(grep("kingdom_Eukaryota;supergroup_Archaeplastida;division_Rhodophyta", parsed_eLSA_int_5m$Y_taxonomy))
parsed_eLSA_int_5m$Y_tax_abr[grep("kingdom_Eukaryota;supergroup_Archaeplastida;division_Rhodophyta", parsed_eLSA_int_5m$Y_taxonomy)]="Archeplastida"

length(grep("kingdom_Eukaryota;supergroup_Archaeplastida;division_Chlorophyta", parsed_eLSA_int_5m$Y_taxonomy))
parsed_eLSA_int_5m$Y_tax_abr[grep("kingdom_Eukaryota;supergroup_Archaeplastida;division_Chlorophyta", parsed_eLSA_int_5m$Y_taxonomy)]="Chlorophyta"

length(grep("kingdom_Eukaryota;supergroup_Archaeplastida;division_Chlorophyta;class_Mamiellophyceae;order_Mamiellales;family_Bathycoccaceae", parsed_eLSA_int_5m$Y_taxonomy))
parsed_eLSA_int_5m$Y_tax_abr[grep("kingdom_Eukaryota;supergroup_Archaeplastida;division_Chlorophyta;class_Mamiellophyceae;order_Mamiellales;family_Bathycoccaceae", parsed_eLSA_int_5m$Y_taxonomy)]="Bathycoccus"

#Dinoflagellates
parsed_eLSA_int_5m$Y_tax_abr[grep(temp[grep("Alveolata", temp)][1], parsed_eLSA_int_5m$Y_taxonomy)]="Alveolata"

length(grep("kingdom_Eukaryota;supergroup_Alveolata;division_Dinoflagellata", parsed_eLSA_int_5m$Y_taxonomy)) #886 #WOWZA! 
parsed_eLSA_int_5m$Y_tax_abr[grep("kingdom_Eukaryota;supergroup_Alveolata;division_Dinoflagellata", parsed_eLSA_int_5m$Y_taxonomy)]="Dinoflagellate"

length(grep("kingdom_Eukaryota;supergroup_Alveolata;division_Dinoflagellata;class_Syndiniales", parsed_eLSA_int_5m$Y_taxonomy)) #496
parsed_eLSA_int_5m$Y_tax_abr[grep("kingdom_Eukaryota;supergroup_Alveolata;division_Dinoflagellata;class_Syndiniales", parsed_eLSA_int_5m$Y_taxonomy)]="Syndiniales"

length(grep("kingdom_Eukaryota;supergroup_Alveolata;division_Dinoflagellata;class_Dinophyta", parsed_eLSA_int_5m$Y_taxonomy)) #20
parsed_eLSA_int_5m$Y_tax_abr[grep("kingdom_Eukaryota;supergroup_Alveolata;division_Dinoflagellata;class_Dinophyta", parsed_eLSA_int_5m$Y_taxonomy)]="Dinophyta"

length(grep("kingdom_Eukaryota;supergroup_Alveolata;division_Dinoflagellata;class_Dinophyceae", parsed_eLSA_int_5m$Y_taxonomy)) #358
parsed_eLSA_int_5m$Y_tax_abr[grep("kingdom_Eukaryota;supergroup_Alveolata;division_Dinoflagellata;class_Dinophyceae", parsed_eLSA_int_5m$Y_taxonomy)]="Dinophyceae"

length(grep("kingdom_Eukaryota;supergroup_Alveolata;division_Dinoflagellata;class_Dinophyceae;order_Gymnodiniales;family_Gymnodiniaceae;genus_Lepidodinium", parsed_eLSA_int_5m$Y_taxonomy)) #5
parsed_eLSA_int_5m$Y_tax_abr[grep("kingdom_Eukaryota;supergroup_Alveolata;division_Dinoflagellata;class_Dinophyceae;order_Gymnodiniales;family_Gymnodiniaceae;genus_Lepidodinium", parsed_eLSA_int_5m$Y_taxonomy)]="Lepidodinium"

length(grep("kingdom_Eukaryota;supergroup_Alveolata;division_Dinoflagellata;class_Dinophyceae;order_Gymnodiniales", parsed_eLSA_int_5m$Y_taxonomy)) #76
parsed_eLSA_int_5m$Y_tax_abr[grep("kingdom_Eukaryota;supergroup_Alveolata;division_Dinoflagellata;class_Dinophyceae;order_Gymnodiniales", parsed_eLSA_int_5m$Y_taxonomy)]="Gymnodinium"

length(grep("kingdom_Eukaryota;supergroup_Alveolata;division_Ciliophora;class_Spirotrichea;order_Tintinnida", parsed_eLSA_int_5m$Y_taxonomy))
parsed_eLSA_int_5m$Y_tax_abr[grep("kingdom_Eukaryota;supergroup_Alveolata;division_Ciliophora;class_Spirotrichea;order_Tintinnida", parsed_eLSA_int_5m$Y_taxonomy)]="Tintinnida"

length(grep(" kingdom_Eukaryota;supergroup_Alveolata;division_Ciliophora", parsed_eLSA_int_5m$Y_taxonomy)) #176
parsed_eLSA_int_5m$Y_tax_abr[grep(" kingdom_Eukaryota;supergroup_Alveolata;division_Ciliophora;class_Spirotrichea;order_Strombidiida", parsed_eLSA_int_5m$Y_taxonomy)]="Ciliophora"

for(i in grep("Ciliophora", temp)){
  parsed_eLSA_int_5m$Y_tax_abr[grep(temp[i], parsed_eLSA_int_5m$Y_taxonomy)]="Ciliophora"
}

length(grep("kingdom_Eukaryota;supergroup_Alveolata;division_Apicomplexa", parsed_eLSA_int_5m$Y_taxonomy)) #21
parsed_eLSA_int_5m$Y_tax_abr[grep("kingdom_Eukaryota;supergroup_Alveolata;division_Apicomplexa", parsed_eLSA_int_5m$Y_tax)]="Ampicomplexa"


#one more!
parsed_eLSA_int_5m$Y_tax_abr[grep(temp[grep("Amoebozoa", temp)], parsed_eLSA_int_5m$Y_taxonomy)]="Amoebozoa"
length(grep("kingdom_Eukaryota;supergroup_Amoebozoa", parsed_eLSA_int_5m$Y_taxonomy))
parsed_eLSA_int_5m$Y_tax_abr[grep("kingdom_Eukaryota;supergroup_Apusozoa", parsed_eLSA_int_5m$Y_taxonomy)]="Amoebozoa"

#OKAY DOUBLE CHECK IT 
length(grep("SAR", parsed_eLSA_int_5m$Y_tax_abr)) #198 #203
temp=sort(unique(parsed_eLSA_int_5m$Y_taxonomy[grep("SAR", parsed_eLSA_int_5m$Y_tax_abr)]))
length(temp) #20
length(grep("Ciliophora", temp)) #11
temp[grep("Ciliophora", temp, invert=T)]

#Write in UCYN-A
parsed_eLSA_int_5m$Y_tax_abr[grep("UCYN-A", parsed_eLSA_int_5m$Y_taxonomy)]=parsed_eLSA_int_5m$Y_taxonomy[grep("UCYN-A", parsed_eLSA_int_5m$Y_taxonomy)]

#FINALLY, all the remainder are just "Eukaryotes" 
parsed_eLSA_int_5m$Y_taxonomy[grep("SAR", parsed_eLSA_int_5m$Y_tax_abr)]
parsed_eLSA_int_5m$Y_tax_abr[grep("SAR", parsed_eLSA_int_5m$Y_tax_abr)]="Eukaryote"
unique(parsed_eLSA_int_5m$Y_tax_abr)

#WRITE THAT SHIT OUT!!! 
dim(parsed_eLSA_int_5m) #2227 33
write.table(x=parsed_eLSA_int_5m, file="eLSA_output/piecewise_eLSA/int_5m/parsed_eLSA_output_int_5m_12.08.2020.tsv", sep="\t", quote=F, row.names=F)

#4. Non-CLR transformed data from DCM----
directory <- getwd()
filenames<-system(paste('ls ', directory, '/eLSA_output/piecewise_eLSA/int_DCM', sep=''), intern=TRUE)
filenames

#Read in the data in a for loop
for(i in c(1:6)){
  a <- strsplit(filenames[i], split="_", fixed=T) 
  print(paste("eLSA_data_int", a[[1]][3], a[[1]][4], a[[1]][5], sep="_"))
  assign(x=paste("eLSA_data_int", a[[1]][3], a[[1]][4], a[[1]][5], sep="_"), value=read.table(file=paste("eLSA_output/piecewise_eLSA/int_DCM/", filenames[i], sep=""), sep="\t", stringsAsFactors = F, header=T))
}

#Rbind em together! 
eLSA_int_DCM_concat=rbind(eLSA_data_int_subs_A_B, eLSA_data_int_subs_A_C, eLSA_data_int_subs_A_int, eLSA_data_int_subs_B_C, eLSA_data_int_subs_B_int, eLSA_data_int_subs_C_int)
dim(eLSA_int_DCM_concat) #268278     27

#Calculate q values------
library(qvalue)

#First, normal Q:
Q=qvalue(eLSA_int_DCM_concat$P)
length(Q$qvalues)==length(eLSA_int_DCM_concat$P)
eLSA_int_DCM_concat$Q=Q$qvalues

#Then Q values for Pearsons's: Qpcc and Qspcc
pearsons_Q=qvalue(eLSA_int_DCM_concat$Ppcc)
length(pearsons_Q$qvalues)==length(eLSA_int_DCM_concat$Ppcc)
eLSA_int_DCM_concat$Qpcc=pearsons_Q$qvalues

pearsons_Q_delay=qvalue(eLSA_int_DCM_concat$Pspcc)
length(pearsons_Q_delay$qvalues)==length(eLSA_int_DCM_concat$Pspcc)
eLSA_int_DCM_concat$Qspcc=pearsons_Q_delay$qvalues

#Next Q values for Spearman's: Qscc and Qsscc
spearmans_Q=qvalue(eLSA_int_DCM_concat$Pscc)
length(spearmans_Q$qvalues)==length(eLSA_int_DCM_concat$Pscc)
eLSA_int_DCM_concat$Qscc=spearmans_Q$qvalues

spearmans_Q_delay=qvalue(eLSA_int_DCM_concat$Psscc)
length(spearmans_Q_delay$qvalues)==length(eLSA_int_DCM_concat$Psscc)
eLSA_int_DCM_concat$Qsscc=spearmans_Q_delay$qvalues

#Write this table out!
write.table(x=eLSA_int_DCM_concat, file="eLSA_output/piecewise_eLSA/int_DCM/concat_eLSA_output_int_DCM_12.18.2020.tsv", sep="\t", quote=F, row.names=F)

#Parse the data: all P<0.005 and Q< 0.01
length(which(eLSA_int_DCM_concat$P<0.005 & eLSA_int_DCM_concat$Q<0.01)) #9326
length(which(eLSA_int_DCM_concat$P<0.005 & eLSA_int_DCM_concat$Q<0.01 & eLSA_int_DCM_concat$Ppcc<0.005 & eLSA_int_DCM_concat$Qpcc<0.01)) #3507
length(which(eLSA_int_DCM_concat$P<0.005 & eLSA_int_DCM_concat$Q<0.01 & eLSA_int_DCM_concat$Ppcc<0.005 & eLSA_int_DCM_concat$Qpcc<0.01 & eLSA_int_DCM_concat$Pscc<0.005 & eLSA_int_DCM_concat$Qscc<0.01)) #3472

parsed_eLSA_int_DCM=eLSA_int_DCM_concat[which(eLSA_int_DCM_concat$P<0.005 & eLSA_int_DCM_concat$Q<0.01 & eLSA_int_DCM_concat$Ppcc<0.005 & eLSA_int_DCM_concat$Qpcc<0.01 & eLSA_int_DCM_concat$Pscc<0.005 & eLSA_int_DCM_concat$Qscc<0.01),]
dim(parsed_eLSA_int_DCM) #3472 27 

#Add in X euks tax and find indices for UCYN-A
ix <- c()
for(i in c(1:nrow(parsed_eLSA_int_DCM))){
  a <- parsed_eLSA_int_DCM$X[i]
  b <- grep(a, euks_tax$Feature.ID)
  if(as.numeric(length(i)<1)){
    print(b)
    ix <- c(ix, i)
  }else{
    parsed_eLSA_int_DCM$X_tax[i]=euks_tax$Taxon[b]
  }
}
ix #empty
head(parsed_eLSA_int_DCM$X_tax)

parsed_eLSA_int_DCM$Y[grep("_AE|_D", parsed_eLSA_int_DCM$Y)] #Crap, UCYN-A only shows up once at the DCM apparently

#Y taxonomy 
ix <- c()
for(i in c(1:nrow(parsed_eLSA_int_DCM))){
  a <- parsed_eLSA_int_DCM$Y[i]
  b <- grep(a, euks_tax$Feature.ID)
  if(as.numeric(length(b)<1)){
    print(b)
    ix <- c(ix, i)
  }else{
    parsed_eLSA_int_DCM$Y_tax[i]=euks_tax$Taxon[b]
  }
}
ix
grep("_AE|_D", parsed_eLSA_int_DCM$Y)
parsed_eLSA_int_DCM$Y_tax[ix]="UCYNA-1_AE"

#Abbreviate X taxonomy-----
#Start by labeling everything "SAR"
parsed_eLSA_int_DCM$X_tax_abr="SAR"

#Manually change taxa, oh god :( #At the "supergroup" level 

##Start with the Stramenopiles: All stramenopiles except Pseudo-nitzschia and Chaetoceros set to "Stramenopiles"

length(grep("kingdom_Eukaryota;supergroup_Stramenopiles", parsed_eLSA_int_DCM$X_tax)) #424 #That's like 10% of the X taxa
parsed_eLSA_int_DCM$X_tax_abr[grep("kingdom_Eukaryota;supergroup_Stramenopiles;", parsed_eLSA_int_DCM$X_tax)]="Stramenopile"
parsed_eLSA_int_DCM$X_tax_abr[grep(temp[grep("Stramenopile", temp)], parsed_eLSA_int_DCM$X_tax)]="Stramenopile"

length(grep("kingdom_Eukaryota;supergroup_Stramenopiles;division_Ochrophyta;class_Bacillariophyta;order_Bacillariophyta_X;family_Raphid-pennate;genus_Pseudo-nitzschia", parsed_eLSA_int_DCM$X_tax)) #28
parsed_eLSA_int_DCM$X_tax_abr[grep("kingdom_Eukaryota;supergroup_Stramenopiles;division_Ochrophyta;class_Bacillariophyta;order_Bacillariophyta_X;family_Raphid-pennate;genus_Pseudo-nitzschia", parsed_eLSA_int_DCM$X_taxonomy)]="Pseudo-nitzchia"

length(grep("kingdom_Eukaryota;supergroup_Stramenopiles;division_Ochrophyta;class_Bacillariophyta;order_Bacillariophyta_X;family_Polar-centric-Mediophyceae;genus_Chaetoceros", parsed_eLSA_int_DCM$X_tax))
parsed_eLSA_int_DCM$X_tax_abr[grep("kingdom_Eukaryota;supergroup_Stramenopiles;division_Ochrophyta;class_Bacillariophyta;order_Bacillariophyta_X;family_Polar-centric-Mediophyceae;genus_Chaetoceros", parsed_eLSA_int_DCM$X_tax)]="Chaetoceros"

length(grep("kingdom_Eukaryota;supergroup_Rhizaria", parsed_eLSA_int_DCM$X_tax))
parsed_eLSA_int_DCM$X_tax_abr[grep("Rhizaria", parsed_eLSA_int_DCM$X_tax)]="Rhizaria"

length(grep("kingdom_Eukaryota;supergroup_Opisthokonta", parsed_eLSA_int_DCM$X_tax))
parsed_eLSA_int_DCM$X_tax_abr[grep("kingdom_Eukaryota;supergroup_Opisthokonta", parsed_eLSA_int_DCM$X_tax)]="Opisthokonta"

##All hacrobia except prymnesiophytes, Haptophytes, and Phaeocystis set to "Hacrobia" 

length(grep("kingdom_Eukaryota;supergroup_Hacrobia", parsed_eLSA_int_DCM$X_tax)) 
parsed_eLSA_int_DCM$X_tax_abr[grep("kingdom_Eukaryota;supergroup_Hacrobia", parsed_eLSA_int_DCM$X_tax)]="Hacrobia"

length(grep("kingdom_Eukaryota;supergroup_Hacrobia;division_Haptophyta;class_Prymnesiophyceae;order_Prymnesiales;family_Prymnesiaceae;genus_Prymnesium;species_Prymnesium_sp.", parsed_eLSA_int_DCM$X_tax))
parsed_eLSA_int_DCM$X_tax_abr[grep("kingdom_Eukaryota;supergroup_Hacrobia;division_Haptophyta;class_Prymnesiophyceae;order_Prymnesiales;family_Prymnesiaceae;genus_Prymnesium;species_Prymnesium_sp.", parsed_eLSA_int_DCM$X_tax)]="Prymnesiophyte"

length(grep("kingdom_Eukaryota;supergroup_Hacrobia;division_Haptophyta;class_Prymnesiophyceae;order_Prymnesiales;family_Chrysochromulinaceae;genus_Chrysochromulina", parsed_eLSA_int_DCM$X_tax)) #87
parsed_eLSA_int_DCM$X_tax_abr[grep("kingdom_Eukaryota;supergroup_Hacrobia;division_Haptophyta;class_Prymnesiophyceae;order_Prymnesiales;family_Chrysochromulinaceae;genus_Chrysochromulina", parsed_eLSA_int_DCM$X_tax)]="Chrysochromulina"

#Phaeocystis
length(grep("kingdom_Eukaryota;supergroup_Hacrobia;division_Haptophyta;class_Prymnesiophyceae;order_Phaeocystales;family_Phaeocystaceae;genus_Phaeocystis", parsed_eLSA_int_DCM$X_tax)) #46
parsed_eLSA_int_DCM$X_tax_abr[grep("kingdom_Eukaryota;supergroup_Hacrobia;division_Haptophyta;class_Prymnesiophyceae;order_Phaeocystales;family_Phaeocystaceae;genus_Phaeocystis", parsed_eLSA_int_DCM$X_tax)]="Phaeocystis"

#Haptophyte
length(grep("kingdom_Eukaryota;supergroup_Hacrobia;division_Haptophyta;class_Haptophyta", parsed_eLSA_int_DCM$X_tax)) #3
parsed_eLSA_int_DCM$X_tax_abr[grep("kingdom_Eukaryota;supergroup_Hacrobia;division_Haptophyta;class_Haptophyta", parsed_eLSA_int_DCM$X_tax)]="Haptophyte"

#Insert Brad ix's
length(grep("kingdom_Eukaryota;supergroup_Hacrobia;division_Haptophyta;class_Prymnesiophyceae;order_Prymnesiophyceae_X;family_Braarudosphaeraceae;genus_Braarudosphaeraceae_X", parsed_eLSA_int_DCM$X_tax)) #24
parsed_eLSA_int_DCM$X[grep("kingdom_Eukaryota;supergroup_Hacrobia;division_Haptophyta;class_Prymnesiophyceae;order_Prymnesiophyceae_X;family_Braarudosphaeraceae;genus_Braarudosphaeraceae_X", parsed_eLSA_int_DCM$X_tax)] #They are all #7

for(i in Brad_ix){
  if(length(grep(i, parsed_eLSA_int_DCM$X))>1){
    print(i)
    print(grep(i, parsed_eLSA_int_DCM$X)) 
  }
} #Hmm ok, only Brad3 and Brad7 

parsed_eLSA_int_DCM$X_tax_abr[grep(Brad_ix[4], parsed_eLSA_int_DCM$X)]="Braarudospharea_bigelowii_3"
parsed_eLSA_int_DCM$X_tax_abr[grep(Brad_ix[7], parsed_eLSA_int_DCM$X)]="Braarudosphaera_sp_7x"


#Moving on!
length(grep("kingdom_Eukaryota;supergroup_Archaeplastida;division_Rhodophyta", parsed_eLSA_int_DCM$X_tax))
parsed_eLSA_int_DCM$X_tax_abr[grep("kingdom_Eukaryota;supergroup_Archaeplastida;division_Rhodophyta", parsed_eLSA_int_DCM$X_tax)]="Archeplastida"

length(grep("kingdom_Eukaryota;supergroup_Archaeplastida;division_Chlorophyta", parsed_eLSA_int_DCM$X_tax))
parsed_eLSA_int_DCM$X_tax_abr[grep("kingdom_Eukaryota;supergroup_Archaeplastida;division_Chlorophyta", parsed_eLSA_int_DCM$X_tax)]="Chlorophyta"

length(grep("kingdom_Eukaryota;supergroup_Archaeplastida;division_Chlorophyta;class_Mamiellophyceae;order_Mamiellales;family_Bathycoccaceae", parsed_eLSA_int_DCM$X_tax))
parsed_eLSA_int_DCM$X_tax_abr[grep("kingdom_Eukaryota;supergroup_Archaeplastida;division_Chlorophyta;class_Mamiellophyceae;order_Mamiellales;family_Bathycoccaceae", parsed_eLSA_int_DCM$X_tax)]="Bathycoccus"

#Dinoflagellates
parsed_eLSA_int_DCM$X_tax_abr[grep(temp[grep("Alveolata", temp)][1], parsed_eLSA_int_DCM$X_tax)]="Alveolata"

length(grep("kingdom_Eukaryota;supergroup_Alveolata;division_Dinoflagellata", parsed_eLSA_int_DCM$X_tax)) #1554 #WOWZA! 
parsed_eLSA_int_DCM$X_tax_abr[grep("kingdom_Eukaryota;supergroup_Alveolata;division_Dinoflagellata", parsed_eLSA_int_DCM$X_tax)]="Dinoflagellate"

length(grep("kingdom_Eukaryota;supergroup_Alveolata;division_Dinoflagellata;class_Syndiniales", parsed_eLSA_int_DCM$X_tax)) #920
parsed_eLSA_int_DCM$X_tax_abr[grep("kingdom_Eukaryota;supergroup_Alveolata;division_Dinoflagellata;class_Syndiniales", parsed_eLSA_int_DCM$X_tax)]="Syndiniales"

length(grep("kingdom_Eukaryota;supergroup_Alveolata;division_Dinoflagellata;class_Dinophyta", parsed_eLSA_int_DCM$X_tax)) #0
#parsed_eLSA_int_5m$X_tax_abr[grep("kingdom_Eukaryota;supergroup_Alveolata;division_Dinoflagellata;class_Dinophyta", parsed_eLSA_int_5m$X_taxonomy)]="Dinophyta"

length(grep("kingdom_Eukaryota;supergroup_Alveolata;division_Dinoflagellata;class_Dinophyceae", parsed_eLSA_int_DCM$X_tax)) #438
parsed_eLSA_int_DCM$X_tax_abr[grep("kingdom_Eukaryota;supergroup_Alveolata;division_Dinoflagellata;class_Dinophyceae", parsed_eLSA_int_DCM$X_tax)]="Dinophyceae"

length(grep("kingdom_Eukaryota;supergroup_Alveolata;division_Dinoflagellata;class_Dinophyceae;order_Gymnodiniales;family_Gymnodiniaceae;genus_Lepidodinium", parsed_eLSA_int_DCM$X_tax)) #39
parsed_eLSA_int_DCM$X_tax_abr[grep("kingdom_Eukaryota;supergroup_Alveolata;division_Dinoflagellata;class_Dinophyceae;order_Gymnodiniales;family_Gymnodiniaceae;genus_Lepidodinium", parsed_eLSA_int_DCM$X_tax)]="Lepidodinium"

length(grep("kingdom_Eukaryota;supergroup_Alveolata;division_Dinoflagellata;class_Dinophyceae;order_Gymnodiniales", parsed_eLSA_int_DCM$X_tax)) #138
parsed_eLSA_int_DCM$X_tax_abr[grep("kingdom_Eukaryota;supergroup_Alveolata;division_Dinoflagellata;class_Dinophyceae;order_Gymnodiniales", parsed_eLSA_int_DCM$X_tax)]="Gymnodinium"

length(grep("kingdom_Eukaryota;supergroup_Alveolata;division_Ciliophora;class_Spirotrichea;order_Tintinnida", parsed_eLSA_int_DCM$X_tax))
parsed_eLSA_int_DCM$X_tax_abr[grep("kingdom_Eukaryota;supergroup_Alveolata;division_Ciliophora;class_Spirotrichea;order_Tintinnida", parsed_eLSA_int_DCM$X_tax)]="Tintinnida"

length(grep(" kingdom_Eukaryota;supergroup_Alveolata;division_Ciliophora", parsed_eLSA_int_DCM$X_tax)) #461
parsed_eLSA_int_DCM$X_tax_abr[grep(" kingdom_Eukaryota;supergroup_Alveolata;division_Ciliophora;class_Spirotrichea;order_Strombidiida", parsed_eLSA_int_DCM$X_tax)]="Ciliophora"

for(i in grep("Ciliophora", temp)){
  parsed_eLSA_int_DCM$X_tax_abr[grep(temp[i], parsed_eLSA_int_DCM$X_taxonomy)]="Ciliophora"
}

length(grep("kingdom_Eukaryota;supergroup_Alveolata;division_Apicomplexa", parsed_eLSA_int_DCM$X_tax)) #21
parsed_eLSA_int_DCM$X_tax_abr[grep("kingdom_Eukaryota;supergroup_Alveolata;division_Apicomplexa", parsed_eLSA_int_DCM$X_tax)]="Ampicomplexa"


#one more!
##parsed_eLSA_int_DCM$X_tax_abr[grep(temp[grep("Amoebozoa", temp)], parsed_eLSA_int_DCM$X_taxonomy)]="Amoebozoa"
length(grep("kingdom_Eukaryota;supergroup_Amoebozoa", parsed_eLSA_int_DCM$X_tax)) #0
parsed_eLSA_int_DCM$X_tax_abr[grep(temp[grep("kingdom_Eukaryota;supergroup_Apusozoa",temp)], parsed_eLSA_int_DCM$X_tax)]="Apusozoa"

#OKAY DOUBLE CHECK IT 
length(grep("SAR", parsed_eLSA_int_DCM$X_tax_abr)) #216  
temp=sort(unique(parsed_eLSA_int_DCM$X_tax[grep("SAR", parsed_eLSA_int_DCM$X_tax_abr)]))
length(temp) #17
length(grep("Ciliophora", temp)) #14
temp[grep("Ciliophora", temp, invert=T)]

#FINALLY, all the remainder are just "Eukaryotes" 
parsed_eLSA_int_DCM$X_tax_abr[grep("SAR", parsed_eLSA_int_DCM$X_tax_abr)]="Eukaryote"
unique(parsed_eLSA_int_DCM$X_tax_abr) #21 long 

#parsed_eLSA_int_DCM$X_tax_abr[grep("UCYN-A", parsed_eLSA_int_DCM$X_tax)]="UCYN-A_1_AE"

#Abbreviate Y taxonomy-----
#Start by labeling everything "SAR"
parsed_eLSA_int_DCM$Y_tax_abr="SAR"

#Manually change taxa, oh god :( #At the "supergroup" level 

##Start with the Stramenopiles: All stramenopiles except Pseudo-nitzschia and Chaetoceros set to "Stramenopiles"

length(grep("kingdom_Eukaryota;supergroup_Stramenopiles", parsed_eLSA_int_DCM$Y_tax)) #608 #That's like 10% of the X taxa
parsed_eLSA_int_DCM$Y_tax_abr[grep("kingdom_Eukaryota;supergroup_Stramenopiles;", parsed_eLSA_int_DCM$Y_tax)]="Stramenopile"
parsed_eLSA_int_DCM$Y_tax_abr[grep(temp[grep("Stramenopile", temp)], parsed_eLSA_int_DCM$Y_tax)]="Stramenopile"

length(grep("kingdom_Eukaryota;supergroup_Stramenopiles;division_Ochrophyta;class_Bacillariophyta;order_Bacillariophyta_X;family_Raphid-pennate;genus_Pseudo-nitzschia", parsed_eLSA_int_DCM$Y_tax)) #20
parsed_eLSA_int_DCM$Y_tax_abr[grep("kingdom_Eukaryota;supergroup_Stramenopiles;division_Ochrophyta;class_Bacillariophyta;order_Bacillariophyta_X;family_Raphid-pennate;genus_Pseudo-nitzschia", parsed_eLSA_int_DCM$Y_taxonomy)]="Pseudo-nitzchia"

length(grep("kingdom_Eukaryota;supergroup_Stramenopiles;division_Ochrophyta;class_Bacillariophyta;order_Bacillariophyta_X;family_Polar-centric-Mediophyceae;genus_Chaetoceros", parsed_eLSA_int_DCM$Y_tax))
parsed_eLSA_int_DCM$Y_tax_abr[grep("kingdom_Eukaryota;supergroup_Stramenopiles;division_Ochrophyta;class_Bacillariophyta;order_Bacillariophyta_X;family_Polar-centric-Mediophyceae;genus_Chaetoceros", parsed_eLSA_int_DCM$Y_tax)]="Chaetoceros"

length(grep("kingdom_Eukaryota;supergroup_Rhizaria", parsed_eLSA_int_DCM$Y_tax))
parsed_eLSA_int_DCM$Y_tax_abr[grep("Rhizaria", parsed_eLSA_int_DCM$Y_tax)]="Rhizaria"

length(grep("kingdom_Eukaryota;supergroup_Opisthokonta", parsed_eLSA_int_DCM$Y_tax))
parsed_eLSA_int_DCM$Y_tax_abr[grep("kingdom_Eukaryota;supergroup_Opisthokonta", parsed_eLSA_int_DCM$Y_tax)]="Opisthokonta"

##All hacrobia except prymnesiophytes, Haptophytes, and Phaeocystis set to "Hacrobia" 

length(grep("kingdom_Eukaryota;supergroup_Hacrobia", parsed_eLSA_int_DCM$Y_tax)) 
parsed_eLSA_int_DCM$Y_tax_abr[grep("kingdom_Eukaryota;supergroup_Hacrobia", parsed_eLSA_int_DCM$Y_tax)]="Hacrobia"

length(grep("kingdom_Eukaryota;supergroup_Hacrobia;division_Haptophyta;class_Prymnesiophyceae;order_Prymnesiales;family_Prymnesiaceae;genus_Prymnesium;species_Prymnesium_sp.", parsed_eLSA_int_DCM$Y_tax))
parsed_eLSA_int_DCM$Y_tax_abr[grep("kingdom_Eukaryota;supergroup_Hacrobia;division_Haptophyta;class_Prymnesiophyceae;order_Prymnesiales;family_Prymnesiaceae;genus_Prymnesium;species_Prymnesium_sp.", parsed_eLSA_int_DCM$Y_tax)]="Prymnesiophyte"

length(grep("kingdom_Eukaryota;supergroup_Hacrobia;division_Haptophyta;class_Prymnesiophyceae;order_Prymnesiales;family_Chrysochromulinaceae;genus_Chrysochromulina", parsed_eLSA_int_DCM$Y_tax)) #124
parsed_eLSA_int_DCM$Y_tax_abr[grep("kingdom_Eukaryota;supergroup_Hacrobia;division_Haptophyta;class_Prymnesiophyceae;order_Prymnesiales;family_Chrysochromulinaceae;genus_Chrysochromulina", parsed_eLSA_int_DCM$Y_tax)]="Chrysochromulina"

#Phaeocystis
length(grep("kingdom_Eukaryota;supergroup_Hacrobia;division_Haptophyta;class_Prymnesiophyceae;order_Phaeocystales;family_Phaeocystaceae;genus_Phaeocystis", parsed_eLSA_int_DCM$Y_tax)) #66
parsed_eLSA_int_DCM$Y_tax_abr[grep("kingdom_Eukaryota;supergroup_Hacrobia;division_Haptophyta;class_Prymnesiophyceae;order_Phaeocystales;family_Phaeocystaceae;genus_Phaeocystis", parsed_eLSA_int_DCM$Y_tax)]="Phaeocystis"

#Haptophyte
length(grep("kingdom_Eukaryota;supergroup_Hacrobia;division_Haptophyta;class_Haptophyta", parsed_eLSA_int_DCM$Y_tax)) #0
#parsed_eLSA_int_DCM$Y_tax_abr[grep("kingdom_Eukaryota;supergroup_Hacrobia;division_Haptophyta;class_Haptophyta", parsed_eLSA_int_DCM$Y_tax)]="Haptophyte"

#Insert Brad ix's
length(grep("kingdom_Eukaryota;supergroup_Hacrobia;division_Haptophyta;class_Prymnesiophyceae;order_Prymnesiophyceae_X;family_Braarudosphaeraceae;genus_Braarudosphaeraceae", parsed_eLSA_int_DCM$Y_tax)) #1
parsed_eLSA_int_DCM$Y[grep("kingdom_Eukaryota;supergroup_Hacrobia;division_Haptophyta;class_Prymnesiophyceae;order_Prymnesiophyceae_X;family_Braarudosphaeraceae;genus_Braarudosphaeraceae_X", parsed_eLSA_int_DCM$Y_tax)] #It's only #7

#for(i in Brad_ix){
  if(length(grep(i, parsed_eLSA_int_DCM$X))>1){
    print(i)
    print(grep(i, parsed_eLSA_int_DCM$X)) 
  }
} #Hmm ok, only Brad3 and Brad7 

#parsed_eLSA_int_DCM$Y_tax_abr[grep(Brad_ix[4], parsed_eLSA_int_DCM$X)]="Braarudospharea_bigelowii_3"
parsed_eLSA_int_DCM$Y_tax_abr[grep(Brad_ix[7], parsed_eLSA_int_DCM$Y)]="Braarudosphaera_sp_7x"


#Moving on!
length(grep("kingdom_Eukaryota;supergroup_Archaeplastida;division_Rhodophyta", parsed_eLSA_int_DCM$Y_tax))
parsed_eLSA_int_DCM$Y_tax_abr[grep("kingdom_Eukaryota;supergroup_Archaeplastida;division_Rhodophyta", parsed_eLSA_int_DCM$Y_tax)]="Archeplastida"

length(grep("kingdom_Eukaryota;supergroup_Archaeplastida;division_Chlorophyta", parsed_eLSA_int_DCM$Y_tax))
parsed_eLSA_int_DCM$Y_tax_abr[grep("kingdom_Eukaryota;supergroup_Archaeplastida;division_Chlorophyta", parsed_eLSA_int_DCM$Y_tax)]="Chlorophyta"

length(grep("kingdom_Eukaryota;supergroup_Archaeplastida;division_Chlorophyta;class_Mamiellophyceae;order_Mamiellales;family_Bathycoccaceae", parsed_eLSA_int_DCM$Y_tax))
parsed_eLSA_int_DCM$Y_tax_abr[grep("kingdom_Eukaryota;supergroup_Archaeplastida;division_Chlorophyta;class_Mamiellophyceae;order_Mamiellales;family_Bathycoccaceae", parsed_eLSA_int_DCM$Y_tax)]="Bathycoccus"

#Dinoflagellates
##parsed_eLSA_int_DCM$Y_tax_abr[grep(temp[grep("Alveolata", temp)][1], parsed_eLSA_int_DCM$Y_tax)]="Alveolata"

length(grep("kingdom_Eukaryota;supergroup_Alveolata;division_Dinoflagellata", parsed_eLSA_int_DCM$Y_tax)) #1307 #WOWZA! 
parsed_eLSA_int_DCM$Y_tax_abr[grep("kingdom_Eukaryota;supergroup_Alveolata;division_Dinoflagellata", parsed_eLSA_int_DCM$Y_tax)]="Dinoflagellate"

length(grep("kingdom_Eukaryota;supergroup_Alveolata;division_Dinoflagellata;class_Syndiniales", parsed_eLSA_int_DCM$Y_tax)) #846
parsed_eLSA_int_DCM$Y_tax_abr[grep("kingdom_Eukaryota;supergroup_Alveolata;division_Dinoflagellata;class_Syndiniales", parsed_eLSA_int_DCM$Y_tax)]="Syndiniales"

length(grep("kingdom_Eukaryota;supergroup_Alveolata;division_Dinoflagellata;class_Dinophyta", parsed_eLSA_int_DCM$Y_tax)) #0
#parsed_eLSA_int_5m$Y_tax_abr[grep("kingdom_Eukaryota;supergroup_Alveolata;division_Dinoflagellata;class_Dinophyta", parsed_eLSA_int_5m$Y_taxonomy)]="Dinophyta"

length(grep("kingdom_Eukaryota;supergroup_Alveolata;division_Dinoflagellata;class_Dinophyceae", parsed_eLSA_int_DCM$Y_tax)) #438
parsed_eLSA_int_DCM$Y_tax_abr[grep("kingdom_Eukaryota;supergroup_Alveolata;division_Dinoflagellata;class_Dinophyceae", parsed_eLSA_int_DCM$Y_tax)]="Dinophyceae"

length(grep("kingdom_Eukaryota;supergroup_Alveolata;division_Dinoflagellata;class_Dinophyceae;order_Gymnodiniales;family_Gymnodiniaceae;genus_Lepidodinium", parsed_eLSA_int_DCM$Y_tax)) #4
parsed_eLSA_int_DCM$Y_tax_abr[grep("kingdom_Eukaryota;supergroup_Alveolata;division_Dinoflagellata;class_Dinophyceae;order_Gymnodiniales;family_Gymnodiniaceae;genus_Lepidodinium", parsed_eLSA_int_DCM$Y_tax)]="Lepidodinium"

length(grep("kingdom_Eukaryota;supergroup_Alveolata;division_Dinoflagellata;class_Dinophyceae;order_Gymnodiniales", parsed_eLSA_int_DCM$Y_tax)) #80
parsed_eLSA_int_DCM$Y_tax_abr[grep("kingdom_Eukaryota;supergroup_Alveolata;division_Dinoflagellata;class_Dinophyceae;order_Gymnodiniales", parsed_eLSA_int_DCM$Y_tax)]="Gymnodinium"

length(grep(" kingdom_Eukaryota;supergroup_Alveolata;division_Ciliophora", parsed_eLSA_int_DCM$Y_tax))
parsed_eLSA_int_DCM$Y_tax_abr[grep(" kingdom_Eukaryota;supergroup_Alveolata;division_Ciliophora", parsed_eLSA_int_DCM$Y_tax)]="Ciliophora"

length(grep("kingdom_Eukaryota;supergroup_Alveolata;division_Ciliophora;class_Spirotrichea;order_Tintinnida", parsed_eLSA_int_DCM$Y_tax))
parsed_eLSA_int_DCM$Y_tax_abr[grep("kingdom_Eukaryota;supergroup_Alveolata;division_Ciliophora;class_Spirotrichea;order_Tintinnida", parsed_eLSA_int_DCM$Y_tax)]="Tintinnida"

length(grep(" kingdom_Eukaryota;supergroup_Alveolata;division_Ciliophora", parsed_eLSA_int_DCM$Y_tax)) #358
parsed_eLSA_int_DCM$Y_tax_abr[grep(" kingdom_Eukaryota;supergroup_Alveolata;division_Ciliophora;class_Spirotrichea;order_Strombidiida", parsed_eLSA_int_DCM$Y_tax)]="Ciliophora"

for(i in grep("Ciliophora", temp)){
  parsed_eLSA_int_DCM$Y_tax_abr[grep(temp[i], parsed_eLSA_int_DCM$Y_taxonomy)]="Ciliophora"
}

length(grep("kingdom_Eukaryota;supergroup_Alveolata;division_Apicomplexa", parsed_eLSA_int_DCM$Y_tax)) #22
parsed_eLSA_int_DCM$Y_tax_abr[grep("kingdom_Eukaryota;supergroup_Alveolata;division_Apicomplexa", parsed_eLSA_int_DCM$Y_tax)]="Ampicomplexa"


#one more!
length(grep(" kingdom_Eukaryota;supergroup_Apusozoa", parsed_eLSA_int_DCM$Y_tax)) #19
parsed_eLSA_int_DCM$Y_tax_abr[grep(" kingdom_Eukaryota;supergroup_Apusozoa", parsed_eLSA_int_DCM$Y_tax)]="Apusozoa"
##parsed_eLSA_int_DCM$Y_tax_abr[grep(temp[grep("Amoebozoa", temp)], parsed_eLSA_int_DCM$Y_taxonomy)]="Amoebozoa"
length(grep("kingdom_Eukaryota;supergroup_Amoebozoa", parsed_eLSA_int_DCM$Y_tax)) #0
#parsed_eLSA_int_DCM$Y_tax_abr[grep(temp[grep("kingdom_Eukaryota;supergroup_Apusozoa",temp)], parsed_eLSA_int_DCM$Y_tax)]="Apusozoa"

#OKAY DOUBLE CHECK IT 
length(grep("SAR", parsed_eLSA_int_DCM$Y_tax_abr)) #296
temp=sort(unique(parsed_eLSA_int_DCM$Y_tax[grep("SAR", parsed_eLSA_int_DCM$Y_tax_abr)]))
length(temp) #12
length(grep("Ciliophora", temp)) #9
temp[grep("Ciliophora", temp, invert=T)]

#Add in UCYN-A
parsed_eLSA_int_DCM$Y_tax_abr[grep("UCYNA", parsed_eLSA_int_DCM$Y_tax)]="UCYNA-1_AE"

#FINALLY, all the remainder are just "Eukaryotes" 
parsed_eLSA_int_DCM$Y_tax_abr[grep("SAR", parsed_eLSA_int_DCM$Y_tax_abr)]="Eukaryote"
unique(parsed_eLSA_int_DCM$Y_tax_abr) #21 long


#Abbreivate X and Y hashes because I forgot :( 
length(unique(parsed_eLSA_int_DCM$X)) #465
length(unique(abbreviate(parsed_eLSA_int_DCM$X, minlength=4))) #465
parsed_eLSA_int_DCM$X_abr=abbreviate(parsed_eLSA_int_DCM$X, minlength=4)
head(unique(parsed_eLSA_int_DCM$X_abr), n=15)

length(unique(parsed_eLSA_int_DCM$Y)) #454
length(unique(abbreviate(parsed_eLSA_int_DCM$Y, minlength=4))) #454
parsed_eLSA_int_DCM$Y_abr=abbreviate(parsed_eLSA_int_DCM$Y, minlength=4)
head(unique(parsed_eLSA_int_DCM$Y_abr), n=15)

#Write it out!  
dim(parsed_eLSA_int_DCM)
names(parsed_eLSA_int_DCM)
write.table(x=parsed_eLSA_int_DCM, file="eLSA_output/piecewise_eLSA/int_DCM/parsed_eLSA_output_int_DCM_01.11.2020.tsv", sep="\t", quote=F, row.names=F)

