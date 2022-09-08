#Plot normalized data of UCYN-A ASVs vs. 18S taxa!!!
#Ver. 12.09.21 #Updated with different colors for lines

rm(list=ls())
getwd()
setwd("Desktop/DataAnalysis/SPOT_16S_18S_eLSA/")

#1. Read in files-----
proks_tax=read.table("ModifiedFiles/Proks_tax_classified_23072020_SILVA132.tsv", sep="\t", header=T, stringsAsFactors = F)
dim(proks_tax) #78953 3
names(proks_tax)
head(proks_tax)

UCYNA_ix=proks_tax$Feature.ID[grep("UCYN-A", proks_tax$Taxon)]
UCYNA_ix

euks_tax=read.table("ModifiedFiles/Euks_tax_classified_19052020_SILVA132.tsv", sep="\t", header=T, stringsAsFactors = F)
dim(euks_tax) #30537 3
names(euks_tax)
head(euks_tax)

Brad_ix=euks_tax$Feature.ID[grep("Braarudosphaera", euks_tax$Taxon)]

data=read.table("ModifiedFiles/9.SPOT_16S_no_chloro_18S_no_metaz_normalized_proportions.tsv", sep="\t", header=T, stringsAsFactors = F, row.names = 1)  #OTU ID = row names
dim(data) #100686 1208
head(colnames(data))
tail(colnames(data))
length(grep("5m", colnames(data))) #313
length(grep("5m.AE", colnames(data))) #125 :) 

#Read in eLSA outputs to grep hashes vs. hashes abbreviations 
eLSA_CLR_5m=read.table("eLSA_output/piecewise_eLSA/CLR_5m/parsed_eLSA_output_CLR_5m_11.24.2020.tsv", sep="\t", header=T, stringsAsFactors = F)
dim(eLSA_CLR_5m) #12313 33
length(unique(eLSA_CLR_5m$X_abr)) #304

#Set up vector of dates
all_dates <- c()
for(i in c(1:ncol(data))){
  a <- strsplit(colnames(data)[i], split=".", fixed=T)
  b <- paste(a[[1]][2],a[[1]][3],a[[1]][4], sep="-")
  all_dates <- c(all_dates, b)
}
head(all_dates)
length(all_dates) #1208
length(unique(all_dates)) #213

all_dates <- as.Date(all_dates, format="%Y-%m-%d")
length(all_dates)==ncol(data) #TRUE
length(all_dates[grep("5m.AE", colnames(data))])

dates_5m_AE=all_dates[grep("5m.AE", colnames(data))]
dates_5m_D=all_dates[grep("5m.D", colnames(data))]
dates_DCM_AE=all_dates[grep("DCM.AE", colnames(data))]
dates_DCM_D=all_dates[grep("DCM.D", colnames(data))]

#Which dates overlap? 
include_ix <- c()
for(i in c(1:length(dates_5m_D))){
  if(length(grep(dates_5m_D[i], dates_5m_AE))<1){
    print(dates_5m_D[i])
  } else{
    include_ix <- c(include_ix, i)
  }
}
include_ix
length(include_ix) #we only have 107 dates with AE and durapore both #Scratch that plan


#2. UCYN-A ASVs------
#UCYNA ASV1----- 
summary(as.numeric(data[grep(UCYNA_ix[1], rownames(data)),grep("5m.AE", colnames(data))]))
#WOW UCYNA1 ASV1 is at most 2.6% of whole community - really?? 
#When? 
which(as.numeric(data[grep(UCYNA_ix[1], rownames(data)),grep("5m.AE", colnames(data))])>0.02) #19
colnames(data)[grep("5m.AE", colnames(data))][19] #"SPOT.2008.07.09.5m.AE"
#Are there many other taxa on this date? 
length(which(data[,grep("SPOT.2008.07.09.5m.AE", colnames(data))]>0)) #There are 1450 taxa present on this date 
#Compare to another random date:
length(which(data[,grep("SPOT.2017.02.15.5m.AE", colnames(data))]>0)) #1011 #ok
length(which(data[,grep("SPOT.2017.07.11.5m.AE", colnames(data))]>0)) #513 

#Plot UCYNA ASV1 at 5m
grep(UCYNA_ix[1], rownames(data)) #430
pdf("Plots/12.2021_Plots_time/UCYNA_ASV1_5m_time_03.14.2022.pdf", width=7, height=3, onefile=F)
plot(as.Date(dates_5m_AE), as.numeric(data[grep(UCYNA_ix[1], rownames(data)), grep("5m.AE", colnames(data))]), type="b", pch=16, col="darkgrey", cex=1.5, lwd=3, xlab="SPOT Sampling Date", ylab="Relative abundance (16S + 18S)", main="Rel. abundance of UCYN-A ASV1 (3d852410f44d21c92c9c55fbbb25187e)/(16S+18S) in 5m AE") #this works fine, but x axis labels are 2005, 2010, 2015, that's it
dev.off()
#Left justified axis, no main title
grep(UCYNA_ix[1], rownames(data)) #430
plot(as.Date(dates_5m_AE), as.numeric(data[grep(UCYNA_ix[1], rownames(data)), grep("5m.AE", colnames(data))]), type="b", pch=16, lwd=5, lty=2, xlab="SPOT Sampling Date", ylab=, axes=F, ann=T)
axis(4) 

#When does A1 peak in the AE size fraction? 
colnames(data)[grep("5m.AE", colnames(data))][which(as.numeric(data[(grep(UCYNA_ix[1], rownames(data))), grep("5m.AE", colnames(data))])>0.005)]

#Plot up UCYNA-1 5m in durapore
grep(UCYNA_ix[1], rownames(data)) #430
pdf("Plots/12.2021_Plots_time/UCYNA_ASV1_5m_time_durapore.pdf", width=7, height=3)
plot(as.Date(dates_5m_D), as.numeric(data[grep(UCYNA_ix[1], rownames(data)), grep("5m.D", colnames(data))]), type="b", pch=16, lwd=5, lty=2, xlab="SPOT Sampling Date", ylab="Relative abundance (16S + 18S)", main="Rel. abundance of UCYN-A ASV1 (3d852410f44d21c92c9c55fbbb25187e)/(16S+18S) in 5m Durapore")
dev.off()

#When does A1 peak in the durapore? What dates after 2015, specifically?
colnames(data)[grep("5m.D", colnames(data))][which(as.numeric(data[(grep(UCYNA_ix[1], rownames(data))), grep("5m.D", colnames(data))])>0)]
colnames(data)[grep("5m.D", colnames(data))][which(as.numeric(data[(grep(UCYNA_ix[1], rownames(data))), grep("5m.D", colnames(data))])>0.0005)]
#Skip L justified axis, because this will probably be an inset or standalone figure 

#Plot up UCYNA-1 at DCM
pdf("Plots/03.16.22_Plots_time/UCYNA_ASV1_DCM_time_03.16.2022.pdf", width=7, height=3, onefile=F)
plot(as.Date(dates_DCM_AE), as.numeric(data[grep(UCYNA_ix[1], rownames(data)), grep("DCM.AE", colnames(data))]), type="b", pch=16, col="darkgrey", cex=1.5, lwd=3, xlab="SPOT Sampling Date", ylab="Relative abundance (16S + 18S)", main="Rel. abundance of UCYN-A ASV1 (3d852410f44d21c92c9c55fbbb25187e)/(16S+18S) at DCM AE")
dev.off()
#L justified axis
plot(as.Date(dates_DCM_AE), as.numeric(data[grep(UCYNA_ix[1], rownames(data)), grep("DCM.AE", colnames(data))]), type="b", pch=16, lwd=5, lty=2, xlab="SPOT Sampling Date", ylab="Relative abundance (16S + 18S)", axes=F, ann=T)
axis(4)

#Plot UCYN-A1 at DCM in Durapore
plot(as.Date(dates_DCM_D), as.numeric(data[grep(UCYNA_ix[1], rownames(data)), grep("DCM.D", colnames(data))]), type="b", pch=16, lwd=5, lty=2, xlab="SPOT Sampling Date", ylab="Relative abundance (16S + 18S)", main="Rel. abundance of UCYN-A ASV1 (3d852410f44d21c92c9c55fbbb25187e)/(16S+18S) at DCM Dura")

#Which dates did UCYNA peak at DCM? 
#AE first
dates_DCM_AE[which(as.numeric(data[grep(UCYNA_ix[1], row.names(data)), grep("DCM.AE", colnames(data))])>0)]
#Durapore second
dates_DCM_D[which(as.numeric(data[grep(UCYNA_ix[1], row.names(data)), grep("DCM.D", colnames(data))])>0)]

#UCYNA ASV3 (6115)-----
UCYNA_ix[3]
plot(as.Date(dates_5m_AE), as.numeric(data[grep(UCYNA_ix[3], rownames(data)), grep("5m.AE", colnames(data))]), type="b", pch=0, lwd=5, lty=2, xlab="SPOT Sampling Date", ylab="Relative abundance (16S + 18S)", main=paste("Rel. abundance of UCYN-A ASV3 (", UCYNA_ix[3], ") at 5m in AE", sep=""))

#DCM
plot(as.Date(dates_DCM_AE), as.numeric(data[grep(UCYNA_ix[3], rownames(data)), grep("DCM.AE", colnames(data))]), type="b", pch=0, lwd=5, lty=2, xlab="SPOT Sampling Date", ylab="Relative abundance (16S + 18S)", main=paste("Rel. abundance of UCYN-A ASV3 (", UCYNA_ix[3], ") at the DCM in AE", sep=""))

#UCYNA ASV4(a6411) at 5m-----
UCYNA_ix[4]
plot(as.Date(dates_5m_AE), as.numeric(data[grep(UCYNA_ix[4], rownames(data)), grep("5m.AE", colnames(data))]), type="b", pch=0, lwd=5, lty=2, xlab="SPOT Sampling Date", ylab="Relative abundance (16S + 18S)", main=paste("Rel. abundance of UCYN-A ASV4 (", UCYNA_ix[4], ") at 5m in AE", sep=""))

UCYNA_ix[4]
plot(as.Date(dates_DCM_AE), as.numeric(data[grep(UCYNA_ix[4], rownames(data)), grep("DCM.AE", colnames(data))]), type="b", pch=0, lwd=5, lty=2, xlab="SPOT Sampling Date", ylab="Relative abundance (16S + 18S)", main=paste("Rel. abundance of UCYN-A ASV4 (", UCYNA_ix[4], ") at the DCM in AE", sep=""))


#UCYNA ASV5-----
pdf("Plots/03.16.22_Plots_time/UCYNA_ASV5_5m_time_03.16.22.pdf", width=7, height=3, onefile=F)
plot(as.Date(dates_5m_AE), as.numeric(data[grep(UCYNA_ix[5], rownames(data)), grep("5m.AE", colnames(data))]), type="b", pch=25, col="darkgrey", cex=1, lwd=5, xlab="SPOT Sampling Date", ylab="Relative abundance (16S+18S)", main="Rel. abundance of UCYN-A ASV 5 (af1bb1f9fb1c3f3d18571e711df407bb) in AE 5m")
dev.off()
#L axis
plot(as.Date(dates_5m_AE), as.numeric(data[grep(UCYNA_ix[5], rownames(data)), grep("5m.AE", colnames(data))]), type="b", pch=25, col="darkgrey", cex=2, lwd=7, xlab="SPOT Sampling Date", ylab="Relative abundance (16S+18S)", axes=F, ann=T)
axis(4)

#Inset
pdf("Plots/12.2021_Plots_time/UCYNA_ASV5_time_5m_inset.pdf", width=2, height=3, onefile=F)
plot(as.Date(dates_5m_AE[grep("2009", dates_5m_AE)[1]:grep("2010", dates_5m_AE)[3]]), as.numeric(data[grep(UCYNA_ix[5], rownames(data)), grep("5m.AE", colnames(data))][grep("2009", dates_5m_AE)[1]:grep("2010", dates_5m_AE)[3]]), type="b", pch=25, col="darkgrey", cex=2, lwd=5, xlab="SPOT Sampling Date", ylab="Relative abundance (16S+18S)", ann=T, axes=T)
dev.off()

#When is that big spike? And how big is it? 
dates_5m_AE[which(as.numeric(data[grep(UCYNA_ix[5], rownames(data)), grep("5m.AE", colnames(data))])>0.003)]
data[grep(UCYNA_ix[5], rownames(data)), grep("2008.06.18.5m.AE|2015.07.14.5m.AE", colnames(data))]

#DCM
pdf("Plots/03.16.22_Plots_time/UCYNA_ASV5_DCM_time_03.16.22.pdf", width=7, height=3, onefile=F)
plot(as.Date(dates_DCM_AE), as.numeric(data[grep(UCYNA_ix[5], rownames(data)), grep("DCM.AE", colnames(data))]), type="b", pch=25,  col="darkgrey", cex=1, lwd=5, xlab="SPOT Sampling Date", ylab="Relative abundance (16S+18S)", main="Rel. abundance of UCYN-A ASV 5 (af1bb1f9fb1c3f3d18571e711df407bb) in AE DCM")
dev.off()
#L justify axis
plot(as.Date(dates_DCM_AE), as.numeric(data[grep(UCYNA_ix[5], rownames(data)), grep("DCM.AE", colnames(data))]), type="b", pch=25, lwd=5, lty=2, xlab="SPOT Sampling Date", ylab="Relative abundance (16S+18S)", axes=F, ann=T)

#UCYNA ASV6-----
pdf("Plots/12.2021_Plots_time/UCYNA_ASV6_time_5m_R.axis.pdf", width=7, height=3)
plot(as.Date(dates_5m_AE), as.numeric(data[grep(UCYNA_ix[6], rownames(data)), grep("5m.AE", colnames(data))]), type="b", pch=23, col="darkgrey", cex=2, lwd=5, lty=2, xlab="SPOT Sampling Date", ylab="Relative abundance (16S+18S)", main=paste("Rel. abundance of UCYNA ASV 6 (", UCYNA_ix[6], ") in AE 5m", sep=""))
dev.off()

plot(as.Date(dates_DCM_AE), as.numeric(data[grep(UCYNA_ix[6], rownames(data)), grep("DCM.AE", colnames(data))]), type="b", pch=23, lwd=5, lty=2, xlab="SPOT Sampling Date", ylab="Relative abundance (16S+18S)", main=paste("Rel. abundance of UCYNA ASV 6 (", UCYNA_ix[6], ") in AE DCM", sep=""))
#L justified
plot(as.Date(dates_DCM_AE), as.numeric(data[grep(UCYNA_ix[6], rownames(data)), grep("DCM.AE", colnames(data))]), type="b", pch=23, lwd=5, lty=2, xlab="SPOT Sampling Date", ylab="Relative abundance (16S+18S)", axes=F, ann=T)

#UCYN-A ASV1 vs. UCYN-A ASV5----
#All depths
pdf("Plots/03.16.22_Plots_time/UCYNA_ASV1_vs_UCYNA_ASV5_all_depths_03.25.2022.pdf", height=4, width=4)
plot(y=as.numeric(data[grep(UCYNA_ix[1], rownames(data)),]), x=as.numeric(data[grep(UCYNA_ix[5], rownames(data)),]), pch=16, ylab="Relative abundance of UCYN-A ASV1", xlab="Relative abundance of UCYN-A ASV5")
abline(LM, col="red")
text(y=0.025, x=0.003, "R^2=0.08")
dev.off()

LM=lm(as.numeric(data[grep(UCYNA_ix[1], rownames(data)),])~as.numeric(data[grep(UCYNA_ix[5], rownames(data)),]))
summary(LM)

#5m depth
plot(x=as.numeric(data[grep(UCYNA_ix[1], rownames(data)), grep("5m.AE", colnames(data))]), y=as.numeric(data[grep(UCYNA_ix[5], rownames(data)), grep("5m.AE", colnames(data))]), pch=16, xlab="Relative abundance of UCYN-A ASV1 (5m)", ylab="Relative abundance of UCYN-A ASV5 (5m)")
abline(LM, col="red")
text(x=0.025, y=0.003, "R^2=0.008")


LM=lm(as.numeric(data[grep(UCYNA_ix[5], rownames(data)),grep("5m.AE", colnames(data))])~as.numeric(data[grep(UCYNA_ix[1], rownames(data)), grep("5m.AE", colnames(data))]))
summary(LM)

#Switch X and Y
plot(y=as.numeric(data[grep(UCYNA_ix[1], rownames(data)), grep("5m.AE", colnames(data))]), x=as.numeric(data[grep(UCYNA_ix[5], rownames(data)), grep("5m.AE", colnames(data))]), pch=16, ylab="Relative abundance of UCYN-A ASV1 (5m)", xlab="Relative abundance of UCYN-A ASV5 (5m)")
abline(LM, col="red")
text(y=0.025, x=0.003, "R^2=0.008")


LM=lm(as.numeric(data[grep(UCYNA_ix[1], rownames(data)),grep("5m.AE", colnames(data))])~as.numeric(data[grep(UCYNA_ix[5], rownames(data)), grep("5m.AE", colnames(data))]))
summary(LM)

#DCM
plot(y=as.numeric(data[grep(UCYNA_ix[1], rownames(data)), grep("DCM.AE", colnames(data))]), x=as.numeric(data[grep(UCYNA_ix[5], rownames(data)), grep("DCM.AE", colnames(data))]), pch=16, ylab="Relative abundance of UCYN-A ASV1 (DCM)", xlab="Relative abundance of UCYN-A ASV5 (DCM)")
abline(LM, col="red")
text(y=0.025, x=0.003, "R^2=0.02")


LM=lm(as.numeric(data[grep(UCYNA_ix[1], rownames(data)),grep("DCM.AE", colnames(data))])~as.numeric(data[grep(UCYNA_ix[5], rownames(data)), grep("DCM.AE", colnames(data))]))
summary(LM)


#3.Braarudospharea and other prymnesiophytes------
#Brad ASV7 (be3cd...)-----
#pdf("Plots/12.2021_Plots_time/Brad_ASV7x_5m_time.pdf", width=7, height=3, onefile=F)
plot(as.Date(dates_5m_AE), as.numeric(data[grep(Brad_ix[7], rownames(data)), grep("5m.AE", colnames(data))]), type="b", lwd=2, lty=2, pch=24, col="forestgreen", xlab="SPOT Sampling Date", ylab="Relative abundance (16S+18S)",main="Rel. abundance of Brad ASV 7x (be3cdecefbceb0d8b25a2e42ed058b50)/ (16S+18S) in 5m AE")
#dev.off()

#Left justified axis
pdf("Plots/12.2021_Plots_time/Brad_ASV7x_5m_time.pdf", width=7, height=3, onefile=F)
plot(as.Date(dates_5m_AE), as.numeric(data[grep(Brad_ix[7], rownames(data)), grep("5m.AE", colnames(data))]), type="b", lwd=2, lty=2, pch=24, col="forestgreen", xlab="SPOT Sampling Date", ylab="Relative abundance (16S+18S)", axes=F, ann=T)
axis(4)
dev.off()

#Which dates was Brad7 peaking in abundance after 2015 in 5m samples?
colnames(data)[grep("5m.AE", colnames(data))][which(as.numeric(data[(grep(Brad_ix[7], rownames(data))), grep("5m.AE", colnames(data))])>0.005)]

#Plot Brad7 at DCM AE
plot(as.Date(dates_DCM_AE), as.numeric(data[grep(Brad_ix[7], rownames(data)), grep("DCM.AE", colnames(data))]), type="b", lwd=2, lty=2, pch=24, col="forestgreen", xlab="SPOT Sampling Date", ylab="Relative abundance (16S+18S)",main="Rel. abundance of Brad ASV 7x (be3cdecefbceb0d8b25a2e42ed058b50)/ (16S+18S) at DCM AE")
#L justify axis
pdf(width=7, height=3, onefile=F, "Plots/03.16.22_Plots_time/Brad_ASV7x_time_DCM_03.16.2022.pdf")
plot(as.Date(dates_DCM_AE), as.numeric(data[grep(Brad_ix[7], rownames(data)), grep("DCM.AE", colnames(data))]), type="b", lwd=5, lty=2, pch=24, col="forestgreen", xlab="SPOT Sampling Date", ylab="Relative abundance (16S+18S)",axes=F, ann=T)
axis(4)
dev.off()

#On which dates did Brad 7 appear in the DCM?
dates_DCM_AE[which(as.numeric(data[grep(Brad_ix[7], rownames(data)), grep("DCM.AE", colnames(data))])>0)]

#Brad ASV3 (technically 4) (70a5...)-----
plot(as.Date(dates_5m_AE), as.numeric(data[grep(Brad_ix[4], rownames(data)), grep("5m.AE", colnames(data))]), type="b", pch=0, col="dodgerblue", lwd=5, lty=2, xlab="SPOT Sampling Date", ylab="Relative abundance (16S+18S)", main=paste("Rel. abundance of Brad ASV 3 (", Brad_ix[4], ") in AE 5m", sep=""))
#Add in a line for Brad ASV1
lines(x=as.Date(dates_5m_AE), y=as.numeric(data[grep(Brad_ix[1], rownames(data)), grep("5m.AE", colnames(data))])) 
#, type="b", pch=16, col="purple", lwd=5, lty=2, xlab="SPOT Sampling Date", ylab="Relative abundance (16S+18S)", main=paste("Rel. abundance of Brad ASV 1 (", Brad_ix[1], ") in AE 5m", sep=""))

#L justify axis
pdf("Plots/12.2021_Plots_time/Brad_3_5m_time_Laxis.pdf", width=7, height=3, onefile=F)
plot(as.Date(dates_5m_AE), as.numeric(data[grep(Brad_ix[4], rownames(data)), grep("5m.AE", colnames(data))]), type="b", pch=0, col="dodgerblue", lwd=2, lty=2, xlab="SPOT Sampling Date", ylab="Relative abundance (16S+18S)", ann=T, axes=F)
axis(4)
dev.off()

#At the DCM

plot(as.Date(dates_DCM_AE), as.numeric(data[grep(Brad_ix[4], rownames(data)), grep("DCM.AE", colnames(data))]), type="b", pch=0, col="dodgerblue", lwd=5, lty=2, xlab="SPOT Sampling Date", ylab="Relative abundance (16S+18S)", main=paste("Rel. abundance of Brad ASV 3 (", Brad_ix[4], ") in AE DCM", sep=""))

#L justify
pdf("Plots/03.16.22_Plots_time/Brad_3_DCM_time_03.16.22.pdf", width=7, height=3, onefile=F)
plot(as.Date(dates_DCM_AE), as.numeric(data[grep(Brad_ix[4], rownames(data)), grep("DCM.AE", colnames(data))]), type="b", pch=0, col="dodgerblue", lwd=5, lty=2, xlab="SPOT Sampling Date", ylab="Relative abundance (16S+18S)", axes=F, ann=T)
axis(4)
dev.off()

#Brad ASV1 (0492)----- 
#not actually in Cytoskape networks, but WTH - showed up previously
plot(as.Date(dates_5m_AE), as.numeric(data[grep(Brad_ix[1], rownames(data)), grep("5m.AE", colnames(data))]), type="b", pch=16, col="purple", lwd=5, lty=2, xlab="SPOT Sampling Date", ylab="Relative abundance (16S+18S)", main=paste("Rel. abundance of Brad ASV 1 (", Brad_ix[1], ") in AE 5m", sep=""))
#L justify axis
pdf("Plots/12.2021_Plots_time/Brad_1_time_5m_L.axis.pdf", width=7, height=3)
plot(as.Date(dates_5m_AE), as.numeric(data[grep(Brad_ix[1], rownames(data)), grep("5m.AE", colnames(data))]), type="b", pch=16, col="purple", lwd=5, lty=2, xlab="SPOT Sampling Date", ylab="Relative abundance (16S+18S)", ann=T, axes=F)
axis(4)
dev.off()

#Inset
pdf("Plots/12.2021_Plots_time/Brad_3_time_5m_inset.pdf", width=2, height=3, onefile=F)
plot(as.Date(dates_5m_AE[grep("2009", dates_5m_AE)[1]:grep("2010", dates_5m_AE)[3]]), as.numeric(data[grep(Brad_ix[1], rownames(data)), grep("5m.AE", colnames(data))][grep("2009", dates_5m_AE)[1]:grep("2010", dates_5m_AE)[3]]), type="b", pch=16, col="purple", lwd=5, lty=2, xlab="SPOT Sampling Date", ylab="Relative abundance (16S+18S)", ann=T, axes=T)
dev.off()

#DCM
plot(as.Date(dates_DCM_AE), as.numeric(data[grep(Brad_ix[1], rownames(data)), grep("DCM.AE", colnames(data))]), type="b", pch=16, col="purple", lwd=5, lty=2, xlab="SPOT Sampling Date", ylab="Relative abundance (16S+18S)", main=paste("Rel. abundance of Brad ASV 1 (", Brad_ix[1], ") in AE DCM", sep=""))

#Brad ASV4 (technically 5) (8c144...)-----
plot(as.Date(dates_5m_AE), as.numeric(data[grep(Brad_ix[5], rownames(data)), grep("5m.AE", colnames(data))]), type="b", pch=23, col="mediumturquoise", lwd=5, lty=2, xlab="SPOT Sampling Date", ylab="Relative abundance (16S+18S)", main=paste("Rel. abundance of Brad ASV 4 (", Brad_ix[5], ") in AE 5m", sep=""))
#L justify axis
pdf("Plots/12.2021_Plots_time/Brad_4_time_5m_L.axis.pdf", width=7, height=3)
plot(as.Date(dates_5m_AE), as.numeric(data[grep(Brad_ix[5], rownames(data)), grep("5m.AE", colnames(data))]), type="b", pch=23, col="mediumturquoise", lwd=5, lty=2, xlab="SPOT Sampling Date", ylab="Relative abundance (16S+18S)", axes=F, ann=T)
axis(4)
dev.off()

#At the DCM
plot(as.Date(dates_DCM_AE), as.numeric(data[grep(Brad_ix[5], rownames(data)), grep("DCM.AE", colnames(data))]), type="b", pch=23, col="mediumturquoise", lwd=5, lty=2, xlab="SPOT Sampling Date", ylab="Relative abundance (16S+18S)", main=paste("Rel. abundance of Brad ASV 4 (", Brad_ix[5], ") in AE DCM", sep="")) #And it's zero...

#Brad ASV 5x (technically 2, 3246...)------
plot(as.Date(dates_5m_AE), as.numeric(data[grep(Brad_ix[2], rownames(data)), grep("5m.AE", colnames(data))]), type="b", pch=23, col="yellow", lwd=5, lty=2, xlab="SPOT Sampling Date", ylab="Relative abundance (16S+18S)", main=paste("Rel. abundance of Brad ASV 5x (", Brad_ix[2], ") in AE 5m", sep=""))
#L justify 
pdf("Plots/12.2021_Plots_time/Brad_5x_time_5m_L.axis.pdf", width=7, height=3)
plot(as.Date(dates_5m_AE), as.numeric(data[grep(Brad_ix[2], rownames(data)), grep("5m.AE", colnames(data))]), type="b", pch=23, col="yellow", lwd=5, lty=2, xlab="SPOT Sampling Date", ylab="Relative abundance (16S+18S)", axes=F, ann=T)
axis(4)
dev.off()

#At the DCM 
plot(as.Date(dates_DCM_AE), as.numeric(data[grep(Brad_ix[2], rownames(data)), grep("DCM.AE", colnames(data))]), type="b", pch=23, col="mediumturquoise", lwd=5, lty=2, xlab="SPOT Sampling Date", ylab="Relative abundance (16S+18S)", main=paste("Rel. abundance of Brad ASV 4 (", Brad_ix[2], ") in AE DCM", sep="")) #And it's zero

#Plot up UCYN-A ASV6 vs. all prymnesiophytes in a for-loop----
dim(euks_tax)
length(grep("class_Prymnesiophyceae", euks_tax$Taxon)) #643
head(euks_tax$Taxon[grep("class_Prymnesiophyceae", euks_tax$Taxon)], n=20)

prym_ix=euks_tax$Feature.ID[grep("class_Prymnesiophyceae", euks_tax$Taxon)]
length(prym_ix)

pdf("Plots/03.16.22_Plots_time/massive_file_UCYNA_ASV6_vs_all_prymnesiophytes_03.24.2022.pdf")
for(i in c(1:length(prym_ix))){
  print(i)
  plot(x=as.numeric(data[grep(prym_ix[i], rownames(data)), grep("5m.AE", colnames(data))]), y=as.numeric(data[grep(UCYNA_ix[6], rownames(data)), grep("5m.AE", colnames(data))]), pch=16, col="black", xlab=paste("Relative abundance of ", prym_ix[i]), ylab="Relative abundance of UCYNA ASV6")
}
dev.off()

#Okay, will probably need to delete this file
#Look through plots and find ASVs that look kinda sorta relationship with UCYN-A ASV6 
#Nothing too inspiring #Ok, no prymnesiophyte host

#4. Dinoflagellates, including Lepidodinium and Syndiniales (predators?)-------
#Lepidodinium, c114
unique(eLSA_CLR_5m$Y[grep("Lepidodinium", eLSA_CLR_5m$Y_tax_abr)])
pdf("Plots/12.2021_Plots_time/Lepidodinium_time_5m_R.axis.pdf", width=7, height=3)
plot(as.Date(dates_5m_AE), as.numeric(data[grep("c114523e0bef5840b096693e46f441a2", rownames(data)), grep("5m.AE", colnames(data))]), pch=18, type="b", col="hotpink", lwd=5, lty=2, xlab="SPOT Sampling Date", ylab="Relative abundance (16S+18S)", main="Rel. abundance of Lepidodinium ASV (c114523e0bef5840b096693e46f441a2) in AE 5m")
dev.off()

#L axis
pdf("Plots/12.2021_Plots_time/Lepidodinium_time_5m.pdf", width=7, height=3)
plot(as.Date(dates_5m_AE), as.numeric(data[grep("c114523e0bef5840b096693e46f441a2", rownames(data)), grep("5m.AE", colnames(data))]), pch=18, type="b", col="hotpink", lwd=5, lty=2, xlab="SPOT Sampling Date", ylab="Relative abundance (16S+18S)", main="Rel. abundance of Lepidodinium ASV (c114523e0bef5840b096693e46f441a2) in AE 5m", axes=F, ann=T)
axis(4)
dev.off()

#Like 15 ASVs of Syndiniales - write into a .PDF 
SYN=unique(eLSA_CLR_5m$X[grep("f143|e315|b310|b238|b220|8090|6944|6821|6603|6246|4372|4077|3650|1251|1082", eLSA_CLR_5m$X_abr)])
length(SYN)
pdf("normalized_plots_12.2020/Syndiniales_ASVS_time_5m_AE_12.15.20.pdf")
for(i in SYN){
  plot(as.Date(dates_5m_AE), as.numeric(data[grep(i, rownames(data)), grep("5m.AE", colnames(data))]), pch=23, type="b", col="red", lwd=5, lty=2, xlab="SPOT Sampling Date", ylab="Relative abundance (16S+18S)", main=paste("Rel. abundance of Syndiniales ASV \n", i, "in AE 5m"))
}
dev.off()
#CHE BELLO! Although the dimensions are not the same as UCYN-A graphs. Exported PDF graphs are each a square. 
#Here are the Syndiniales ASVs that look particularly similar to UCYN-A 1 ASV patterns: eefc..., 12516 (kinda), 6821, 80d9 - plot those up, not in a PDF 

plot(as.Date(dates_5m_AE), as.numeric(data[grep("eefc3152825b60051b3f78c504aca2a9", rownames(data)), grep("5m.AE", colnames(data))]), pch=23, type="b", col="red", lwd=5, lty=2, xlab="SPOT Sampling Date", ylab="Relative abundance (16S+18S)", main="Rel. abundance of Syndiniales ASV (eefc3152825b60051b3f78c504aca2a9) in AE 5m")

unique(eLSA_CLR_5m$X[grep("12516", eLSA_CLR_5m$X)])
plot(as.Date(dates_5m_AE), as.numeric(data[grep("12516469ae16dbe47afd1383f1a38d29", rownames(data)), grep("5m.AE", colnames(data))]), pch=23, type="b", col="red", lwd=5, lty=2, xlab="SPOT Sampling Date", ylab="Relative abundance (16S+18S)", main="Rel. abundance of Syndiniales ASV (12516469ae16dbe47afd1383f1a38d29) in AE 5m")

unique(eLSA_CLR_5m$X[grep("6821", eLSA_CLR_5m$X)])
plot(as.Date(dates_5m_AE), as.numeric(data[grep("68216c70c4cc7ae340a7f17cb6973f5b", rownames(data)), grep("5m.AE", colnames(data))]), pch=23, type="b", col="red", lwd=5, lty=2, xlab="SPOT Sampling Date", ylab="Relative abundance (16S+18S)", main="Rel. abundance of Syndiniales ASV (68216c70c4cc7ae340a7f17cb6973f5b) in AE 5m")

unique(eLSA_CLR_5m$Y[grep("80d9", eLSA_CLR_5m$Y)])
plot(as.Date(dates_5m_AE), as.numeric(data[grep("80d903e0cb7271e69b2ec8bd1f87f45c", rownames(data)), grep("5m.AE", colnames(data))]), pch=23, type="b", col="red", lwd=5, lty=2, xlab="SPOT Sampling Date", ylab="Relative abundance (16S+18S)", main="Rel. abundance of Syndiniales ASV (80d903e0cb7271e69b2ec8bd1f87f45c) in AE 5m")

#Plot ratio of UCYNA-1 and all these associated taxa vs. each other, send it to a .PDF
taxa=c(Brad_ix[7], "d55992e6da65321a9b3c0ce3426e73ac", "c114523e0bef5840b096693e46f441a2", "ffa92571884bcbe0193ce4a1e4e843be", "eefc3152825b60051b3f78c504aca2a9", "12516469ae16dbe47afd1383f1a38d29", "68216c70c4cc7ae340a7f17cb6973f5b", "80d903e0cb7271e69b2ec8bd1f87f45c")

#Now plot up Syn vs. Lepidinodinium into a .PDF
unique(eLSA_CLR_5m$X[grep("c114", eLSA_CLR_5m$X_abr)])
pdf("normalized_plots_12.2020/Syndiniales_ASVs_vs_Lepidodinium_5m_AE_12.17.20.pdf")
for(i in SYN){
  plot(as.numeric(data[grep(i, rownames(data)), grep("5m.AE", colnames(data))]), as.numeric(data[grep("c114523e0bef5840b096693e46f441a2", rownames(data)), grep("5m.AE",colnames(data))]), pch=8, lwd=5, lty=2, xlab=paste("Rel. abundance of Syndiniales ASV \n", i, "in AE 5m"), ylab="Relative abundance of Lepidodinium ASV \n (c114523e0bef5840b096693e46f441a2) in AE 5m")
}
dev.off()



#5. UCYNA ASVs vs. host ASVs------

#Best linear relationship is UCYN-A ASV1 with Brad7
#Plot across all depths
pdf("Plots/UCYNA_ASV1_Vs_Brad_7x_w.trendline_02.09.2022.pdf", width=4, height=4)
plot(y=as.numeric(data[grep(UCYNA_ix[1], rownames(data)),]), x=as.numeric(data[grep(Brad_ix[7], rownames(data)),]), ylab="Relative abundance UCYN-A ASV1", xlab=paste("Relative abundance Braarudospharea ASV7x"), pch=16, cex=1.5)
#Add in a trendline 
LM <- lm(as.numeric(data[grep(Brad_ix[7], rownames(data)),])~as.numeric(data[grep(UCYNA_ix[1], rownames(data)),]))
abline(LM, col="red", lwd=2)
text(x=0.002, y=0.025, labels="R^2=0.691")
dev.off()

#Plot at 5m only
pdf("Plots/03.16.22_Plots_time/UCYNA_ASV1_Vs_Brad_7x_5m_w.trendline_03.16.2022.pdf", width=4, height=4)
plot(y=as.numeric(data[grep(UCYNA_ix[1], rownames(data)),grep("5m.AE", colnames(data))]), x=as.numeric(data[grep(Brad_ix[7], rownames(data)),grep("5m.AE", colnames(data))]), ylab="Relative abundance UCYN-A ASV1", xlab=paste("Relative abundance Braarudospharea ASV7x"), pch=16, cex=1.5)
abline(LM, col="red", lwd=2)
text(x=0.002, y=0.025, labels="R^2=0.699")
dev.off()

LM <- lm(as.numeric(data[grep(Brad_ix[7], rownames(data)),grep("5m.AE", colnames(data))])~as.numeric(data[grep(UCYNA_ix[1], rownames(data)), grep("5m.AE", colnames(data))]))
summary(LM)

#At the DCM only
pdf("Plots/03.16.22_Plots_time/UCYNA_ASV1_Vs_Brad_7x_DCM_w.trendline_03.16.2022.pdf", width=4, height=4)
plot(y=as.numeric(data[grep(UCYNA_ix[1], rownames(data)),grep("DCM.AE", colnames(data))]), x=as.numeric(data[grep(Brad_ix[7], rownames(data)),grep("DCM.AE", colnames(data))]), ylab="Relative abundance UCYN-A ASV1", xlab=paste("Relative abundance Braarudospharea ASV7x"), pch=16, cex=1.5)
abline(LM, col="red", lwd=2)
text(x=0.001, y=0.0015, labels="R^2=0.347")
dev.off()

LM <- lm(as.numeric(data[grep(Brad_ix[7], rownames(data)),grep("DCM.AE", colnames(data))])~as.numeric(data[grep(UCYNA_ix[1], rownames(data)), grep("DCM.AE", colnames(data))]))
summary(LM)

#UCYN-A ASV4 and Brad 4
#All depths 
pdf("Plots/UCYNA_ASV5_Vs_Brad_4_w.trendline_02.10.2021.pdf", width=4, height=4)
plot(y=as.numeric(data[grep(UCYNA_ix[5], rownames(data)),]), x=as.numeric(data[grep(Brad_ix[4], rownames(data)),]), ylab="Relative abundance \n UCYN-A ASV5", xlab=paste("Relative abundance \n Braarudospharea ASV4"), pch=25, cex=1.5, bg="black")
#Add in a trendline 
LM <- lm(as.numeric(data[grep(Brad_ix[4], rownames(data)),])~as.numeric(data[grep(UCYNA_ix[5], rownames(data)),]))
abline(LM, col="red", lwd=2)
text(x=0.0031, y=0.0032, labels="R^2=0.420")
dev.off()

#5m only
pdf("Plots/03.16.22_Plots_time/UCYNA_ASV5_Vs_Brad_4_5m_w.trendline_03.16.2021.pdf", width=4, height=4)
plot(y=as.numeric(data[grep(UCYNA_ix[5], rownames(data)),grep("5m.AE", colnames(data))]), x=as.numeric(data[grep(Brad_ix[4], rownames(data)),grep("5m.AE", colnames(data))]), ylab="Relative abundance UCYN-A ASV5", xlab=paste("Relative abundance Braarudospharea ASV3"), pch=25, cex=1.5, bg="black")
#Add in a trendline 
abline(LM, col="red", lwd=2)
text(x=0.0031, y=0.0032, labels="R^2=0.402")
dev.off()

LM <- lm(as.numeric(data[grep(Brad_ix[4], rownames(data)),grep("5m.AE", colnames(data))])~as.numeric(data[grep(UCYNA_ix[5], rownames(data)),grep("5m.AE", colnames(data))]))
summary(LM)

#DCM only
pdf("Plots/03.16.22_Plots_time/UCYNA_ASV5_Vs_Brad_4_DCM_w.trendline_03.16.2021.pdf", width=4, height=4)
plot(y=as.numeric(data[grep(UCYNA_ix[5], rownames(data)),grep("DCM.AE", colnames(data))]), x=as.numeric(data[grep(Brad_ix[4], rownames(data)),grep("DCM.AE", colnames(data))]), ylab="Relative abundance UCYN-A ASV5", xlab=paste("Relative abundance Braarudospharea ASV4"), pch=25, cex=1.5, bg="black")
#Add in a trendline 
abline(LM, col="red", lwd=2)
text(x=0.0002, y=0.0015, labels="R^2=0.814")
dev.off()

LM <- lm(as.numeric(data[grep(Brad_ix[4], rownames(data)),grep("DCM.AE", colnames(data))])~as.numeric(data[grep(UCYNA_ix[5], rownames(data)),grep("DCM.AE", colnames(data))]))
summary(LM)

#Plot all ASVs of UCYN-A vs all ASVs of Brad with the package Liv recommended, ggpairs
install.packages("GGally")
library(GGally)

#Set up data
include_rows=c()
for(i in c(UCYNA_ix, Brad_ix)){
  a=grep(i, rownames(data))
  include_rows=c(include_rows, a)
}

#5m
data_5m=data[include_rows, grep("5m.AE", colnames(data))]
dim(data_5m) #13 125
rownames(data_5m)
head(colnames(data_5m), n=10)
data_5m_t=t(data_5m)
colnames(data_5m_t)
colnames(data_5m_t)=c(paste(rep("UCYNA", times=6), c(1:6), sep="_"),paste(rep("Brad", times=7), c("1", "5x", "2", "3", "4", "6x", "7x"), sep="_"))

#DCM
data_DCM=data[include_rows, grep("DCM.AE", colnames(data))]
dim(data_DCM) #13 125
rownames(data_DCM)
head(colnames(data_DCM), n=10)
data_DCM_t=t(data_DCM)
colnames(data_DCM_t)
colnames(data_DCM_t)=c(paste(rep("UCYNA", times=6), c(1:6), sep="_"),paste(rep("Brad", times=7), c("1", "5x", "2", "3", "4", "6x", "7x"), sep="_"))

#Now plot! 
pdf("Plots/03.16.22_Plots_time/All_UCYNA_Brad_ASVs_pairwise_03.16.2022.pdf", width=14, height=14)
pairs(as.data.frame(data_5m_t), labels=colnames(data_5m_t), pch=16) #Spruce this up a little in Affinity design
dev.off()


#6. Chrysochromulina ASVs-----
#0f7d43a484e7828d4d276ffcf4368daa
#unique(eLSA_CLR_5m$Y[grep("Lepidodinium", eLSA_CLR_5m$Y_tax_abr)])
pdf("Plots/12.2021_Plots_time/Chrysochromulina_1_time_5m_R.axis.pdf", width=7, height=3)
plot(as.Date(dates_5m_AE), as.numeric(data[grep("0f7d43a484e7828d4d276ffcf4368daa", rownames(data)), grep("5m.AE", colnames(data))]), pch=18, type="b", col="blue", lwd=5, lty=2, xlab="SPOT Sampling Date", ylab="Relative abundance (16S+18S)", main="Rel. abundance of Chrysochromulina ASV1 (0f7d43a484e7828d4d276ffcf4368daa) in AE 5m")
dev.off()

pdf("Plots/12.2021_Plots_time/Chrysochromulina_2_time_5m_R.axis.pdf", width=7, height=3)
plot(as.Date(dates_5m_AE), as.numeric(data[grep("d55992e6da65321a9b3c0ce3426e73ac", rownames(data)), grep("5m.AE", colnames(data))]), pch=18, type="b", col="blue", lwd=5, lty=2, xlab="SPOT Sampling Date", ylab="Relative abundance (16S+18S)", main="Rel. abundance of Chrysochromulina ASV2 (d55992e6da65321a9b3c0ce3426e73ac) in AE 5m")
dev.off()

pdf("Plots/12.2021_Plots_time/Chrysochromulina_3_time_5m_R.axis.pdf", width=7, height=3)
plot(as.Date(dates_5m_AE), as.numeric(data[grep("eaaf40a3c970e0ec2167de48c4b001eb", rownames(data)), grep("5m.AE", colnames(data))]), pch=18, type="b", col="blue", lwd=5, lty=2, xlab="SPOT Sampling Date", ylab="Relative abundance (16S+18S)", main="Rel. abundance of Chrysochromulina ASV3 (eaaf40a3c970e0ec2167de48c4b001eb) in AE 5m")
dev.off()


#Chrysochromulina
#Find out the right hash
unique(eLSA_CLR_5m$X[grep("d559", eLSA_CLR_5m$X_abr)])
unique(eLSA_CLR_5m$Y[grep("d559", eLSA_CLR_5m$Y_abr)])
plot(as.Date(dates_5m_AE), as.numeric(data[grep("d55992e6da65321a9b3c0ce3426e73ac", rownames(data)), grep("5m.AE", colnames(data))]), pch=18, type="b", col="blue3", lwd=5, lty=2, xlab="SPOT Sampling Date", ylab="Relative abundance (16S+18S)", main="Rel. abundance of Chrysochromulina ASV (d55992e6da65321a9b3c0ce3426e73ac) in AE 5m")

#That one prymnesiophyte, f925
unique(eLSA_CLR_5m$X[grep("f925", eLSA_CLR_5m$X_abr)])
plot(as.Date(dates_5m_AE), as.numeric(data[grep("ffa92571884bcbe0193ce4a1e4e843be", rownames(data)), grep("5m.AE", colnames(data))]), pch=18, type="b", col="maroon", lwd=5, lty=2, xlab="SPOT Sampling Date", ylab="Relative abundance (16S+18S)", main="Rel. abundance of Prymnesiophyte ASV (ffa92571884bcbe0193ce4a1e4e843be) in AE 5m")
#Doesn't look like it co-occurs with UCYN-A at all! 

#Chrysochromulina e403
unique(eLSA_CLR_5m$X[grep("e403", eLSA_CLR_5m$X_abr)])
plot(as.Date(dates_5m_AE), as.numeric(data[grep("eaaf40a3c970e0ec2167de48c4b001eb", rownames(data)), grep("5m.AE", colnames(data))]), type="b", pch=18, col="blue3", lwd=5, lty=2, xlab="SPOT Sampling Date", ylab="Relative abundance (16S+18S)", main="Rel. abundance of Chrysochromulina (eaaf40a3c970e0ec2167de48c4b001eb) in AE 5m")

#7. Compare absolute counts of UCYN-A and Brads-----
#Read in data 
proks_data <- read.table("ModifiedFiles/1.SPOT_16S_w.chloro_counts.tsv", sep=" ", header=T, stringsAsFactors = F, row.names=1) #Why is the sep not "\t"??? 
dim(proks_data)
length(grep("5m.AE", colnames(proks_data))) #125
length(grep("DCM.AE", colnames(proks_data))) #118

euks_data <- read.table("ModifiedFiles/6.SPOT_18S_no_metaz_counts.tsv", sep="\t", header=T, stringsAsFactors = F, row.names=1)
dim(euks_data)
length(grep("5m.AE", colnames(euks_data))) #124
length(grep("DCM.AE", colnames(euks_data))) #116

#Missing a few dates <- remove them?

#5m-----
colnames(proks_data)[grep("5m.AE", colnames(proks_data))]==colnames(euks_data)[grep("5m.AE", colnames(euks_data))]
colnames(proks_data)[grep("5m.AE", colnames(proks_data))][28:31]
colnames(euks_data)[grep("5m.AE", colnames(euks_data))][28:31]
#Missing 9.24.2009 in the euks data
proks_data[UCYNA_ix, grep("SPOT.2009.09.24.5m.AE", colnames(proks_data))] #Fortunately UCYN-A is absent on this date
colnames(proks_data)[grep("5m.AE", colnames(proks_data))][-29]==colnames(euks_data)[grep("5m.AE", colnames(euks_data))] #All TRUE! :) 
grep("FALSE", colnames(proks_data)[grep("5m.AE", colnames(proks_data))][-29]==colnames(euks_data)[grep("5m.AE", colnames(euks_data))])
which(grepl("FALSE", colnames(proks_data)[grep("5m.AE", colnames(proks_data))][-29]==colnames(euks_data)[grep("5m.AE", colnames(euks_data))])==TRUE) 
#Ok so for 5m, just ignore column #478 which is #29 on the grep("5m.AE", colnames())

#Create a spreadsheet
counts_5m=as.data.frame(matrix(ncol=(length(UCYNA_ix)+length(Brad_ix)+1), nrow=(length(grep("5m.AE", colnames(proks_data))[-29]))))
dim(counts_5m) #124 14
colnames(counts_5m)=c("SampleID", paste("UCYNA", c(1:6), sep="_"), "Brad_1", "Brad_5x", "Brad_2", "Brad_3", "Brad_4", "Brad_6x", "Brad_7x")

#Write in data
counts_5m$SampleID=colnames(proks_data)[grep("5m.AE", colnames(proks_data))][-29]
#Write in a column for dates
counts_dates <- c()
for(i in counts_5m$SampleID){
  a <- strsplit(i, split=".", fixed=T)
  b <- paste(a[[1]][2], a[[1]][3], a[[1]][4], sep="-")
  counts_dates <- c(counts_dates, b)
}
counts_dates <- as.Date(counts_dates, format="%Y-%m-%d")
print(head(counts_dates, n=10))
length(counts_dates)
counts_5m$Dates <- counts_dates
head(counts_5m$Dates)

#Add in numeric data
head(as.numeric(proks_data[grep(UCYNA_ix[1], row.names(proks_data)), grep("5m.AE", colnames(proks_data))]))
length(as.numeric(proks_data[grep(UCYNA_ix[1], row.names(proks_data)), grep("5m.AE", colnames(proks_data))])) #125
length(as.numeric(proks_data[grep(UCYNA_ix[1], row.names(proks_data)), grep("5m.AE", colnames(proks_data))[-29]])) #124

#UCYN-A
counts_5m$UCYNA_1=as.numeric(proks_data[grep(UCYNA_ix[1], row.names(proks_data)), grep("5m.AE", colnames(proks_data))[-29]])
counts_5m$UCYNA_2=as.numeric(proks_data[grep(UCYNA_ix[2], row.names(proks_data)), grep("5m.AE", colnames(proks_data))[-29]])
counts_5m$UCYNA_3=as.numeric(proks_data[grep(UCYNA_ix[3], row.names(proks_data)), grep("5m.AE", colnames(proks_data))[-29]])
counts_5m$UCYNA_4=as.numeric(proks_data[grep(UCYNA_ix[4], row.names(proks_data)), grep("5m.AE", colnames(proks_data))[-29]])
counts_5m$UCYNA_5=as.numeric(proks_data[grep(UCYNA_ix[5], row.names(proks_data)), grep("5m.AE", colnames(proks_data))[-29]])
counts_5m$UCYNA_6=as.numeric(proks_data[grep(UCYNA_ix[6], row.names(proks_data)), grep("5m.AE", colnames(proks_data))[-29]])

#Brad
counts_5m$Brad_1=as.numeric(euks_data[grep(Brad_ix[1], rownames(euks_data)), grep("5m.AE", colnames(euks_data))])
counts_5m$Brad_5x=as.numeric(euks_data[grep(Brad_ix[2], rownames(euks_data)), grep("5m.AE", colnames(euks_data))])
counts_5m$Brad_2=as.numeric(euks_data[grep(Brad_ix[3], rownames(euks_data)), grep("5m.AE", colnames(euks_data))])
counts_5m$Brad_3=as.numeric(euks_data[grep(Brad_ix[4], rownames(euks_data)), grep("5m.AE", colnames(euks_data))])
counts_5m$Brad_4=as.numeric(euks_data[grep(Brad_ix[5], rownames(euks_data)), grep("5m.AE", colnames(euks_data))])
counts_5m$Brad_6x=as.numeric(euks_data[grep(Brad_ix[6], rownames(euks_data)), grep("5m.AE", colnames(euks_data))])
counts_5m$Brad_7x=as.numeric(euks_data[grep(Brad_ix[7], rownames(euks_data)), grep("5m.AE", colnames(euks_data))])

#Calculate ratios 
counts_5m$A1_Brad7=counts_5m$UCYNA_1/counts_5m$Brad_7x
counts_5m$A5_Brad3=counts_5m$UCYNA_5/counts_5m$Brad_3

#Write it out!
#write.table(x=counts_5m, file="counts_UCYNA_Brad_ASVs_5m_02.06.2021.tsv", quote=F, sep="\t", row.names=F)

#Now plot ratios over time 
#5m: A1 and Brad7
#First, change ratios slightly 
#Both absent (true NA) -> change to be 8.00
length(which(counts_5m$UCYNA_1==0 & counts_5m$Brad_7x==0))
counts_5m$A1_Brad7[which(counts_5m$UCYNA_1==0 & counts_5m$Brad_7x==0)]=8.00
#Change "Inf" to be 9.00
length(which(counts_5m$UCYNA_1>0 & counts_5m$Brad_7x==0))
counts_5m$A1_Brad7[which(counts_5m$UCYNA_1>0 & counts_5m$Brad_7x==0)]=9.00
#Also change that one date where ratio was 100 for some reason to be 9.50 #This date was 3/12/2015
#counts_5m$A1_Brad7[which(counts_5m$A1_Brad7>100)]=9.50
counts_5m$A1_Brad7[grep("2015-03-12", counts_5m$Dates)]=9.50

#Now plot!
pdf("Plots/12.2021_Plots_time/Ratio_A1_Brad7_time_5m.pdf", width=7, height=3)
plot(as.Date(counts_5m$Dates), counts_5m$A1_Brad7,xlab="SPOT Sampling Date", ylab="Ratio of UCYN-A ASV1: \n Braarudospharea ASV7x", main="Ratio of UCYNA ASV1: Braarudosphaera ASV7x / time, 5m AE filters", pch=16)
abline(h=2.176)
abline(h=8.00, col="red")
abline(h=9.00, col="red")
abline(h=2.176+std.error(counts_5m$A1_Brad7), lty=2) #add in standard error #0.297
abline(h=2.176-std.error(counts_5m$A1_Brad7), lty=2) 
dev.off()

##std.error(counts_5m$A1_Brad7[which(counts_5m$A1_Brad7<7 & counts_5m$A1_Brad7>0)])

#Change ratios of A5/B3 and plot #Change with intention of editing in Affinity Design
#True NAs (both absent) -> 7.00 #In AD, change these dots to be red, and move down to zero
counts_5m$A5_Brad3[which(counts_5m$UCYNA_5==0 & counts_5m$Brad_3==0)]=7.00
#"Inf" -> change to be 8.00 #Delete these data points in AD
counts_5m$A5_Brad3[which(counts_5m$UCYNA_5>0 & counts_5m$Brad_3==0)]=8.00
#Ok now plot!
pdf("Plots/12.2021_Plots_time/Ratio_A5_Brad3_time_5m.pdf", height=3, width=7)
plot(as.Date(counts_5m$Dates), counts_5m$A5_Brad3, xlab="SPOT Sampling Date", ylab="Ratio of UCYN-A ASV\5: \n Braarudospharea ASV3", main="Ratio of UCYNA ASV5: Braarudosphaera ASV3 / time, 5m AE filters", pch=16)
abline(h=2.729)
abline(h=7.00, col="red")
abline(h=8.00, col="red")
abline(h=2.729+std.error(counts_5m$A5_Brad3), lty=2) #0.2438
abline(h=2.729-std.error(counts_5m$A5_Brad3), lty=2)
dev.off()


#DCM-----
length(grep("DCM.AE", colnames(proks_data)))
length(grep("DCM.AE", colnames(euks_data))) #116
which(colnames(proks_data)[grep("DCM.AE", colnames(proks_data))]==colnames(euks_data)[grep("DCM.AE", colnames(euks_data))]) #1:28 match, just as in 5m data - conicidence? It's 6.18.2009
colnames(proks_data)[grep("DCM.AE", colnames(proks_data))][28:35]
colnames(euks_data)[grep("DCM.AE", colnames(euks_data))][28:35]
colnames(proks_data)[grep("DCM.AE", colnames(proks_data))][30:length(grep("DCM.AE", colnames(proks_data)))]==colnames(euks_data)[grep("DCM.AE", colnames(euks_data))][29:length(grep("DCM.AE", colnames(euks_data)))] #F's up again 2 names later, really? #32, which is SPOT.2009.09.24.DCM.AE
colnames(proks_data)[grep("DCM.AE", colnames(proks_data))][-c(29, 32)]==colnames(euks_data)[grep("DCM.AE", colnames(euks_data))] #Looks good! 
grep("FALSE", colnames(proks_data)[grep("DCM.AE", colnames(proks_data))][-c(29, 32)]==colnames(euks_data)[grep("DCM.AE", colnames(euks_data))]) #Integer(0)! Figured out which one's to get rid of

#Set up a data frame
counts_DCM=as.data.frame(matrix(nrow=length(grep("DCM.AE", colnames(euks_data))), ncol=(length(UCYNA_ix)+ length(Brad_ix)+1)))
dim(counts_DCM) #116 14 - that's right! 
colnames(counts_DCM)=c("SampleID", paste("UCYNA", c(1:6), sep="_"), "Brad_1", "Brad_5x", "Brad_2", "Brad_3", "Brad_4", "Brad_6x", "Brad_7x")
counts_DCM$SampleID=colnames(proks_data)[grep("DCM.AE", colnames(proks_data))][-c(29, 32)]

#Add in a column for dates  
counts_dates <- c()
for(i in counts_DCM$SampleID){
  a <- strsplit(i, split=".", fixed=T)
  b <- paste(a[[1]][2], a[[1]][3], a[[1]][4], sep="-")
  counts_dates <- c(counts_dates, b)
}
length(counts_dates)
counts_dates <- as.Date(counts_dates, format="%Y-%m-%d")
counts_DCM$Dates=counts_dates
head(counts_DCM$Dates)

#Ok now write in numerical data
#UCYNA
counts_DCM$UCYNA_1=as.numeric(proks_data[grep(UCYNA_ix[1], row.names(proks_data)), grep("DCM.AE", colnames(proks_data))[-c(29, 32)]])
counts_DCM$UCYNA_2=as.numeric(proks_data[grep(UCYNA_ix[2], row.names(proks_data)), grep("DCM.AE", colnames(proks_data))[-c(29, 32)]])
counts_DCM$UCYNA_3=as.numeric(proks_data[grep(UCYNA_ix[3], row.names(proks_data)), grep("DCM.AE", colnames(proks_data))[-c(29, 32)]])
counts_DCM$UCYNA_4=as.numeric(proks_data[grep(UCYNA_ix[4], row.names(proks_data)), grep("DCM.AE", colnames(proks_data))[-c(29, 32)]])
counts_DCM$UCYNA_5=as.numeric(proks_data[grep(UCYNA_ix[5], row.names(proks_data)), grep("DCM.AE", colnames(proks_data))[-c(29, 32)]])
counts_DCM$UCYNA_6=as.numeric(proks_data[grep(UCYNA_ix[6], row.names(proks_data)), grep("DCM.AE", colnames(proks_data))[-c(29, 32)]])

#Brad
counts_DCM$Brad_1=as.numeric(euks_data[grep(Brad_ix[1], rownames(euks_data)), grep("DCM.AE", colnames(euks_data))])
counts_DCM$Brad_5x=as.numeric(euks_data[grep(Brad_ix[2], rownames(euks_data)), grep("DCM.AE", colnames(euks_data))])
counts_DCM$Brad_2=as.numeric(euks_data[grep(Brad_ix[3], rownames(euks_data)), grep("DCM.AE", colnames(euks_data))])
counts_DCM$Brad_3=as.numeric(euks_data[grep(Brad_ix[4], rownames(euks_data)), grep("DCM.AE", colnames(euks_data))])
counts_DCM$Brad_4=as.numeric(euks_data[grep(Brad_ix[5], rownames(euks_data)), grep("DCM.AE", colnames(euks_data))])
counts_DCM$Brad_6x=as.numeric(euks_data[grep(Brad_ix[6], rownames(euks_data)), grep("DCM.AE", colnames(euks_data))])
counts_DCM$Brad_7x=as.numeric(euks_data[grep(Brad_ix[7], rownames(euks_data)), grep("DCM.AE", colnames(euks_data))])

#Calculate ratios 
counts_DCM$A1_Brad7=counts_DCM$UCYNA_1/counts_DCM$Brad_7x
counts_DCM$A5_Brad3=counts_DCM$UCYNA_5/counts_DCM$Brad_3

#Write it out!
write.table(x=counts_DCM, file="counts_UCYNA_brad_ASVs_DCM_02.11.2021.tsv", sep="\t", quote=F, row.names=F)

#Change ratios and graph 
#A1 and Brad7
#Set NAs to be 3.00 #Again, these are true zeros, in Affinity Design, dye red and move down
length(which(counts_DCM$UCYNA_1==0 & counts_DCM$Brad_7x==0)) #95, or 81.90% of dates both absent
counts_DCM$A1_Brad7[which(counts_DCM$UCYNA_1==0 & counts_DCM$Brad_7x==0)]=3.00
#Set "infs" to be 2.50 #In AD, delete these points
length(which(counts_DCM$UCYNA_1>0 & counts_DCM$Brad_7x==0)) #9
counts_DCM$A1_Brad7[which(counts_DCM$UCYNA_1>0 & counts_DCM$Brad_7x==0)]=2.50

#Graph!
plot(as.Date(counts_DCM$Dates), as.numeric(counts_DCM$A1_Brad7), pch=16, xlab="SPOT Sampling Date", ylab="Ratio of UCYN-A ASV1: \n Braarudospharea ASV7x", main="Ratio of UCYNA ASV1: Braarudosphaera ASV7x / time, DCM AE filters")
abline(h=2.232)
abline(h=3.00, col="red")
abline(h=2.50, col="red")

#Moving on to A5 and Brad3
#Set Nas to be 1.500
length(which(counts_DCM$UCYNA_5==0 & counts_DCM$Brad_3==0)) #Both absent 108 times, or ~94% of days
counts_DCM$A5_Brad3[which(counts_DCM$UCYNA_5==0 & counts_DCM$Brad_3==0)]=1.50
#Set "Infs" to be 1.00
which(counts_DCM$UCYNA_5>0 & counts_DCM$Brad_3==0)
counts_DCM$A5_Brad3[which(counts_DCM$UCYNA_5>0 & counts_DCM$Brad_3==0)]=1.00

#NOW PLOT
plot(as.Date(counts_DCM$Date), as.numeric(counts_DCM$A5_Brad3), pch=16, xlab="SPOT Sampling Date", ylab="Ratio of UCYN-A ASV5: \n Braarudospharea ASV3", main="Ratio of UCYNA ASV5: Braarudosphaera ASV3 / time, DCM AE filters")
abline(h=1.407)
abline(h=1.50, col="red")
abline(h=1.00, col="red")

#Boxplots of average ratio of UCYN-A: Brad ASVs-----
#5m
#Read table back in 
ratio_5m=read.table("counts_UCYNA_Brad_ASVs_5m_02.06.2021.tsv", header=T, stringsAsFactors = F, sep="\t")
head(ratio_5m)
colnames(ratio_5m)
head(ratio_5m$A1_Brad7, n=20)

#Plot A1 vs. Brad7x, excluding infinities and zeros (NAs)
ratio_5m$A1_Brad7[which(ratio_5m$A1_Brad7>100)] #Why is this one ratio >100? #Look into that, exclude this date for now
pdf("Plots/03.16.22_Plots_time/boxplot_A1_vs_Brad_7x_5m_03.21.2022.pdf", height=3, width=2)
boxplot(x=ratio_5m$A1_Brad7[-which(ratio_5m$A1_Brad7>100)], ylab="Ratio of UCYN-A \n ASV1 16S:Braarudospharea\n ASV7x 18S", xlab="5m depth")
abline(h=1.991, col="red")
dev.off()
#Subset to exclude ratios >10
pdf("Plots/03.16.22_Plots_time/boxplot_A1_vs_Brad_7x_5m_subset_03.21.2022.pdf", height=3, width=2)
boxplot(x=ratio_5m$A1_Brad7[-which(ratio_5m$A1_Brad7>10)], ylab="Ratio of UCYN-A \n ASV1 16S:Braarudospharea \n ASV7x 18S", xlab="5m depth")
abline(h=1.991, col="red")
dev.off()

#Wait, why is this average apparently 1? The average ratio was 2? 
boxplot(ratio_5m$A1_Brad7[-c(which(ratio_5m$A1_Brad7>100), which(is.na(ratio_5m$A1_Brad7)==TRUE))]) #Ok, NAs not included 
mean(ratio_5m$A1_Brad7[-c(which(ratio_5m$A1_Brad7>100), which(is.na(ratio_5m$A1_Brad7)==TRUE))]) #Mean=1.99 #Add in an abline with this value

#Plot A5 vs. Brad 3
pdf("Plots/03.16.22_Plots_time/boxplot_A5_vs_Brad_3_5m_03.21.2022.pdf", height=3, width=2)
boxplot(ratio_5m$A5_Brad3[-c(which(ratio_5m$A5_Brad3>100), which(is.na(ratio_5m$A5_Brad3)==TRUE))], xlab="5m depth", ylab="Ratio of UCYN-A \n ASV5 16S:Braarudospharea ASV3 18S")
abline(h=2.139, col="red")
dev.off()
#Subset to exclude ratios >8
pdf("Plots/03.16.22_Plots_time/boxplot_A5_vs_Brad_3_subsetted_5m_03.21.2022.pdf", height=3, width=2)
boxplot(ratio_5m$A5_Brad3[-c(which(ratio_5m$A5_Brad3>8), which(is.na(ratio_5m$A5_Brad3)==TRUE))], xlab="5m depth", ylab="Ratio of UCYN-A \n ASV5 16S:Braarudospharea ASV3 18S")
abline(h=2.139, col="red")
dev.off()

#Figure out the average
mean(ratio_5m$A5_Brad3[-c(which(ratio_5m$A5_Brad3>100), which(is.na(ratio_5m$A5_Brad3)==TRUE))]) #2.139

#DCM----
#Read in data
ratio_DCM=read.table("counts_UCYNA_brad_ASVs_DCM_02.11.2021.tsv", header=T, stringsAsFactors = F, sep="\t")
dim(ratio_DCM) #116 17
colnames(ratio_DCM)

#A1 vs Brad 7
pdf("Plots/03.16.22_Plots_time/boxplot_A1_brad_7x_DCM_03.21.2022.pdf", width=2, height=3)
boxplot(ratio_DCM$A1_Brad7[-c(which(ratio_DCM$A1_Brad7>100),which(is.na(ratio_DCM$A1_Brad7)==TRUE))], xlab="DCM", ylab="Ratio of UCYN-A ASV1 16S: \n Braarodosphaera ASV 7x 18S")
abline(h=mean(ratio_DCM$A1_Brad7[-c(which(ratio_DCM$A1_Brad7>100),which(is.na(ratio_DCM$A1_Brad7)==TRUE))]), col="red") #1.674
dev.off()

#A5 vs Brad 3
pdf("Plots/03.16.22_Plots_time/boxplot_A5_brad_3_DCM_03.21.2022.pdf", width=2, height=3)
boxplot(ratio_DCM$A5_Brad3[-c(which(ratio_DCM$A5_Brad3>100),which(is.na(ratio_DCM$A5_Brad3)==TRUE))], xlab="DCM", ylab="Ratio of UCYN-A ASV5 16S: \n Braarodosphaera ASV 3 18S")
abline(h=mean(ratio_DCM$A5_Brad3[-c(which(ratio_DCM$A5_Brad3>100),which(is.na(ratio_DCM$A5_Brad3)==TRUE))]), col="red")
dev.off()






