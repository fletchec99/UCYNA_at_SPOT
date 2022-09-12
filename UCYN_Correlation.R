## 20220601: http://www.sthda.com/english/wiki/correlation-matrix-a-quick-start-guide-to-analyze-format-and-visualize-a-correlation-matrix-using-r-software
library(corrplot)
library(Hmisc)
ucyn <- read.table("/Users/yubinraut/Documents/Fuhrman_Lab/UCYN/Final/Raw_Data/int_env_data_UCYNA_relabun_pres_abs_05.26.2022.tsv", sep = " ", header = T)
ucyn.clr <- read.table("/Users/yubinraut/Documents/Fuhrman_Lab/UCYN/Final/Raw_Data/corr.df_07.14.2022.tsv", sep = " ", header = T)

## brad abundance
brad <- read.table("/Users/yubinraut/Documents/Fuhrman_Lab/UCYN/Final/Raw_Data/relabun_UCYNA_Brad_ASVs_01.12.2021.tsv", sep = "\t", header = T)
ucyn$Year <- sub("_.*", "", ucyn$Year_Month)
ucyn$Month <- sprintf("%02s",sub("^[^_]*_", "", ucyn$Year_Month))
ucyn$Month <- sprintf("%02s", ucyn$Month)
ucyn$Date <- paste(ucyn$Year, ucyn$Month, sep = "-")
brad$Date2 <- sub("(.*-.*)-(.*)", "\\1", brad$Date)
ucyn$Brad3_abundance <- brad$Brad_3[match(ucyn$Date, brad$Date2)]
ucyn$Brad7_abundance <- brad$Brad_7x[match(ucyn$Date, brad$Date2)]
## select variables to do correlation against
cor.df <- ucyn[,c("NO2_NO3","PO4","BactProd_Leu",
   "BactProd_Thy","MODIS_Chl","MODIS_SST","CUTI_ix",
   "BEUTI_ix","MEI_ix","UCYN_ASV1_relabun","UCYNA_ASV5_relabun", "Brad3_abundance", "Brad7_abundance")]
colnames(cor.df) <- c("In-situ NOx", "In-situ PO4", "Bacterial Production (Leu)",
                        "Bacterial Production (Thy)", "MODIS Chl-a", "MODIS SST",
                         "CUTI", "BEUTI", "MEI",
                        "UCYN-A (I) ASV-1 Abundance", "UCYN-A (II) ASV-5 Abundance", "Braarudosphaera ASV-3", "Braarudosphaera ASV-7" )
## New df Colette provided with clr transformed abundances and same environmental data
cor.df2 <- ucyn.clr[,c("X.NO2.NO3...uM.","X.PO4...uM.",
                       "Bacterial.production..Leucine...cells.mL.day.",
                       "MODIS..Chl...mg.m.3.","MODIS.SST..ºC.","CUTI",
                       "BEUTI","MEI","UCYN.A1.Abundance",
                       "UCYN.A2.Abundance","B..bigelowii.ASV3.Abundance",
                       "Braarudosphaera.7x.Abundance")]
colnames(cor.df2) <- c("In-situ NOx", "In-situ PO4", "Bacterial Production (Leu)",
                      "MODIS Chl-a", "MODIS SST","CUTI", "BEUTI", "MEI",
                      "UCYN-A (I) ASV-1 Abundance", "UCYN-A (II) ASV-5 Abundance", 
                      "Braarudosphaera ASV-3", "Braarudosphaera ASV-7" )
cor.df3 <- ucyn.clr[,c("X.NO2.NO3...uM.","X.PO4...uM.",
                       "Bacterial.production..Leucine...cells.mL.day.",
                       "MODIS..Chl...mg.m.3.","MODIS.SST..ºC.","CUTI",
                       "BEUTI","MEI","Ekman.Transport.N.S..kg.m.s.",
                       "Ekman.Transport.E.W..kg.m.s.","Surface.wind..m.s.","UCYN.A1.Abundance",
                       "UCYN.A2.Abundance","B..bigelowii.ASV3.Abundance",
                       "Braarudosphaera.7x.Abundance")]
colnames(cor.df3) <- c("In-situ NOx", "In-situ PO4", "Bacterial Production (Leu)",
                       "MODIS Chl-a", "MODIS SST","CUTI", "BEUTI", "MEI",
                       "Ekman NS", "Ekman EW", "Surface Wind",
                       "UCYN-A (I) ASV-1 Abundance", "UCYN-A (II) ASV-5 Abundance", 
                       "Braarudosphaera ASV-3", "Braarudosphaera ASV-7" )

## reorder variables so it's alpha and remove buoy SST and estimated nitrate

## We are using Spearman's rank correlation - nonparametric measure of the correlation that uses rank of observations
## in its calculations rather than the original numeric values. This measures the
## MONOTONIC relationship between two variables: + = if Y tends to increase as X increases
## - is if Y tends to decrease as X decreases. *** Makes no assumptions about distribution of
## data..

## Computing significance levels and correlation coefficients
ucyn.cor <- rcorr(as.matrix(cor.df), type = "spearman")
ucyn.cor2 <- rcorr(as.matrix(cor.df2), type = "spearman")
ucyn.cor3 <- rcorr(as.matrix(cor.df3), type = "spearman")
## Plot with only correlations with p-value < 0.05 displayed
setwd("/Users/yubinraut/Documents/Fuhrman_Lab/UCYN/Final/Plots/")
pdf(file = "UCYN_Clr_Bonus_Correlation.pdf")
corrplot(ucyn.cor3$r, type = "upper", order = "alphabet",tl.col = "black",
         p.mat = ucyn.cor3$P, sig.level = 0.05, insig = "blank")
dev.off()
pdf(file = "UCYN_Clr_Correlation.pdf")
corrplot(ucyn.cor2$r, type = "upper", order = "alphabet",tl.col = "black",
         p.mat = ucyn.cor2$P, sig.level = 0.05, insig = "blank")
dev.off()
pdf(file = "UCYN_Correlation.pdf")
corrplot(ucyn.cor$r, type = "upper", order = "alphabet",tl.col = "black",
         p.mat = ucyn.cor$P, sig.level = 0.05, insig = "blank")
dev.off()

## Table of coefficients and p-values
flattenCorrMatrix <- function(cormat, pmat) {
  ut <- upper.tri(cormat)
  data.frame(
    row = rownames(cormat)[row(cormat)[ut]],
    column = rownames(cormat)[col(cormat)[ut]],
    cor  =(cormat)[ut],
    p = pmat[ut]
  )
}
all.cor.p <- flattenCorrMatrix(ucyn.cor$r, ucyn.cor$P)
all.cor.p$significance <- ifelse(all.cor.p$p < 0.05, "P-value < 0.05", "Not Significant")
write.csv(all.cor.p, "/Users/yubinraut/Documents/Fuhrman_Lab/UCYN/Final/UCYN_Correlation_Coefficients_Pvalue.csv", row.names = F)

# Spearman's rank-order correlation, a non-parametric alternative to Pearson's correlation, was used to determine the monotonic relationship between environmental variables and UCYN-A relative abundances.
# This was done using the rcorr() function with type = "spearman" from the "Hmisc" package in R. There's a dataframe of the correlation coefficients alongside the associated p-values.
# The correlogram was produced using the function corrplot() from the "corrplot" package in R with correlation coefficient values excluded if the p-values < 0.05.