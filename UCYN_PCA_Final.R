## 20220828
library(ggfortify)
library(cowplot)
## read in data provided by Colette for presence of UCYN-A and host
ucyn <- read.csv("/Users/yubinraut/Documents/Fuhrman_Lab/UCYN/Final/Raw_Data/20220425_int_env_data_UCYNA_Brad_pres_abs_cutoff_0.0001_04.18.2022.csv")
## CUTI (2), BEUTI (3), MEI (4), MODIS_Chl (14), MODIS_SST (15)
## Model_NO2 (16), Model_SST (17)
## BactProd_Leu (21), BactProd_Thy (22), PO4 (23),NO2+NO3 (24)
ucyn.transform <- na.omit(ucyn[,c(2:4,14:17,21:24)])
pca <- prcomp(ucyn.transform, center = T, scale. = T)
## UCYN-A (I) ASV-1 and Host (B7)
pca.a1b7 <- autoplot(pca, data = na.omit(ucyn[,c(2:4,14:17,21:30)]), 
         size = 3, colour = "black", fill = "A1_abund", shape = "Brad7_pres", loadings = T, 
         loadings.label = F, loadings.colour = "black")+
  theme_classic()+
  theme(legend.position = c(0.88,0.2))+
  scale_fill_manual(values = c("#ffffff","#969696","#000000"), 
                    name = "UCYN-A (I) ASV-1 Abundance",  labels = c("Absent", "Low", "High"))+
  scale_shape_manual(values = c(22,24), name = expression(paste(italic("Braarudosphaera "), "ASV-7 Abundance")))+
  guides(fill = guide_legend(override.aes = list(shape = c(21,21,21))))+
  annotate("text", label = "SST", x = 0.3, y = 0.04, size = 4)+
  annotate("text", label = "MEI", x = 0.127, y = 0.15, size = 4)+
  annotate("text", label = "CUTI", x = -0.09, y = 0.337, size = 4)+
  annotate("text", label = "BEUTI", x = -0.225, y = 0.24, size = 4)+
  annotate("text", label = "Bacterial\n Production", x = -0.21, y = 0.1, size = 4)+
  annotate("text", label = "Estimated\n Nitrate", x = -0.225, y = -0.145, size = 4)+
  annotate("text", label = "Chl-a", x = -0.13, y = -0.25, size = 4)+
  annotate("text", label = "In-situ\n NOx", x = -0.08, y = -0.205, size = 4)+
  annotate("text", label = "In-situ\n Phosphate", x = -0.025, y = -0.25, size = 4)
ggsave("A1_B7_Abundance_Environmental_PCA_Scale.pdf", plot=pca.a1b7, path="/Users/yubinraut/Documents/Fuhrman_Lab/UCYN/Final/Plots/", width=12, height=6)

## UCYN-A (II) ASV-5 & Host (B3)
pca.a2b3 <- autoplot(pca, data = na.omit(ucyn[,c(2:4,14:17,21:30)]), 
                     size = 3, colour = "black", fill = "A2_abund", shape = "Brad3_pres", loadings = T, 
                     loadings.label = F, loadings.colour = "black")+
  theme_classic()+
  theme(legend.position = c(0.88,0.2))+
  scale_fill_manual(values = c("#ffffff","#969696","#000000"), 
                    name = "UCYN-A (II) ASV-5 Abundance",  labels = c("Absent", "Low", "High"))+
  scale_shape_manual(values = c(22,24), name = expression(paste(italic("Braarudosphaera "), "ASV-3 Abundance")))+
  guides(fill = guide_legend(override.aes = list(shape = c(21,21,21))))+
  annotate("text", label = "SST", x = 0.3, y = 0.04, size = 4)+
  annotate("text", label = "MEI", x = 0.127, y = 0.15, size = 4)+
  annotate("text", label = "CUTI", x = -0.09, y = 0.337, size = 4)+
  annotate("text", label = "BEUTI", x = -0.225, y = 0.24, size = 4)+
  annotate("text", label = "Bacterial\n Production", x = -0.21, y = 0.1, size = 4)+
  annotate("text", label = "Estimated\n Nitrate", x = -0.225, y = -0.145, size = 4)+
  annotate("text", label = "Chl-a", x = -0.13, y = -0.25, size = 4)+
  annotate("text", label = "In-situ\n NOx", x = -0.08, y = -0.205, size = 4)+
  annotate("text", label = "In-situ\n Phosphate", x = -0.025, y = -0.25, size = 4)
ggsave("A2_B3_Abundance_Environmental_PCA_Scale.pdf", plot=pca.a2b3, path="/Users/yubinraut/Documents/Fuhrman_Lab/UCYN/Final/Plots/", width=12, height=6)

## Plot both side by side
p1 <- autoplot(pca, data = na.omit(ucyn[,c(2:4,14:17,21:30)]), 
               size = 3, colour = "black", fill = "A1_abund", shape = "Brad7_pres", loadings = T, 
               loadings.label = F, loadings.colour = "black")+
  theme_classic()+
  theme(legend.position = c(0.88,0.15))+
  scale_fill_manual(values = c("#ffffff","#969696","#000000"), 
                    name = "UCYN-A Abundance",  labels = c("Absent", "Low", "High"))+
  scale_shape_manual(values = c(22,24), guide = "none")+
  guides(fill = guide_legend(override.aes = list(shape = c(21,21,21))))+
  annotate("text", label = "SST", x = 0.31, y = 0.04, size = 3)+
  annotate("text", label = "MEI", x = 0.135, y = 0.155, size = 3)+
  annotate("text", label = "CUTI", x = -0.09, y = 0.337, size = 3)+
  annotate("text", label = "BEUTI", x = -0.22, y = 0.25, size = 3)+
  annotate("text", label = "Bacterial\n Production", x = -0.23, y = 0.1, size = 3)+
  annotate("text", label = "Estimated\n Nitrate", x = -0.235, y = -0.145, size = 3)+
  annotate("text", label = "Chl-a", x = -0.13, y = -0.25, size = 3)+
  annotate("text", label = "In-situ\n NOx", x = -0.09, y = -0.20, size = 3)+
  annotate("text", label = "In-situ\n PO4", x = -0.025, y = -0.24, size = 3)+
  annotate("text", label = expression(paste("UCYN-A (I) ASV-1 ", italic("Braarudosphaera "), "ASV-7")), 
           x = 0.05, y = 0.4, size = 4)
## 
p2 <- autoplot(pca, data = na.omit(ucyn[,c(2:4,14:17,21:30)]), 
               size = 3, colour = "black", fill = "A2_abund", shape = "Brad3_pres", loadings = T, 
               loadings.label = F, loadings.colour = "black")+
  theme_classic()+
  theme(legend.position = c(0.82,0.15))+
  scale_fill_manual(values = c("#ffffff","#969696","#000000"), 
                    guide = "none")+
  scale_shape_manual(values = c(22,24), name = expression(paste(italic("Braarudosphaera "), "Abundance")))+
  annotate("text", label = "SST", x = 0.31, y = 0.04, size = 3)+
  annotate("text", label = "MEI", x = 0.135, y = 0.155, size = 3)+
  annotate("text", label = "CUTI", x = -0.09, y = 0.337, size = 3)+
  annotate("text", label = "BEUTI", x = -0.22, y = 0.25, size = 3)+
  annotate("text", label = "Bacterial\n Production", x = -0.23, y = 0.1, size = 3)+
  annotate("text", label = "Estimated\n Nitrate", x = -0.235, y = -0.145, size = 3)+
  annotate("text", label = "Chl-a", x = -0.13, y = -0.25, size = 3)+
  annotate("text", label = "In-situ\n NOx", x = -0.09, y = -0.20, size = 3)+
  annotate("text", label = "In-situ\n PO4", x = -0.025, y = -0.24, size = 3)+
  annotate("text", label = expression(paste("UCYN-A (II) ASV-5 ", italic("Braarudosphaera "), "ASV-3")), 
           x = 0.05, y = 0.4, size = 4)
pca.both <- plot_grid(p1, p2)
ggsave("UCYN-A_Host_Both_Abundance_Environmental_PCA_Scale.pdf", plot=pca.both, path="/Users/yubinraut/Documents/Fuhrman_Lab/UCYN/Final/Plots/", width=12, height=6)

