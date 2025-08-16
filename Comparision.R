library(gplots)
library(dplyr)
#BiocManager::install("ggVennDiagram")
library(ggVennDiagram)
library(ggplot2)


## This script is to compare between DEG sets to find what is shared and what is not

#Starting within Pairwise Groups

### Brain High-Low-DMSO ####

#In Brain High vs DMSO we had one significant DEG, this was ttc6, which was upregulated 5.2x in the High DEHP exposure.
head(Br.res_s_HD)

#In Brain Low vs DMSO we had one significant DEG, this was hyal3, which was upregulated 2.6x in the Low DEHP exposure
head(Br.res_s_LD)

### Brain Microplastics ####

#In Brain Microplastics control vs DMSO we had 2 significant DEGs, these were histh1l and si:dkey-32n7.4.
#histh1l was downregulated 1.2x in Microplastics
# si:dkey-32n7.4. was upregulated 2.8x in Microplastics
Br.res_s_MP.DMSO


#In Brain Microplastics + DEHP vs DMSO control we had 3 significant DEGs, these were si:dkey-32n7.4, krt5 and krt5_1
#si:dkey-32n7.4 which was upregualted 2.8x  in Microplastics + DEHP
# krt5 and krt5_1 are the same feature/gene, it is downregulated 2.0x in Microplastics +DEHP
Br.res_s_MPD.DMSO


### Gonad High_Low-DMSO ####

#In Gonad High DEHP vs Low DEHP we had 5032 significant DEGs.

Go.HD_sig_ordered <- Go.res_HD[Go.res_HD$padj <= 0.05 & abs(Go.res_HD$log2FoldChange) > 1, ]
Go.HD_sig_ordered <- Go.res_s_HD[order(-abs(Go.res_s_HD$log2FoldChange), Go.res_s_HD$log2FoldChange), ]
View(Go.HD_sig_ordered)
Go.HD_sig_ordered <- as.data.frame(Go.HD_sig_ordered)

#In Gonad High DEHP vs DMSO control we had 5184 significant DEGs

#Because there are so many genes we will make a venn diagram to see the overlap
Go_HLD_SL <- list(HighvLow=rownames(Go.res_s_HL), HighvDMSO=rownames(Go.res_s_HD), LowvDMSO=rownames(Go.res_s_LD)) 

ggVennDiagram(Go_HLD_SL, label_alpha = 0, label = "count") +
  scale_fill_gradient(low = "white", high = "pink")+
  theme(legend.position = "none") + 
  theme(text = element_text(size = 10)) 
  
  
  
#Between the HighvsLow and HighvsDMSO we found that 3935 DEG were differentially expressed in both
#1249 DEGs unique to HighvsDMSO
#1097 DEGs unique to HighvsLow


shared_HL <- Go.res_s_HL[rownames(Go.res_s_HL) %in% genes_shared, ]
shared_HD <- Go.res_s_HD[rownames(Go.res_s_HD) %in% genes_shared, ]

# Unique to High vs Low
unique_HL <- Go.res_s_HL[rownames(Go.res_s_HL) %in% genes_unique_HL, ]

# Unique to High vs DMSO
unique_HD <- Go.res_s_HD[rownames(Go.res_s_HD) %in% genes_unique_HD, ]

#In Gonad Low DEHP vs DMSO control we had 3 DEGs. These were acbd5b, kifap3b and ttpa
#acbd5b was upregulated 1.5x in High
#kifap3b was downregulated 0.84x in High
#ttpa was downregulated 2.17x in High
Go.res_s_LD




### Gonad Microplastics ####

#In Gonad Microplastics control vs Microplastic+DEHP there were 1517 DEGs

#In Gonad Microplastics control vs DMSO control there were 5758 DEGs

#In Gonad Microplastics+DEHP vs DMSO control there were 9 DEGs
Go.res_s_MPD.DMSO
#apof was upregulated 6.8x in in Microplastics+DEHP
#fam126b was upregulated 6.2x 
#ddc was upregulated 6.7x
#prox1a was upregulated 7.5x
#arhgef10lb was upregulated 5.4x
#ptdss2 was upregulated 3.8x
#apof_1 was upregulated 6.8x
#got1_1 was upregulated 6.6x
#arhgef10lb_1 was upregulated 5.4x


Go_MP_SL <- list(MPvMPDEHP=rownames(Go.res_s_MP.MPD), MPvDMSO=rownames(Go.res_s_MP.DMSO), MPDEHPvDMSO=rownames(Go.res_s_MPD.DMSO)) 
ggVennDiagram(Go_MP_SL, label_alpha = 0, label = "count") +
  scale_fill_gradient(low = "white", high = "pink")+
  theme(legend.position = "none") + 
  theme(text = element_text(size = 10)) 


#95 DEGs unique to MPvMPD - may be related to the DEHP side - is there overlap between HIGH DEHP DEGs
#1422 in both MPvDMSO and MPvMPD - shared response to microplastics
#4327 unique to MPvDMSO - ????
#9 Shared between MPvDMSO and MPDvDMSO - ????


### Liver High-Low-DMSO ####

#In Liver High vs Low there were 4 significant DEGs
#uchl1 was upregulated 6.9x in High
#si:dkey-14d8.6 was downregulated 8.5x in High
#scg3 and scg3_1 are the same feature/genes and was upregulated 7.7x

#In Liver High vs DMSO there were 2 significant DEGs
#si:dkey-1h4.4 was upregulated 4.3x in High
#LOC110439129 was upregulated 4.4x in High

#In liver Low vs DMSO there were 415 significant DEGs

Li.LD_sig_ordered <- Li.res_s_LD[order(-abs(Li.res_s_LD$log2FoldChange), Li.res_s_LD$log2FoldChange), ]
Li.LD_sig_ordered <- as.data.frame(Li.LD_sig_ordered)
View(Li.LD_sig_ordered)

Li_HLD_SL <- list(HighvLow=rownames(Li.res_s_HL), HighvDMSO=rownames(Li.res_s_HD), LowvDMSO=rownames(Li.res_s_LD))

ggVennDiagram(Li_HLD_SL, label_alpha = 0, label = "count") +
  scale_fill_gradient(low = "white", high = "pink")+
  theme(legend.position = "none") + 
  theme(text = element_text(size = 10)) 


#410 unique DEGs in Low vs DMSO
#


### Brain F1 vs Brain F0 ####


Br_F1_SL <- list(F1=rownames(Br.res_s_F1), F0=rownames(Br.res_s_HD))

ggVennDiagram(Br_F1_SL, label_alpha = 0, label = "count") +
  scale_fill_gradient(low = "white", high = "gray")+
  theme(legend.position = "none") + 
  theme(text = element_text(size = 10))

### Gonad F1 vs Gonad F0 ####


Go_F1_SL <- list(F1=rownames(Go.res_s_F1), F0=rownames(Go.res_s_HD))

ggVennDiagram(Go_F1_SL, label_alpha = 0, label = "count") +
  scale_fill_gradient(low = "white", high = "pink")+
  theme(legend.position = "none") + 
  theme(text = element_text(size = 10)) 

#Overlap is 9 genes
F1F0_shared <- intersect(rownames(Go.res_s_F1), rownames(Go.res_s_HD))

GoF1F0_overlap <- Go.res_s_F1[rownames(Go.res_s_F1) %in% F1F0_shared, ]

head(GoF1F0_overlap, 10)

### Liver F1 vs Liver F0 ####


Li_F1_SL <- list(F1=rownames(Li.res_s_F1), F0=rownames(Li.res_s_HD))

ggVennDiagram(Li_F1_SL, label_alpha = 0, label = "count") +
  scale_fill_gradient(low = "white", high = "pink")+
  theme(legend.position = "none") + 
  theme(text = element_text(size = 10)) 

#Out of curiosity 
### Low vs DMSO compared to F1

Test_Li_F1_SL <- list(F1=rownames(Li.res_s_F1), F0=rownames(Li.res_s_LD))
venn(Test_Li_F1_SL)

#How many genes are upregulated

# BF1
upregulated_BF1 <- Br.res_s_F1[Br.res_s_F1$log2FoldChange > 0, ]
r_upregulated_BF1 <- upregulated_BF1
# Calculate ranking metric
r_upregulated_BF1$ranking_metric <- -log10(r_upregulated_BF1$padj + 1e-10) * sign(r_upregulated_BF1$log2FoldChange)
# Sort by ranking
r_upregulated_BF1 <- r_upregulated_BF1[order(-r_upregulated_BF1$ranking_metric), ]
# Extract gene IDs
r_upregulated_BF1 <- rownames(r_upregulated_BF1)

write.table(r_upregulated_BF1, file = "Significant genes/Upregulated/r_upregulated_BF1.txt", 
            quote = FALSE, row.names = FALSE, col.names = FALSE)

Downregulated_BF1 <- Br.res_s_F1[Br.res_s_F1$log2FoldChange < 0, ]
r_Downregulated_BF1 <- Downregulated_BF1
# Calculate ranking metric
r_Downregulated_BF1$ranking_metric <- -log10(r_Downregulated_BF1$padj + 1e-10) * sign(r_Downregulated_BF1$log2FoldChange)
# Sort by ranking
r_Downregulated_BF1 <- r_Downregulated_BF1[order(-r_Downregulated_BF1$ranking_metric), ]
# Extract gene IDs
r_Downregulated_BF1 <- rownames(r_Downregulated_BF1)

write.table(r_Downregulated_BF1, file = "Significant genes/Downregulated/r_Downregulated_BF1.txt", 
            quote = FALSE, row.names = FALSE, col.names = FALSE)


# HD
upregulated_HD <- Go.res_s_HD[Go.res_s_HD$log2FoldChange > 0, ]
r_upregulated_HD <-upregulated_HD

# Calculate ranking metric
r_upregulated_HD$ranking_metric <- -log10(r_upregulated_HD$padj + 1e-10) * sign(r_upregulated_HD$log2FoldChange)
# Sort by ranking
r_upregulated_HD <- r_upregulated_HD[order(-r_upregulated_HD$ranking_metric), ]
# Extract gene IDs
r_upregulated_HD <- rownames(r_upregulated_HD)

write.table(r_upregulated_HD, file = "Significant genes/Upregulated/r_upregulated_HD.txt", 
            quote = FALSE, row.names = FALSE, col.names = FALSE)


Downregulated_HD <- Go.res_s_HD[Go.res_s_HD$log2FoldChange < 0, ]
r_Downregulated_HD <- Downregulated_HD

# Calculate ranking metric
r_Downregulated_HD$ranking_metric <- -log10(r_Downregulated_HD$padj + 1e-10) * sign(r_Downregulated_HD$log2FoldChange)
# Sort by ranking
r_Downregulated_HD <- r_Downregulated_HD[order(-r_Downregulated_HD$ranking_metric), ]
# Extract gene IDs
r_Downregulated_HD <- rownames(r_Downregulated_HD)

write.table(r_Downregulated_HD, file = "Significant genes/Downregulated/r_Downregulated_HD.txt", 
            quote = FALSE, row.names = FALSE, col.names = FALSE)


# GF1
upregulated_GF1 <- Go.res_s_F1[Go.res_s_F1$log2FoldChange > 0, ]
r_upregulated_GF1 <- upregulated_GF1

# Calculate ranking metric
r_upregulated_GF1$ranking_metric <- -log10(r_upregulated_GF1$padj + 1e-10) * sign(r_upregulated_GF1$log2FoldChange)
# Sort by ranking
r_upregulated_GF1 <- r_upregulated_GF1[order(-r_upregulated_GF1$ranking_metric), ]
# Extract gene IDs
r_upregulated_GF1 <- rownames(r_upregulated_GF1)

write.table(r_upregulated_GF1, file = "Significant genes/Upregulated/r_upregulated_GF1.txt", 
            quote = FALSE, row.names = FALSE, col.names = FALSE)


Downregulated_GF1 <- Go.res_s_F1[Go.res_s_F1$log2FoldChange < 0, ]
r_Downregulated_GF1 <- Downregulated_GF1

# Calculate ranking metric
r_Downregulated_GF1$ranking_metric <- -log10(r_Downregulated_GF1$padj + 1e-10) * sign(r_Downregulated_GF1$log2FoldChange)
# Sort by ranking
r_Downregulated_GF1 <- r_Downregulated_GF1[order(-r_Downregulated_GF1$ranking_metric), ]
# Extract gene IDs
r_Downregulated_GF1 <- rownames(r_Downregulated_GF1)

write.table(r_Downregulated_GF1, file = "Significant genes/Downregulated/r_Downregulated_GF1.txt", 
            quote = FALSE, row.names = FALSE, col.names = FALSE)

# LLD
upregulated_LLD <- Li.res_s_LD[Li.res_s_LD$log2FoldChange > 0, ]
r_upregulated_LLD <- upregulated_LLD

# Calculate ranking metric
r_upregulated_LLD$ranking_metric <- -log10(r_upregulated_LLD$padj + 1e-10) * sign(r_upregulated_LLD$log2FoldChange)
# Sort by ranking
r_upregulated_LLD <- r_upregulated_LLD[order(-r_upregulated_LLD$ranking_metric), ]
# Extract gene IDs
r_upregulated_LLD <- rownames(r_upregulated_LLD)

write.table(r_upregulated_LLD, file = "Significant genes/Upregulated/r_upregulated_LLD.txt", 
            quote = FALSE, row.names = FALSE, col.names = FALSE)


Downregulated_LLD <- Li.res_s_LD[Li.res_s_LD$log2FoldChange < 0, ]
r_Downregulated_LLD <- Downregulated_LLD

# Calculate ranking metric
r_Downregulated_LLD$ranking_metric <- -log10(r_Downregulated_LLD$padj + 1e-10) * sign(r_Downregulated_LLD$log2FoldChange)
# Sort by ranking
r_Downregulated_LLD <- r_Downregulated_LLD[order(-r_Downregulated_LLD$ranking_metric), ]
# Extract gene IDs
r_Downregulated_LLD <- rownames(r_Downregulated_LLD)

write.table(r_Downregulated_LLD, file = "Significant genes/Downregulated/r_Downregulated_LLD.txt", 
            quote = FALSE, row.names = FALSE, col.names = FALSE)

# LF1 
upregulated_LF1 <- Li.res_s_F1[Li.res_s_F1$log2FoldChange > 0, ]
r_upregulated_LF1 <- upregulated_LF1

# Calculate ranking metric
r_upregulated_LF1$ranking_metric <- -log10(r_upregulated_LF1$padj + 1e-10) * sign(r_upregulated_LF1$log2FoldChange)
# Sort by ranking
r_upregulated_LF1 <- r_upregulated_LF1[order(-r_upregulated_LF1$ranking_metric), ]
# Extract gene IDs
r_upregulated_LF1 <- rownames(r_upregulated_LF1)

write.table(r_upregulated_LF1, file = "Significant genes/Upregulated/r_upregulated_LF1.txt", 
            quote = FALSE, row.names = FALSE, col.names = FALSE)

Downregulated_LF1 <- Li.res_s_F1[Li.res_s_F1$log2FoldChange < 0, ]
r_Downregulated_LF1 <- Downregulated_LF1

# Calculate ranking metric
r_Downregulated_LF1$ranking_metric <- -log10(r_Downregulated_LF1$padj + 1e-10) * sign(r_Downregulated_LF1$log2FoldChange)
# Sort by ranking
r_Downregulated_LF1 <- r_Downregulated_LF1[order(-r_Downregulated_LF1$ranking_metric), ]
# Extract gene IDs
r_Downregulated_LF1 <- rownames(r_Downregulated_LF1)

write.table(r_Downregulated_LF1, file = "Significant genes/Downregulated/r_Downregulated_LF1.txt", 
            quote = FALSE, row.names = FALSE, col.names = FALSE)

# GMPDM
upregulated_GMPDM <- Go.res_s_MP.DMSO[Go.res_s_MP.DMSO$log2FoldChange > 0, ]
r_upregulated_GMPDM <- upregulated_GMPDM

# Calculate ranking metric
r_upregulated_GMPDM$ranking_metric <- -log10(r_upregulated_GMPDM$padj + 1e-10) * sign(r_upregulated_GMPDM$log2FoldChange)
# Sort by ranking
r_upregulated_GMPDM <- r_upregulated_GMPDM[order(-r_upregulated_GMPDM$ranking_metric), ]
# Extract gene IDs
r_upregulated_GMPDM <- rownames(r_upregulated_GMPDM)

write.table(r_upregulated_GMPDM, file = "Significant genes/Upregulated/r_upregulated_GMPDM.txt", 
            quote = FALSE, row.names = FALSE, col.names = FALSE)

Downregulated_GMPDM <- Go.res_s_MP.DMSO[Go.res_s_MP.DMSO$log2FoldChange < 0, ]
r_Downregulated_GMPDM <- Downregulated_GMPDM

# Calculate ranking metric
r_Downregulated_GMPDM$ranking_metric <- -log10(r_Downregulated_GMPDM$padj + 1e-10) * sign(r_Downregulated_GMPDM$log2FoldChange)
# Sort by ranking
r_Downregulated_GMPDM <- r_Downregulated_GMPDM[order(-r_Downregulated_GMPDM$ranking_metric), ]
# Extract gene IDs
r_Downregulated_GMPDM <- rownames(r_Downregulated_GMPDM)

write.table(r_Downregulated_GMPDM, file = "Significant genes/Downregulated/r_Downregulated_GMPDM.txt", 
            quote = FALSE, row.names = FALSE, col.names = FALSE)


