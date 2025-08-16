#install.packages("BiocManager")
#BiocManager::install("Rsubread")
library(dplyr)
library(DESeq2)
library(edgeR)
library(limma)

#The "Script" R script makes the following files 

#All_Counts <- read.table("All_Counts.txt", sep="\t", header = TRUE,  row.names = 1)
#Brain_Counts <- read.table("Brain_Counts.txt", sep="\t", header=TRUE, row.names =1)
#Gonad_Counts <- read.table("Gonad_Counts.txt", sep="\t", header=TRUE, row.names = 1)
#Liver_Counts <- read.table("Liver_Counts.txt", sep="\t", header=TRUE, row.names = 1)

### DESeq Analysis ####

#### Separate into conditions ####

#HLD Brain
HLD.Brain.cols <- c("D_F0_M1_B", "L_F0_M1_B", "D_F0_M2_B", "D_F0_M3_B", 
                  "L_F0_M2_B", "L_F0_M3_B", "H_F0_M2_B", "H_F0_M1_B", 
                  "H_F0_M3_B", "H_F0_M4_B", "D_F0_M4_B")
H.L.D.Brain <- Brain_Counts[, HLD.Brain.cols]
rownames(H.L.D.Brain) <- rownames(Brain_Counts)

#Mciroplastics Brain
MP.MPD.D.Brain.cols <- c("MP_F0_M1_B", "D_F0_M1_B", "D_F0_M2_B", "MPD_F0_M1_B", "D_F0_M3_B", "MP_F0_M3_B", "MP_F0_M2_B", "MPD_F0_M2_B", "MPD_F0_M3_B", "MP_F0_M4_B", "D_F0_M4_B")
MP.MPD.D.Brain <- Brain_Counts[, MP.MPD.D.Brain.cols]
rownames(MP.MPD.D.Brain) <- rownames(Brain_Counts)

#Brain F1- High and DMSO
F1.Brain.cols <- c("H_F1_M3_B", "H_F1_M3_B.2.", "D_F1_M3_B", "D_F1_M1_B", "D_F1_M1_B.2.", "H_F1_M2_B")
F1.Brain <- Brain_Counts[,F1.Brain.cols]
rownames(F1.Brain) <-rownames(Brain_Counts)

#HLD Gonad
H.L.D.Gonad.cols <- c("D_F0_M2_G", "L_F0_M2_G", "L_F0_M3_G", "H_F0_M2_G", "H_F0_M1_G", "H_F0_M3_G", "D_F0_M3_G", "D_F0_M1_G", "L_F0_M1_G", "H_F0_M4_G", "D_F0_M4_G")
H.L.D.Gonad <-Gonad_Counts[,H.L.D.Gonad.cols]
rownames(H.L.D.Gonad) <- rownames(Gonad_Counts)

#Microplastics Gonad
MP.MPD.D.Gonad.cols <- c("D_F0_M2_G", "MP_F0_M1_G", "MPD_F0_M1_G", "MP_F0_M3_G", "MP_F0_M2_G", "MPD_F0_M2_G", "MPD_F0_M3_G", "D_F0_M3_G", "D_F0_M1_G", "MP_F0_M4_G", "D_F0_M4_G")
MP.MPD.D.Gonad <- Gonad_Counts[, MP.MPD.D.Gonad.cols]
rownames(MP.MPD.D.Gonad) <-rownames(Gonad_Counts)

#F1 Gonads 
F1.Gonad.cols <- c("D_F1_M1_G", "D_F1_M1_G.2.", "D_F1_M3_G", "H_F1_M3_G", "H_F1_M3_G.2.", "H_F1_M2_G")
F1.Gonad <- Gonad_Counts[, F1.Gonad.cols]
rownames(F1.Gonad) <- rownames(Gonad_Counts)

#Liver
#HLD Liver
H.L.D.Liver.cols <- c("L_F0_M1_L", "L_F0_M2_L", "L_F0_M3_L", "H_F0_M1_L", "D_F0_M3_L", "D_F0_M2_L", "D_F0_M1_L", "H_F0_M2_L", "H_F0_M3_L", "H_F0_M4_L", "D_F0_M4_L")
H.L.D.Liver <- Liver_Counts[, H.L.D.Liver.cols]
rownames(H.L.D.Liver) <- rownames(Liver_Counts)

#Microplastics Liver
MP.MPD.D.Liver.cols <- c("MP_F0_M2_L", "MP_F0_M1_L", "MPD_F0_M1_L", "D_F0_M3_L", "D_F0_M2_L", "D_F0_M1_L", "MP_F0_M3_L", "MPD_F0_M2_L", "MPD_F0_M3_L", "MP_F0_M4_L", "D_F0_M4_L")
MP.MPD.D.Liver <- Liver_Counts[, MP.MPD.D.Liver.cols]
rownames(MP.MPD.D.Liver) <- rownames(Liver_Counts)

#F1 Liver
F1.Liver.cols <- c("D_F1_M3_L", "H_F1_M3_L", "H_F1_M3_L.2.", "H_F1_M2_L", "D_F1_M1_L", "D_F1_M1_L.2.")
F1.Liver <- Liver_Counts[, F1.Liver.cols]
rownames(F1.Liver) <- rownames(Liver_Counts)

#### Brain High-Low- DMSO ####

#install.packages("BiocManager")  # if not already installed
#BiocManager::install(version = "3.21")  # or latest version for your R
#BiocManager::install()  # update all Bioconductor packages
#BiocManager::install("DESeq2")
library(DESeq2)

#HLD Brain 
HLD.Br_conds <- c("DMSO", "Low", "DMSO", "DMSO", "Low", "Low", "High", "High", "High", "High", "DMSO")

HLD_metadata <- data.frame(row.names = colnames(H.L.D.Brain),
                           conds = factor(HLD.Br_conds))
HLD.Br_conds_dds <- DESeqDataSetFromMatrix(countData = as.matrix(H.L.D.Brain),
                                           colData = HLD_metadata,
                                           design = ~ conds)

keep <- rowSums(cpm(counts(HLD.Br_conds_dds)) >= 1) >= 6
HLD.Br_conds_dds <- HLD.Br_conds_dds[keep,]

HLD.Br_conds_dds <- DESeq(HLD.Br_conds_dds)

#Pairwise
# High vs Low
Br.res_HL <- results(HLD.Br_conds_dds, contrast = c("conds", "High", "Low"))
# High vs DMSO
Br.res_HD <- results(HLD.Br_conds_dds, contrast = c("conds", "High", "DMSO"))
# Low vs DMSO
Br.res_LD <- results(HLD.Br_conds_dds, contrast = c("conds", "Low", "DMSO"))

#remove na values
Br.res_HL <- na.omit(Br.res_HL)
Br.res_HD <- na.omit(Br.res_HD)
Br.res_LD <- na.omit(Br.res_LD)
#significant p-values
Br.res_s_HL <- Br.res_HL[
  Br.res_HL$padj <= 0.05 & 
    (Br.res_HL$log2FoldChange > 1 | Br.res_HL$log2FoldChange < -1), ]
dim(Br.res_s_HL) #0 significant
Br.res_s_HD <- Br.res_HD[
  Br.res_HD$padj <= 0.05 & 
    (Br.res_HD$log2FoldChange > 1 | Br.res_HD$log2FoldChange < -1), ]
dim(Br.res_s_HD) #1 significant
Br.res_s_HD
Br.res_s_LD <- Br.res_LD[
  Br.res_LD$padj <= 0.05 & 
    (Br.res_LD$log2FoldChange > 1 | Br.res_LD$log2FoldChange < -1), ]
dim(Br.res_s_LD) # 1 significant
Br.res_s_LD

#PCA
plotPCA(rlog(HLD.Br_conds_dds, blind=TRUE), intgroup="conds")


#### Brain Microplastics ####

#conditions
MP.Br_conds <- c("MP", "DMSO", "DMSO", "MPD", "DMSO", "MP", "MP", "MPD", "MPD", "MP", "DMSO")

MP.Br_metadata <- data.frame(row.names = colnames(MP.MPD.D.Brain),
                           conds = factor(MP.Br_conds))
MP.Br_conds_dds <- DESeqDataSetFromMatrix(countData = as.matrix(MP.MPD.D.Brain),
                                           colData = MP.Br_metadata,
                                           design = ~ conds)

keep <- rowSums(cpm(counts(MP.Br_conds_dds)) >= 1) >= 6
MP.Br_conds_dds <- MP.Br_conds_dds[keep,]

MP.Br_conds_dds <- DESeq(MP.Br_conds_dds)

#Pairwise
#MP vs MPDEHP
Br.res_MP.MPD <- results(MP.Br_conds_dds, contrast = c("conds", "MP", "MPD"))
# MP vs DMSO
Br.res_MP.DMSO <- results(MP.Br_conds_dds, contrast = c("conds", "MP", "DMSO"))
# MPDEHP vs DMSO
Br.res_MPD.DMSO <- results(MP.Br_conds_dds, contrast = c("conds", "MPD", "DMSO"))

#remove na values
Br.res_MP.MPD <- na.omit(Br.res_MP.MPD)
Br.res_MP.DMSO <- na.omit(Br.res_MP.DMSO)
Br.res_MPD.DMSO <- na.omit(Br.res_MPD.DMSO)
#Significant p-values
Br.res_s_MP.MPD <- Br.res_MP.MPD[Br.res_MP.MPD$padj <= 0.05 & (Br.res_MP.MPD$log2FoldChange > 1 | Br.res_MP.MPD$log2FoldChange < -1), ]
dim(Br.res_s_MP.MPD) #0 significant genes
Br.res_s_MP.DMSO <- Br.res_MP.DMSO[Br.res_MP.DMSO$padj <= 0.05 & (Br.res_MP.DMSO$log2FoldChange > 1 | Br.res_MP.DMSO$log2FoldChange < -1), ]
dim(Br.res_s_MP.DMSO) #2 significant genes
Br.res_s_MP.DMSO
Br.res_s_MPD.DMSO <- Br.res_MPD.DMSO[Br.res_MPD.DMSO$padj <= 0.05 & (Br.res_MPD.DMSO$log2FoldChange >1 | Br.res_MPD.DMSO$log2FoldChange < -1), ]
dim(Br.res_s_MPD.DMSO) # 3 significant gene
Br.res_s_MPD.DMSO

#PCA
plotPCA(rlog(MP.Br_conds_dds, blind=TRUE), intgroup="conds")

#### Brain F1 ####

#conditions
F1.Br_conds <- c("High", "High", "DMSO", "DMSO", "DMSO", "High")

F1.Br_metadata <- data.frame(row.names = colnames(F1.Brain),
                             conds = factor(F1.Br_conds))
F1.Br_conds_dds <- DESeqDataSetFromMatrix(countData = as.matrix(F1.Brain),
                                          colData = F1.Br_metadata,
                                          design = ~ conds)

keep <- rowSums(cpm(counts(F1.Br_conds_dds)) >= 1) >= 6
F1.Br_conds_dds <- F1.Br_conds_dds[keep,]

F1.Br_conds_dds <- DESeq(F1.Br_conds_dds)

#results
Br.res_F1 <- results(F1.Br_conds_dds, contrast = c("conds", "High", "DMSO"))
#remove na values
Br.res_F1 <- na.omit(Br.res_F1)
#Significant p-values
Br.res_s_F1 <- Br.res_F1[
  Br.res_F1$padj <= 0.05 & 
    (Br.res_F1$log2FoldChange > 1 | Br.res_F1$log2FoldChange < -1), 
]
dim(Br.res_s_F1) #59 significant genes
#106 genes with a p-value less than 0.05
as.data.frame(tail(Br.res_s_F1, 50))


#PCA
plotPCA(rlog(F1.Br_conds_dds, blind=TRUE), intgroup="conds")

#### Gonad High-Low-DMSO ####

HLD.Go_conds <- c("DMSO", "Low", "Low", "High", "High", "High", "DMSO", "DMSO", "Low", "High", "DMSO")

HLD.Go_metadata <- data.frame(row.names = colnames(H.L.D.Gonad),
                             conds = factor(HLD.Go_conds))
HLD.Go_conds_dds <- DESeqDataSetFromMatrix(countData = as.matrix(H.L.D.Gonad),
                                          colData = HLD.Go_metadata,
                                          design = ~ conds)

keep <- rowSums(cpm(counts(HLD.Go_conds_dds)) >= 1) >= 6
HLD.Go_conds_dds <- HLD.Go_conds_dds[keep,]

HLD.Go_conds_dds <- DESeq(HLD.Go_conds_dds)

#Pairwise
#High vs Low
Go.res_HL <- results(HLD.Go_conds_dds, contrast = c("conds", "High", "Low"))
# High vs DMSO
Go.res_HD <- results(HLD.Go_conds_dds, contrast = c("conds", "High", "DMSO"))
# Low vs DMSO
Go.res_LD <- results(HLD.Go_conds_dds, contrast = c("conds", "Low", "DMSO"))

#remove na values
Go.res_HL <- na.omit(Go.res_HL)
Go.res_HD <- na.omit(Go.res_HD)
Go.res_LD <- na.omit(Go.res_LD)
#Significant p-values
Go.res_s_HL <- Go.res_HL[Go.res_HL$padj <= 0.05 & (Go.res_HL$log2FoldChange > 1 | Go.res_HL$log2FoldChange < -1),]
dim(Go.res_s_HL) #4507 significant genes
#5032 genes with a p-value less than 0.05
Go.res_s_HD <- Go.res_HD[Go.res_HD$padj <= 0.05 & (Go.res_HD$log2FoldChange > 1 | Go.res_HD$log2FoldChange < -1),]
dim(Go.res_s_HD) #4820 significant genes
#5184 genes with a p-value less than 0.05
Go.res_s_LD <- Go.res_LD[Go.res_LD$padj <= 0.05 & (Go.res_LD$log2FoldChange >1 | Go.res_LD$log2FoldChange < -1),]
dim(Go.res_s_LD) # 2 significant gene
Go.res_s_LD
# 3 genes with a p-value above 0.05

#PCA
plotPCA(rlog(HLD.Go_conds_dds, blind=TRUE), intgroup="conds")

#### Gonad Microplastics ####

#conditions
MP.Go_conds <- c("DMSO", "MP", "MPD", "MP", "MP", "MPD", "MPD", "DMSO", "DMSO", "MP", "DMSO")

MP.Go_metadata <- data.frame(row.names = colnames(MP.MPD.D.Gonad),
                             conds = factor(MP.Go_conds))
MP.Go_conds_dds <- DESeqDataSetFromMatrix(countData = as.matrix(MP.MPD.D.Gonad),
                                          colData = MP.Go_metadata,
                                          design = ~ conds)

keep <- rowSums(cpm(counts(MP.Go_conds_dds)) >= 1) >= 6
MP.Go_conds_dds <- MP.Go_conds_dds[keep,]

MP.Go_conds_dds <- DESeq(MP.Go_conds_dds)

#Pairwise
#MP vs MPDEHP
Go.res_MP.MPD <- results(MP.Go_conds_dds, contrast = c("conds", "MP", "MPD"))
# MP vs DMSO
Go.res_MP.DMSO <- results(MP.Go_conds_dds, contrast = c("conds", "MP", "DMSO"))
# MPDEHP vs DMSO
Go.res_MPD.DMSO <- results(MP.Go_conds_dds, contrast = c("conds", "MPD", "DMSO"))

#remove na values
Go.res_MP.MPD <- na.omit(Go.res_MP.MPD)
Go.res_MP.DMSO <- na.omit(Go.res_MP.DMSO)
Go.res_MPD.DMSO <- na.omit(Go.res_MPD.DMSO)
#Significant p-values
Go.res_s_MP.MPD <- Go.res_MP.MPD[Go.res_MP.MPD$padj <= 0.05 & (Go.res_MP.MPD$log2FoldChange > 1 | Go.res_MP.MPD$log2FoldChange < -1),]
dim(Go.res_s_MP.MPD) #1478 significant genes
#1517 genes with a p-value less than 0.05
Go.res_s_MP.DMSO <-Go.res_MP.DMSO[Go.res_MP.DMSO$padj <= 0.05 & (Go.res_MP.DMSO$log2FoldChange >1 | Go.res_MP.DMSO$log2FoldChange < -1), ]
dim(Go.res_s_MP.DMSO) #5388 significant genes
#5758 genes with a p-value less than 0.05
Go.res_s_MPD.DMSO <- Go.res_MPD.DMSO[Go.res_MPD.DMSO$padj <= 0.05 & (Go.res_MPD.DMSO$log2FoldChange > 1| Go.res_MPD.DMSO$log2FoldChange < -1), ]
dim(Go.res_s_MPD.DMSO) # 9 significant genes
Go.res_s_MPD.DMSO

#PCA
plotPCA(rlog(MP.Go_conds_dds, blind=TRUE), intgroup="conds")

#### Gonad F1 ####

#conditions
F1.Go_conds <- c("DMSO", "DMSO", "DMSO", "High", "High", "High")

F1.Go_metadata <- data.frame(row.names = colnames(F1.Gonad),
                             conds = factor(F1.Go_conds))
F1.Go_conds_dds <- DESeqDataSetFromMatrix(countData = as.matrix(F1.Gonad),
                                          colData = F1.Go_metadata,
                                          design = ~ conds)

keep <- rowSums(cpm(counts(F1.Go_conds_dds)) >= 1) >= 6
F1.Go_conds_dds <- F1.Go_conds_dds[keep,]

F1.Go_conds_dds <- DESeq(F1.Go_conds_dds)



#results
Go.res_F1 <- results(F1.Go_conds_dds, contrast = c("conds", "High", "DMSO"))
#remove na values
Go.res_F1 <- na.omit(Go.res_F1)
#Significant p-values
Go.res_s_F1 <- Go.res_F1[Go.res_F1$padj <= 0.05 & (Go.res_F1$log2FoldChange > 1 | Go.res_F1$log2FoldChange < -1), ]
dim(Go.res_s_F1) #60 significant genes
#71 genes with a p-value less than 0.05


#PCA
plotPCA(rlog(F1.Go_conds_dds, blind=TRUE), intgroup="conds")

#### Liver High-Low-DMSO ####

HLD.Li_conds <- c("Low", "Low", "Low", "High", "DMSO", "DMSO", "DMSO", "High", "High", "High", "DMSO")

HLD.Li_metadata <- data.frame(row.names = colnames(H.L.D.Liver),
                              conds = factor(HLD.Li_conds))
HLD.Li_conds_dds <- DESeqDataSetFromMatrix(countData = as.matrix(H.L.D.Liver),
                                           colData = HLD.Li_metadata,
                                           design = ~ conds)

keep <- rowSums(cpm(counts(HLD.Li_conds_dds)) >= 1) >= 6
HLD.Li_conds_dds <- HLD.Li_conds_dds[keep,]

HLD.Li_conds_dds <- DESeq(HLD.Li_conds_dds)

#Pairwise
#High vs Low
Li.res_HL <- results(HLD.Li_conds_dds, contrast = c("conds", "High", "Low"))
# High vs DMSO
Li.res_HD <- results(HLD.Li_conds_dds, contrast = c("conds", "High", "DMSO"))
# Low vs DMSO
Li.res_LD <- results(HLD.Li_conds_dds, contrast = c("conds", "Low", "DMSO"))

#remove na values
Li.res_HL <- na.omit(Li.res_HL)
Li.res_HD <- na.omit(Li.res_HD)
Li.res_LD <- na.omit(Li.res_LD)
#Significant p-values
Li.res_s_HL <- Li.res_HL[Li.res_HL$padj <= 0.05 & (Li.res_HL$log2FoldChange > 1 | Li.res_HL$log2FoldChange < -1), ]
dim(Li.res_s_HL) #4 significant genes
Li.res_s_HD <-Li.res_HD[Li.res_HD$padj <= 0.05 & (Li.res_HD$log2FoldChange > 1 | Li.res_HD$log2FoldChange < -1), ]
dim(Li.res_s_HD) #2 significant genes
Li.res_s_HD
Li.res_s_LD <- Li.res_LD[Li.res_LD$padj <= 0.05 & (Li.res_LD$log2FoldChange > 1 | Li.res_LD$log2FoldChange < -1), ]
dim(Li.res_s_LD) #414 significant genes
# 415 genes with a p-value less than 0.05


#PCA
plotPCA(rlog(HLD.Li_conds_dds, blind=TRUE), intgroup="conds")

#### Liver Microplastics ####

MP.Li_conds <- c("MP", "MP", "MPD", "DMSO", "DMSO", "DMSO", "MP", "MPD", "MPD", "MP", "DMSO")

MP.Li_metadata <- data.frame(row.names = colnames(MP.MPD.D.Liver),
                             conds = factor(MP.Li_conds))
MP.Li_conds_dds <- DESeqDataSetFromMatrix(countData = as.matrix(MP.MPD.D.Liver),
                                          colData = MP.Li_metadata,
                                          design = ~ conds)

keep <- rowSums(cpm(counts(MP.Li_conds_dds)) >= 1) >= 6
MP.Li_conds_dds <- MP.Li_conds_dds[keep,]

MP.Li_conds_dds <- DESeq(MP.Li_conds_dds)

#Pairwise
#MP vs MPDEHP
Li.res_MP.MPD <- results(MP.Li_conds_dds, contrast = c("conds", "MP", "MPD"))
# MP vs DMSO
Li.res_MP.DMSO <- results(MP.Li_conds_dds, contrast = c("conds", "MP", "DMSO"))
# MPDEHP vs DMSO
Li.res_MPD.DMSO <- results(MP.Li_conds_dds, contrast = c("conds", "MPD", "DMSO"))

#remove na values
Li.res_MP.MPD <- na.omit(Li.res_MP.MPD)
Li.res_MP.DMSO <- na.omit(Li.res_MP.DMSO)
Li.res_MPD.DMSO <- na.omit(Li.res_MPD.DMSO)
#Significant p-values
Li.res_s_MP.MPD <- Li.res_MP.MPD[Li.res_MP.MPD$padj <= 0.05 & (Li.res_MP.MPD$log2FoldChange > 1 | Li.res_MP.MPD$log2FoldChange < -1), ]
dim(Li.res_s_MP.MPD) #0 significant genes
Li.res_s_MP.DMSO <-Li.res_MP.DMSO[Li.res_MP.DMSO$padj <= 0.05 & (Li.res_MP.DMSO$log2FoldChange > 1 | Li.res_MP.DMSO$log2FoldChange < -1), ]
dim(Li.res_s_MP.DMSO) #0 significant genes
Li.res_s_MPD.DMSO <- Li.res_MPD.DMSO[Li.res_MPD.DMSO$padj <= 0.05 & (Li.res_MPD.DMSO$log2FoldChange >1 | Li.res_MPD.DMSO$log2FoldChange < -1), ]
dim(Li.res_s_MPD.DMSO) # 0 significant gene

#PCA
plotPCA(rlog(MP.Li_conds_dds, blind=TRUE), intgroup="conds")

#### F1 Liver ####

#conditions

F1.Li_conds <- c("DMSO", "High", "High", "High", "DMSO", "DMSO")

F1.Li_metadata <- data.frame(row.names = colnames(F1.Liver),
                             conds = factor(F1.Li_conds))
F1.Li_conds_dds <- DESeqDataSetFromMatrix(countData = as.matrix(F1.Liver),
                                          colData = F1.Li_metadata,
                                          design = ~ conds)

keep <- rowSums(cpm(counts(F1.Li_conds_dds)) >= 1) >= 6
F1.Li_conds_dds <- F1.Li_conds_dds[keep,]

F1.Li_conds_dds <- DESeq(F1.Li_conds_dds)

#results
Li.res_F1 <- results(F1.Li_conds_dds, contrast = c("conds", "High", "DMSO"))
#remove na values
Li.res_F1 <- na.omit(Li.res_F1)
#Significant p-values
Li.res_s_F1 <- Li.res_F1[Li.res_F1$padj <= 0.05 & (Li.res_F1$log2FoldChange > 1 | Li.res_F1$log2FoldChange < -1), ]
dim(Li.res_s_F1) #52 significant genes

#PCA
plotPCA(rlog(F1.Li_conds_dds, blind=TRUE), intgroup="conds")

#Remove metadata
rm(HLD.Li_metadata, HLD_metadata, HLD.Go_metadata, MP.Br_metadata, MP.Go_metadata, MP.Li_metadata, F1.Br_metadata, F1.Go_metadata, F1.Li_metadata)

#What are our genes of interest
head(Br.res_s_HD)
head(Br.res_s_LD)

head(Br.res_s_MP.DMSO)
head(Br.res_s_MPD.DMSO)

head(Br.res_s_F1, 10)


### Only keep results - Free space
#rm(HLD.Br_conds_dds, MP.Br_conds_dds, F1.Br_conds_dds, HLD.Go_conds_dds,  MP.Go_conds_dds, F1.Go_conds_dds, HLD.Li_conds_dds, MP.Li_conds_dds, F1.Li_conds_dds, keep)

### Save the significant genes 
write.csv(as.data.frame(Br.res_s_HD), file = "Significant genes/Br_res_s_HD.csv")
write.csv(as.data.frame(Br.res_s_LD), file = "Significant genes/Br_res_s_LD.csv")
write.csv(as.data.frame(Br.res_s_HL), file = "Significant genes/Br_res_s_HL.csv")
write.csv(as.data.frame(Br.res_s_MP.DMSO), file = "Significant genes/Br.res_s_MP.DMSO.csv")
write.csv(as.data.frame(Br.res_s_MP.MPD), file = "Significant genes/Br.res_s_MP.MPD.csv")
write.csv(as.data.frame(Br.res_s_MPD.DMSO), file = "Significant genes/Br.res_s_MPD.DMSO.csv")
write.csv(as.data.frame(Br.res_s_F1), file = "Significant genes/Br.res_s_F1.csv")

write.csv(as.data.frame(Go.res_s_HD), file = "Significant genes/Go.res_s_HD.csv")
write.csv(as.data.frame(Go.res_s_LD), file = "Significant genes/Go.res_s_LD.csv")
write.csv(as.data.frame(Go.res_s_HL), file = "Significant genes/Go.res_s_HL.csv")
write.csv(as.data.frame(Go.res_s_MP.DMSO), file = "Significant genes/Go.res_s_MP.DMSO.csv")
write.csv(as.data.frame(Go.res_s_MP.MPD), file = "Significant genes/Go.res_s_MP.MPD.csv")
write.csv(as.data.frame(Go.res_s_MPD.DMSO), file = "Significant genes/Go.res_s_MPD.DMSO.csv")
write.csv(as.data.frame(Go.res_s_F1), file = "Significant genes/Go.res_s_F1.csv")

write.csv(as.data.frame(Li.res_s_HD), file = "Significant genes/Li.res_s_HD.csv")
write.csv(as.data.frame(Li.res_s_LD), file = "Significant genes/Li.res_s_LD.csv")
write.csv(as.data.frame(Li.res_s_HL), file = "Significant genes/Li.res_s_HL.csv")
write.csv(as.data.frame(Li.res_s_MP.DMSO), file = "Significant genes/Li.res_s_MP.DMSO.csv")
write.csv(as.data.frame(Li.res_s_MP.MPD), file = "Significant genes/Li.res_s_MP.MPD.csv")
write.csv(as.data.frame(Li.res_s_MPD.DMSO), file = "Significant genes/Li.res_s_MPD.DMSO.csv")
write.csv(as.data.frame(Li.res_s_F1), file = "Significant genes/Li.res_s_F1.csv")

#All genes

# Brain data
write.csv(as.data.frame(Br.res_HD), file = "All genes/Br.res_HD.csv")
write.csv(as.data.frame(Br.res_LD), file = "All genes/Br_res_LD.csv")
write.csv(as.data.frame(Br.res_HL), file = "All genes/Br.res_HL.csv")
write.csv(as.data.frame(Br.res_MP.DMSO), file = "All genes/Br.res_MP.DMSO.csv")
write.csv(as.data.frame(Br.res_MP.MPD), file = "All genes/Br.res_MP.MPD.csv")
write.csv(as.data.frame(Br.res_MPD.DMSO), file = "All genes/Br.res_MPD.DMSO.csv")
write.csv(as.data.frame(Br.res_F1), file = "All genes/Br.res_F1.csv")

# Gonad data
write.csv(as.data.frame(Go.res_HD), file = "All genes/Go.res_HD.csv")
write.csv(as.data.frame(Go.res_LD), file = "All genes/Go.res_LD.csv")
write.csv(as.data.frame(Go.res_HL), file = "All genes/Go.res_HL.csv")
write.csv(as.data.frame(Go.res_MP.DMSO), file = "All genes/Go.res_MP.DMSO.csv")
write.csv(as.data.frame(Go.res_MP.MPD), file = "All genes/Go.res_MP.MPD.csv")
write.csv(as.data.frame(Go.res_MPD.DMSO), file = "All genes/Go.res_MPD.DMSO.csv")
write.csv(as.data.frame(Go.res_F1), file = "All genes/Go.res_F1.csv")

# Liver data
write.csv(as.data.frame(Li.res_HD), file = "All genes/Li.res_HD.csv")
write.csv(as.data.frame(Li.res_LD), file = "All genes/Li.res_LD.csv")
write.csv(as.data.frame(Li.res_HL), file = "All genes/Li.res_HL.csv")
write.csv(as.data.frame(Li.res_MP.DMSO), file = "All genes/Li.res_MP.DMSO.csv")
write.csv(as.data.frame(Li.res_MP.MPD), file = "All genes/Li.res_MP.MPD.csv")
write.csv(as.data.frame(Li.res_MPD.DMSO), file = "All genes/Li.res_MPD.DMSO.csv")
write.csv(as.data.frame(Li.res_F1), file = "All genes/Li.res_F1.csv")

#UGH

write.table(rownames(Br.res_F1), file = "All genes text/Br.res_F1.txt",quote = FALSE, row.names = FALSE, col.names = FALSE)

write.table(rownames(Go.res_HD), file = "All genes text/Go.res_HD.txt",quote = FALSE, row.names = FALSE, col.names = FALSE)

write.table(rownames(Go.res_F1), file = "All genes text/Go.res_F1.txt",quote = FALSE, row.names = FALSE, col.names = FALSE)

write.table(rownames(Li.res_LD), file = "All genes text/Li.res_LD.txt",quote = FALSE, row.names = FALSE, col.names = FALSE)

write.table(rownames(Go.res_MP.DMSO), file = "All genes text/Go.res_MP.DMSO.txt",quote = FALSE, row.names = FALSE, col.names = FALSE)
