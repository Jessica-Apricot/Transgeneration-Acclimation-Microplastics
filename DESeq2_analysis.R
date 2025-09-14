#install.packages("BiocManager")
#BiocManager::install("Rsubread")
#install.packages("ggfortify")
library(dplyr)
library(DESeq2)
library(edgeR)
library(limma)
library(ggplot2)
library(plotly)
library(ggfortify)


#The "Script" R script makes the following files 

#All_Counts <- read.table("All_Counts.txt", sep="\t", header = TRUE,  row.names = 1)
#Brain_Counts <- read.table("Brain_Counts.txt", sep="\t", header=TRUE, row.names =1)
#Gonad_Counts <- read.table("Gonad_Counts.txt", sep="\t", header=TRUE, row.names = 1)
#Liver_Counts <- read.table("Liver_Counts.txt", sep="\t", header=TRUE, row.names = 1)

### DESeq Analysis ####

#### Separate into conditions ####

#HLD Brain
HLD_Brain_cols <- c("D_F0_M1_B", "L_F0_M1_B", "D_F0_M2_B", "D_F0_M3_B", 
                  "L_F0_M2_B", "L_F0_M3_B", "H_F0_M2_B", "H_F0_M1_B", 
                  "H_F0_M3_B", "H_F0_M4_B", "D_F0_M4_B")
H_L_D_Brain <- Brain_Counts[, HLD_Brain_cols]
rownames(H_L_D_Brain) <- rownames(Brain_Counts)

#Mciroplastics Brain
MP_MPD_D_Brain_cols <- c("MP_F0_M1_B", "D_F0_M1_B", "D_F0_M2_B", "MPD_F0_M1_B", "D_F0_M3_B", "MP_F0_M3_B", "MP_F0_M2_B", "MPD_F0_M2_B", "MPD_F0_M3_B", "MP_F0_M4_B", "D_F0_M4_B")
MP_MPD_D_Brain <- Brain_Counts[, MP_MPD_D_Brain_cols]
rownames(MP_MPD_D_Brain) <- rownames(Brain_Counts)

#Brain F1- High and DMSO
F1_Brain_cols <- c("H_F1_M3_B", "H_F1_M3_B.2.", "D_F1_M3_B", "D_F1_M1_B", "D_F1_M1_B.2.", "H_F1_M2_B")
F1_Brain <- Brain_Counts[,F1_Brain_cols]
rownames(F1_Brain) <-rownames(Brain_Counts)

#HLD Gonad
H_L_D_Gonad_cols <- c("D_F0_M2_G", "L_F0_M2_G", "L_F0_M3_G", "H_F0_M2_G", "H_F0_M1_G", "H_F0_M3_G", "D_F0_M3_G", "D_F0_M1_G", "L_F0_M1_G", "H_F0_M4_G", "D_F0_M4_G")
H_L_D_Gonad <-Gonad_Counts[,H_L_D_Gonad_cols]
rownames(H_L_D_Gonad) <- rownames(Gonad_Counts)

#Microplastics Gonad
MP_MPD_D_Gonad_cols <- c("D_F0_M2_G", "MP_F0_M1_G", "MPD_F0_M1_G", "MP_F0_M3_G", "MP_F0_M2_G", "MPD_F0_M2_G", "MPD_F0_M3_G", "D_F0_M3_G", "D_F0_M1_G", "MP_F0_M4_G", "D_F0_M4_G")
MP_MPD_D_Gonad <- Gonad_Counts[, MP_MPD_D_Gonad_cols]
rownames(MP_MPD_D_Gonad) <-rownames(Gonad_Counts)

#F1 Gonads 
F1_Gonad_cols <- c("D_F1_M1_G", "D_F1_M1_G.2.", "D_F1_M3_G", "H_F1_M3_G", "H_F1_M3_G.2.", "H_F1_M2_G")
F1_Gonad <- Gonad_Counts[, F1_Gonad_cols]
rownames(F1_Gonad) <- rownames(Gonad_Counts)

#Liver
#HLD Liver
H_L_D_Liver_cols <- c("L_F0_M1_L", "L_F0_M2_L", "L_F0_M3_L", "H_F0_M1_L", "D_F0_M3_L", "D_F0_M2_L", "D_F0_M1_L", "H_F0_M2_L", "H_F0_M3_L", "H_F0_M4_L", "D_F0_M4_L")
H_L_D_Liver <- Liver_Counts[, H_L_D_Liver_cols]
rownames(H_L_D_Liver) <- rownames(Liver_Counts)

#Microplastics Liver
MP_MPD_D_Liver_cols <- c("MP_F0_M2_L", "MP_F0_M1_L", "MPD_F0_M1_L", "D_F0_M3_L", "D_F0_M2_L", "D_F0_M1_L", "MP_F0_M3_L", "MPD_F0_M2_L", "MPD_F0_M3_L", "MP_F0_M4_L", "D_F0_M4_L")
MP_MPD_D_Liver <- Liver_Counts[, MP_MPD_D_Liver_cols]
rownames(MP_MPD_D_Liver) <- rownames(Liver_Counts)

#F1 Liver
F1_Liver_cols <- c("D_F1_M3_L", "H_F1_M3_L", "H_F1_M3_L.2.", "H_F1_M2_L", "D_F1_M1_L", "D_F1_M1_L.2.")
F1_Liver <- Liver_Counts[, F1_Liver_cols]
rownames(F1_Liver) <- rownames(Liver_Counts)

#### Brain High-Low- DMSO ####

#install.packages("BiocManager")  # if not already installed
#BiocManager::install(version = "3.21")  # or latest version for your R
#BiocManager::install()  # update all Bioconductor packages
#BiocManager::install("DESeq2")
library(DESeq2)

#HLD Brain 
HLD_Br_conds <- c("DMSO", "Low", "DMSO", "DMSO", "Low", "Low", "High", "High", "High", "High", "DMSO")

HLD_metadata <- data.frame(row.names = colnames(H_L_D_Brain),
                           conds = factor(HLD_Br_conds))
HLD_Br_conds_dds <- DESeqDataSetFromMatrix(countData = as.matrix(H_L_D_Brain),
                                           colData = HLD_metadata,
                                           design = ~ conds)

keep <- rowSums(cpm(counts(HLD_Br_conds_dds)) >= 1) >= 6
HLD_Br_conds_dds <- HLD_Br_conds_dds[keep,]

HLD_Br_conds_dds <- DESeq(HLD_Br_conds_dds)

#Pairwise
# High vs Low
Br_res_HL <- results(HLD_Br_conds_dds, contrast = c("conds", "High", "Low"))
# High vs DMSO
Br_res_HD <- results(HLD_Br_conds_dds, contrast = c("conds", "High", "DMSO"))
# Low vs DMSO
Br_res_LD <- results(HLD_Br_conds_dds, contrast = c("conds", "Low", "DMSO"))

#remove na values
Br_res_HL <- na.omit(Br_res_HL)
Br_res_HD <- na.omit(Br_res_HD)
Br_res_LD <- na.omit(Br_res_LD)
#significant p-values
Br_res_s_HL <- Br_res_HL[
  Br_res_HL$padj <= 0.05 & 
    (Br_res_HL$log2FoldChange > 1 | Br_res_HL$log2FoldChange < -1), ]
dim(Br_res_s_HL) #0 significant
Br_res_s_HD <- Br_res_HD[
  Br_res_HD$padj <= 0.05 & 
    (Br_res_HD$log2FoldChange > 1 | Br_res_HD$log2FoldChange < -1), ]
dim(Br_res_s_HD) #1 significant
Br_res_s_HD
Br_res_s_LD <- Br_res_LD[
  Br_res_LD$padj <= 0.05 & 
    (Br_res_LD$log2FoldChange > 1 | Br_res_LD$log2FoldChange < -1), ]
dim(Br_res_s_LD) # 1 significant
Br_res_s_LD

#PCA
#plotPCA(rlog(HLD_Br_conds_dds, blind=TRUE), intgroup="conds")

HLD_Br_PCA <- rlog(HLD_Br_conds_dds, blind = TRUE) 
HLD_Br_rv <- rowVars(assay(HLD_Br_PCA)) 
HLD_Br_top_500 <- order(HLD_Br_rv, decreasing = TRUE)[1:500] 

HLD_Br_top_matrix <- assay(HLD_Br_PCA)[HLD_Br_top_500, ] 

HLD_Br_pca_res <- prcomp(t(HLD_Br_top_matrix), scale. = TRUE)

HLD_Br_metadata <- as.data.frame(colData(HLD_Br_conds_dds)) 

HLD_Br_pca_df <- as.data.frame(HLD_Br_pca_res$x)

HLD_Br_pca_df$conds <- HLD_Br_metadata$conds 

percentVar <- round(100 * (HLD_Br_pca_res$sdev)^2 / sum(HLD_Br_pca_res$sdev^2), 1)

#with jitter 
set.seed(9)
HLD_Br_pca_df$jitter_PC1 <- HLD_Br_pca_df$PC1 + runif(nrow(HLD_Br_pca_df), -0.3, 0.3)
HLD_Br_pca_df$jitter_PC2 <- HLD_Br_pca_df$PC2 + runif(nrow(HLD_Br_pca_df), -0.3, 0.3)

ggplot(HLD_Br_pca_df, aes(x = jitter_PC1, y = jitter_PC2, color = conds)) +
  geom_point(size = 3, alpha = 0.8) +
  theme_minimal(base_size = 14) +
  labs(title = "",
       x = paste0("PC1: ", percentVar[1], "%"),
       y = paste0("PC2: ", percentVar[2], "%"))



#Without jitter
ggplot(HLD_Br_pca_df, aes(x = PC1, y = PC2, color = conds)) +
  geom_point(size = 3, alpha = 0.8) +
  theme_minimal(base_size = 14) +
  labs(title = "",
       x = paste0("PC1: ", percentVar[1], "%"),
       y = paste0("PC2: ", percentVar[2], "%")) 



#### Brain Microplastics ####

#conditions
MP_Br_conds <- c("MP", "DMSO", "DMSO", "MPD", "DMSO", "MP", "MP", "MPD", "MPD", "MP", "DMSO")

MP_Br_metadata <- data.frame(row.names = colnames(MP_MPD_D_Brain),
                           conds = factor(MP_Br_conds))
MP_Br_conds_dds <- DESeqDataSetFromMatrix(countData = as.matrix(MP_MPD_D_Brain),
                                           colData = MP_Br_metadata,
                                           design = ~ conds)

keep <- rowSums(cpm(counts(MP_Br_conds_dds)) >= 1) >= 6
MP_Br_conds_dds <- MP_Br_conds_dds[keep,]

MP_Br_conds_dds <- DESeq(MP_Br_conds_dds)

#Pairwise
#MP vs MPDEHP
Br_res_MP_MPD <- results(MP_Br_conds_dds, contrast = c("conds", "MP", "MPD"))
# MP vs DMSO
Br_res_MP_DMSO <- results(MP_Br_conds_dds, contrast = c("conds", "MP", "DMSO"))
# MPDEHP vs DMSO
Br_res_MPD_DMSO <- results(MP_Br_conds_dds, contrast = c("conds", "MPD", "DMSO"))

#remove na values
Br_res_MP_MPD <- na.omit(Br_res_MP_MPD)
Br_res_MP_DMSO <- na.omit(Br_res_MP_DMSO)
Br_res_MPD_DMSO <- na.omit(Br_res_MPD_DMSO)
#Significant p-values
Br_res_s_MP_MPD <- Br_res_MP_MPD[Br_res_MP_MPD$padj <= 0.05 & (Br_res_MP_MPD$log2FoldChange > 1 | Br_res_MP_MPD$log2FoldChange < -1), ]
dim(Br_res_s_MP_MPD) #0 significant genes
Br_res_s_MP_DMSO <- Br_res_MP_DMSO[Br_res_MP_DMSO$padj <= 0.05 & (Br_res_MP_DMSO$log2FoldChange > 1 | Br_res_MP_DMSO$log2FoldChange < -1), ]
dim(Br_res_s_MP_DMSO) #2 significant genes
Br_res_s_MP_DMSO
Br_res_s_MPD_DMSO <- Br_res_MPD_DMSO[Br_res_MPD_DMSO$padj <= 0.05 & (Br_res_MPD_DMSO$log2FoldChange >1 | Br_res_MPD_DMSO$log2FoldChange < -1), ]
dim(Br_res_s_MPD_DMSO) # 3 significant gene
Br_res_s_MPD_DMSO

#PCA
#plotPCA(rlog(MP_Br_conds_dds, blind=TRUE), intgroup="conds")

MP_Br_PCA <- rlog(MP_Br_conds_dds, blind = TRUE)

# Step 2: Extract top 500 most variable genes
MP_Br_rv <- rowVars(assay(MP_Br_PCA))
MP_Br_top_500 <- order(MP_Br_rv, decreasing = TRUE)[1:500]
MP_Br_top_matrix <- assay(MP_Br_PCA)[MP_Br_top_500, ]

MP_Br_pca_res <- prcomp(t(MP_Br_top_matrix), scale. = TRUE)


MP_Br_metadata <- as.data.frame(colData(MP_Br_conds_dds))

MP_Br_pca_df <- as.data.frame(MP_Br_pca_res$x)

MP_Br_pca_df$conds <- MP_Br_metadata$conds 

percentVar <- round(100 * (MP_Br_pca_res$sdev)^2 / sum(MP_Br_pca_res$sdev^2), 1)

#with jitter 

set.seed(9)
MP_Br_pca_df$jitter_PC1 <- MP_Br_pca_df$PC1 + runif(nrow(MP_Br_pca_df), -0.3, 0.3)
MP_Br_pca_df$jitter_PC2 <- MP_Br_pca_df$PC2 + runif(nrow(MP_Br_pca_df), -0.3, 0.3)

ggplot(MP_Br_pca_df, aes(x = jitter_PC1, y = jitter_PC2, color = conds)) +
  geom_point(size = 3, alpha = 0.8) +
  theme_minimal(base_size = 14) +
  labs(title = "",
       x = paste0("PC1: ", percentVar[1], "%"),
       y = paste0("PC2: ", percentVar[2], "%"))


#Without jitter
ggplot(MP_Br_pca_df, aes(x = PC1, y = PC2, color = conds)) +
  geom_point(size = 3, alpha = 0.8) +
  theme_minimal(base_size = 14) +
  labs(title = "",
       x = paste0("PC1: ", percentVar[1], "%"),
       y = paste0("PC2: ", percentVar[2], "%")) 


#### Brain F1 ####

#conditions
F1_Br_conds <- c("High", "High", "DMSO", "DMSO", "DMSO", "High")

F1_Br_metadata <- data.frame(row.names = colnames(F1_Brain),
                             conds = factor(F1_Br_conds))
F1_Br_conds_dds <- DESeqDataSetFromMatrix(countData = as.matrix(F1_Brain),
                                          colData = F1_Br_metadata,
                                          design = ~ conds)

keep <- rowSums(cpm(counts(F1_Br_conds_dds)) >= 1) >= 6
F1_Br_conds_dds <- F1_Br_conds_dds[keep,]

F1_Br_conds_dds <- DESeq(F1_Br_conds_dds)

#results
Br_res_F1 <- results(F1_Br_conds_dds, contrast = c("conds", "High", "DMSO"))
#remove na values
Br_res_F1 <- na.omit(Br_res_F1)
#Significant p-values
Br_res_s_F1 <- Br_res_F1[
  Br_res_F1$padj <= 0.05 & 
    (Br_res_F1$log2FoldChange > 1 | Br_res_F1$log2FoldChange < -1), 
]
dim(Br_res_s_F1) #59 significant genes
#106 genes with a p-value less than 0.05
as.data.frame(tail(Br_res_s_F1, 50))


#PCA
#plotPCA(rlog(F1_Br_conds_dds, blind=TRUE), intgroup="conds")

F1_Br_PCA <- rlog(F1_Br_conds_dds, blind = TRUE)

# Step 2: Extract top 500 most variable genes
F1_Br_rv <- rowVars(assay(F1_Br_PCA))
F1_Br_top_500 <- order(F1_Br_rv, decreasing = TRUE)[1:500]
F1_Br_top_matrix <- assay(F1_Br_PCA)[F1_Br_top_500, ]

F1_Br_pca_res <- prcomp(t(F1_Br_top_matrix), scale. = TRUE)

F1_Br_metadata <- as.data.frame(colData(F1_Br_conds_dds))

F1_Br_pca_df <- as.data.frame(F1_Br_pca_res$x)

F1_Br_pca_df$conds <- F1_Br_metadata$conds 

percentVar <- round(100 * (F1_Br_pca_res$sdev)^2 / sum(F1_Br_pca_res$sdev^2), 1)

#with jitter 
set.seed(9)
F1_Br_pca_df$jitter_PC1 <- F1_Br_pca_df$PC1 + runif(nrow(F1_Br_pca_df), -0.3, 0.3)
F1_Br_pca_df$jitter_PC2 <- F1_Br_pca_df$PC2 + runif(nrow(F1_Br_pca_df), -0.3, 0.3)

ggplot(F1_Br_pca_df, aes(x = jitter_PC1, y = jitter_PC2, color = conds))  + geom_point(size = 3, alpha = 0.8)  + theme_minimal(base_size = 14) + labs(title = "",
    x = paste0("PC1: ", percentVar[1], "%"),
     y = paste0("PC2: ", percentVar[2], "%"))

#Without jitter
ggplot(F1_Br_pca_df, aes(x = PC1, y = PC2, color = conds)) +
  geom_point(size = 3, alpha = 0.8) +
  theme_minimal(base_size = 14) +
  labs(title = "",
       x = paste0("PC1: ", percentVar[1], "%"),
       y = paste0("PC2: ", percentVar[2], "%")) 


#### Gonad High-Low-DMSO ####

HLD_Go_conds <- c("DMSO", "Low", "Low", "High", "High", "High", "DMSO", "DMSO", "Low", "High", "DMSO")

HLD_Go_metadata <- data.frame(row.names = colnames(H_L_D_Gonad),
                             conds = factor(HLD_Go_conds))
HLD_Go_conds_dds <- DESeqDataSetFromMatrix(countData = as.matrix(H_L_D_Gonad),
                                          colData = HLD_Go_metadata,
                                          design = ~ conds)

keep <- rowSums(cpm(counts(HLD_Go_conds_dds)) >= 1) >= 6
HLD_Go_conds_dds <- HLD_Go_conds_dds[keep,]

HLD_Go_conds_dds <- DESeq(HLD_Go_conds_dds)

#Pairwise
#High vs Low
Go_res_HL <- results(HLD_Go_conds_dds, contrast = c("conds", "High", "Low"))
# High vs DMSO
Go_res_HD <- results(HLD_Go_conds_dds, contrast = c("conds", "High", "DMSO"))
# Low vs DMSO
Go_res_LD <- results(HLD_Go_conds_dds, contrast = c("conds", "Low", "DMSO"))

#remove na values
Go_res_HL <- na.omit(Go_res_HL)
Go_res_HD <- na.omit(Go_res_HD)
Go_res_LD <- na.omit(Go_res_LD)
#Significant p-values
Go_res_s_HL <- Go_res_HL[Go_res_HL$padj <= 0.05 & (Go_res_HL$log2FoldChange > 1 | Go_res_HL$log2FoldChange < -1),]
dim(Go_res_s_HL) #4507 significant genes
#5032 genes with a p-value less than 0.05
Go_res_s_HD <- Go_res_HD[Go_res_HD$padj <= 0.05 & (Go_res_HD$log2FoldChange > 1 | Go_res_HD$log2FoldChange < -1),]
dim(Go_res_s_HD) #4820 significant genes
#5184 genes with a p-value less than 0.05
Go_res_s_LD <- Go_res_LD[Go_res_LD$padj <= 0.05 & (Go_res_LD$log2FoldChange >1 | Go_res_LD$log2FoldChange < -1),]
dim(Go_res_s_LD) # 2 significant gene
Go_res_s_LD
# 3 genes with a p-value above 0.05

#PCA

#plotPCA(rlog(HLD_Go_conds_dds, blind=TRUE), intgroup="conds")

HLD_Go_PCA <- rlog(HLD_Go_conds_dds, blind = TRUE)

HLD_Go_rv <- rowVars(assay(HLD_Go_PCA))
HLD_Go_top_500 <- order(HLD_Go_rv, decreasing = TRUE)[1:500]
HLD_Go_top_matrix <- assay(HLD_Go_PCA)[HLD_Go_top_500, ]

HLD_Go_pca_res <- prcomp(t(HLD_Go_top_matrix), scale. = TRUE)

HLD_Go_metadata <- as.data.frame(colData(HLD_Go_conds_dds))


HLD_Go_pca_df <- as.data.frame(HLD_Go_pca_res$x)

HLD_Go_pca_df$conds <- HLD_Go_metadata$conds 

percentVar <- round(100 * (HLD_Go_pca_res$sdev)^2 / sum(HLD_Go_pca_res$sdev^2), 1)

#with jitter 
set.seed(9)
HLD_Go_pca_df$jitter_PC1 <- HLD_Go_pca_df$PC1 + runif(nrow(HLD_Go_pca_df), -0.3, 0.3)
HLD_Go_pca_df$jitter_PC2 <- HLD_Go_pca_df$PC2 + runif(nrow(HLD_Go_pca_df), -0.3, 0.3)

ggplot(HLD_Go_pca_df, aes(x = jitter_PC1, y = jitter_PC2, color = conds))  + geom_point(size = 3, alpha = 0.8)  + theme_minimal(base_size = 14) + labs(title = "",
   x = paste0("PC1: ", percentVar[1], "%"),
   y = paste0("PC2: ", percentVar[2], "%"))

#Without jitter
ggplot(HLD_Go_pca_df, aes(x = PC1, y = PC2, color = conds)) +
  geom_point(size = 3, alpha = 0.8) +
  theme_minimal(base_size = 14) +
  labs(title = "",
       x = paste0("PC1: ", percentVar[1], "%"),
       y = paste0("PC2: ", percentVar[2], "%")) 


#### Gonad Microplastics ####

#conditions
MP_Go_conds <- c("DMSO", "MP", "MPD", "MP", "MP", "MPD", "MPD", "DMSO", "DMSO", "MP", "DMSO")

MP_Go_metadata <- data.frame(row.names = colnames(MP_MPD_D_Gonad),
                             conds = factor(MP_Go_conds))
MP_Go_conds_dds <- DESeqDataSetFromMatrix(countData = as.matrix(MP_MPD_D_Gonad),
                                          colData = MP_Go_metadata,
                                          design = ~ conds)

keep <- rowSums(cpm(counts(MP_Go_conds_dds)) >= 1) >= 6
MP_Go_conds_dds <- MP_Go_conds_dds[keep,]

MP_Go_conds_dds <- DESeq(MP_Go_conds_dds)

#Pairwise
#MP vs MPDEHP
Go_res_MP_MPD <- results(MP_Go_conds_dds, contrast = c("conds", "MP", "MPD"))
# MP vs DMSO
Go_res_MP_DMSO <- results(MP_Go_conds_dds, contrast = c("conds", "MP", "DMSO"))
# MPDEHP vs DMSO
Go_res_MPD_DMSO <- results(MP_Go_conds_dds, contrast = c("conds", "MPD", "DMSO"))

#remove na values
Go_res_MP_MPD <- na.omit(Go_res_MP_MPD)
Go_res_MP_DMSO <- na.omit(Go_res_MP_DMSO)
Go_res_MPD_DMSO <- na.omit(Go_res_MPD_DMSO)
#Significant p-values
Go_res_s_MP_MPD <- Go_res_MP_MPD[Go_res_MP_MPD$padj <= 0.05 & (Go_res_MP_MPD$log2FoldChange > 1 | Go_res_MP_MPD$log2FoldChange < -1),]
dim(Go_res_s_MP_MPD) #1478 significant genes
#1517 genes with a p-value less than 0.05
Go_res_s_MP_DMSO <-Go_res_MP_DMSO[Go_res_MP_DMSO$padj <= 0.05 & (Go_res_MP_DMSO$log2FoldChange >1 | Go_res_MP_DMSO$log2FoldChange < -1), ]
dim(Go_res_s_MP_DMSO) #5388 significant genes
#5758 genes with a p-value less than 0.05
Go_res_s_MPD_DMSO <- Go_res_MPD_DMSO[Go_res_MPD_DMSO$padj <= 0.05 & (Go_res_MPD_DMSO$log2FoldChange > 1| Go_res_MPD_DMSO$log2FoldChange < -1), ]
dim(Go_res_s_MPD_DMSO) # 9 significant genes
Go_res_s_MPD_DMSO

#PCA
#plotPCA(rlog(MP_Go_conds_dds, blind=TRUE), intgroup="conds")

MP_Go_PCA <- rlog(MP_Go_conds_dds, blind = TRUE)

MP_Go_rv <- rowVars(assay(MP_Go_PCA))
MP_Go_top_500 <- order(MP_Go_rv, decreasing = TRUE)[1:500]
MP_Go_top_matrix <- assay(MP_Go_PCA)[MP_Go_top_500, ]

MP_Go_pca_res <- prcomp(t(MP_Go_top_matrix), scale. = TRUE)

MP_Go_metadata <- as.data.frame(colData(MP_Go_conds_dds))
MP_Go_pca_df <- as.data.frame(MP_Go_pca_res$x)

MP_Go_pca_df$conds <- MP_Go_metadata$conds 

percentVar <- round(100 * (MP_Go_pca_res$sdev)^2 / sum(MP_Go_pca_res$sdev^2), 1)

#with jitter 
set.seed(9)
MP_Go_pca_df$jitter_PC1 <- MP_Go_pca_df$PC1 + runif(nrow(MP_Go_pca_df), -0.3, 0.3)
MP_Go_pca_df$jitter_PC2 <- MP_Go_pca_df$PC2 + runif(nrow(MP_Go_pca_df), -0.3, 0.3)

ggplot(MP_Go_pca_df, aes(x = jitter_PC1, y = jitter_PC2, color = conds))  + geom_point(size = 3, alpha = 0.8)  + theme_minimal(base_size = 14) + labs(title = "",
    x = paste0("PC1: ", percentVar[1], "%"),
    y = paste0("PC2: ", percentVar[2], "%"))

#Without jitter
ggplot(MP_Go_pca_df, aes(x = PC1, y = PC2, color = conds)) +
  geom_point(size = 3, alpha = 0.8) +
  theme_minimal(base_size = 14) +
  labs(title = "",
       x = paste0("PC1: ", percentVar[1], "%"),
       y = paste0("PC2: ", percentVar[2], "%")) 

#### Gonad F1 ####

#conditions
F1_Go_conds <- c("DMSO", "DMSO", "DMSO", "High", "High", "High")

F1_Go_metadata <- data.frame(row.names = colnames(F1_Gonad),
                             conds = factor(F1_Go_conds))
F1_Go_conds_dds <- DESeqDataSetFromMatrix(countData = as.matrix(F1_Gonad),
                                          colData = F1_Go_metadata,
                                          design = ~ conds)

keep <- rowSums(cpm(counts(F1_Go_conds_dds)) >= 1) >= 6
F1_Go_conds_dds <- F1_Go_conds_dds[keep,]

F1_Go_conds_dds <- DESeq(F1_Go_conds_dds)



#results
Go_res_F1 <- results(F1_Go_conds_dds, contrast = c("conds", "High", "DMSO"))
#remove na values
Go_res_F1 <- na.omit(Go_res_F1)
#Significant p-values
Go_res_s_F1 <- Go_res_F1[Go_res_F1$padj <= 0.05 & (Go_res_F1$log2FoldChange > 1 | Go_res_F1$log2FoldChange < -1), ]
dim(Go_res_s_F1) #60 significant genes
#71 genes with a p-value less than 0.05


#PCA
#plotPCA(rlog(F1_Go_conds_dds, blind=TRUE), intgroup="conds")

F1_Go_PCA <- rlog(F1_Go_conds_dds, blind = TRUE)

F1_Go_rv <- rowVars(assay(F1_Go_PCA))
F1_Go_top_500 <- order(F1_Go_rv, decreasing = TRUE)[1:500]
F1_Go_top_matrix <- assay(F1_Go_PCA)[F1_Go_top_500, ]

F1_Go_pca_res <- prcomp(t(F1_Go_top_matrix), scale. = TRUE)

F1_Go_metadata <- as.data.frame(colData(F1_Go_conds_dds))
F1_Go_pca_df <- as.data.frame(F1_Go_pca_res$x)

F1_Go_pca_df$conds <- F1_Go_metadata$conds 

percentVar <- round(100 * (F1_Go_pca_res$sdev)^2 / sum(F1_Go_pca_res$sdev^2), 1)

#with jitter 
set.seed(9)
F1_Go_pca_df$jitter_PC1 <- F1_Go_pca_df$PC1 + runif(nrow(F1_Go_pca_df), -0.3, 0.3)
F1_Go_pca_df$jitter_PC2 <- F1_Go_pca_df$PC2 + runif(nrow(F1_Go_pca_df), -0.3, 0.3)

ggplot(F1_Go_pca_df, aes(x = jitter_PC1, y = jitter_PC2, color = conds))  + geom_point(size = 3, alpha = 0.8)  + theme_minimal(base_size = 14) + labs(title = "",
      x = paste0("PC1: ", percentVar[1], "%"),
      y = paste0("PC2: ", percentVar[2], "%"))

#Without jitter
ggplot(F1_Go_pca_df, aes(x = PC1, y = PC2, color = conds)) +
  geom_point(size = 3, alpha = 0.8) +
  theme_minimal(base_size = 14) +
  labs(title = "",
       x = paste0("PC1: ", percentVar[1], "%"),
       y = paste0("PC2: ", percentVar[2], "%")) 


#### Liver High-Low-DMSO ####

HLD_Li_conds <- c("Low", "Low", "Low", "High", "DMSO", "DMSO", "DMSO", "High", "High", "High", "DMSO")

HLD_Li_metadata <- data.frame(row.names = colnames(H_L_D_Liver),
                              conds = factor(HLD_Li_conds))
HLD_Li_conds_dds <- DESeqDataSetFromMatrix(countData = as.matrix(H_L_D_Liver),
                                           colData = HLD_Li_metadata,
                                           design = ~ conds)

keep <- rowSums(cpm(counts(HLD_Li_conds_dds)) >= 1) >= 6
HLD_Li_conds_dds <- HLD_Li_conds_dds[keep,]

HLD_Li_conds_dds <- DESeq(HLD_Li_conds_dds)

#Pairwise
#High vs Low
Li_res_HL <- results(HLD_Li_conds_dds, contrast = c("conds", "High", "Low"))
# High vs DMSO
Li_res_HD <- results(HLD_Li_conds_dds, contrast = c("conds", "High", "DMSO"))
# Low vs DMSO
Li_res_LD <- results(HLD_Li_conds_dds, contrast = c("conds", "Low", "DMSO"))

#remove na values
Li_res_HL <- na.omit(Li_res_HL)
Li_res_HD <- na.omit(Li_res_HD)
Li_res_LD <- na.omit(Li_res_LD)
#Significant p-values
Li_res_s_HL <- Li_res_HL[Li_res_HL$padj <= 0.05 & (Li_res_HL$log2FoldChange > 1 | Li_res_HL$log2FoldChange < -1), ]
dim(Li_res_s_HL) #4 significant genes
Li_res_s_HD <-Li_res_HD[Li_res_HD$padj <= 0.05 & (Li_res_HD$log2FoldChange > 1 | Li_res_HD$log2FoldChange < -1), ]
dim(Li_res_s_HD) #2 significant genes
Li_res_s_HD
Li_res_s_LD <- Li_res_LD[Li_res_LD$padj <= 0.05 & (Li_res_LD$log2FoldChange > 1 | Li_res_LD$log2FoldChange < -1), ]
dim(Li_res_s_LD) #414 significant genes
# 415 genes with a p-value less than 0.05


#PCA
plotPCA(rlog(HLD_Li_conds_dds, blind=TRUE), intgroup="conds")


HLD_Li_PCA <- rlog(HLD_Li_conds_dds, blind = TRUE)

# Step 2: Extract top 500 most variable genes
HLD_Li_rv <- rowVars(assay(HLD_Li_PCA))
HLD_Li_top_500 <- order(HLD_Li_rv, decreasing = TRUE)[1:500]
HLD_Li_top_matrix <- assay(HLD_Li_PCA)[HLD_Li_top_500, ]

HLD_Li_pca_res <- prcomp(t(HLD_Li_top_matrix), scale. = TRUE)

HLD_Li_metadata <- as.data.frame(colData(HLD_Li_conds_dds))

HLD_Li_pca_df <- as.data.frame(HLD_Li_pca_res$x)

HLD_Li_pca_df$conds <- HLD_Li_metadata$conds 

percentVar <- round(100 * (HLD_Li_pca_res$sdev)^2 / sum(HLD_Li_pca_res$sdev^2), 1)

#with jitter 
set.seed(9)
HLD_Li_pca_df$jitter_PC1 <- HLD_Li_pca_df$PC1 + runif(nrow(HLD_Li_pca_df), -0.3, 0.3)
HLD_Li_pca_df$jitter_PC2 <- HLD_Li_pca_df$PC2 + runif(nrow(HLD_Li_pca_df), -0.3, 0.3)

ggplot(HLD_Li_pca_df, aes(x = jitter_PC1, y = jitter_PC2, color = conds))  + geom_point(size = 3, alpha = 0.8)  + theme_minimal(base_size = 14) + labs(title = "",
     x = paste0("PC1: ", percentVar[1], "%"),
     y = paste0("PC2: ", percentVar[2], "%"))

#Without jitter
ggplot(HLD_Li_pca_df, aes(x = PC1, y = PC2, color = conds)) +
  geom_point(size = 3, alpha = 0.8) +
  theme_minimal(base_size = 14) +
  labs(title = "",
       x = paste0("PC1: ", percentVar[1], "%"),
       y = paste0("PC2: ", percentVar[2], "%")) 

#### Liver Microplastics ####

MP_Li_conds <- c("MP", "MP", "MPD", "DMSO", "DMSO", "DMSO", "MP", "MPD", "MPD", "MP", "DMSO")

MP_Li_metadata <- data.frame(row.names = colnames(MP_MPD_D_Liver),
                             conds = factor(MP_Li_conds))
MP_Li_conds_dds <- DESeqDataSetFromMatrix(countData = as.matrix(MP_MPD_D_Liver),
                                          colData = MP_Li_metadata,
                                          design = ~ conds)

keep <- rowSums(cpm(counts(MP_Li_conds_dds)) >= 1) >= 6
MP_Li_conds_dds <- MP_Li_conds_dds[keep,]

MP_Li_conds_dds <- DESeq(MP_Li_conds_dds)

#Pairwise
#MP vs MPDEHP
Li_res_MP_MPD <- results(MP_Li_conds_dds, contrast = c("conds", "MP", "MPD"))
# MP vs DMSO
Li_res_MP_DMSO <- results(MP_Li_conds_dds, contrast = c("conds", "MP", "DMSO"))
# MPDEHP vs DMSO
Li_res_MPD_DMSO <- results(MP_Li_conds_dds, contrast = c("conds", "MPD", "DMSO"))

#remove na values
Li_res_MP_MPD <- na.omit(Li_res_MP_MPD)
Li_res_MP_DMSO <- na.omit(Li_res_MP_DMSO)
Li_res_MPD_DMSO <- na.omit(Li_res_MPD_DMSO)
#Significant p-values
Li_res_s_MP_MPD <- Li_res_MP_MPD[Li_res_MP_MPD$padj <= 0.05 & (Li_res_MP_MPD$log2FoldChange > 1 | Li_res_MP_MPD$log2FoldChange < -1), ]
dim(Li_res_s_MP_MPD) #0 significant genes
Li_res_s_MP_DMSO <-Li_res_MP_DMSO[Li_res_MP_DMSO$padj <= 0.05 & (Li_res_MP_DMSO$log2FoldChange > 1 | Li_res_MP_DMSO$log2FoldChange < -1), ]
dim(Li_res_s_MP_DMSO) #0 significant genes
Li_res_s_MPD_DMSO <- Li_res_MPD_DMSO[Li_res_MPD_DMSO$padj <= 0.05 & (Li_res_MPD_DMSO$log2FoldChange >1 | Li_res_MPD_DMSO$log2FoldChange < -1), ]
dim(Li_res_s_MPD_DMSO) # 0 significant gene

#PCA
#plotPCA(rlog(MP_Li_conds_dds, blind=TRUE), intgroup="conds")

MP_Li_PCA <- rlog(MP_Li_conds_dds, blind = TRUE)

# Step 2: Extract top 500 most variable genes
MP_Li_rv <- rowVars(assay(MP_Li_PCA))
MP_Li_top_500 <- order(MP_Li_rv, decreasing = TRUE)[1:500]
MP_Li_top_matrix <- assay(MP_Li_PCA)[MP_Li_top_500, ]

MP_Li_pca_res <- prcomp(t(MP_Li_top_matrix), scale. = TRUE)

MP_Li_metadata <- as.data.frame(colData(MP_Li_conds_dds))

MP_Li_pca_df <- as.data.frame(MP_Li_pca_res$x)

MP_Li_pca_df$conds <- MP_Li_metadata$conds 

percentVar <- round(100 * (MP_Li_pca_res$sdev)^2 / sum(MP_Li_pca_res$sdev^2), 1)

#with jitter 
set.seed(9)
MP_Li_pca_df$jitter_PC1 <- MP_Li_pca_df$PC1 + runif(nrow(MP_Li_pca_df), -0.3, 0.3)
MP_Li_pca_df$jitter_PC2 <- MP_Li_pca_df$PC2 + runif(nrow(MP_Li_pca_df), -0.3, 0.3)

ggplot(MP_Li_pca_df, aes(x = jitter_PC1, y = jitter_PC2, color = conds))  + geom_point(size = 3, alpha = 0.8)  + theme_minimal(base_size = 14) + labs(title = "",
       x = paste0("PC1: ", percentVar[1], "%"),
       y = paste0("PC2: ", percentVar[2], "%"))

#Without jitter
ggplot(MP_Li_pca_df, aes(x = PC1, y = PC2, color = conds)) +
  geom_point(size = 3, alpha = 0.8) +
  theme_minimal(base_size = 14) +
  labs(title = "",
       x = paste0("PC1: ", percentVar[1], "%"),
       y = paste0("PC2: ", percentVar[2], "%")) 

#### F1 Liver ####

#conditions

F1_Li_conds <- c("DMSO", "High", "High", "High", "DMSO", "DMSO")

F1_Li_metadata <- data.frame(row.names = colnames(F1_Liver),
                             conds = factor(F1_Li_conds))
F1_Li_conds_dds <- DESeqDataSetFromMatrix(countData = as.matrix(F1_Liver),
                                          colData = F1_Li_metadata,
                                          design = ~ conds)

keep <- rowSums(cpm(counts(F1_Li_conds_dds)) >= 1) >= 6
F1_Li_conds_dds <- F1_Li_conds_dds[keep,]

F1_Li_conds_dds <- DESeq(F1_Li_conds_dds)

#results
Li_res_F1 <- results(F1_Li_conds_dds, contrast = c("conds", "High", "DMSO"))
#remove na values
Li_res_F1 <- na.omit(Li_res_F1)
#Significant p-values
Li_res_s_F1 <- Li_res_F1[Li_res_F1$padj <= 0.05 & (Li_res_F1$log2FoldChange > 1 | Li_res_F1$log2FoldChange < -1), ]
dim(Li_res_s_F1) #52 significant genes

#PCA
#plotPCA(rlog(F1_Li_conds_dds, blind=TRUE), intgroup="conds")

F1_Li_PCA <- rlog(F1_Li_conds_dds, blind = TRUE)

F1_Li_rv <- rowVars(assay(F1_Li_PCA))
F1_Li_top_500 <- order(F1_Li_rv, decreasing = TRUE)[1:500]
F1_Li_top_matrix <- assay(F1_Li_PCA)[F1_Li_top_500, ]

F1_Li_pca_res <- prcomp(t(F1_Li_top_matrix), scale. = TRUE)

F1_Li_metadata <- as.data.frame(colData(F1_Li_conds_dds))


F1_Li_pca_df <- as.data.frame(F1_Li_pca_res$x)

F1_Li_pca_df$conds <- F1_Li_metadata$conds 

percentVar <- round(100 * (F1_Li_pca_res$sdev)^2 / sum(F1_Li_pca_res$sdev^2), 1)

#with jitter 
set.seed(9)
F1_Li_pca_df$jitter_PC1 <- F1_Li_pca_df$PC1 + runif(nrow(F1_Li_pca_df), -0.3, 0.3)
F1_Li_pca_df$jitter_PC2 <- F1_Li_pca_df$PC2 + runif(nrow(F1_Li_pca_df), -0.3, 0.3)

ggplot(F1_Li_pca_df, aes(x = jitter_PC1, y = jitter_PC2, color = conds))  + geom_point(size = 3, alpha = 0.8)  + theme_minimal(base_size = 14) + labs(title = "",
      x = paste0("PC1: ", percentVar[1], "%"),
      y = paste0("PC2: ", percentVar[2], "%"))

#Without jitter
ggplot(F1_Li_pca_df, aes(x = PC1, y = PC2, color = conds)) +
  geom_point(size = 3, alpha = 0.8) +
  theme_minimal(base_size = 14) +
  labs(title = "",
       x = paste0("PC1: ", percentVar[1], "%"),
       y = paste0("PC2: ", percentVar[2], "%")) 

#Remove metadata
rm(HLD_Li_metadata, HLD_metadata, HLD_Go_metadata, MP_Br_metadata, MP_Go_metadata, MP_Li_metadata, F1_Br_metadata, F1_Go_metadata, F1_Li_metadata)

#What are our genes of interest
head(Br_res_s_HD)
head(Br_res_s_LD)

head(Br_res_s_MP_DMSO)
head(Br_res_s_MPD_DMSO)

head(Br_res_s_F1, 10)


### Only keep results - Free space
#rm(HLD_Br_conds_dds, MP_Br_conds_dds, F1_Br_conds_dds, HLD_Go_conds_dds,  MP_Go_conds_dds, F1_Go_conds_dds, HLD_Li_conds_dds, MP_Li_conds_dds, F1_Li_conds_dds, keep)

### Save the significant genes 
write.csv(as.data.frame(Br_res_s_HD), file = "Significant_genes/Br_res_s_HD.csv")
write.csv(as.data.frame(Br_res_s_LD), file = "Significant_genes/Br_res_s_LD.csv")
write.csv(as.data.frame(Br_res_s_HL), file = "Significant_genes/Br_res_s_HL.csv")
write.csv(as.data.frame(Br_res_s_MP_DMSO), file = "Significant_genes/Br_res_s_MP_DMSO.csv")
write.csv(as.data.frame(Br_res_s_MP_MPD), file = "Significant_genes/Br_res_s_MP_MPD.csv")
write.csv(as.data.frame(Br_res_s_MPD_DMSO), file = "Significant_genes/Br_res_s_MPD_DMSO.csv")
write.csv(as.data.frame(Br_res_s_F1), file = "Significant_genes/Br_res_s_F1.csv")

write.csv(as.data.frame(Go_res_s_HD), file = "Significant_genes/Go_res_s_HD.csv")
write.csv(as.data.frame(Go_res_s_LD), file = "Significant_genes/Go_res_s_LD.csv")
write.csv(as.data.frame(Go_res_s_HL), file = "Significant_genes/Go_res_s_HL.csv")
write.csv(as.data.frame(Go_res_s_MP_DMSO), file = "Significant_genes/Go_res_s_MP_DMSO.csv")
write.csv(as.data.frame(Go_res_s_MP_MPD), file = "Significant_genes/Go_res_s_MP_MPD.csv")
write.csv(as.data.frame(Go_res_s_MPD_DMSO), file = "Significant_genes/Go_res_s_MPD_DMSO.csv")
write.csv(as.data.frame(Go_res_s_F1), file = "Significant_genes/Go_res_s_F1.csv")

write.csv(as.data.frame(Li_res_s_HD), file = "Significant_genes/Li_res_s_HD.csv")
write.csv(as.data.frame(Li_res_s_LD), file = "Significant_genes/Li_res_s_LD.csv")
write.csv(as.data.frame(Li_res_s_HL), file = "Significant_genes/Li_res_s_HL.csv")
write.csv(as.data.frame(Li_res_s_MP_DMSO), file = "Significant_genes/Li_res_s_MP_DMSO.csv")
write.csv(as.data.frame(Li_res_s_MP_MPD), file = "Significant_genes/Li_res_s_MP_MPD.csv")
write.csv(as.data.frame(Li_res_s_MPD_DMSO), file = "Significant_genes/Li_res_s_MPD_DMSO.csv")
write.csv(as.data.frame(Li_res_s_F1), file = "Significant_genes/Li_res_s_F1.csv")

#All genes

# Brain data
write.csv(as.data.frame(Br_res_HD), file = "All genes/Br_res_HD.csv")
write.csv(as.data.frame(Br_res_LD), file = "All genes/Br_res_LD.csv")
write.csv(as.data.frame(Br_res_HL), file = "All genes/Br_res_HL.csv")
write.csv(as.data.frame(Br_res_MP_DMSO), file = "All genes/Br_res_MP_DMSO.csv")
write.csv(as.data.frame(Br_res_MP_MPD), file = "All genes/Br_res_MP_MPD.csv")
write.csv(as.data.frame(Br_res_MPD_DMSO), file = "All genes/Br_res_MPD_DMSO.csv")
write.csv(as.data.frame(Br_res_F1), file = "All genes/Br_res_F1.csv")

# Gonad data
write.csv(as.data.frame(Go_res_HD), file = "All genes/Go_res_HD.csv")
write.csv(as.data.frame(Go_res_LD), file = "All genes/Go_res_LD.csv")
write.csv(as.data.frame(Go_res_HL), file = "All genes/Go_res_HL.csv")
write.csv(as.data.frame(Go_res_MP_DMSO), file = "All genes/Go_res_MP_DMSO.csv")
write.csv(as.data.frame(Go_res_MP_MPD), file = "All genes/Go_res_MP_MPD.csv")
write.csv(as.data.frame(Go_res_MPD_DMSO), file = "All genes/Go_res_MPD_DMSO.csv")
write.csv(as.data.frame(Go_res_F1), file = "All genes/Go_res_F1.csv")

# Liver data
write.csv(as.data.frame(Li_res_HD), file = "All genes/Li_res_HD.csv")
write.csv(as.data.frame(Li_res_LD), file = "All genes/Li_res_LD.csv")
write.csv(as.data.frame(Li_res_HL), file = "All genes/Li_res_HL.csv")
write.csv(as.data.frame(Li_res_MP_DMSO), file = "All genes/Li_res_MP_DMSO.csv")
write.csv(as.data.frame(Li_res_MP_MPD), file = "All genes/Li_res_MP_MPD.csv")
write.csv(as.data.frame(Li_res_MPD_DMSO), file = "All genes/Li_res_MPD_DMSO.csv")
write.csv(as.data.frame(Li_res_F1), file = "All genes/Li_res_F1.csv")

#UGH

write.table(rownames(Br_res_F1), file = "All genes text/Br_res_F1.txt",quote = FALSE, row.names = FALSE, col.names = FALSE)

write.table(rownames(Go_res_HD), file = "All genes text/Go_res_HD.txt",quote = FALSE, row.names = FALSE, col.names = FALSE)

write.table(rownames(Go_res_F1), file = "All genes text/Go_res_F1.txt",quote = FALSE, row.names = FALSE, col.names = FALSE)

write.table(rownames(Li_res_LD), file = "All genes text/Li_res_LD.txt",quote = FALSE, row.names = FALSE, col.names = FALSE)

write.table(rownames(Go_res_MP_DMSO), file = "All genes text/Go_res_MP_DMSO.txt",quote = FALSE, row.names = FALSE, col.names = FALSE)
